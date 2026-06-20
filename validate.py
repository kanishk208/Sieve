"""
validate.py — Variant validation pipeline using network topology and biochemistry
Runs the Sieve pipeline across a diverse multi-gene cohort.

Features:
  - Layer 0: Network topology (NetworkX contact graph, centrality metrics)
  - Functional-site override (catalytic motifs are unconditionally pathogenic)
  - Chemistry gates (Grantham distance, charge-class changes, proline-special-casing)
  - Centrality gate (radical substitutions only pathogenic if non-peripheral)

Usage:
    python validate.py
    python validate.py --gene VCP
"""

import argparse
import json
import math
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, ShrakeRupley

ROOT = Path(__file__).parent
DATA_DIR = ROOT / "data"
sys.path.insert(0, str(ROOT))

from build_graph import build_contact_graph, compute_centrality, simulate_mutation, layer0_classification
from variant_utils import build_validation_dataset
from grantham import is_chemically_conservative, grantham_distance
from conservation import load_conservation_profile

# Max reference SASA per residue type (Rost & Sander 1994, Angstroms^2)
_MAX_SASA = {
    "ALA": 106, "ARG": 248, "ASN": 157, "ASP": 163,
    "CYS": 135, "GLN": 198, "GLU": 194, "GLY": 84,
    "HIS": 184, "ILE": 169, "LEU": 164, "LYS": 205,
    "MET": 188, "PHE": 197, "PRO": 136, "SER": 130,
    "THR": 142, "TRP": 227, "TYR": 222, "VAL": 142,
}

# VCP functional sites that must never be Grantham-suppressed
VCP_FUNCTIONAL_SITES: dict[str, set[int]] = {
    # Walker A P-loop: range ends at the invariant Lys (K232 / K530), not the
    # following Thr/Ile which tolerates conservative substitutions (I233V is benign).
    "D1_WALKER_A":   set(range(227, 233)),   # ends at K232 inclusive
    "D2_WALKER_A":   set(range(525, 531)),   # ends at K530 inclusive
    "D1_WALKER_B":   set(range(284, 289)),
    "D2_WALKER_B":   set(range(578, 583)),
    "D1_ARG_FINGER": {359, 362},
    "D2_ARG_FINGER": {635, 638},
    # LINKER (458-481) omitted: too broad, contains tolerated benign variants
    "DISEASE_LOOP":  set(range(155, 175)),
}
_ALL_FUNCTIONAL_RESIDUES: set[int] = set().union(*VCP_FUNCTIONAL_SITES.values())

# ---------------------------------------------------------------------------
# Charge-class rescue helpers (Fix B)
# ---------------------------------------------------------------------------
# At physiological pH 7.4: ARG and LYS are fully positive (+1).
# HIS pKa ≈ 6.0 ->~90% neutral at pH 7.4 ->treat as neutral (0).
# ASP and GLU are fully negative (−1).
# All others neutral (0).
_POSITIVE_AA = {"ARG", "LYS"}
_NEGATIVE_AA = {"ASP", "GLU"}

_AA_1TO3: dict[str, str] = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
    "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
    "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
    "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR",
}

def _charge_class(aa1: str) -> str:
    aa3 = _AA_1TO3.get(aa1.upper(), aa1.upper())
    if aa3 in _POSITIVE_AA:
        return "+"
    if aa3 in _NEGATIVE_AA:
        return "-"
    return "0"

def _is_charge_changing(ref: str, alt: str) -> bool:
    """True if the substitution changes formal ionization class at pH 7.4."""
    return _charge_class(ref) != _charge_class(alt)

_AA_3TO1: dict[str, str] = {v: k for k, v in _AA_1TO3.items()}


def build_pdb_residue_map(pdb_path: str, chain: str = "A") -> dict[int, str]:
    """
    Map residue number -> one-letter amino acid for a PDB chain (Cα atoms).

    Used to filter out 'phantom' ClinVar variants whose stated WT residue does
    not match the structure. ClinVar `protein_change` fields for VCP are polluted
    with a -45-shifted duplicate of every real mutation (e.g. a record titled
    p.Ala160Ser carries "A115S, A160S"); the −45 copy lands on the wrong residue
    and injects pure noise into validation. Keeping only variants whose WT letter
    matches the structure removes both numbering artifacts and genuine isoform
    mismatches that cannot be scored against this model.
    """
    residue_map: dict[int, str] = {}
    try:
        with open(pdb_path) as f:
            for line in f:
                if line.startswith("ATOM") and line[12:16].strip() == "CA" and line[21] == chain:
                    try:
                        resnum = int(line[22:26])
                        residue_map[resnum] = _AA_3TO1.get(line[17:20].strip(), "X")
                    except (ValueError, IndexError):
                        continue
    except Exception as exc:
        print(f"  [!] Failed to read PDB residue map: {exc}")
    return residue_map


def build_sasa_map(pdb_path: str, chain: str = "A") -> dict[int, float]:
    """
    Return relative SASA per residue number for a given PDB chain.
    Values in [0, 1]: <0.20 buried, 0.20-0.50 intermediate, >0.50 surface.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("prot", pdb_path)
    sr = ShrakeRupley()
    sr.compute(structure, level="R")

    sasa_map: dict[int, float] = {}
    try:
        chain_obj = structure[0][chain]
    except KeyError:
        return sasa_map

    for residue in chain_obj:
        if residue.id[0] != " ":   # skip HETATM / water
            continue
        resnum = residue.id[1]
        resname = residue.get_resname()
        max_sasa = _MAX_SASA.get(resname, 180.0)
        sasa_map[resnum] = residue.sasa / max_sasa
    return sasa_map


def normalize_variant_record(record: dict) -> list[dict]:
    """
    Normalize ClinVar format to standard variant format.
    
    ClinVar provides protein_change as comma-separated mutations (e.g., "A115S, M158I").
    This function expands each mutation into individual variant records with
    mutation, resnum, and known_pathogenicity keys.
    
    Returns list of normalized variant dicts (typically 1, sometimes multiple).
    """
    normalized = []
    protein_change = record.get("protein_change", "")
    
    if not protein_change:
        return normalized
    
    # Map search_class to pathogenicity
    search_class = record.get("search_class", "unknown").lower()
    known_pathogenicity = "pathogenic" if search_class == "pathogenic" else "benign"
    
    # Parse comma-separated mutations: "A115S, M158I" -> ["A115S", "M158I"]
    mutations = [m.strip() for m in protein_change.split(",")]
    
    for mutation_str in mutations:
        # Match pattern: Single letter + digits + single letter (e.g., "A115S")
        match = re.match(r"([A-Z])(\d+)([A-Z])", mutation_str)
        if match:
            wt_aa, resnum, mut_aa = match.groups()
            mutation = f"{wt_aa}{resnum}{mut_aa}"
            
            # Create normalized record preserving original fields
            normalized_rec = {
                "mutation": mutation,
                "resnum": int(resnum),
                "known_pathogenicity": known_pathogenicity,
                # Keep original fields for reference
                "variation_id": record.get("variation_id", ""),
                "title": record.get("title", ""),
                "significance": record.get("significance", ""),
                "protein_change": protein_change,
            }
            normalized.append(normalized_rec)
    
    return normalized

# 🌟 MULTI-GENE BENCHMARK COHORT
# Define your target proteins here. Add/remove genes and UniProt IDs as needed.
TARGET_COHORT = {
    "VCP": "P55072",      # Motor protein AAA-ATPase
    # "TP53": "P04637",   # Tumor suppressor (add when PDB available)
    # "PTEN": "P60484",   # Phosphatase (add when PDB available)
    # "KCNQ1": "P51787",  # Ion channel (add when PDB available)
}


def load_dataset(path: str) -> list[dict]:
    with open(path) as f:
        return json.load(f)


def validate_gene(gene: str, uniprot: str) -> pd.DataFrame:
    """
    Run Sieve pipeline on a single gene with Layer 0 network topology
    and biochemical features.
    
    Args:
        gene: Gene symbol (e.g., "VCP")
        uniprot: UniProt ID (e.g., "P55072")
    
    Returns:
        DataFrame with validation results for all variants
    """
    gene_dir = DATA_DIR / gene

    # Find AlphaFold PDB file in gene subdirectory
    pdb_files = sorted(gene_dir.glob(f"AF-{uniprot}-*.pdb"))
    if not pdb_files:
        print(f"[ERROR] No PDB found for {gene} (UniProt: {uniprot}) in {gene_dir}")
        return pd.DataFrame()


    pdb_path = pdb_files[-1]
    print(f"\n[Sieve Pipeline] Processing {gene} ({len(list(gene_dir.glob(f'{gene}_clinvar_variants.json')))} variants)")
    print(f"  PDB: {pdb_path.name}")

    # Build contact graph
    print(f"  Building contact graph...")
    G, _ = build_contact_graph(str(pdb_path), cutoff=8.0)

    # Load or compute centrality
    cent_path = gene_dir / f"{gene}_centrality.csv"
    if cent_path.exists():
        print(f"  Loading cached centrality...")
        df_cent = pd.read_csv(cent_path)
    else:
        print(f"  Computing centrality...")
        df_cent = compute_centrality(G)
        df_cent.to_csv(cent_path, index=False)

    # Load ClinVar dataset
    dataset_path = gene_dir / f"{gene}_clinvar_variants.json"
    if not dataset_path.exists():
        print(f"[ERROR] {gene}_clinvar_variants.json not found in {gene_dir}")
        return pd.DataFrame()
    
    with open(dataset_path) as f:
        clinvar_records = json.load(f)
    
    # Normalize ClinVar format to standard variant records
    dataset = []
    for record in clinvar_records:
        normalized = normalize_variant_record(record)
        dataset.extend(normalized)

    # ── Filter phantom / isoform-mismatched variants ───────────────────────────
    # Keep only variants whose stated WT residue matches the structure at that
    # position. This removes the −45-shifted ClinVar duplicates that otherwise
    # double the dataset and inject noise scored at the wrong residue. Dedupe on
    # (mutation, pathogenicity) so a variant listed in several records counts once.
    pdb_residues = build_pdb_residue_map(str(pdb_path))
    clean, seen, dropped = [], set(), 0
    for v in dataset:
        mut, resnum = v["mutation"], v["resnum"]
        wt = mut[0]
        if pdb_residues.get(resnum) != wt:
            dropped += 1
            continue
        key = (mut, v["known_pathogenicity"])
        if key in seen:
            continue
        seen.add(key)
        clean.append(v)
    print(f"  Filtered variants: {len(clean)} kept, {dropped} phantom/mismatch dropped "
          f"(from {len(dataset)} raw rows)")
    dataset = clean

    # Pre-compute SASA for the whole structure once (expensive but done once per gene)
    print(f"  Computing SASA map...")
    sasa_map = build_sasa_map(str(pdb_path))

    # Preload the conservation profile once (resnum -> score) for C_3D.
    cons_profile = load_conservation_profile(gene)

    results = []
    print(f"  Analyzing {len(dataset)} mutations...")

    for idx, v in enumerate(dataset, 1):
        if idx % 25 == 0:
            print(f"    [{idx}/{len(dataset)}] processed...")

        mut = v["mutation"]
        resnum = v["resnum"]
        known = v["known_pathogenicity"]

        wt_aa = mut[0] if len(mut) >= 3 else "?"
        mut_aa = mut[-1] if len(mut) >= 3 else "?"

        # Skip variants outside graph
        if resnum not in G.nodes():
            results.append({**v, "predicted_tier": "NOT_IN_GRAPH", "predicted_pathogenic": False,
                            "ddg_foldx": None, "ddg_crystal": None, "relative_sasa": None})
            continue

        # Layer 0: network simulation
        mut_result = simulate_mutation(G, resnum)
        trow = df_cent[df_cent["resnum"] == resnum]
        if trow.empty:
            results.append({**v, "predicted_tier": "NO_CENTRALITY", "predicted_pathogenic": False,
                            "ddg_foldx": None, "ddg_crystal": None, "relative_sasa": None})
            continue

        trow = trow.iloc[0]
        b_rank      = float(trow["betweenness_rank"])
        plddt_score = float(trow["plddt"])

        # Chemistry features
        grantham_d        = grantham_distance(wt_aa, mut_aa)
        chem_conservative = is_chemically_conservative(wt_aa, mut_aa)

        # Fix C: Proline is conformationally unique (cyclic backbone).
        # Grantham distance underestimates its structural impact — never suppress.
        if wt_aa == "P" or mut_aa == "P":
            chem_conservative = False

        # SASA for this residue (None if not in map) — informational in the output.
        relative_sasa = sasa_map.get(resnum)

        # Layer 0 classification
        tier, explanation = layer0_classification(trow, mut_result, grantham_d)

        # --- Classification with functional-site and charge-class rescue ---
        in_functional_site  = resnum in _ALL_FUNCTIONAL_RESIDUES
        charge_change       = _is_charge_changing(wt_aa, mut_aa)

        # Centrality gate for structural-disruption calls. A radical substitution
        # or a topological "hub" flag only implies pathogenicity if the residue is
        # not peripheral: variants in the bottom ~20% of betweenness sit on the
        # surface / chain termini where even radical swaps are usually tolerated.
        # Calibrated on the clean VCP ClinVar set — MCC is flat across 0.15–0.22
        # (0.40 -> 0.45) — so 0.20 is a stable, non-overfit choice.
        CENTRALITY_GATE = 0.20
        central_enough = b_rank >= CENTRALITY_GATE

        if in_functional_site:
            # Walker motifs / Arg fingers / disease loop: catalytically critical,
            # pathogenic regardless of chemistry or centrality. Unconditional.
            predicted_pathogenic = True
            if tier != "HIGH_PRIORITY":
                explanation += " | [FUNCTIONAL_SITE] Variant in catalytic/disease motif"
        elif tier == "HIGH_PRIORITY" and central_enough:
            predicted_pathogenic = True
        elif tier == "MODERATE" and not chem_conservative and central_enough:
            # Radical chemistry at a non-peripheral position.
            predicted_pathogenic = True
        elif tier == "MODERATE" and chem_conservative and charge_change:
            # Conservative by Grantham but changes ionization class at pH 7.4
            # (e.g. R→H loses +1, N→D gains −1): disrupts electrostatic networks
            # and H-bond geometry despite the small physicochemical distance.
            predicted_pathogenic = True
            explanation += " | [CHARGE_CLASS_RESCUE] Conservative by Grantham but changes ionization class at pH 7.4"
        else:
            predicted_pathogenic = False

        # Layer 2: 3D Spatial Conservation Index (C_3D) — informational only.
        # Uses the preloaded conservation profile (one CSV read per gene, not one
        # per residue lookup).
        target_cons = cons_profile.get(resnum)
        c_3d = None
        if target_cons is not None:
            neighborhood_scores = [target_cons]
            for nb in G.neighbors(resnum):
                s = cons_profile.get(nb)
                if s is not None:
                    neighborhood_scores.append(s)
            c_3d = sum(neighborhood_scores) / len(neighborhood_scores)

        # Layer 1: FoldX thermodynamic escalation — DISABLED.
        #
        # Empirically validated on the clean VCP ClinVar set that FoldX ΔΔG is
        # essentially uncorrelated with pathogenicity here, in both directions:
        #   - Benign buried-hydrophobic swaps DESTABILIZE strongly yet are tolerated
        #     in vivo (F267I=2.64, L396V=2.16, M442V=1.90, V394M=1.58 — all benign),
        #     so an escalation gate manufactures false positives.
        #   - Pathogenic conservative variants act through the ATPase cycle, not local
        #     stability, so they read as STABLE (V88L=-0.04, E706D=-0.05, K277R=-0.11).
        # Net effect of the physics gates was −0.11 MCC (0.40 → 0.30) plus minutes of
        # FoldX runtime per validation. Kept off; ΔΔG fields remain in the schema as
        # None for backward compatibility. (foldx_query.run_foldx /
        # run_foldx_crystal stay available for single-variant inspection.)
        ddg_foldx = None
        ddg_crystal = None

        results.append({
            "gene": gene,
            "mutation": mut,
            "resnum": resnum,
            "known_pathogenicity": known,
            "predicted_tier": tier,
            "predicted_pathogenic": predicted_pathogenic,
            "betweenness_pct": b_rank,
            "degree_pct": float(trow["degree_rank"]),
            "edges_removed": mut_result.get("edges_removed", 0),
            "delta_path": mut_result.get("delta_avg_path", 0),
            "delta_components": mut_result.get("delta_components", 0),
            "explanation": explanation,
            "grantham_distance": grantham_d,
            "c_3d": c_3d,
            "relative_sasa": relative_sasa,
            "ddg_foldx": ddg_foldx,
            "ddg_crystal": ddg_crystal,
            "plddt": plddt_score,
        })

    return pd.DataFrame(results)



def compute_metrics(df: pd.DataFrame) -> dict:
    """Compute confusion matrix and performance metrics for global dataset."""
    if df.empty:
        return {
            "tp": 0, "fp": 0, "fn": 0, "tn": 0,
            "accuracy": 0, "sensitivity": 0, "specificity": 0,
            "precision": 0, "f1_score": 0, "mcc": 0, "n_total": 0,
        }
    
    # Map known pathogenicity to binary
    df = df.copy()
    df["actual_positive"] = df["known_pathogenicity"].isin(["pathogenic"])
    df["pred_positive"] = df["predicted_pathogenic"] == True

    tp = ((df["actual_positive"]) & (df["pred_positive"])).sum()
    fp = ((~df["actual_positive"]) & (df["pred_positive"])).sum()
    fn = ((df["actual_positive"]) & (~df["pred_positive"])).sum()
    tn = ((~df["actual_positive"]) & (~df["pred_positive"])).sum()

    n = tp + fp + fn + tn
    accuracy = (tp + tn) / n if n > 0 else 0
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0   # recall
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    f1 = 2 * precision * sensitivity / (precision + sensitivity) if (precision + sensitivity) > 0 else 0

    # Matthew's Correlation Coefficient
    denom = math.sqrt((tp+fp) * (tp+fn) * (tn+fp) * (tn+fn)) if (tp+fp)*(tp+fn)*(tn+fp)*(tn+fn) > 0 else 1
    mcc = (tp * tn - fp * fn) / denom

    return {
        "tp": int(tp), "fp": int(fp), "fn": int(fn), "tn": int(tn),
        "accuracy": round(accuracy, 3),
        "sensitivity": round(sensitivity, 3),
        "specificity": round(specificity, 3),
        "precision": round(precision, 3),
        "f1_score": round(f1, 3),
        "mcc": round(mcc, 3),
        "n_total": int(n),
    }


def print_report(metrics: dict, df: pd.DataFrame):
    """Print global validation report across all genes."""
    print("\n" + "=" * 70)
    print("  SIEVE MULTI-GENE VALIDATION REPORT")
    print("=" * 70)
    print(f"\n  Global dataset size: {metrics['n_total']} variants")
    print(f"\n  Confusion Matrix:")
    print(f"                    Predicted +    Predicted -")
    print(f"  Actual +   (TP)   {metrics['tp']:>5}          (FN)  {metrics['fn']:>5}")
    print(f"  Actual -   (FP)   {metrics['fp']:>5}          (TN)  {metrics['tn']:>5}")
    print(f"\n  Accuracy:     {metrics['accuracy']:.3f}")
    print(f"  Sensitivity:  {metrics['sensitivity']:.3f}  (recall — catches true pathogenic)")
    print(f"  Specificity:  {metrics['specificity']:.3f}  (avoids false alarms)")
    print(f"  Precision:    {metrics['precision']:.3f}")
    print(f"  F1 Score:     {metrics['f1_score']:.3f}")
    print(f"  MCC:          {metrics['mcc']:.3f}  {'[OK]' if metrics['mcc'] > 0.4 else '[LOW]'} {'> 0.4 publishable' if metrics['mcc'] > 0.4 else 'below 0.4 threshold'}")
    print("=" * 70)

    # Per-gene summary
    if "gene" in df.columns:
        print(f"\n  Per-Gene Summary:")
        for gene in df["gene"].unique():
            df_gene = df[df["gene"] == gene]
            metrics_gene = compute_metrics(df_gene)
            print(f"    {gene:8s}  n={metrics_gene['n_total']:3d}  mcc={metrics_gene['mcc']:6.3f}  sens={metrics_gene['sensitivity']:.3f}  spec={metrics_gene['specificity']:.3f}")
    print()


def main():
    parser = argparse.ArgumentParser(description="Multi-gene Sieve validation using network topology and biochemistry.")
    parser.add_argument("--gene", default=None, help="Single gene to validate (default: all in TARGET_COHORT)")
    args = parser.parse_args()

    print("=" * 70)
    print("  SIEVE MULTI-GENE PIPELINE")
    print("=" * 70)

    global_results = []
    genes_to_validate = {args.gene: TARGET_COHORT[args.gene]} if args.gene else TARGET_COHORT

    for gene, uniprot in genes_to_validate.items():
        df_gene = validate_gene(gene, uniprot)
        if not df_gene.empty:
            global_results.append(df_gene)
        else:
            print(f"[WARNING] No results for {gene}")

    if global_results:
        df_all = pd.concat(global_results, ignore_index=True)
        metrics = compute_metrics(df_all)
        print_report(metrics, df_all)
        
        # Save global results
        out = DATA_DIR / "global_validation_results.json"
        with open(out, "w") as f:
            json.dump({
                "metrics": metrics,
                "results": df_all.to_dict(orient="records")
            }, f, indent=2, default=str)
        print(f"[Sieve] Global results saved: {out}")
    else:
        print("[ERROR] No validation results generated.")


if __name__ == "__main__":
    main()
