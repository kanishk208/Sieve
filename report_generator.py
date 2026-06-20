"""
report_generator.py — Phase 3, Step 15: Clinical Report Assembler
Reads the outputs from all prior pipeline stages (pLDDT audit, Layer 0 graph,
Layer 1 DynaMut2) and assembles a structured clinical interpretation report.

The report is produced in two formats:
  1. A human-readable Markdown file for the clinician
  2. A machine-readable JSON file for downstream integration

Output classification follows the Five-Tier system from Section 6.1 of the
architecture document.

Usage:
    python report_generator.py --gene VCP --mutation R155H
    python report_generator.py --gene VCP --mutation R155H --format both
"""

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from typing import Optional

DATA_DIR    = os.path.join(os.path.dirname(__file__), "data")
REPORTS_DIR = os.path.join(os.path.dirname(__file__), "reports")
os.makedirs(REPORTS_DIR, exist_ok=True)


# ─── Five-Tier Classification Engine ─────────────────────────────────────────
# Implements Section 6.1 of the architecture document exactly.

TIERS = {
    1: ("Likely pathogenic — structural",   "🔴"),
    2: ("Likely pathogenic — functional",   "🔴"),
    3: ("Uncertain significance — destabilising signal", "🟡"),
    4: ("Uncertain significance — no signal",            "🟡"),
    5: ("Insufficient data",                             "⚪"),
}

def classify_tier(
    *,
    plddt: Optional[float],
    layer0_tier: Optional[str],
    ddg: Optional[float],
    betweenness_pct: Optional[float],
    idr_evidence: bool = False,
) -> tuple[int, str]:
    """
    Phase 3 five-tier logic with automated ΔΔG triage (Technical Reference §6.1).
    """
    plddt = plddt if plddt is not None else 0.0
    bet = betweenness_pct or 0.0

    # ── Tier 2: IDR / functional pathogenic (Strategy C) ───────────────────────
    if idr_evidence:
        return 2, (
            "Mutation in intrinsically disordered / very low pLDDT region (pLDDT < 50). "
            "DynaMut2 blocked. ELM / UniProt PTM analysis indicates functional motif or "
            "modification-site disruption (IDR pathogenesis track)."
        )

    # ── Tier 5: insufficient structural confidence ───────────────────────────
    if plddt < 50:
        return 5, (
            f"pLDDT = {plddt:.0f} (< 50). Cannot run reliable graph or ΔΔG analysis. "
            "Seek experimental structure or functional assays."
        )

    # ── Phase 3 ΔΔG triage (structured track, pLDDT >= 70) ───────────────────
    if ddg is not None and plddt >= 70:
        if ddg > -1.0:
            return 4, (
                f"ΔΔG = {ddg:+.2f} kcal/mol (neutral noise, 0 to −1.0 kcal/mol). "
                "No thermodynamic destabilisation — likely benign / polymorphism by stability criteria."
            )
        if ddg > -2.0:
            return 3, (
                f"ΔΔG = {ddg:+.2f} kcal/mol (mildly destabilising, −1.0 to −2.0). "
                "Uncertain significance — combine with graph and literature."
            )
        if ddg <= -2.0 and bet >= 0.90:
            return 1, (
                f"ΔΔG = {ddg:+.2f} kcal/mol (≤ −2.0) AND betweenness "
                f"{100*bet:.0f}th percentile (≥ 90th). Likely pathogenic — structural."
            )
        if ddg <= -2.0:
            return 3, (
                f"ΔΔG = {ddg:+.2f} kcal/mol (destabilising) but hub centrality below 90th percentile. "
                "Escalate with functional evidence."
            )

    # ── Graph-only fallback when ddG missing ───────────────────────────────────
    if layer0_tier == "HIGH_PRIORITY" and plddt >= 70:
        return 3, (
            "Strong graph disruption (HIGH_PRIORITY) but ΔΔG not available. "
            "Re-run Layer 1 or submit to DynaMut2 manually."
        )

    if layer0_tier in ("HIGH_PRIORITY", "MODERATE") and plddt >= 70:
        return 3, (
            f"Graph tier {layer0_tier}; ΔΔG not calculated. Moderate structural concern."
        )

    if plddt >= 50:
        ddg_str = f"ΔΔG = {ddg:+.2f} kcal/mol" if ddg is not None else "ΔΔG not calculated"
        return 4, (
            f"No convergent pathogenic signal. {ddg_str}. "
            "Variant likely benign by structural criteria."
        )

    return 5, "Insufficient data for classification."


# ─── Load Pipeline Outputs ────────────────────────────────────────────────────

def load_json(path: str) -> Optional[dict]:
    if os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)
    return None


def gather_evidence(gene: str, mutation: str, uniprot: str) -> dict:
    """
    Collect all available pipeline outputs for a given variant.
    Returns a unified evidence dict — missing pieces are None.
    """
    # pLDDT per-residue CSV
    import re
    import pandas as pd
    m = re.match(r"[A-Z](\d+)[A-Z]", mutation.upper())
    resnum = int(m.group(1)) if m else None

    plddt_csv = os.path.join(DATA_DIR, f"{gene}_plddt_per_residue.csv")
    plddt_val = None
    if resnum and os.path.exists(plddt_csv):
        df = pd.read_csv(plddt_csv)
        row = df[df["resnum"] == resnum]
        if not row.empty:
            plddt_val = float(row.iloc[0]["plddt"])

    # Layer 0 graph output
    l0 = load_json(os.path.join(DATA_DIR, f"{gene}_{mutation}_layer0.json"))

    # Layer 1 DynaMut2 output
    l1 = load_json(os.path.join(DATA_DIR, f"{gene}_{mutation}_dynamut2.json"))

    # Phase 3 IDR track
    idr = load_json(os.path.join(DATA_DIR, f"{gene}_{mutation}_idr.json"))

    # Phase 3 Literature Engine
    lit = load_json(os.path.join(DATA_DIR, f"{gene}_{mutation}_literature.json"))

    # UniProt annotations
    annot = load_json(os.path.join(DATA_DIR, f"{uniprot}_annotations.json"))

    # ClinVar variants
    clinvar_path = os.path.join(DATA_DIR, f"{gene}_clinvar_variants.json")
    clinvar_match = None
    if os.path.exists(clinvar_path) and resnum:
        import re as re2
        with open(clinvar_path) as fh:
            cv_list = json.load(fh)
        for v in cv_list:
            pc = v.get("protein_change", "")
            mn = re2.search(r"\d+", pc)
            if mn and int(mn.group()) == resnum:
                clinvar_match = v
                break

    # Domain context
    containing_domains = []
    if annot and resnum:
        for dom in annot.get("domains", []):
            if dom["start"] <= resnum <= dom["end"]:
                containing_domains.append(dom)

    return {
        "gene":              gene,
        "mutation":          mutation,
        "uniprot":           uniprot,
        "resnum":            resnum,
        "plddt":             plddt_val,
        "protein_name":      annot.get("protein_name") if annot else None,
        "protein_length":    annot.get("length") if annot else None,
        "containing_domains": containing_domains,
        "layer0":            l0,
        "layer1":            l1,
        "idr":               idr,
        "literature":        lit,
        "clinvar_match":     clinvar_match,
    }


# ─── Markdown Report ──────────────────────────────────────────────────────────

def format_markdown_report(ev: dict, tier: int, rationale: str) -> str:
    gene       = ev["gene"]
    mutation   = ev["mutation"]
    plddt      = ev["plddt"]
    l0         = ev["layer0"] or {}
    l1         = ev["layer1"] or {}
    clinvar    = ev["clinvar_match"] or {}
    domains    = ev["containing_domains"]
    tier_label, tier_icon = TIERS[tier]
    now        = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    # ddG display
    ddg_val  = l1.get("ddg_kcal_mol")
    ddg_tool = l1.get("tool", "Not calculated")
    ddg_tier = l1.get("tier", "N/A")
    ddg_str  = f"{ddg_val:+.3f} kcal/mol ({ddg_tool}, {ddg_tier})" if ddg_val is not None else "Not calculated"

    # Graph metrics display
    bet_pct = l0.get("betweenness_pct", 0.0) * 100 if l0 else None
    deg_pct = l0.get("degree_pct", 0.0) * 100 if l0 else None

    # Domain context
    if domains:
        domain_str = "; ".join(
            f"{d['desc']} [{d['type']}] (residues {d['start']}–{d['end']})"
            for d in domains
        )
    else:
        domain_str = "No annotated domain at this position"

    # pLDDT context
    if plddt is None:
        plddt_str = "Not determined (pLDDT audit not run)"
    else:
        from plddt_audit import plddt_tier, plddt_routing  # type: ignore
        plddt_str = f"{plddt:.1f} — {plddt_tier(plddt)} ({plddt_routing(plddt).split('.')[0]})"

    # ClinVar match
    if clinvar:
        cv_str = (
            f"**ClinVar variation:** {clinvar.get('title', 'Unknown')}\n"
            f"**Clinical significance:** {clinvar.get('significance', 'Unknown')}\n"
            f"**Associated conditions:** {', '.join(clinvar.get('conditions', ['Unknown']))}"
        )
    else:
        cv_str = "*No direct ClinVar match found for this exact residue position.*"

    doc = f"""# Digital Patient Twin — Clinical Report

**Generated:** {now}
**Variant:** `{gene}:{mutation}`
**Protein:** {ev.get('protein_name', 'Unknown')} (UniProt: {ev.get('uniprot', 'Unknown')})
**Protein length:** {ev.get('protein_length', 'Unknown')} residues

---

## {tier_icon} Classification: Tier {tier} — {tier_label}

> **Rationale:** {rationale}

---

## 1 · Variant Context

| Field | Value |
|---|---|
| Gene | {gene} |
| Mutation | {mutation} |
| Residue position | {ev.get('resnum', 'Unknown')} |
| Structural domain | {domain_str} |
| ClinVar significance | {clinvar.get('significance', 'Not found')} |

{cv_str}

---

## 2 · Structural Confidence

| Metric | Value |
|---|---|
| AlphaFold pLDDT at target residue | {plddt_str} |
| Structure source | AlphaFold DB v4 |
| Analysis basis | {'Structural analysis (pLDDT ≥ 70)' if plddt and plddt >= 70 else 'Limited — low confidence region'} |

> ⚠️ pLDDT is a **predicted** confidence score, not an experimental measurement.
> Structures with pLDDT < 70 should be interpreted with additional caution.

---

## 3 · Layer 0 — Graph-Based Structural Importance

{'*Graph analysis not available — run `build_graph.py` first.*' if not l0 else f"""
| Metric | Value | Percentile |
|---|---|---|
| Degree centrality | {l0.get("degree", 0):.4f} | {deg_pct:.1f}th |
| Betweenness centrality | {l0.get("betweenness", 0):.6f} | {bet_pct:.1f}th |
| Clustering coefficient | {l0.get("clustering", {}):.4f} | — |
| Edges removed by mutation | {l0.get("disruption", {}).get("edges_removed", 0)} | — |
| Δ Avg path length | {l0.get("disruption", {}).get("delta_avg_path", 0):+.4f} | — |
| Δ Connected components | {l0.get("disruption", {}).get("delta_components", 0):+d} | — |

**Layer 0 classification:** `{l0.get("layer0_tier", "N/A")}`

{l0.get("layer0_note", "")}
"""}

> **Interpretation:** Residues in the top 10th percentile of betweenness centrality
> are structural hubs whose removal significantly disrupts long-range protein communication.
> Degree percentile reflects the number of local contacts broken.

---

## 4 · Layer 1 — ΔΔG Stability Prediction (DynaMut2)

{'*DynaMut2 analysis not available — run `dynamut2_query.py` first.*' if not l1 else f"""
| Metric | Value |
|---|---|
| ΔΔG (kcal/mol) | **{ddg_val:+.3f}** |
| Prediction tool | {ddg_tool} |
| Stability tier | `{ddg_tier}` |
| Clinical signal | {l1.get("clinical_signal", "Unknown")} |
| Recommended action | {l1.get("action", "Unknown")} |

**Calibration reference (Section 2.2):**
- ΔΔG > 0: stabilising — typically benign
- ΔΔG −1.0 to 0: neutral/noise
- ΔΔG −1.0 to −2.0: mildly destabilising — uncertain significance
- ΔΔG −2.0 to −4.0: destabilising — **likely pathogenic**
- ΔΔG < −4.0: severely destabilising — **strong pathogenic signal**
"""}

---

## 5 · Literature Intelligence (Gemini AI)

{'*Literature analysis not available.*' if not ev.get('literature') else f"""
| Metric | Value |
|---|---|
| Concordant sources | {ev['literature'].get('concordant_count', 0)} |
| Discordant sources | {ev['literature'].get('discordant_count', 0)} |
| Experimental Methods | {", ".join(ev['literature'].get('experimental_methods', [])) or "None reported"} |

**AI Summary:**
> {ev['literature'].get('summary', 'No summary available.')}
"""}

---

## 6 · Layer 0 ↔ Layer 1 Concordance

{'*Concordance check requires both Layer 0 and Layer 1 results.*' if not l0 or not l1 else f"""
| | Layer 0 (graph) | Layer 1 (ΔΔG) |
|---|---|---|
| Tier | `{l0.get("layer0_tier", "N/A")}` | `{ddg_tier}` |
| Concordant? | {'**YES ✓**' if (l0.get("layer0_tier") in ("HIGH_PRIORITY","MODERATE") and ddg_val and ddg_val <= -1.5) or (l0.get("layer0_tier") == "LOW" and ddg_val and ddg_val > -1.5) else '**NO — manual review recommended**'} | |

> When Layer 0 and Layer 1 conflict, do not escalate to higher tiers without additional evidence.
> Proceed to Layer 2 (PyRosetta) for disambiguation.
"""}

---

## 7 · Mandatory Caveats

> **⚠️ This report is a computational prediction tool, not a diagnostic device.**

1. **Computational predictions only.** Results have not been validated in vitro or in vivo
   for this specific variant. Experimental confirmation is required before clinical action.

2. **pLDDT confidence.** The AlphaFold structure has a pLDDT of **{f'{plddt:.0f}' if plddt is not None else 'Unknown'}** at
   the target residue. {'This is in the high-confidence range (≥70).' if plddt is not None and plddt >= 70 else 'This is BELOW the reliable confidence threshold. Structural conclusions are tentative.'}

3. **ΔΔG uncertainty.** The DynaMut2 stability estimate has a reported correlation
   of ~0.68 with experimental measurements. For clinical reporting, Layer 2 (PyRosetta)
   or Layer 3 (FoldX) confirmation is recommended.

4. **Intended use.** This report is decision support for a qualified clinician.
   It is not a standalone diagnosis. Genetic counselling should accompany any
   patient-facing communication of this result.

5. **Pipeline version.** Digital Patient Twin v0.1 (research prototype).
   Not FDA-cleared. Not CE-marked.

---

## 8 · Evidence Audit Trail

| Stage | Status | Tool | Output file |
|---|---|---|---|
| Structure fetch | {'✓' if os.path.exists(os.path.join(DATA_DIR, f"AF-{ev['uniprot']}-F1-model_v4.pdb")) else '✗ Not run'} | AlphaFold DB | `AF-{ev['uniprot']}-F1-model_v4.pdb` |
| pLDDT audit | {'✓' if plddt else '✗ Not run'} | Biopython | `{gene}_plddt_per_residue.csv` |
| Graph build | {'✓' if l0 else '✗ Not run'} | NetworkX | `{gene}_centrality.csv` |
| Layer 0 sim | {'✓' if l0 else '✗ Not run'} | This pipeline | `{gene}_{mutation}_layer0.json` |
| Layer 1 ddG | {'✓' if l1 else '✗ Not run'} | DynaMut2 | `{gene}_{mutation}_dynamut2.json` |
| Literature  | {'✓' if ev.get('literature') else '✗ Not run'} | Gemini AI | `{gene}_{mutation}_literature.json` |
| Layer 2 ddG | ✗ Not run | PyRosetta | — |
| Layer 3 ddG | ✗ Not run | FoldX | — |

*To escalate: run `python pyrosetta_ddg.py --gene {gene} --mutation {mutation}`*

---

*Digital Patient Twin — Mechanistic AI for Rare Disease*
*Open Science Stack · Zero-budget research prototype*
*For questions contact the generating researcher, not the tool.*
"""
    return doc


# ─── JSON Report ──────────────────────────────────────────────────────────────

def format_json_report(ev: dict, tier: int, rationale: str) -> dict:
    l0  = ev["layer0"] or {}
    l1  = ev["layer1"] or {}
    now = datetime.now(timezone.utc).isoformat()
    tier_label, _ = TIERS[tier]

    return {
        "schema_version":   "1.0",
        "generated_at":     now,
        "tool":             "Digital Patient Twin v0.1",
        "variant": {
            "gene":          ev["gene"],
            "mutation":      ev["mutation"],
            "uniprot":       ev["uniprot"],
            "residue":       ev["resnum"],
            "protein_name":  ev["protein_name"],
        },
        "classification": {
            "tier":          tier,
            "label":         tier_label,
            "rationale":     rationale,
        },
        "structural_confidence": {
            "plddt":         ev["plddt"],
            "source":        "AlphaFold DB v4",
        },
        "domain_context":        ev["containing_domains"],
        "layer0_graph": {
            "tier":             l0.get("layer0_tier"),
            "note":             l0.get("layer0_note"),
            "degree_pct":       l0.get("degree_pct"),
            "betweenness_pct":  l0.get("betweenness_pct"),
            "clustering":       l0.get("clustering"),
            "edges_removed":    l0.get("disruption", {}).get("edges_removed"),
            "delta_path":       l0.get("disruption", {}).get("delta_avg_path"),
            "delta_components": l0.get("disruption", {}).get("delta_components"),
            "c_3d":             l0.get("c_3d"),
            "ddg_foldx":        l0.get("ddg_foldx"),
        },
        "layer1_ddg": {
            "ddg_kcal_mol":     l1.get("ddg_kcal_mol"),
            "ddg_foldx":        l0.get("ddg_foldx"),
            "tool":             l1.get("tool"),
            "tier":             l1.get("tier"),
            "clinical_signal":  l1.get("clinical_signal"),
        },
        "literature": ev.get("literature"),
        "clinvar": ev["clinvar_match"],
        "caveats": [
            "Computational prediction only — not validated in vitro for this variant.",
            "Not a diagnostic device. For clinician decision support only.",
            "Digital Patient Twin v0.1 — research prototype, not FDA-cleared.",
        ],
    }


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Generate clinical interpretation report.")
    parser.add_argument("--gene",     default="VCP",    help="Gene symbol")
    parser.add_argument("--mutation", default="R155H",  help="Mutation, e.g. R155H")
    parser.add_argument("--uniprot",  default="P55072", help="UniProt accession")
    parser.add_argument("--format",   default="both",
                        choices=["markdown", "json", "both"],
                        help="Output format (default: both)")
    args = parser.parse_args()

    print("=" * 65)
    print("  Digital Patient Twin — Phase 3: Report Generator")
    print(f"  Variant: {args.gene}:{args.mutation}")
    print("=" * 65)

    # Collect all evidence
    ev = gather_evidence(args.gene, args.mutation, args.uniprot)

    # Extract key signals for tier classification
    l0  = ev["layer0"] or {}
    l1  = ev["layer1"] or {}
    ddg = l1.get("ddg_kcal_mol")
    bet = l0.get("betweenness_pct")
    l0_tier = l0.get("layer0_tier")

    # Classify
    tier, rationale = classify_tier(
        plddt           = ev["plddt"],
        layer0_tier     = l0_tier,
        ddg             = ddg,
        betweenness_pct = bet,
        idr_evidence    = False,     # set True when IDR analysis is implemented
    )

    tier_label, tier_icon = TIERS[tier]
    print(f"\n  {tier_icon}  Tier {tier}: {tier_label}")
    print(f"  Rationale: {rationale[:80]}...")

    # Generate outputs
    if args.format in ("markdown", "both"):
        md_content  = format_markdown_report(ev, tier, rationale)
        md_path     = os.path.join(REPORTS_DIR, f"{args.gene}_{args.mutation}_report.md")
        with open(md_path, "w", encoding="utf-8") as fh:
            fh.write(md_content)
        print(f"\n  [Report] Markdown saved to: {md_path}")

    if args.format in ("json", "both"):
        json_content = format_json_report(ev, tier, rationale)
        json_path    = os.path.join(REPORTS_DIR, f"{args.gene}_{args.mutation}_report.json")
        with open(json_path, "w", encoding="utf-8") as fh:
            json.dump(json_content, fh, indent=2)
        print(f"  [Report] JSON saved to: {json_path}")

    print("\n" + "═" * 65)
    print("  PIPELINE COMPLETE")
    print(f"  Variant:    {args.gene}:{args.mutation}")
    print(f"  Tier:       {tier} ({tier_label})")
    print(f"  pLDDT:      {ev['plddt']:.1f if ev['plddt'] else 'N/A'}")
    if ddg:
        print(f"  ΔΔG:        {ddg:+.3f} kcal/mol ({l1.get('tier', 'N/A')})")
    print("═" * 65)


if __name__ == "__main__":
    main()
