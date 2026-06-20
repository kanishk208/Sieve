"""
plddt_audit.py — Phase 1, Step 3 (pLDDT Confidence Audit)
Reads a downloaded AlphaFold PDB file and produces:
  - Per-residue pLDDT table (saved as CSV)
  - Summary statistics printed to console
  - A colour-coded confidence plot saved as PNG
  - The gatekeeper routing decision for each ClinVar variant

Usage:
    python plddt_audit.py
    python plddt_audit.py --pdb data/AF-P55072-F1-model_v4.pdb --variants data/VCP_clinvar_variants.json

Requires: biopython numpy pandas matplotlib
"""

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

try:
    from Bio import PDB
except ImportError:
    sys.exit("ERROR: biopython not installed.  Run: pip install biopython")

try:
    import matplotlib
    matplotlib.use("Agg")          # headless — no display needed
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    MATPLOTLIB_OK = True
except ImportError:
    MATPLOTLIB_OK = False
    print("WARNING: matplotlib not installed. Plots will be skipped.")

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
os.makedirs(DATA_DIR, exist_ok=True)


# ─── pLDDT Colour Scheme (mirrors AlphaFoldDB convention) ────────────────────
PLDDT_BINS = [
    (90, 100, "#0053D6", "Very high (90-100)"),
    (70,  90, "#65CBF3", "Confident (70-90)"),
    (50,  70, "#FFDB13", "Low (50-70)"),
    ( 0,  50, "#FF7D45", "Very low (<50) — disordered"),
]

def plddt_colour(score: float) -> str:
    for lo, hi, colour, _ in PLDDT_BINS:
        if lo <= score < hi:
            return colour
    return "#FF7D45"

def plddt_tier(score: float) -> str:
    if score >= 90:
        return "VERY_HIGH"
    elif score >= 70:
        return "CONFIDENT"
    elif score >= 50:
        return "LOW"
    else:
        return "DISORDERED"

def plddt_routing(score: float) -> str:
    """Return the analysis track recommended by Section 3.4 of the architecture document."""
    if score >= 90:
        return "Full analysis — Graph + ddG. High-confidence reporting."
    elif score >= 70:
        return "Full analysis with moderate-confidence caveat."
    elif score >= 50:
        return "Check RCSB for experimental structure. If none, use ESMFold. Flag as low-confidence."
    else:
        return "DISORDERED REGION — do NOT run ddG. Switch to IDR track: ELM / PhosphoSitePlus / catGRANULE."


# ─── Parse PDB for pLDDT ─────────────────────────────────────────────────────

def parse_plddt(pdb_path: str) -> pd.DataFrame:
    """
    Extract per-residue pLDDT scores from AlphaFold PDB.
    AlphaFold stores pLDDT in the B-factor column of CA atoms.
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)

    records = []
    for model in structure:
        for chain in model:
            chain_id = chain.get_id()
            for residue in chain:
                het, resnum, icode = residue.get_id()
                if het.strip():          # skip HETATM records (water, ligands)
                    continue
                resname = residue.get_resname()
                for atom in residue:
                    if atom.get_name() == "CA":   # alpha-carbon only
                        plddt = atom.get_bfactor()
                        records.append({
                            "chain":   chain_id,
                            "resnum":  resnum,
                            "resname": resname,
                            "plddt":   plddt,
                            "tier":    plddt_tier(plddt),
                            "colour":  plddt_colour(plddt),
                            "routing": plddt_routing(plddt),
                        })
                        break  # only one CA per residue

    df = pd.DataFrame(records)
    return df


# ─── Summary Report ───────────────────────────────────────────────────────────

def print_summary(df: pd.DataFrame, pdb_path: str):
    n = len(df)
    mean_plddt  = df["plddt"].mean()
    median_plddt= df["plddt"].median()
    very_high   = (df["plddt"] >= 90).sum()
    confident   = ((df["plddt"] >= 70) & (df["plddt"] < 90)).sum()
    low_conf    = ((df["plddt"] >= 50) & (df["plddt"] < 70)).sum()
    disordered  = (df["plddt"] < 50).sum()

    print("\n" + "═" * 65)
    print("  pLDDT CONFIDENCE AUDIT")
    print(f"  File: {os.path.basename(pdb_path)}")
    print("═" * 65)
    print(f"  Total residues analysed : {n}")
    print(f"  Mean pLDDT              : {mean_plddt:.1f}")
    print(f"  Median pLDDT            : {median_plddt:.1f}")
    print()
    print(f"  Very high  (90-100)  : {very_high:>4}  ({100*very_high/n:.1f}%)  → Full analysis safe")
    print(f"  Confident  (70-90)   : {confident:>4}  ({100*confident/n:.1f}%)  → Full analysis + caveat")
    print(f"  Low conf   (50-70)   : {low_conf:>4}  ({100*low_conf/n:.1f}%)  → Seek experimental structure")
    print(f"  Disordered (<50)     : {disordered:>4}  ({100*disordered/n:.1f}%)  → IDR track only")
    print("═" * 65)


# ─── Variant Routing ─────────────────────────────────────────────────────────

def route_clinvar_variants(df: pd.DataFrame, variants_path: str) -> None:
    """
    For each ClinVar variant, look up pLDDT at the affected residue position
    and print the recommended analysis routing.
    """
    if not os.path.exists(variants_path):
        print(f"\n[Routing] Variants file not found: {variants_path}  — skipping routing step.")
        return

    with open(variants_path) as fh:
        variants = json.load(fh)

    print(f"\n{'─'*65}")
    print("  VARIANT ROUTING DECISIONS")
    print(f"{'─'*65}")
    print(f"  {'Variant':<25} {'Residue':>7}  {'pLDDT':>6}  Tier           Routing")
    print(f"  {'─'*25} {'─'*7}  {'─'*6}  {'─'*13}  {'─'*20}")

    routed = []
    for v in variants:
        protein_change = v.get("protein_change", "")
        title = v.get("title", "")

        # Try to extract residue number from protein_change (e.g. "R155H" → 155)
        res_num = None
        import re
        # protein_change field can look like "R155H" or "p.Arg155His"
        m = re.search(r"(\d+)", protein_change)
        if m:
            res_num = int(m.group(1))

        if res_num is not None:
            row = df[df["resnum"] == res_num]
            if not row.empty:
                plddt = row.iloc[0]["plddt"]
                tier  = row.iloc[0]["tier"]
                route = plddt_routing(plddt)
                short_route = route.split("—")[0].strip() if "—" in route else route[:30]
            else:
                plddt, tier, short_route = float("nan"), "NOT_FOUND", "Residue not in structure"
        else:
            plddt, tier, short_route = float("nan"), "UNKNOWN", "Cannot parse residue"

        label = protein_change or title[:24]
        plddt_str = f"{plddt:.0f}" if not np.isnan(plddt) else "  —"
        print(f"  {label:<25} {res_num if res_num else '?':>7}  {plddt_str:>6}  {tier:<13}  {short_route}")
        routed.append({**v, "resnum": res_num, "plddt": plddt, "tier": tier, "routing": short_route if res_num else "unresolvable"})

    out_path = os.path.join(DATA_DIR, "variant_routing.json")
    with open(out_path, "w") as fh:
        json.dump(routed, fh, indent=2, default=str)
    print(f"\n  Routing table saved to: {out_path}")


# ─── Plot ─────────────────────────────────────────────────────────────────────

def plot_plddt(df: pd.DataFrame, out_png: str, gene: str = "VCP"):
    if not MATPLOTLIB_OK:
        return

    fig, axes = plt.subplots(2, 1, figsize=(16, 8),
                             gridspec_kw={"height_ratios": [3, 1]})
    fig.patch.set_facecolor("#0D1117")

    # ── Top panel: pLDDT line plot ──────────────────────────────────────────
    ax = axes[0]
    ax.set_facecolor("#0D1117")

    # Coloured fill bands
    ax.axhspan(90, 100, color="#0053D6", alpha=0.10)
    ax.axhspan(70,  90, color="#65CBF3", alpha=0.10)
    ax.axhspan(50,  70, color="#FFDB13", alpha=0.10)
    ax.axhspan( 0,  50, color="#FF7D45", alpha=0.10)

    # Threshold lines
    for thresh, col in [(90, "#0053D6"), (70, "#65CBF3"), (50, "#FFDB13")]:
        ax.axhline(thresh, color=col, linewidth=0.6, linestyle="--", alpha=0.6)

    # pLDDT line — colour each segment
    x = df["resnum"].values
    y = df["plddt"].values
    colours = df["colour"].values

    for i in range(len(x) - 1):
        ax.plot(x[i:i+2], y[i:i+2], color=colours[i], linewidth=1.0, alpha=0.9)

    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(0, 100)
    ax.set_ylabel("pLDDT", color="white", fontsize=11)
    ax.set_title(f"{gene} — AlphaFold pLDDT Confidence Profile", color="white", fontsize=13, pad=10)
    ax.tick_params(colors="white")
    for spine in ax.spines.values():
        spine.set_edgecolor("#3D4450")

    # Legend
    patches = [mpatches.Patch(color=c, label=lbl) for _, _, c, lbl in PLDDT_BINS]
    ax.legend(handles=patches, loc="lower right", framealpha=0.3,
              labelcolor="white", fontsize=8, facecolor="#161B22")

    # ── Bottom panel: tier colour strip ────────────────────────────────────
    ax2 = axes[1]
    ax2.set_facecolor("#0D1117")
    colours_strip = df["colour"].values
    ax2.bar(x, [1]*len(x), color=colours_strip, width=1.0, align="center")
    ax2.set_xlim(x[0], x[-1])
    ax2.set_ylim(0, 1)
    ax2.set_yticks([])
    ax2.set_xlabel("Residue Number", color="white", fontsize=10)
    ax2.tick_params(colors="white", axis="x")
    ax2.set_title("Confidence Tier Strip", color="#8B949E", fontsize=9)
    for spine in ax2.spines.values():
        spine.set_edgecolor("#3D4450")

    plt.tight_layout(pad=1.5)
    plt.savefig(out_png, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close()
    print(f"\n[Plot] pLDDT profile saved to: {out_png}")


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Audit AlphaFold structure confidence (pLDDT).")
    parser.add_argument("--pdb",      default=None, help="Path to AlphaFold PDB file")
    parser.add_argument("--variants", default=None, help="Path to ClinVar variants JSON")
    parser.add_argument("--gene",     default="VCP", help="Gene symbol (for plot title)")
    parser.add_argument("--uniprot",  default="P55072", help="UniProt ID (for auto file-finding)")
    args = parser.parse_args()

    # Auto-find PDB if not specified
    if args.pdb is None:
        candidates = [f for f in os.listdir(DATA_DIR)
                      if f.startswith(f"AF-{args.uniprot}") and f.endswith(".pdb")]
        if not candidates:
            sys.exit(
                f"ERROR: No PDB file found in data/ for {args.uniprot}.\n"
                "       Run  python fetch_structure.py  first."
            )
        args.pdb = os.path.join(DATA_DIR, sorted(candidates)[-1])

    # Auto-find variants JSON
    if args.variants is None:
        candidates = [f for f in os.listdir(DATA_DIR)
                      if f.startswith(args.gene) and f.endswith("_clinvar_variants.json")]
        if candidates:
            args.variants = os.path.join(DATA_DIR, candidates[0])

    print(f"[Audit] Reading structure: {args.pdb}")
    df = parse_plddt(args.pdb)

    if df.empty:
        sys.exit("ERROR: No residues parsed from PDB file. Check file path and format.")

    # Save per-residue table
    csv_path = os.path.join(DATA_DIR, f"{args.gene}_plddt_per_residue.csv")
    df.to_csv(csv_path, index=False)
    print(f"[Audit] Per-residue pLDDT saved to: {csv_path}")

    # Print summary
    print_summary(df, args.pdb)

    # Route variants
    if args.variants:
        route_clinvar_variants(df, args.variants)

    # Plot
    png_path = os.path.join(DATA_DIR, f"{args.gene}_plddt_profile.png")
    plot_plddt(df, png_path, gene=args.gene)

    print("\n[Audit] ✓ Complete. Next step: run  python build_graph.py")


if __name__ == "__main__":
    main()
