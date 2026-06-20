"""
run_pipeline.py — End-to-end pipeline orchestrator
Runs all stages in sequence for a given gene/variant:
  Stage 1: fetch_structure      → AlphaFold PDB + UniProt annotations + ClinVar variants
  Stage 2: plddt_audit          → per-residue confidence table + routing + plot
  Stage 3: build_graph          → contact network + centrality + Layer 0 disruption
  Stage 4: dynamut2_query       → Layer 1 ΔΔG (skipped if server unreachable)
  Stage 5: report_generator     → Markdown + JSON clinical report

Usage:
    python run_pipeline.py
    python run_pipeline.py --gene VCP --mutation R155H --uniprot P55072
    python run_pipeline.py --gene VCP --mutation R155H --skip-dynamut2
"""

import argparse
import importlib
import subprocess
import sys
import os
import time
from pathlib import Path

ROOT = Path(__file__).parent


def run_stage(label: str, module_name: str, extra_args: list[str]) -> bool:
    """
    Run a pipeline stage as a subprocess (preserves clean stdout).
    Returns True on success, False on failure.
    """
    cmd = [sys.executable, str(ROOT / f"{module_name}.py")] + extra_args
    print(f"\n{'━'*65}")
    print(f"  ▶  {label}")
    print(f"     {' '.join(cmd[1:])}")
    print(f"{'━'*65}")

    t0     = time.time()
    result = subprocess.run(cmd, cwd=str(ROOT))
    elapsed = time.time() - t0

    if result.returncode == 0:
        print(f"\n  ✓  {label} completed in {elapsed:.1f}s")
        return True
    else:
        print(f"\n  ✗  {label} FAILED (exit code {result.returncode})")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Digital Patient Twin — full pipeline runner",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_pipeline.py
  python run_pipeline.py --gene ATPB1 --mutation G185R --uniprot P13637
  python run_pipeline.py --gene VCP   --mutation R191Q --skip-dynamut2
        """
    )
    parser.add_argument("--gene",          default="VCP",    help="Gene symbol (default: VCP)")
    parser.add_argument("--mutation",      default="R155H",  help="Mutation string (default: R155H)")
    parser.add_argument("--uniprot",       default="P55072", help="UniProt accession (default: P55072)")
    parser.add_argument("--chain",         default="A",      help="PDB chain ID (default: A)")
    parser.add_argument("--cutoff",        default="8.0",    help="Calpha contact cutoff Å (default: 8.0)")
    parser.add_argument("--skip-dynamut2", action="store_true",
                        help="Skip DynaMut2 query (use when offline or server is unreachable)")
    parser.add_argument("--skip-fetch",    action="store_true",
                        help="Skip data fetch (if PDB already downloaded)")
    args = parser.parse_args()

    print()
    print("╔" + "═"*63 + "╗")
    print("║  DIGITAL PATIENT TWIN — Pipeline Runner" + " "*23 + "║")
    print("║  Mechanistic AI for Rare Disease" + " "*30 + "║")
    print("╠" + "═"*63 + "╣")
    print(f"║  Gene:     {args.gene:<51} ║")
    print(f"║  Mutation: {args.mutation:<51} ║")
    print(f"║  UniProt:  {args.uniprot:<51} ║")
    print("╚" + "═"*63 + "╝")

    stages_run   = []
    stages_failed= []

    base_args = [
        "--gene",     args.gene,
        "--mutation", args.mutation,
        "--uniprot",  args.uniprot,
    ]

    # ── Stage 1: Data fetch ───────────────────────────────────────────────────
    if not args.skip_fetch:
        ok = run_stage(
            "Stage 1 — Data Foundation (AlphaFold + UniProt + ClinVar)",
            "fetch_structure",
            ["--gene", args.gene, "--uniprot", args.uniprot],
        )
        (stages_run if ok else stages_failed).append("fetch_structure")
        if not ok:
            print("\n  Pipeline aborted: cannot proceed without structure data.")
            sys.exit(1)
    else:
        print("\n  Stage 1 skipped (--skip-fetch).")

    # ── Stage 2: pLDDT audit ──────────────────────────────────────────────────
    ok = run_stage(
        "Stage 2 — pLDDT Confidence Audit",
        "plddt_audit",
        ["--gene", args.gene, "--uniprot", args.uniprot],
    )
    (stages_run if ok else stages_failed).append("plddt_audit")

    # ── Stage 3: Graph engine ─────────────────────────────────────────────────
    ok = run_stage(
        "Stage 3 — Protein Contact Graph + Layer 0 Disruption",
        "build_graph",
        ["--gene", args.gene, "--mutation", args.mutation,
         "--uniprot", args.uniprot, "--cutoff", args.cutoff],
    )
    (stages_run if ok else stages_failed).append("build_graph")

    # ── Stage 4: DynaMut2 ─────────────────────────────────────────────────────
    if not args.skip_dynamut2:
        ok = run_stage(
            "Stage 4 — Layer 1 ΔΔG (DynaMut2)",
            "dynamut2_query",
            ["--gene", args.gene, "--mutation", args.mutation,
             "--uniprot", args.uniprot, "--chain", args.chain],
        )
        (stages_run if ok else stages_failed).append("dynamut2_query")
        if not ok:
            print("\n  DynaMut2 unavailable — generating report without Layer 1 ddG.")
            print("  Re-run later with:  python dynamut2_query.py " + " ".join(base_args))
    else:
        print("\n  Stage 4 skipped (--skip-dynamut2).")

    # ── Stage 5: Report ───────────────────────────────────────────────────────
    ok = run_stage(
        "Stage 5 — Clinical Report (Markdown + JSON)",
        "report_generator",
        base_args + ["--format", "both"],
    )
    (stages_run if ok else stages_failed).append("report_generator")

    # ── Summary ───────────────────────────────────────────────────────────────
    print()
    print("╔" + "═"*63 + "╗")
    print("║  PIPELINE SUMMARY" + " "*45 + "║")
    print("╠" + "═"*63 + "╣")
    for s in stages_run:
        print(f"║  ✓  {s:<57} ║")
    for s in stages_failed:
        print(f"║  ✗  {s:<57} ║")
    print("╠" + "═"*63 + "╣")

    report_md   = ROOT / "reports" / f"{args.gene}_{args.mutation}_report.md"
    report_json = ROOT / "reports" / f"{args.gene}_{args.mutation}_report.json"
    if report_md.exists():
        print(f"║  📄 Report: reports/{args.gene}_{args.mutation}_report.md" +
              " "*(63 - len(f"  📄 Report: reports/{args.gene}_{args.mutation}_report.md") - 1) + "║")
    if report_json.exists():
        print(f"║  📦 JSON:   reports/{args.gene}_{args.mutation}_report.json" +
              " "*(63 - len(f"  📦 JSON:   reports/{args.gene}_{args.mutation}_report.json") - 1) + "║")
    print("╚" + "═"*63 + "╝")

    if stages_failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
