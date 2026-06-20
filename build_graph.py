"""
build_graph.py — Phase 2: Protein Contact Network
Constructs the protein contact graph from Calpha coordinates,
computes centrality metrics for every residue, simulates a target
mutation by removing its edges, and reports the disruption signal.

Usage:
    python build_graph.py
    python build_graph.py --pdb data/AF-P55072-F1-model_v4.pdb --mutation R155H --gene VCP

Requires: biopython networkx numpy pandas matplotlib
"""

import argparse
import json
import os
import re
import sys
from typing import Optional

import numpy as np
import pandas as pd

try:
    from Bio import PDB
except ImportError:
    sys.exit("ERROR: biopython not installed.  Run: pip install biopython")

try:
    import networkx as nx
except ImportError:
    sys.exit("ERROR: networkx not installed.  Run: pip install networkx")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    MATPLOTLIB_OK = True
except ImportError:
    MATPLOTLIB_OK = False

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
os.makedirs(DATA_DIR, exist_ok=True)

# Standard Calpha-Calpha contact cutoff (Angstroms)
CONTACT_CUTOFF_A = 8.0


# ─── Build Graph ─────────────────────────────────────────────────────────────

def build_contact_graph(pdb_path: str, cutoff: float = CONTACT_CUTOFF_A) -> tuple[nx.Graph, pd.DataFrame]:
    """
    Parse Calpha atoms from PDB. Build an undirected contact graph where
    an edge exists between residues i and j iff dist(CA_i, CA_j) <= cutoff Å.
    Also returns a DataFrame of residue metadata (pLDDT, resname, position).
    """
    print(f"\n[Graph] Parsing structure: {pdb_path}")
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("prot", pdb_path)

    residues = []
    for model in structure:
        for chain in model:
            for residue in chain:
                het, resnum, _ = residue.get_id()
                if het.strip():
                    continue
                for atom in residue:
                    if atom.get_name() == "CA":
                        residues.append({
                            "resnum":  resnum,
                            "resname": residue.get_resname(),
                            "x": atom.get_vector()[0],
                            "y": atom.get_vector()[1],
                            "z": atom.get_vector()[2],
                            "plddt": atom.get_bfactor(),
                        })
                        break

    df = pd.DataFrame(residues)
    n  = len(df)
    print(f"[Graph] {n} Calpha atoms loaded.")

    # Build graph
    G = nx.Graph()
    for _, row in df.iterrows():
        G.add_node(int(row["resnum"]),
                   resname=row["resname"],
                   plddt=row["plddt"])

    coords = df[["x", "y", "z"]].values
    resnums = df["resnum"].values

    edge_count = 0
    for i in range(n):
        for j in range(i + 4, n):   # skip i±1,2,3 (sequential, non-informative)
            dist = np.linalg.norm(coords[i] - coords[j])
            if dist <= cutoff:
                G.add_edge(int(resnums[i]), int(resnums[j]), weight=1.0 / dist)
                edge_count += 1

    print(f"[Graph] Nodes: {G.number_of_nodes()} | Edges: {edge_count} (cutoff = {cutoff} Å, skip ±3)")
    return G, df


# ─── Centrality Metrics ───────────────────────────────────────────────────────

def compute_centrality(G: nx.Graph) -> pd.DataFrame:
    """
    Compute per-residue structural importance metrics.
    Returns a DataFrame indexed by residue number.
    """
    print("[Graph] Computing centrality metrics (this may take 1-3 minutes for large proteins)...")

    degree      = nx.degree_centrality(G)
    clustering  = nx.clustering(G)

    # Betweenness is the expensive one — use k-approximation for speed
    n_nodes = G.number_of_nodes()
    k_sample = min(500, n_nodes)     # sample up to 500 nodes for betweenness
    betweenness = nx.betweenness_centrality(G, k=k_sample, normalized=True, seed=42)

    records = []
    for node in G.nodes():
        records.append({
            "resnum":       node,
            "resname":      G.nodes[node].get("resname", "UNK"),
            "plddt":        G.nodes[node].get("plddt", 0.0),
            "degree":       degree[node],
            "betweenness":  betweenness[node],
            "clustering":   clustering[node],
            "degree_rank":  None,      # filled below
            "betweenness_rank": None,
        })

    df = pd.DataFrame(records).sort_values("resnum").reset_index(drop=True)
    df["degree_rank"]      = df["degree"].rank(pct=True)
    df["betweenness_rank"] = df["betweenness"].rank(pct=True)

    print(f"[Graph] Centrality computed for {len(df)} residues.")
    return df


# ─── Mutation Simulation ──────────────────────────────────────────────────────

def simulate_mutation(G: nx.Graph, resnum: int) -> dict:
    """
    Simulate removing all edges of a residue (the Layer 0 structural disruption model).
    Report: how much does the global network change?
    """
    if resnum not in G.nodes():
        return {"error": f"Residue {resnum} not in graph."}

    # Before metrics
    degree_before = G.degree(resnum)
    neigh_before  = list(G.neighbors(resnum))

    # The wild-type baseline (avg path length, component count, sampling set) is
    # identical for every variant on the same graph, so compute it once and cache
    # it on the graph. On a batch of dozens of variants this removes a full
    # connected-components pass and a 200-source shortest-path sweep per variant.
    baseline = G.graph.get("_wt_baseline")
    if baseline is None:
        lcc_before = max(nx.connected_components(G), key=len)
        G_before   = G.subgraph(lcc_before).copy()
        sample     = list(G_before.nodes())[:200]   # sample for approx avg path
        lengths_before = []
        for src in sample:
            lengths = nx.single_source_shortest_path_length(G_before, src, cutoff=20)
            lengths_before.extend(lengths.values())
        baseline = {
            "sample":            sample,
            "avg_path_before":   float(np.mean(lengths_before)),
            "components_before": nx.number_connected_components(G),
        }
        G.graph["_wt_baseline"] = baseline

    sample          = baseline["sample"]
    avg_path_before = baseline["avg_path_before"]
    components_before = baseline["components_before"]

    # Create mutant graph by removing the target residue's edges. Removing edges
    # (rather than copying the whole graph) and restoring them afterwards avoids a
    # full G.copy() per variant — the dominant per-variant cost on large proteins.
    removed_edges = list(G.edges(resnum, data=True))
    G.remove_edges_from([(u, v) for u, v, _ in removed_edges])
    try:
        lcc_after = max(nx.connected_components(G), key=len)
        G_after   = G.subgraph(lcc_after)
        lengths_after = []
        for src in sample:
            if src in G_after:
                lengths = nx.single_source_shortest_path_length(G_after, src, cutoff=20)
                lengths_after.extend(lengths.values())
        avg_path_after   = np.mean(lengths_after) if lengths_after else float("inf")
        components_after = nx.number_connected_components(G)
    finally:
        # Restore the graph to wild-type so the cached baseline stays valid.
        G.add_edges_from(removed_edges)

    delta_path  = avg_path_after - avg_path_before
    delta_comps = components_after - components_before

    return {
        "resnum":             resnum,
        "edges_removed":      degree_before,
        "neighbours":         neigh_before[:10],   # first 10 for readability
        "avg_path_before":    round(avg_path_before, 4),
        "avg_path_after":     round(avg_path_after, 4),
        "delta_avg_path":     round(delta_path, 4),
        "components_before":  components_before,
        "components_after":   components_after,
        "delta_components":   delta_comps,
    }


# ─── Layer 0 Classification ───────────────────────────────────────────────────

def layer0_classification(centrality_row: pd.Series, mutation_result: dict, grantham_d: Optional[int] = None) -> tuple[str, str]:
    """
    Classify a variant based purely on graph metrics (Layer 0).
    Returns (tier_label, explanation).
    """
    bet_rank  = centrality_row["betweenness_rank"]
    deg_rank  = centrality_row["degree_rank"]
    delta_path= mutation_result.get("delta_avg_path", 0)
    delta_comp= mutation_result.get("delta_components", 0)
    edges     = mutation_result.get("edges_removed", 0)

    reasons = []

    if bet_rank >= 0.90:
        reasons.append(f"top-10% betweenness hub (rank={bet_rank:.2f})")
    if deg_rank >= 0.90:
        reasons.append(f"top-10% degree hub (rank={deg_rank:.2f})")
    if delta_path > 0.10:
        reasons.append(f"path length increased by {delta_path:.3f} steps")
    if delta_comp > 0:
        reasons.append(f"mutation disconnects {delta_comp} additional component(s)")
    if edges >= 10:
        reasons.append(f"{edges} contacts broken")

    score = (
        (1 if bet_rank >= 0.90 else 0) +
        (1 if deg_rank >= 0.90 else 0) +
        (1 if delta_path > 0.10 else 0) +
        (2 if delta_comp > 0 else 0) +
        (1 if edges >= 10 else 0)
    )

    if score >= 3:
        if grantham_d is not None and grantham_d < 50 and bet_rank >= 0.90:
            tier  = "MODERATE"
            label = "[CHEM_OVERRIDE] Conservative substitution in dense hub."
        else:
            tier  = "HIGH_PRIORITY"
            label = "Strong structural hub — escalate to Layer 1"
    elif score >= 1:
        tier  = "MODERATE"
        label = "Moderate structural signal — flag for literature review"
    else:
        tier  = "LOW"
        label = "Weak structural signal — likely benign by graph metrics"

    explanation = "; ".join(reasons) if reasons else "No significant graph disruption detected."
    return tier, f"{label}. Evidence: {explanation}"


# ─── Plot ─────────────────────────────────────────────────────────────────────

def plot_centrality(df: pd.DataFrame, target_resnum: Optional[int], out_png: str, gene: str):
    if not MATPLOTLIB_OK:
        return

    fig, axes = plt.subplots(1, 2, figsize=(16, 5))
    fig.patch.set_facecolor("#0D1117")

    for ax in axes:
        ax.set_facecolor("#0D1117")
        ax.tick_params(colors="white")
        for sp in ax.spines.values():
            sp.set_edgecolor("#3D4450")

    x = df["resnum"].values

    # ── Betweenness centrality ──
    ax = axes[0]
    c  = cm.plasma(df["betweenness_rank"].values)
    ax.scatter(x, df["betweenness"].values, c=c, s=2, alpha=0.7)
    if target_resnum is not None and target_resnum in df["resnum"].values:
        trow = df[df["resnum"] == target_resnum].iloc[0]
        ax.scatter(target_resnum, trow["betweenness"], s=80, c="white",
                   zorder=5, marker="*", label=f"Res {target_resnum}")
        ax.annotate(f"  {target_resnum}", (target_resnum, trow["betweenness"]),
                    color="white", fontsize=8)
    ax.set_title(f"{gene} — Betweenness Centrality", color="white", fontsize=11)
    ax.set_xlabel("Residue", color="white")
    ax.set_ylabel("Betweenness", color="white")

    # ── Degree centrality ──
    ax = axes[1]
    c  = cm.viridis(df["degree_rank"].values)
    ax.scatter(x, df["degree"].values, c=c, s=2, alpha=0.7)
    if target_resnum is not None and target_resnum in df["resnum"].values:
        trow = df[df["resnum"] == target_resnum].iloc[0]
        ax.scatter(target_resnum, trow["degree"], s=80, c="white",
                   zorder=5, marker="*")
        ax.annotate(f"  {target_resnum}", (target_resnum, trow["degree"]),
                    color="white", fontsize=8)
    ax.set_title(f"{gene} — Degree Centrality", color="white", fontsize=11)
    ax.set_xlabel("Residue", color="white")
    ax.set_ylabel("Degree", color="white")

    plt.tight_layout(pad=1.5)
    plt.savefig(out_png, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close()
    print(f"[Plot] Centrality plot saved to: {out_png}")


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Build protein contact graph and simulate mutation.")
    parser.add_argument("--pdb",      default=None,    help="AlphaFold PDB path")
    parser.add_argument("--mutation", default="R155H", help="Mutation in format [WT-AA][ResNum][Mut-AA], e.g. R155H")
    parser.add_argument("--gene",     default="VCP",   help="Gene symbol")
    parser.add_argument("--uniprot",  default="P55072",help="UniProt ID for auto file-finding")
    parser.add_argument("--cutoff",   default=8.0, type=float, help="Calpha-Calpha contact cutoff (Å)")
    args = parser.parse_args()

    # Auto-find PDB
    if args.pdb is None:
        candidates = [f for f in os.listdir(DATA_DIR)
                      if f.startswith(f"AF-{args.uniprot}") and f.endswith(".pdb")]
        if not candidates:
            sys.exit("ERROR: No PDB found. Run  python fetch_structure.py  first.")
        args.pdb = os.path.join(DATA_DIR, sorted(candidates)[-1])

    # Parse mutation residue number
    m = re.match(r"([A-Z])(\d+)([A-Z])", args.mutation.upper())
    if not m:
        sys.exit(f"ERROR: Cannot parse mutation '{args.mutation}'. Expected format: R155H")
    wt_aa, resnum_str, mut_aa = m.group(1), m.group(2), m.group(3)
    target_resnum = int(resnum_str)

    print("=" * 65)
    print(f"  Digital Patient Twin — Phase 2: Graph Engine")
    print(f"  Gene: {args.gene} | Mutation: {args.mutation} (residue {target_resnum})")
    print("=" * 65)

    # Build graph
    G, df_struct = build_contact_graph(args.pdb, cutoff=args.cutoff)

    # Compute centrality
    df_cent = compute_centrality(G)

    # Save centrality table
    cent_path = os.path.join(DATA_DIR, f"{args.gene}_centrality.csv")
    df_cent.to_csv(cent_path, index=False)
    print(f"[Graph] Centrality table saved to: {cent_path}")

    # ── Target residue stats ──────────────────────────────────────────────────
    if target_resnum in df_cent["resnum"].values:
        trow = df_cent[df_cent["resnum"] == target_resnum].iloc[0]
        print(f"\n{'─'*65}")
        print(f"  RESIDUE {target_resnum} ({trow['resname']}) — {args.mutation}")
        print(f"{'─'*65}")
        print(f"  pLDDT:               {trow['plddt']:.1f}")
        print(f"  Degree centrality:   {trow['degree']:.4f}  (percentile: {100*trow['degree_rank']:.1f}th)")
        print(f"  Betweenness:         {trow['betweenness']:.6f}  (percentile: {100*trow['betweenness_rank']:.1f}th)")
        print(f"  Clustering coeff:    {trow['clustering']:.4f}")
    else:
        print(f"\nWARNING: Residue {target_resnum} not found in graph nodes.")
        trow = None

    # ── Mutation simulation ───────────────────────────────────────────────────
    print(f"\n[Sim] Simulating {args.mutation} — removing edges of residue {target_resnum}...")
    mut_result = simulate_mutation(G, target_resnum)

    print(f"\n{'─'*65}")
    print(f"  LAYER 0 MUTATION DISRUPTION — {args.mutation}")
    print(f"{'─'*65}")
    if "error" in mut_result:
        print(f"  ERROR: {mut_result['error']}")
    else:
        print(f"  Edges removed (degree):         {mut_result['edges_removed']}")
        print(f"  Avg path length (WT):           {mut_result['avg_path_before']:.4f}")
        print(f"  Avg path length (mutant):       {mut_result['avg_path_after']:.4f}")
        print(f"  Δ avg path length:              {mut_result['delta_avg_path']:+.4f}")
        print(f"  Connected components (WT):      {mut_result['components_before']}")
        print(f"  Connected components (mutant):  {mut_result['components_after']}")
        print(f"  Δ components:                   {mut_result['delta_components']:+d}")

    # Layer 0 classification
    if trow is not None and "error" not in mut_result:
        tier, explanation = layer0_classification(trow, mut_result)
        print(f"\n  LAYER 0 CLASSIFICATION: {tier}")
        print(f"  {explanation}")

        result_summary = {
            "gene":         args.gene,
            "mutation":     args.mutation,
            "resnum":       target_resnum,
            "plddt":        float(trow["plddt"]),
            "degree":       float(trow["degree"]),
            "degree_pct":   float(trow["degree_rank"]),
            "betweenness":  float(trow["betweenness"]),
            "betweenness_pct": float(trow["betweenness_rank"]),
            "clustering":   float(trow["clustering"]),
            "disruption":   mut_result,
            "layer0_tier":  tier,
            "layer0_note":  explanation,
        }
        out_json = os.path.join(DATA_DIR, f"{args.gene}_{args.mutation}_layer0.json")
        with open(out_json, "w") as fh:
            json.dump(result_summary, fh, indent=2)
        print(f"\n[Graph] Layer 0 result saved to: {out_json}")

    # ── Plots ─────────────────────────────────────────────────────────────────
    cent_png = os.path.join(DATA_DIR, f"{args.gene}_centrality.png")
    plot_centrality(df_cent, target_resnum, cent_png, args.gene)

    print("\n[Graph] ✓ Complete.")


if __name__ == "__main__":
    main()
