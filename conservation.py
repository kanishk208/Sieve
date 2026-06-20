"""
conservation.py — Layer 2: Evolutionary Conservation Engine
Computes per-residue conservation scores from ortholog multiple sequence alignment.

Strategy:
  1. Fetch ortholog protein sequences from UniProt REST API
  2. Align each ortholog against human reference using Biopython PairwiseAligner
  3. Score each human residue position: fraction of orthologs with identical AA
  4. Cache results as {gene}_conservation.csv

Conservation Score Interpretation:
  1.00 = Completely invariant across all orthologs (strict evolutionary constraint)
  0.50 = Half the species have a different amino acid (moderate tolerance)
  0.00 = Every species has a different amino acid (highly tolerant / disordered)

Usage:
    from conservation import get_conservation_score, compute_conservation_profile
    score = get_conservation_score("VCP", 155)  # returns float 0.0-1.0 or None
"""

import csv
import json
import time
from pathlib import Path
from typing import Optional

import requests

try:
    from Bio import Align
    BIOPYTHON_OK = True
except ImportError:
    BIOPYTHON_OK = False

ROOT = Path(__file__).parent
DATA_DIR = ROOT / "data"

# ──────────────────────────────────────────────────────────────────────────────
# ORTHOLOG SPECIES SET
# Covers deep eukaryotic time (~1.5 billion years) without entering prokaryotes.
# Species are ordered by evolutionary distance from human.
# ──────────────────────────────────────────────────────────────────────────────
ORTHOLOG_SPECIES = {
    # Mammals (diverged ~90-180 Mya)
    "MOUSE":   {"query": "gene:VCP AND organism_id:10090 AND reviewed:true",  "label": "Mus musculus"},
    "RAT":     {"query": "gene:VCP AND organism_id:10116 AND reviewed:true",  "label": "Rattus norvegicus"},
    "DOG":     {"query": "gene:VCP AND organism_id:9615 AND reviewed:true",   "label": "Canis lupus familiaris"},
    "BOVIN":   {"query": "gene:VCP AND organism_id:9913 AND reviewed:true",   "label": "Bos taurus"},
    "PIG":     {"query": "gene:VCP AND organism_id:9823 AND reviewed:true",   "label": "Sus scrofa"},
    # Birds (diverged ~320 Mya)
    "CHICK":   {"query": "gene:VCP AND organism_id:9031 AND reviewed:true",   "label": "Gallus gallus"},
    # Amphibians (diverged ~350 Mya)
    "XENLA":   {"query": "gene:vcp AND organism_id:8355 AND reviewed:true",   "label": "Xenopus laevis"},
    # Fish (diverged ~450 Mya)
    "DANRE":   {"query": "gene:vcp AND organism_id:7955 AND reviewed:true",   "label": "Danio rerio"},
    # Invertebrates (diverged ~600-800 Mya)
    "DROME":   {"query": "gene:TER94 AND organism_id:7227 AND reviewed:true", "label": "Drosophila melanogaster"},
    "CAEEL":   {"query": "gene:cdc-48.1 AND organism_id:6239 AND reviewed:true", "label": "C. elegans"},
    # Fungi (diverged ~1 Bya)
    "YEAST":   {"query": "gene:CDC48 AND organism_id:559292 AND reviewed:true", "label": "S. cerevisiae"},
    "SCHPO":   {"query": "gene:cdc48 AND organism_id:284812 AND reviewed:true", "label": "S. pombe"},
    # Plants (diverged ~1.5 Bya)
    "ARATH":   {"query": "gene:CDC48A AND organism_id:3702 AND reviewed:true", "label": "Arabidopsis thaliana"},
}

# Human VCP reference
HUMAN_VCP_UNIPROT = "P55072"


def _fetch_uniprot_sequence(accession: str) -> Optional[str]:
    """Fetch protein sequence from UniProt by accession ID."""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    try:
        resp = requests.get(url, timeout=15)
        if resp.status_code == 200:
            lines = resp.text.strip().split("\n")
            seq = "".join(line.strip() for line in lines if not line.startswith(">"))
            return seq
    except Exception as e:
        print(f"    [!] Failed to fetch {accession}: {e}")
    return None


def _search_uniprot_sequence(query: str) -> Optional[str]:
    """Search UniProt for a sequence matching a query string."""
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": query,
        "format": "fasta",
        "size": 1,
    }
    try:
        resp = requests.get(url, params=params, timeout=15)
        if resp.status_code == 200 and resp.text.strip():
            lines = resp.text.strip().split("\n")
            seq = "".join(line.strip() for line in lines if not line.startswith(">"))
            if len(seq) > 100:  # Sanity check: VCP orthologs are ~800 AA
                return seq
    except Exception as e:
        print(f"    [!] UniProt search failed for '{query}': {e}")
    return None


def fetch_orthologs(gene: str) -> dict:
    """
    Fetch ortholog sequences for a gene from UniProt.
    Returns dict mapping species_code -> amino_acid_sequence.
    """
    print(f"  [Conservation] Fetching ortholog sequences for {gene}...")
    orthologs = {}

    for species_code, info in ORTHOLOG_SPECIES.items():
        seq = _search_uniprot_sequence(info["query"])
        if seq:
            orthologs[species_code] = seq
            print(f"    [+] {info['label']:30s} ({species_code}) -> {len(seq)} AA")
        else:
            print(f"    [-] {info['label']:30s} ({species_code}) -> not found")
        time.sleep(0.5)  # Rate limit: be polite to UniProt

    print(f"  [Conservation] Retrieved {len(orthologs)}/{len(ORTHOLOG_SPECIES)} orthologs")
    return orthologs


def align_and_score(human_seq: str, orthologs: dict) -> list:
    """
    Align each ortholog to the human sequence and compute per-residue conservation.

    Returns a list of dicts, one per human residue position:
        [{"resnum": 1, "aa": "M", "conservation": 0.85, "n_orthologs": 13, "n_identical": 11}, ...]
    """
    if not BIOPYTHON_OK:
        raise ImportError("Biopython is required. Install with: pip install biopython")

    n_orthologs = len(orthologs)
    if n_orthologs == 0:
        return []

    human_len = len(human_seq)
    # Track how many orthologs match at each human position
    identical_counts = [0] * human_len
    aligned_counts = [0] * human_len  # how many orthologs have coverage at this position

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -0.5

    for species_code, orth_seq in orthologs.items():
        try:
            alignments = aligner.align(human_seq, orth_seq)
            best = alignments[0]

            # Walk through the alignment to map ortholog residues to human positions
            # best.aligned returns pairs of (target_intervals, query_intervals)
            human_aligned = best.aligned[0]  # intervals in human (target)
            orth_aligned = best.aligned[1]    # intervals in ortholog (query)

            # Build position-level mapping from the alignment blocks
            human_positions = []
            orth_positions = []
            for (h_start, h_end), (o_start, o_end) in zip(human_aligned, orth_aligned):
                human_positions.extend(range(h_start, h_end))
                orth_positions.extend(range(o_start, o_end))

            for h_pos, o_pos in zip(human_positions, orth_positions):
                if h_pos < human_len and o_pos < len(orth_seq):
                    aligned_counts[h_pos] += 1
                    if human_seq[h_pos] == orth_seq[o_pos]:
                        identical_counts[h_pos] += 1

        except Exception as e:
            print(f"    [!] Alignment failed for {species_code}: {e}")

    # Compute conservation scores
    results = []
    for i in range(human_len):
        n_aligned = aligned_counts[i]
        n_identical = identical_counts[i]
        # Conservation = fraction of aligned orthologs with identical AA
        conservation = n_identical / n_aligned if n_aligned > 0 else 0.0

        results.append({
            "resnum": i + 1,  # 1-indexed
            "aa": human_seq[i],
            "conservation": round(conservation, 4),
            "n_orthologs_aligned": n_aligned,
            "n_identical": n_identical,
        })

    return results


def compute_conservation_profile(gene: str, uniprot_id: str, force: bool = False) -> Path:
    """
    Full pipeline: fetch orthologs, align, score, and cache results.

    Args:
        gene: Gene symbol (e.g., "VCP")
        uniprot_id: UniProt accession (e.g., "P55072")
        force: If True, recompute even if cached file exists

    Returns:
        Path to the conservation CSV file
    """
    gene_dir = DATA_DIR / gene
    gene_dir.mkdir(parents=True, exist_ok=True)
    cache_path = gene_dir / f"{gene}_conservation.csv"

    if cache_path.exists() and not force:
        print(f"  [Conservation] Using cached profile: {cache_path}")
        return cache_path

    print(f"\n{'='*60}")
    print(f"  LAYER 2: EVOLUTIONARY CONSERVATION ENGINE")
    print(f"  Gene: {gene} | UniProt: {uniprot_id}")
    print(f"{'='*60}")

    # Step 1: Fetch human reference
    print(f"  [1/3] Fetching human {gene} sequence...")
    human_seq = _fetch_uniprot_sequence(uniprot_id)
    if not human_seq:
        print(f"  [!] Could not fetch human sequence for {uniprot_id}")
        return cache_path
    print(f"    Human {gene}: {len(human_seq)} residues")

    # Step 2: Fetch orthologs
    print(f"  [2/3] Fetching ortholog sequences...")
    orthologs = fetch_orthologs(gene)

    if not orthologs:
        print(f"  [!] No orthologs found for {gene}")
        return cache_path

    # Step 3: Align and score
    print(f"  [3/3] Computing pairwise alignments and conservation scores...")
    profile = align_and_score(human_seq, orthologs)

    # Save to CSV
    with open(cache_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["resnum", "aa", "conservation", "n_orthologs_aligned", "n_identical"])
        writer.writeheader()
        writer.writerows(profile)

    print(f"  [Conservation] Profile saved: {cache_path}")
    print(f"  [Conservation] {len(profile)} residues scored across {len(orthologs)} orthologs")

    # Quick stats
    scores = [r["conservation"] for r in profile]
    avg = sum(scores) / len(scores) if scores else 0
    highly_conserved = sum(1 for s in scores if s >= 0.85)
    tolerant = sum(1 for s in scores if s <= 0.40)
    print(f"  [Conservation] Mean score: {avg:.3f}")
    print(f"  [Conservation] Highly conserved (>=0.85): {highly_conserved} residues ({100*highly_conserved/len(scores):.1f}%)")
    print(f"  [Conservation] Tolerant (<=0.40): {tolerant} residues ({100*tolerant/len(scores):.1f}%)")

    return cache_path


def get_conservation_score(gene: str, resnum: int) -> Optional[float]:
    """
    Look up the conservation score for a specific residue.

    Args:
        gene: Gene symbol (e.g., "VCP")
        resnum: 1-indexed residue number

    Returns:
        Conservation score (0.0-1.0), or None if not available
    """
    gene_dir = DATA_DIR / gene
    cache_path = gene_dir / f"{gene}_conservation.csv"

    if not cache_path.exists():
        return None

    try:
        with open(cache_path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                if int(row["resnum"]) == resnum:
                    return float(row["conservation"])
    except Exception:
        pass

    return None


def load_conservation_profile(gene: str) -> dict:
    """
    Load the full conservation profile as a dict mapping resnum -> score.
    More efficient than repeated get_conservation_score calls.
    """
    gene_dir = DATA_DIR / gene
    cache_path = gene_dir / f"{gene}_conservation.csv"
    profile = {}

    if cache_path.exists():
        try:
            with open(cache_path) as f:
                reader = csv.DictReader(f)
                for row in reader:
                    profile[int(row["resnum"])] = float(row["conservation"])
        except Exception:
            pass

    return profile


# ──────────────────────────────────────────────────────────────────────────────
# CLI entry point
# ──────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Compute evolutionary conservation profile.")
    parser.add_argument("--gene", default="VCP", help="Gene symbol")
    parser.add_argument("--uniprot", default="P55072", help="UniProt accession")
    parser.add_argument("--force", action="store_true", help="Force recomputation")
    args = parser.parse_args()

    compute_conservation_profile(args.gene, args.uniprot, force=args.force)
