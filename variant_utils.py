"""
variant_utils.py — Parse ClinVar records and build dynamic validation datasets.
"""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Optional

import networkx as nx

AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*", "Stop": "*",
}

HGVS_PROTEIN_RE = re.compile(
    r"p\.(?:Ter|Stop)?([A-Za-z]{3})(\d+)([A-Za-z]{3}|Ter|Stop|\*)",
    re.IGNORECASE,
)
ONE_LETTER_MUT_RE = re.compile(r"^([A-Z\*])(\d+)([A-Z\*])$")
PROTEIN_CHANGE_RE = re.compile(r"([A-Z\*])(\d+)([A-Z\*])")


def normalize_significance(sig: str) -> str:
    """Map ClinVar significance text to pathogenic | benign | uncertain."""
    s = (sig or "").lower()
    if "conflict" in s:
        return "uncertain"
    if "likely benign" in s or "likely_benign" in s or s == "likely_benign":
        return "benign"
    if "likely pathogenic" in s or "likely_pathogenic" in s:
        return "pathogenic"
    if "pathogenic" in s:
        return "pathogenic"
    if "benign" in s:
        return "benign"
    return "uncertain"


def hgvs_to_mutation(wt: str, pos: str, mut: str) -> Optional[str]:
    wt1 = AA3_TO_1.get(wt.capitalize() if len(wt) == 3 else wt.upper(), wt.upper()[:1])
    mut_upper = mut.capitalize()
    if mut_upper in ("Ter", "Stop", "*"):
        mut1 = "*"
    else:
        mut1 = AA3_TO_1.get(mut_upper, mut_upper[:1].upper())
    if len(wt1) != 1 or len(mut1) != 1:
        return None
    return f"{wt1}{pos}{mut1}"


def parse_mutation_string(text: str) -> Optional[tuple[str, int]]:
    """Return (mutation, resnum) from HGVS protein or 1-letter notation."""
    if not text:
        return None

    m = HGVS_PROTEIN_RE.search(text)
    if m:
        mutation = hgvs_to_mutation(m.group(1), m.group(2), m.group(3))
        if mutation:
            return mutation, int(m.group(2))

    m = ONE_LETTER_MUT_RE.match(text.strip().upper())
    if m:
        return f"{m.group(1)}{m.group(2)}{m.group(3)}", int(m.group(2))

    for chunk in text.replace(" ", "").split(","):
        m = PROTEIN_CHANGE_RE.search(chunk.strip())
        if m:
            return f"{m.group(1)}{m.group(2)}{m.group(3)}", int(m.group(2))

    return None


def parse_clinvar_record(record: dict, gene: str) -> Optional[dict]:
    """Extract a validation-ready variant dict from a ClinVar summary record."""
    gene_field = (record.get("gene_id") or record.get("gene_sort") or "").upper()
    if gene_field and gene_field != gene.upper():
        return None

    if record.get("search_class") in ("pathogenic", "benign"):
        sig = record["search_class"]
    else:
        sig = normalize_significance(record.get("significance", ""))
    if sig == "uncertain":
        return None

    parsed = None
    title = record.get("title", "")
    parsed = parse_mutation_string(title)
    if not parsed:
        parsed = parse_mutation_string(record.get("protein_change", ""))

    if not parsed:
        return None

    mutation, resnum = parsed
    return {
        "mutation": mutation,
        "resnum": resnum,
        "known_pathogenicity": sig,
        "source": f"ClinVar {record.get('variation_id', '')}".strip(),
    }


def load_clinvar_json(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with open(path) as f:
        data = json.load(f)
    return data if isinstance(data, list) else []


def load_curated_benchmark(gene: str, data_dir: Path) -> list[dict]:
    """Optional expert-curated benchmark (e.g. validation_dataset.json for VCP)."""
    gene_path = data_dir / f"{gene}_validation_benchmark.json"
    legacy = data_dir / "validation_dataset.json"
    path = gene_path if gene_path.exists() else (legacy if gene.upper() == "VCP" and legacy.exists() else None)
    if not path:
        return []
    with open(path) as f:
        return json.load(f)


def build_validation_dataset(
    gene: str,
    G: nx.Graph,
    clinvar_path: Path,
    *,
    data_dir: Optional[Path] = None,
    max_per_class: int = 25,
) -> list[dict]:
    """
    Build a balanced validation set from ClinVar variants that lie in the contact graph.
    Merges optional curated benchmarks when ClinVar benign missense coverage is sparse.
    Deduplicates by residue number (keeps first seen per class).
    """
    records = load_clinvar_json(clinvar_path)
    pathogenic: list[dict] = []
    benign: list[dict] = []
    seen_p: set[int] = set()
    seen_b: set[int] = set()

    for rec in records:
        entry = parse_clinvar_record(rec, gene)
        if not entry:
            continue
        resnum = entry["resnum"]
        if resnum not in G.nodes():
            continue

        if entry["known_pathogenicity"] == "pathogenic" and resnum not in seen_p:
            pathogenic.append(entry)
            seen_p.add(resnum)
        elif entry["known_pathogenicity"] == "benign" and resnum not in seen_b:
            benign.append(entry)
            seen_b.add(resnum)

    if data_dir and len(benign) < max(5, max_per_class // 2):
        for entry in load_curated_benchmark(gene, data_dir):
            sig = normalize_significance(entry.get("known_pathogenicity", ""))
            if sig == "uncertain":
                continue
            resnum = entry.get("resnum")
            if resnum not in G.nodes():
                continue
            bucket = pathogenic if sig == "pathogenic" else benign
            seen = seen_p if sig == "pathogenic" else seen_b
            if resnum in seen:
                continue
            bucket.append({
                "mutation": entry["mutation"],
                "resnum": resnum,
                "known_pathogenicity": sig,
                "source": entry.get("source", "curated benchmark"),
            })
            seen.add(resnum)

    pathogenic = pathogenic[:max_per_class]
    benign = benign[:max_per_class]
    return pathogenic + benign
