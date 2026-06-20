"""
fetch_structure.py — Phase 1, Step 3 & 4
Downloads the AlphaFold structure for a target protein and fetches
pathogenic/likely-pathogenic variants from ClinVar via the NCBI E-utilities API.

Usage:
    python fetch_structure.py
    python fetch_structure.py --uniprot P55072 --gene VCP

Defaults to VCP (Valosin-Containing Protein), UniProt P55072.
"""

import argparse
import json
import os
import time
import requests

# ─── Configuration ────────────────────────────────────────────────────────────

ALPHAFOLD_BASE = "https://alphafold.ebi.ac.uk/api/prediction"
NCBI_ESEARCH   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_ESUMMARY  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
NCBI_EFETCH    = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
UNIPROT_BASE   = "https://rest.uniprot.org/uniprotkb"

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
os.makedirs(DATA_DIR, exist_ok=True)


# ─── AlphaFold Fetch ──────────────────────────────────────────────────────────

def fetch_alphafold_structure(uniprot_id: str) -> str:
    """
    Download the latest AlphaFold PDB for a given UniProt accession.
    Returns the local path to the saved PDB file.
    """
    print(f"\n[AlphaFold] Querying prediction metadata for {uniprot_id}...")
    url = f"{ALPHAFOLD_BASE}/{uniprot_id}"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()

    entries = resp.json()
    if not entries:
        raise ValueError(f"No AlphaFold entry found for UniProt ID: {uniprot_id}")

    entry   = entries[0]
    pdb_url = entry["pdbUrl"]
    version = entry.get("latestVersion", "v4")

    print(f"[AlphaFold] Found entry: {entry['entryId']}  (version {version})")
    print(f"[AlphaFold] Downloading PDB from: {pdb_url}")

    pdb_resp = requests.get(pdb_url, timeout=60)
    pdb_resp.raise_for_status()

    local_path = os.path.join(DATA_DIR, f"AF-{uniprot_id}-F1-model_{version}.pdb")
    with open(local_path, "w") as fh:
        fh.write(pdb_resp.text)

    size_kb = os.path.getsize(local_path) / 1024
    print(f"[AlphaFold] Saved to: {local_path}  ({size_kb:.1f} KB)")
    return local_path


# ─── UniProt Domain Annotations ───────────────────────────────────────────────

def fetch_uniprot_annotations(uniprot_id: str) -> dict:
    """
    Fetch sequence, length, and domain annotations from UniProt REST API.
    Returns a dict with keys: sequence, length, domains, gene_name, protein_name.
    """
    print(f"\n[UniProt] Fetching annotations for {uniprot_id}...")
    url = f"{UNIPROT_BASE}/{uniprot_id}.json"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()

    data = resp.json()

    # Extract sequence
    seq_info  = data.get("sequence", {})
    sequence  = seq_info.get("value", "")
    length    = seq_info.get("length", 0)

    # Gene name
    gene_names = data.get("genes", [{}])
    gene_name  = gene_names[0].get("geneName", {}).get("value", "Unknown") if gene_names else "Unknown"

    # Protein name
    prot_desc  = data.get("proteinDescription", {})
    prot_name  = (
        prot_desc.get("recommendedName", {})
                 .get("fullName", {})
                 .get("value", "Unknown")
    )

    # Domain/feature annotations
    features = data.get("features", [])
    domains  = [
        {
            "type":  f["type"],
            "desc":  f.get("description", ""),
            "start": f["location"]["start"]["value"],
            "end":   f["location"]["end"]["value"],
        }
        for f in features
        if f["type"] in ("Domain", "Region", "Active site", "Binding site", "Motif")
    ]

    result = {
        "uniprot_id":   uniprot_id,
        "gene_name":    gene_name,
        "protein_name": prot_name,
        "length":       length,
        "sequence":     sequence,
        "domains":      domains,
    }

    print(f"[UniProt] Gene: {gene_name} | Protein: {prot_name}")
    print(f"[UniProt] Length: {length} residues | Annotated features: {len(domains)}")

    # Save to disk
    out_path = os.path.join(DATA_DIR, f"{uniprot_id}_annotations.json")
    with open(out_path, "w") as fh:
        json.dump(result, fh, indent=2)
    print(f"[UniProt] Saved annotations to: {out_path}")

    return result


# ─── ClinVar Fetch ────────────────────────────────────────────────────────────

def _clinvar_search_ids(search_term: str, max_results: int) -> list[str]:
    search_resp = requests.get(NCBI_ESEARCH, params={
        "db": "clinvar",
        "term": search_term,
        "retmax": max_results,
        "retmode": "json",
    }, timeout=30)
    search_resp.raise_for_status()
    return search_resp.json()["esearchresult"].get("idlist", [])


def _clinvar_fetch_summaries(id_list: list[str]) -> dict:
    if not id_list:
        return {}
    summary_resp = requests.get(NCBI_ESUMMARY, params={
        "db": "clinvar",
        "id": ",".join(id_list),
        "retmode": "json",
    }, timeout=60)
    summary_resp.raise_for_status()
    return summary_resp.json().get("result", {})


def fetch_clinvar_variants(gene_symbol: str, max_results: int = 150) -> list[dict]:
    """
    Fetch pathogenic and benign ClinVar missense variants for a gene.
    Saves a merged JSON file used by the validation suite.
    """
    gene = gene_symbol.upper()
    print(f"\n[ClinVar] Fetching labeled variants for gene: {gene}...")

    searches = [
        (
            "pathogenic",
            f"{gene}[gene] AND (pathogenic[CLNSIG] OR likely_pathogenic[CLNSIG]) AND missense variant",
        ),
        (
            "benign",
            f"{gene}[gene] AND (benign[CLNSIG] OR likely_benign[CLNSIG]) AND missense variant",
        ),
    ]

    id_to_class: dict[str, str] = {}
    for label, term in searches:
        ids = _clinvar_search_ids(term, max_results)
        print(f"[ClinVar] {label}: {len(ids)} missense records (up to {max_results})")
        for vid in ids:
            id_to_class.setdefault(vid, label)
        time.sleep(0.34)

    id_list = list(id_to_class.keys())
    if not id_list:
        print("[ClinVar] No variants found.")
        return []

    result_set = _clinvar_fetch_summaries(id_list)
    variants = []

    for var_id in id_list:
        if var_id not in result_set or var_id == "uids":
            continue
        v = result_set[var_id]

        if (v.get("gene_sort") or "").upper() != gene:
            continue

        germline = v.get("germline_classification") or {}
        sig_desc = germline.get("description", "Unknown")
        traits = germline.get("trait_set") or v.get("trait_set") or []
        conditions = [t.get("trait_name", "") for t in traits if t.get("trait_name")]

        variants.append({
            "variation_id": var_id,
            "title": v.get("title", ""),
            "significance": sig_desc,
            "search_class": id_to_class[var_id],
            "conditions": conditions,
            "molecular_consequence": v.get("molecular_consequence_list", []),
            "protein_change": v.get("protein_change", ""),
            "gene_id": v.get("gene_sort", gene),
            "location": v.get("location_sort", ""),
        })

    patho_n = sum(1 for v in variants if "pathogenic" in v["significance"].lower())
    benign_n = sum(1 for v in variants if "benign" in v["significance"].lower())
    print(f"[ClinVar] Parsed {len(variants)} {gene} variants ({patho_n} pathogenic, {benign_n} benign).")

    out_path = os.path.join(DATA_DIR, f"{gene}_clinvar_variants.json")
    with open(out_path, "w") as fh:
        json.dump(variants, fh, indent=2)
    print(f"[ClinVar] Saved to: {out_path}")

    return variants


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Fetch AlphaFold structure and ClinVar variants.")
    parser.add_argument("--uniprot", default="P55072", help="UniProt accession (default: P55072 = VCP)")
    parser.add_argument("--gene",    default="VCP",    help="Gene symbol for ClinVar search (default: VCP)")
    args = parser.parse_args()

    print("=" * 65)
    print("  Digital Patient Twin — Phase 1: Data Foundation")
    print(f"  Target: {args.gene} ({args.uniprot})")
    print("=" * 65)

    # 1. Fetch AlphaFold structure
    pdb_path = fetch_alphafold_structure(args.uniprot)

    # 2. Fetch UniProt annotations
    annotations = fetch_uniprot_annotations(args.uniprot)

    # 3. Fetch ClinVar variants
    variants = fetch_clinvar_variants(args.gene)

    print("\n" + "=" * 65)
    print("  Phase 1 Data Fetch Complete")
    print(f"  PDB file:    {os.path.basename(pdb_path)}")
    print(f"  Protein:     {annotations['protein_name']}")
    print(f"  Length:      {annotations['length']} residues")
    print(f"  Domains:     {len(annotations['domains'])} annotated features")
    print(f"  ClinVar VUS: {len(variants)} pathogenic/likely-pathogenic variants")
    print("=" * 65)
    print("\nNext step: run  python plddt_audit.py")


if __name__ == "__main__":
    main()
