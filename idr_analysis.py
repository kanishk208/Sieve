"""
idr_analysis.py — Phase 3 Strategy C: IDR / low-confidence functional track.
Queries ELM (motifs) and UniProt (PTM / disorder) when pLDDT < 50.
"""
from __future__ import annotations

import json
import os
import re
from pathlib import Path
from typing import Optional

import requests

DATA_DIR = Path(__file__).parent / "data"
ELM_SEARCH_URL = "http://elm.eu.org/start_search/{uniprot}.tsv"
UNIPROT_JSON = "https://rest.uniprot.org/uniprotkb/{uniprot}.json"


def _motif_overlaps_residue(start: int, end: int, resnum: int) -> bool:
    return start <= resnum <= end


def query_elm_motifs(uniprot: str, resnum: int, timeout: int = 25) -> list[dict]:
    """Return ELM motif instances overlapping the target residue."""
    url = ELM_SEARCH_URL.format(uniprot=uniprot)
    try:
        resp = requests.get(url, timeout=timeout)
        if resp.status_code == 429:
            return [{"elm_id": "RATE_LIMIT", "note": "ELM API rate limit — retry in 3 minutes"}]
        if resp.status_code != 200:
            return []
    except requests.RequestException:
        return []

    motifs = []
    lines = resp.text.strip().splitlines()
    if len(lines) < 2:
        return motifs

    header = lines[0].split("\t")
    for line in lines[1:]:
        cols = line.split("\t")
        if len(cols) < len(header):
            continue
        row = dict(zip(header, cols))
        try:
            start = int(row.get("start", row.get("Start", 0)))
            end = int(row.get("stop", row.get("Stop", row.get("end", 0))))
        except (TypeError, ValueError):
            continue
        if not _motif_overlaps_residue(start, end, resnum):
            continue
        motifs.append({
            "elm_id": row.get("elm_identifier", row.get("ELMIdentifier", "unknown")),
            "start": start,
            "end": end,
            "annotated": row.get("is_annotated", ""),
            "sequence": row.get("sequence", ""),
        })
    return motifs


def query_uniprot_ptm(uniprot: str, resnum: int, timeout: int = 25) -> dict:
    """PTM and disorder features from UniProt REST (PhosphoSitePlus proxy via annotations)."""
    out = {"ptm_sites": [], "disorder_regions": [], "motif_features": []}
    try:
        resp = requests.get(UNIPROT_JSON.format(uniprot=uniprot), timeout=timeout)
        resp.raise_for_status()
        data = resp.json()
    except requests.RequestException:
        return out

    for feat in data.get("features", []):
        loc = feat.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")
        if start is None or end is None:
            continue
        ftype = feat.get("type", "")
        desc = feat.get("description", "") or feat.get("ligand", {}).get("name", "")
        entry = {"type": ftype, "desc": desc, "start": start, "end": end}

        if ftype in ("Modified residue", "Cross-link", "Glycosylation", "Disulfide bond"):
            if _motif_overlaps_residue(start, end, resnum):
                out["ptm_sites"].append(entry)
        elif ftype == "Region" and "disorder" in (desc or "").lower():
            out["disorder_regions"].append(entry)
        elif ftype in ("Motif", "Binding site", "Site") and _motif_overlaps_residue(start, end, resnum):
            out["motif_features"].append(entry)

    return out


def analyze_idr_track(
    gene: str,
    mutation: str,
    uniprot: str,
    resnum: int,
    plddt: float,
) -> dict:
    """
    Run disordered-region functional analysis; skip DynaMut2.
    """
    elm = query_elm_motifs(uniprot, resnum)
    ptm = query_uniprot_ptm(uniprot, resnum)

    slm_disrupted = len(elm) > 0 and not any(m.get("elm_id") == "RATE_LIMIT" for m in elm)
    ptm_disrupted = len(ptm["ptm_sites"]) > 0 or len(ptm["motif_features"]) > 0
    idr_pathogenic_signal = slm_disrupted or ptm_disrupted

    notes = []
    if plddt < 50:
        notes.append(
            f"WARNING: Target pLDDT = {plddt:.1f} (< 50). AlphaFold confidence is too low for "
            "reliable graph or ΔΔG analysis. DynaMut2 was not run."
        )
    if slm_disrupted:
        notes.append(f"ELM: {len(elm)} linear motif instance(s) overlap residue {resnum}.")
    if ptm_disrupted:
        notes.append(
            f"UniProt: {len(ptm['ptm_sites'])} PTM/modification site(s) and "
            f"{len(ptm['motif_features'])} motif/binding feature(s) at this position."
        )
    if not notes:
        notes.append("No ELM or PTM overlap detected; IDR pathogenesis cannot be confirmed computationally.")

    result = {
        "gene": gene,
        "mutation": mutation,
        "uniprot": uniprot,
        "resnum": resnum,
        "plddt": plddt,
        "track": "disordered",
        "elm_motifs": elm,
        "ptm_analysis": ptm,
        "idr_pathogenic_signal": idr_pathogenic_signal,
        "layer0_skipped": True,
        "dynamut2_skipped": True,
        "clinical_note": " ".join(notes),
    }

    out_path = DATA_DIR / f"{gene}_{mutation}_idr.json"
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    return result
