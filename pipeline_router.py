"""
pipeline_router.py — Phase 3 pLDDT gatekeeper and analysis track routing.
"""
from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import Optional

import pandas as pd

from plddt_audit import plddt_routing, plddt_tier


class AnalysisTrack(str, Enum):
    STRUCTURED = "structured"           # pLDDT >= 70
    LOW_CONFIDENCE = "low_confidence"   # 50 <= pLDDT < 70
    DISORDERED = "disordered"           # pLDDT < 50


def gatekeeper(plddt: Optional[float]) -> tuple[AnalysisTrack, str, bool]:
    """
    Returns (track, routing_message, run_dynamut2).
    """
    if plddt is None:
        return AnalysisTrack.DISORDERED, "pLDDT unknown — insufficient data for structural ddG.", False
    if plddt >= 70:
        return AnalysisTrack.STRUCTURED, plddt_routing(plddt), True
    if plddt >= 50:
        return AnalysisTrack.LOW_CONFIDENCE, plddt_routing(plddt), False
    return AnalysisTrack.DISORDERED, plddt_routing(plddt), False


def get_residue_plddt(gene: str, resnum: int, data_dir: Path) -> Optional[float]:
    path = data_dir / f"{gene}_plddt_per_residue.csv"
    if not path.exists():
        return None
    df = pd.read_csv(path)
    row = df[df["resnum"] == resnum]
    if row.empty:
        return None
    return float(row.iloc[0]["plddt"])


def evaluate_track_with_window(gene: str, resnum: int, data_dir: Path, fallback_plddt: Optional[float] = None) -> tuple[AnalysisTrack, str, bool, float]:
    """
    Evaluates a localized 3-residue window to prevent lone terminal residues
    from prematurely short-circuiting into the Disordered track.
    Returns (track, routing_msg, run_dynamut2, effective_plddt)
    """
    plddt_current = get_residue_plddt(gene, resnum, data_dir)
    if plddt_current is None:
        plddt_current = fallback_plddt

    if plddt_current is None:
        return gatekeeper(None) + (0.0,)

    plddt_prev = get_residue_plddt(gene, resnum - 1, data_dir)
    plddt_next = get_residue_plddt(gene, resnum + 1, data_dir)

    vals = [v for v in [plddt_prev, plddt_current, plddt_next] if v is not None]
    local_window_avg = sum(vals) / len(vals) if vals else plddt_current

    # Smooth the drop-off if the local neighborhood is structured
    effective_plddt = plddt_current
    if plddt_current < 50 and local_window_avg >= 50:
        effective_plddt = local_window_avg

    track, routing_msg, run_ddg = gatekeeper(effective_plddt)
    return track, routing_msg, run_ddg, plddt_current
