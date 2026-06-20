"""
clinical_brief.py — Phase 3 printer-friendly clinical consultation brief (HTML).
"""
from __future__ import annotations

from datetime import datetime, timezone

from report_generator import TIERS, gather_evidence, classify_tier


def build_clinical_brief_html(gene: str, mutation: str, uniprot: str) -> str:
    ev = gather_evidence(gene, mutation, uniprot)
    l0 = ev.get("layer0") or {}
    l1 = ev.get("layer1") or {}
    idr = ev.get("idr") or {}

    idr_evidence = bool(idr.get("idr_pathogenic_signal"))
    tier, rationale = classify_tier(
        plddt=ev.get("plddt"),
        layer0_tier=l0.get("layer0_tier"),
        ddg=l1.get("ddg_kcal_mol"),
        betweenness_pct=l0.get("betweenness_pct"),
        idr_evidence=idr_evidence,
    )
    tier_label, _ = TIERS[tier]
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    ddg = l1.get("ddg_kcal_mol")
    ddg_display = f"{ddg:+.3f} kcal/mol" if ddg is not None else "Not calculated (blocked or unavailable)"

    idr_block = ""
    if idr:
        idr_block = f"""
        <h2>IDR / Functional Track</h2>
        <p>{idr.get('clinical_note', '')}</p>
        <ul>
          <li>ELM motifs at site: {len(idr.get('elm_motifs', []))}</li>
          <li>PTM features at site: {len((idr.get('ptm_analysis') or {}).get('ptm_sites', []))}</li>
        </ul>
        """

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>Clinical Brief — {gene} {mutation}</title>
  <style>
    @media print {{ @page {{ margin: 1.5cm; }} }}
    body {{ font-family: Georgia, serif; max-width: 800px; margin: 40px auto; color: #222; line-height: 1.5; }}
    h1 {{ font-size: 22px; border-bottom: 2px solid #C67B4E; padding-bottom: 8px; }}
    h2 {{ font-size: 16px; color: #5C4D3C; margin-top: 24px; }}
    table {{ width: 100%; border-collapse: collapse; margin: 12px 0; }}
    th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; font-size: 13px; }}
    th {{ background: #f5f0ea; }}
    .tier-box {{ background: #fff8f0; border-left: 4px solid #C67B4E; padding: 16px; margin: 16px 0; }}
    .disclaimer {{ font-size: 11px; color: #666; margin-top: 32px; border-top: 1px solid #ddd; padding-top: 12px; }}
    .meta {{ font-size: 12px; color: #666; }}
  </style>
</head>
<body>
  <p class="meta">Digital Patient Twin · Clinical Consultation Brief · {now}</p>
  <h1>{gene}:{mutation} — {ev.get('protein_name', 'Protein')}</h1>
  <p>UniProt: {uniprot} · Residue: {ev.get('resnum', '—')}</p>

  <div class="tier-box">
    <strong>Classification: Tier {tier}</strong> — {tier_label}<br>
    {rationale}
  </div>

  <h2>Structural Confidence</h2>
  <table>
    <tr><th>Target pLDDT</th><td>{ev.get('plddt', 'N/A')}</td></tr>
    <tr><th>Source</th><td>AlphaFold</td></tr>
  </table>

  <h2>Graph Analytics (Layer 0)</h2>
  <table>
    <tr><th>Tier</th><td>{l0.get('layer0_tier', 'N/A')}</td></tr>
    <tr><th>Betweenness percentile</th><td>{(l0.get('betweenness_pct') or 0) * 100:.0f}th</td></tr>
    <tr><th>Edges removed</th><td>{(l0.get('disruption') or {}).get('edges_removed', '—')}</td></tr>
  </table>
  <p>{l0.get('layer0_note', '')}</p>

  <h2>Thermodynamic Stability (Layer 1)</h2>
  <table>
    <tr><th>ΔΔG</th><td>{ddg_display}</td></tr>
    <tr><th>Tool</th><td>{l1.get('tool', '—')}</td></tr>
    <tr><th>Signal</th><td>{l1.get('clinical_signal', '—')}</td></tr>
  </table>

  {idr_block}

  <p class="disclaimer">
    <strong>Disclaimer:</strong> Computational decision support only. Not a diagnostic device.
    Results require experimental validation before clinical action.
  </p>
  <script>if (window.location.search.includes('print=1')) window.print();</script>
</body>
</html>"""
