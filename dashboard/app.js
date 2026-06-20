/**
 * SIEVE Dashboard — Mechanistic Variant Interpretation Engine
 *
 * Wires the 7-panel index.html to the FastAPI backend:
 *   - Runs the real pipeline (POST /api/analyze or fast POST /api/variant)
 *   - 3D structure (NGL), confidence profile (Chart.js), contact network (D3)
 *   - Layer-0 assessment, five-tier clinical report, validation suite
 */

const API = '';
const $ = (sel) => document.querySelector(sel);
const $$ = (sel) => document.querySelectorAll(sel);

const STAGE_LABELS = ['Data Foundation', 'pLDDT Audit', 'Graph Engine', 'DynaMut2', 'Clinical Report'];

let state = {
  gene: 'VCP',
  mutation: 'R155H',
  uniprot: 'P55072',
  nglStage: null,
  nglComponent: null,
  reportData: null,
  layer0Data: null,      // freshest layer0 from POST (has c_3d, layer0_note)
  centralityMap: null,   // resnum -> centrality row
  plddtChart: null,
  currentRepr: 'cartoon',
  colorScheme: 'bfactor',
};

let serverOnline = false;
let analyzeRunning = false;
let debounceTimer = null;

// ─── utilities ──────────────────────────────────────────────────────
function setApiStatus(text, kind = 'ok') {
  const el = $('#api-status-text');
  if (!el) return;
  el.textContent = text;
  el.className = kind === 'err' ? 'api-err' : kind === 'busy' ? 'api-busy' : 'api-ok';
}

function apiFetch(path, opts = {}) {
  const bust = `_=${Date.now()}`;
  const url = path.includes('?') ? `${API}${path}&${bust}` : `${API}${path}?${bust}`;
  return fetch(url, {
    cache: 'no-store',
    headers: { 'Cache-Control': 'no-cache', ...(opts.headers || {}) },
    ...opts,
  });
}

function postJSON(path, body) {
  return apiFetch(path, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body),
  });
}

function fmtPct(x, digits = 0) {
  return (x === null || x === undefined || isNaN(x)) ? '—' : `${(x * 100).toFixed(digits)}%`;
}
function fmtNum(x, digits = 2) {
  return (x === null || x === undefined || isNaN(x)) ? '—' : Number(x).toFixed(digits);
}
function residueOf(mut) {
  const m = (mut || '').match(/[A-Za-z](\d+)/);
  return m ? parseInt(m[1], 10) : null;
}
function esc(s) {
  return String(s ?? '').replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
}

function updateActiveVariantChip() {
  const chip = $('#active-variant-chip');
  if (!chip) return;
  if (state.gene && state.mutation && state.uniprot) {
    chip.textContent = `${state.gene} · ${state.mutation} · ${state.uniprot}`;
    chip.style.display = 'inline-flex';
  } else {
    chip.style.display = 'none';
  }
}

function readInputs() {
  return {
    gene: $('#input-gene').value.trim().toUpperCase(),
    mutation: $('#input-mutation').value.trim().toUpperCase(),
    uniprot: $('#input-uniprot').value.trim().toUpperCase(),
  };
}
function syncStateFromInputs() {
  Object.assign(state, readInputs());
  updateActiveVariantChip();
}

// ─── pipeline tracker ───────────────────────────────────────────────
function resetTracker() {
  $$('.tracker-step').forEach((s) => s.classList.remove('active', 'done', 'error'));
}
function setTrackerActive(idx) {
  $$('.tracker-step').forEach((s) => {
    const i = Number(s.dataset.stage);
    s.classList.remove('active', 'done', 'error');
    if (i < idx) s.classList.add('done');
    else if (i === idx) s.classList.add('active');
  });
}
function applyTrackerStages(stages) {
  if (!Array.isArray(stages)) { markTrackerAllDone(); return; }
  stages.forEach((stg) => {
    const idx = STAGE_LABELS.indexOf(stg.name);
    if (idx < 0) return;
    const step = $(`.tracker-step[data-stage="${idx}"]`);
    if (!step) return;
    step.classList.remove('active', 'done', 'error');
    if (stg.status === 'error') step.classList.add('error');
    else step.classList.add('done');
  });
}
function markTrackerAllDone() {
  $$('.tracker-step').forEach((s) => { s.classList.remove('active', 'error'); s.classList.add('done'); });
}

// ─── boot ───────────────────────────────────────────────────────────
document.addEventListener('DOMContentLoaded', async () => {
  if (window.location.protocol === 'file:') {
    const b = $('#offline-banner');
    b.style.display = 'block';
    b.innerHTML = '<strong>Open via the server:</strong> run <code>python server.py</code> and visit <a href="http://localhost:8000">http://localhost:8000</a>';
    setApiStatus('Opened as a local file — start the backend', 'err');
    return;
  }

  $('#input-gene').value = state.gene;
  $('#input-mutation').value = state.mutation;
  $('#input-uniprot').value = state.uniprot;
  updateActiveVariantChip();

  $('#btn-analyze').addEventListener('click', () => runAnalysis());
  ['input-gene', 'input-mutation', 'input-uniprot'].forEach((id) => {
    const el = $(`#${id}`);
    el.addEventListener('input', scheduleAnalysis);
    el.addEventListener('keydown', (e) => {
      if (e.key === 'Enter') { e.preventDefault(); runAnalysis(); }
    });
  });

  // structure representation buttons
  $$('.ctrl-btn[data-repr]').forEach((btn) => {
    btn.addEventListener('click', () => {
      $$('.ctrl-btn[data-repr]').forEach((b) => b.classList.remove('active'));
      btn.classList.add('active');
      state.currentRepr = btn.dataset.repr;
      updateStructureRepr();
    });
  });
  $('#color-scheme').addEventListener('change', (e) => {
    state.colorScheme = e.target.value;
    updateStructureRepr();
  });
  $('#btn-focus-mutation').addEventListener('click', focusMutation);

  // residue detail close
  $('#detail-close').addEventListener('click', () => { $('#residue-detail').style.display = 'none'; });

  // network radius
  $('#network-radius').addEventListener('input', (e) => { $('#radius-value').textContent = e.target.value; });
  $('#network-radius').addEventListener('change', (e) => {
    if (serverOnline && state.gene && state.mutation && state.uniprot) loadGraph(e.target.value);
  });

  // exports + validation
  $('#btn-export-json').addEventListener('click', exportJSON);
  $('#btn-export-md').addEventListener('click', exportMarkdown);
  $('#btn-export-brief').addEventListener('click', exportBrief);
  $('#btn-run-validation').addEventListener('click', runValidation);

  await checkHealth();
  if (serverOnline) await runAnalysis();
});

async function checkHealth() {
  try {
    const r = await apiFetch('/api/health');
    serverOnline = r.ok;
    if (r.ok) {
      $('#offline-banner').style.display = 'none';
      setApiStatus('Connected to backend', 'ok');
    } else throw new Error('health failed');
  } catch {
    serverOnline = false;
    $('#offline-banner').style.display = 'block';
    setApiStatus('Backend offline — run: python server.py', 'err');
  }
}

function scheduleAnalysis() {
  syncStateFromInputs();
  clearTimeout(debounceTimer);
  setApiStatus('Waiting…', 'busy');
  debounceTimer = setTimeout(() => runAnalysis(), 800);
}

// ─── main analysis flow ─────────────────────────────────────────────
async function runAnalysis() {
  if (!serverOnline || analyzeRunning) return;
  syncStateFromInputs();
  if (!state.gene || !state.mutation || !state.uniprot) {
    setApiStatus('Enter gene, variant and UniProt', 'err');
    return;
  }

  analyzeRunning = true;
  $('#btn-analyze').classList.add('running');
  resetTracker();
  state.centralityMap = null;

  try {
    // Decide fast path (mutation-only) vs full pipeline.
    let ready = { centrality: false, pdb: false };
    try {
      const rr = await apiFetch(`/api/ready/${state.gene}/${state.uniprot}`);
      if (rr.ok) ready = await rr.json();
    } catch { /* ignore */ }

    if (ready.centrality && ready.pdb) {
      setApiStatus('Recomputing variant…', 'busy');
      setTrackerActive(2);
      const r = await postJSON('/api/variant', { gene: state.gene, mutation: state.mutation, uniprot: state.uniprot });
      const vdata = await r.json().catch(() => ({}));
      if (!r.ok) throw new Error(vdata.detail || vdata.error || 'Variant analysis failed');
      if (vdata.layer0) state.layer0Data = vdata.layer0;
      markTrackerAllDone();
    } else {
      setApiStatus('Running full pipeline (this can take 1–3 min)…', 'busy');
      setTrackerActive(0);
      const r = await postJSON('/api/analyze', { gene: state.gene, mutation: state.mutation, uniprot: state.uniprot });
      const data = await r.json().catch(() => ({}));
      if (!r.ok || data.error) throw new Error(data.error || 'Pipeline failed');
      applyTrackerStages(data.stages);
    }

    setApiStatus('Loading results…', 'busy');
    await loadCentrality();
    await Promise.allSettled([
      loadStructure(),
      loadReport(),
      loadPlddt(),
      loadGraph($('#network-radius').value || 2),
    ]);
    setApiStatus(`Analysis complete — ${state.gene} ${state.mutation}`, 'ok');
  } catch (err) {
    console.error('Analysis error:', err);
    const step = [...$$('.tracker-step')].find((s) => s.classList.contains('active'));
    if (step) { step.classList.remove('active'); step.classList.add('error'); }
    setApiStatus(`Error: ${err.message}`, 'err');
  } finally {
    analyzeRunning = false;
    $('#btn-analyze').classList.remove('running');
  }
}

async function safeError(r) {
  try { const j = await r.json(); return j.detail || j.error; } catch { return null; }
}

// ─── centrality cache (for residue picking) ─────────────────────────
async function loadCentrality() {
  try {
    const r = await apiFetch(`/api/centrality/${state.gene}`);
    if (!r.ok) return;
    const rows = await r.json();
    const map = {};
    rows.forEach((row) => { map[row.resnum] = row; });
    state.centralityMap = map;
  } catch { /* non-fatal */ }
}

// ─── 3D structure ───────────────────────────────────────────────────
async function loadStructure() {
  try {
    if (!state.nglStage) {
      state.nglStage = new NGL.Stage('ngl-viewport', { backgroundColor: 'white' });
      window.addEventListener('resize', () => state.nglStage.handleResize());
      state.nglStage.signals.clicked.add(onResiduePick);
    }
    if (state.nglComponent) {
      state.nglStage.removeComponent(state.nglComponent);
      state.nglComponent = null;
    }

    const url = `${API}/api/structure/${state.uniprot}`;
    state.nglComponent = await state.nglStage.loadFile(url, { ext: 'pdb' });

    $('#ngl-placeholder').style.display = 'none';
    $('#structure-controls').style.display = 'flex';
    updateStructureRepr();
  } catch (err) {
    console.warn('Structure loading:', err.message);
  }
}

function updateStructureRepr() {
  const c = state.nglComponent;
  if (!c) return;
  c.removeAllRepresentations();
  c.addRepresentation(state.currentRepr, { colorScheme: state.colorScheme });

  // highlight mutation residue
  const resnum = residueOf(state.mutation);
  if (resnum != null) {
    c.addRepresentation('ball+stick', { sele: `${resnum}`, color: '#EF4444', aspectRatio: 2.5, scale: 1.6 });
    c.addRepresentation('label', {
      sele: `${resnum} and .CA`, labelType: 'format', labelFormat: '%(resname)s%(resno)s',
      color: '#0A0A0A', showBackground: true, backgroundColor: 'white', backgroundOpacity: 0.6, fontSize: 1.2,
    });
    c.autoView(`${resnum}`, 1500);
  } else {
    c.autoView();
  }
}

function focusMutation() {
  const resnum = residueOf(state.mutation);
  if (state.nglComponent && resnum != null) state.nglComponent.autoView(`${resnum}`, 1000);
}

function onResiduePick(pp) {
  if (!pp || !pp.atom) { return; }
  const resno = pp.atom.resno;
  const resname = pp.atom.resname;
  const cent = state.centralityMap ? state.centralityMap[resno] : null;
  $('#detail-resname').textContent = `${resname} ${resno}`;
  const rows = [];
  if (cent) {
    rows.push(['pLDDT', fmtNum(cent.plddt, 1)]);
    rows.push(['Degree', fmtNum(cent.degree, 0)]);
    rows.push(['Degree pct', fmtPct(cent.degree_rank, 0)]);
    rows.push(['Betweenness pct', fmtPct(cent.betweenness_rank, 1)]);
    rows.push(['Clustering', fmtNum(cent.clustering, 3)]);
  } else {
    rows.push(['Info', 'No centrality data']);
  }
  $('#detail-body').innerHTML = rows.map(([k, v]) =>
    `<div class="data-row"><span class="data-label">${esc(k)}</span><span class="data-value">${esc(v)}</span></div>`).join('');
  $('#residue-detail').style.display = 'flex';
}

// ─── confidence profile ─────────────────────────────────────────────
async function loadPlddt() {
  try {
    const r = await apiFetch(`/api/plddt/${state.gene}`);
    if (!r.ok) throw new Error('pLDDT not found');
    const rows = await r.json();
    if (!rows.length) throw new Error('empty pLDDT');

    const values = rows.map((d) => d.plddt);
    const mean = values.reduce((a, b) => a + b, 0) / values.length;
    const min = Math.min(...values);
    const max = Math.max(...values);
    const lowFrac = values.filter((v) => v < 70).length / values.length;

    $('#plddt-mean-chip').textContent = `mean ${mean.toFixed(1)}`;
    $('#plddt-stats').innerHTML = `
      <div class="stat-box"><div class="stat-value">${mean.toFixed(1)}</div><div class="stat-label">Mean pLDDT</div></div>
      <div class="stat-box"><div class="stat-value">${min.toFixed(0)}</div><div class="stat-label">Min</div></div>
      <div class="stat-box"><div class="stat-value">${max.toFixed(0)}</div><div class="stat-label">Max</div></div>
      <div class="stat-box"><div class="stat-value">${(lowFrac * 100).toFixed(0)}%</div><div class="stat-label">&lt; 70</div></div>`;

    const mutRes = residueOf(state.mutation);
    renderPlddtChart(rows, mutRes);
  } catch (err) {
    console.warn('pLDDT loading:', err.message);
  }
}

function renderPlddtChart(rows, mutRes) {
  const host = $('#plddt-chart');
  host.innerHTML = '<canvas></canvas>';
  const ctx = host.querySelector('canvas');

  const labels = rows.map((d) => d.resnum);
  const data = rows.map((d) => d.plddt);
  const colors = rows.map((d) =>
    d.plddt >= 90 ? '#3B82F6' : d.plddt >= 70 ? '#10B981' : d.plddt >= 50 ? '#F59E0B' : '#EF4444');

  const mutPoint = rows.map((d) => (mutRes != null && d.resnum === mutRes ? d.plddt : null));

  if (state.plddtChart) state.plddtChart.destroy();
  state.plddtChart = new Chart(ctx, {
    type: 'line',
    data: {
      labels,
      datasets: [
        {
          label: 'pLDDT', data, borderColor: '#0A0A0A', borderWidth: 1,
          pointRadius: 0, fill: true,
          segment: { borderColor: (c) => colors[c.p0DataIndex] },
          backgroundColor: 'rgba(10,10,10,0.04)', tension: 0.2,
        },
        {
          label: `Mutation (${state.mutation})`, data: mutPoint,
          borderColor: '#EF4444', backgroundColor: '#EF4444',
          pointRadius: 6, pointHoverRadius: 8, showLine: false,
        },
      ],
    },
    options: {
      responsive: true, maintainAspectRatio: false,
      plugins: { legend: { display: false }, tooltip: { callbacks: { title: (i) => `Residue ${i[0].label}` } } },
      scales: {
        x: { title: { display: true, text: 'Residue position' }, ticks: { maxTicksLimit: 12 } },
        y: { min: 0, max: 100, title: { display: true, text: 'pLDDT' } },
      },
    },
  });
}

// ─── contact network ────────────────────────────────────────────────
async function loadGraph(radius = 2) {
  try {
    const url = `/api/graph/${state.gene}/${state.mutation}?radius=${radius}&uniprot=${state.uniprot}`;
    const r = await apiFetch(url);
    if (!r.ok) throw new Error((await safeError(r)) || 'Graph not found');
    const data = await r.json();
    renderD3Graph(data);
    $('#radius-control').style.display = 'flex';
    $('#network-legend').style.display = 'flex';
  } catch (err) {
    console.warn('Graph loading:', err.message);
    $('#network-graph').innerHTML = `<div class="placeholder"><span>Network unavailable</span><small>${esc(err.message)}</small></div>`;
  }
}

function renderD3Graph(data) {
  const container = $('#network-graph');
  container.innerHTML = '<svg class="network-svg"></svg>';
  if (!data.nodes || !data.nodes.length) {
    container.innerHTML = '<div class="placeholder"><span>No network data</span></div>';
    return;
  }

  const width = container.clientWidth || 400;
  const height = container.clientHeight || 320;
  const tooltip = $('#tooltip');

  const svg = d3.select(container).select('svg').attr('width', width).attr('height', height);
  const links = (data.edges || []).map((e) => ({ source: e.source, target: e.target }));

  const betw = data.nodes.map((n) => n.betweenness);
  const betwMax = Math.max(...betw, 1e-9);

  const sim = d3.forceSimulation(data.nodes)
    .force('link', d3.forceLink(links).id((d) => d.id).distance(34))
    .force('charge', d3.forceManyBody().strength(-160))
    .force('center', d3.forceCenter(width / 2, height / 2))
    .force('collision', d3.forceCollide(14));

  const link = svg.append('g').selectAll('line').data(links).enter().append('line')
    .attr('stroke', '#D7D7D7').attr('stroke-width', 1).attr('opacity', 0.7);

  const node = svg.append('g').selectAll('circle').data(data.nodes).enter().append('circle')
    .attr('r', (d) => (d.isTarget ? 11 : 6))
    .attr('fill', (d) => (d.isTarget ? '#EF4444' : d.betweenness / betwMax > 0.4 ? '#3B82F6' : '#10B981'))
    .attr('stroke', '#FFF').attr('stroke-width', 1.5)
    .style('cursor', 'pointer')
    .call(d3.drag()
      .on('start', (e, d) => { if (!e.active) sim.alphaTarget(0.3).restart(); d.fx = d.x; d.fy = d.y; })
      .on('drag', (e, d) => { d.fx = e.x; d.fy = e.y; })
      .on('end', (e, d) => { if (!e.active) sim.alphaTarget(0); d.fx = null; d.fy = null; }));

  node.on('mouseover', (e, d) => {
    tooltip.style.display = 'block';
    tooltip.innerHTML = `${d.resname || ''} ${d.id}<br>betw ${fmtNum(d.betweenness, 4)} · deg ${fmtNum(d.degree, 0)} · pLDDT ${fmtNum(d.plddt, 0)}`;
  }).on('mousemove', (e) => {
    tooltip.style.left = `${e.clientX + 12}px`; tooltip.style.top = `${e.clientY + 12}px`;
  }).on('mouseout', () => { tooltip.style.display = 'none'; });

  sim.on('tick', () => {
    link.attr('x1', (d) => d.source.x).attr('y1', (d) => d.source.y)
        .attr('x2', (d) => d.target.x).attr('y2', (d) => d.target.y);
    node.attr('cx', (d) => d.x).attr('cy', (d) => d.y);
  });
}

// ─── report + assessment ────────────────────────────────────────────
async function loadReport() {
  try {
    const r = await apiFetch(`/api/report/${state.gene}/${state.mutation}`);
    if (!r.ok) throw new Error('Report not found');
    const data = await r.json();
    state.reportData = data;
    renderAssessment(data);
    renderReport(data);
  } catch (err) {
    console.warn('Report loading:', err.message);
  }
}

function renderAssessment(data) {
  // layer0_graph from the report JSON; merge freshest live data if available
  const l0rep = data.layer0_graph || {};
  const l0live = state.layer0Data || {};
  // prefer live data for fields that may be richer
  const tier     = l0live.layer0_tier  || l0rep.tier  || '—';
  const noteRaw  = l0live.layer0_note  || l0rep.note  || '';
  const c3d      = (l0live.c_3d != null) ? l0live.c_3d : l0rep.c_3d;

  // Pull pLDDT from either source
  const plddt = (l0live.plddt != null) ? l0live.plddt : data.structural_confidence?.plddt;

  // Parse FoldX ΔΔG out of the note string if DynaMut2 was skipped
  const fxMatch = noteRaw.match(/FoldX.*?ΔΔG\s*=\s*([-+]?\d+\.?\d*)\s*kcal/i);
  const fxDdg   = fxMatch ? parseFloat(fxMatch[1]) : null;

  // Tier pill colour
  const tierClass = tier.includes('HIGH') || tier.includes('PATHOGENIC') ? 'high'
                  : tier === 'MODERATE' ? 'mod' : 'low';

  // C3D block
  let c3dBlock = '';
  if (c3d != null) {
    const c3dLabel = c3d >= 0.92 ? 'Invariant core' : c3d >= 0.70 ? 'Conserved shell' : 'Tolerant region';
    const c3dClass = c3d >= 0.92 ? '#EF4444' : c3d >= 0.70 ? '#F59E0B' : '#10B981';
    c3dBlock = `
      <div class="c3d-block" style="margin-bottom:14px;padding:12px;border:1px solid var(--border);background:var(--bg-alt);">
        <div style="font-size:11px;font-weight:700;text-transform:uppercase;letter-spacing:.05em;color:var(--text-sec);margin-bottom:8px;">
          Layer 2 · Evolutionary Conservation (C3D)
        </div>
        <div class="data-row" style="border:none;padding:4px 0;background:none;">
          <span class="data-label">3D neighborhood score</span>
          <span class="data-value" style="color:${c3dClass}">${fmtNum(c3d, 3)}</span>
        </div>
        <div class="data-row" style="border:none;padding:4px 0;background:none;">
          <span class="data-label">Environment</span>
          <span class="data-value" style="color:${c3dClass}">${c3dLabel}</span>
        </div>
      </div>`;
  }

  // Note block — strip the pipe-delimited segments into readable lines
  const noteParts = noteRaw.split('|').map(s => s.trim()).filter(Boolean);
  const noteBlock = noteParts.length ? `
    <div style="margin-top:12px;padding:10px 12px;background:var(--peach-dim);border:1px solid var(--peach);font-size:11px;line-height:1.6;color:var(--text-sec);">
      ${noteParts.map(p => `<div>${esc(p)}</div>`).join('')}
    </div>` : '';

  const rows = [
    ['pLDDT confidence', plddt != null ? fmtNum(plddt, 1) : '—'],
    ['Betweenness pct', fmtPct(l0rep.betweenness_pct, 1)],
    ['Degree pct', fmtPct(l0rep.degree_pct, 1)],
    ['Contacts broken', l0rep.edges_removed ?? '—'],
    ['Δ avg path length', fmtNum(l0rep.delta_path, 4)],
    ['Components split', l0rep.delta_components ?? '—'],
    ['Clustering coeff.', fmtNum(l0rep.clustering, 3)],
  ];
  if (fxDdg != null) rows.splice(1, 0, ['FoldX ΔΔG', `${fxDdg > 0 ? '+' : ''}${fmtNum(fxDdg, 2)} kcal/mol`]);

  $('#assessment-content').innerHTML = `
    <div class="assess-tier">
      <span class="tier-pill tier-pill--${tierClass}">${esc(tier)}</span>
      <span class="tier-note">${esc(data.pipeline?.routing || 'Layer-0 contact-network assessment')}</span>
    </div>
    ${c3dBlock}
    <div class="data-list">
      ${rows.map(([k, v]) => `<div class="data-row"><span class="data-label">${esc(k)}</span><span class="data-value">${esc(String(v))}</span></div>`).join('')}
    </div>
    ${noteBlock}`;
}

function renderReport(data) {
  const cls = data.classification || {};
  const chip = $('#report-tier-chip');
  if (cls.tier != null) {
    chip.textContent = `Tier ${cls.tier}`;
    chip.dataset.tier = cls.tier;
    chip.style.display = 'inline-flex';
  } else {
    chip.style.display = 'none';
  }

  const cv  = data.clinvar;
  const lit = data.literature;
  const sigClass = cv ? (cv.search_class || '').toLowerCase() : '';
  const sigCss = sigClass.includes('patho') ? 'pathogenic' : sigClass.includes('benign') ? 'benign' : 'uncertain';

  // ΔΔG: prefer DynaMut2; then FoldX from layer0 direct field; then FoldX parsed from note string
  const dm2ddg      = data.layer1_ddg?.ddg_kcal_mol;
  const dm2tool     = data.layer1_ddg?.tool || 'DynaMut2';
  const fxDdgDirect = data.layer0_graph?.ddg_foldx ?? data.layer1_ddg?.ddg_foldx ?? null;
  const noteRaw     = state.layer0Data?.layer0_note || data.layer0_graph?.note || '';
  const fxMatch     = noteRaw.match(/FoldX.*?ΔΔG\s*=\s*([-+]?\d+\.?\d*)\s*kcal/i);
  const fxDdg       = fxDdgDirect != null ? fxDdgDirect : (fxMatch ? parseFloat(fxMatch[1]) : null);

  // Derive predictive pathway label and structural confidence badge from layer0 tier
  const l0tier   = data.layer0_graph?.tier || '';
  const physPath = l0tier.includes('THERMO_OVERTURN') ? 'PATHOGENIC (DYNAMIC_THERMO_OVERTURN)'
                 : l0tier.includes('DYNAMIC_THERMO')  ? 'PATHOGENIC (DYNAMIC_THERMO)'
                 : null;
  const structPlddt = typeof data.structural_confidence?.plddt === 'number'
                        ? data.structural_confidence.plddt : null;
  const plddtBadge = structPlddt != null
    ? `<div style="margin-top:6px"><strong>Confidence Profile:</strong> <span style="color:var(--success);font-family:var(--mono)">${fmtNum(structPlddt, 1)} pLDDT</span> <span style="font-size:11px;color:var(--text-muted)">(High Structural Fidelity)</span></div>`
    : '';

  let ddgHtml;
  if (dm2ddg != null) {
    const destab = dm2ddg <= -1.5;
    ddgHtml = `<span style="font-family:var(--mono);font-weight:700;color:${destab ? 'var(--error)' : 'var(--success)'}">${dm2ddg > 0 ? '+' : ''}${fmtNum(dm2ddg, 2)} kcal/mol</span>
               <span style="font-size:11px;color:var(--text-muted);margin-left:6px">${dm2tool}</span>
               ${destab ? '<div style="font-size:11px;color:var(--error);margin-top:4px">⚠ Destabilising (≤ −1.5 kcal/mol)</div>' : ''}`;
  } else if (fxDdg != null) {
    const destab = fxDdg >= 1.5;
    ddgHtml = `<span style="font-family:var(--mono);font-weight:700;color:${destab ? 'var(--error)' : 'var(--success)'}">${fxDdg > 0 ? '+' : ''}${fmtNum(fxDdg, 2)} kcal/mol</span>
               <span style="font-size:11px;color:var(--text-muted);margin-left:6px">FoldX (physics escalation)</span>
               ${destab ? '<div style="font-size:11px;color:var(--error);margin-top:4px">⚠ Destabilising (≥ +1.5 kcal/mol)</div>' : ''}`;
  } else {
    ddgHtml = '<span style="color:var(--text-muted);font-size:12px">Not calculated — tool unavailable or analysis skipped</span>';
  }

  $('#report-content').innerHTML = `
    <div class="report-summary">
      <div class="summary-badge"><span class="big">${esc(cls.tier ?? '—')}</span><span class="small">Tier</span></div>
      <div class="summary-text">
        <div class="label">${esc(cls.label || 'Pending')}</div>
        <div class="rationale">${esc(cls.rationale || '')}</div>
      </div>
    </div>
    <div class="evidence-grid">
      <div class="evidence-card">
        <h4>ClinVar</h4>
        <div class="ev-body">
          ${cv ? `<span class="sig-tag ${sigCss}">${esc(cv.significance || '—')}</span>
                  <div>${esc(cv.title || '')}</div>
                  ${(cv.conditions && cv.conditions.length) ? `<div class="muted" style="margin-top:6px">${esc(cv.conditions.slice(0, 2).join('; '))}</div>` : ''}`
               : '<span class="muted">No ClinVar entry found.</span>'}
        </div>
      </div>
      <div class="evidence-card">
        <h4>Stability &amp; Literature</h4>
        <div class="ev-body">
          <div><strong>ΔΔG:</strong></div>
          <div style="margin-top:4px">${ddgHtml}</div>
          ${physPath ? `<div style="margin-top:6px"><strong>Pathway:</strong> <span style="font-family:var(--mono);font-size:11px;color:var(--error)">${esc(physPath)}</span></div>` : ''}
          ${plddtBadge}
          <div style="margin-top:10px"><strong>Literature:</strong> ${lit ? `${lit.concordant_count ?? 0} concordant / ${lit.discordant_count ?? 0} discordant` : '<span class="muted">none</span>'}</div>
          ${lit && lit.summary ? `<div class="muted" style="margin-top:4px;font-size:11px">${esc(lit.summary)}</div>` : ''}
        </div>
      </div>
    </div>`;
  $('#report-actions').style.display = 'flex';
}

// ─── validation suite ───────────────────────────────────────────────
async function runValidation() {
  if (!serverOnline) return;
  const btn = $('#btn-run-validation');
  btn.classList.add('running');
  btn.disabled = true;
  $('#validation-content').innerHTML = '<div class="placeholder"><span>Running validation…</span><small>Scoring ClinVar variants through Layer 0. FoldX escalation on borderline cases can take several minutes — leave this running.</small></div>';
  setApiStatus('Running validation (may take several minutes)…', 'busy');
  try {
    const r = await apiFetch(`/api/validate/${state.gene}?uniprot=${state.uniprot}`);
    if (!r.ok) throw new Error((await safeError(r)) || 'Validation failed');
    const data = await r.json();
    renderValidation(data);
    setApiStatus(`Validation complete — ${data.n_variants} variants`, 'ok');
  } catch (err) {
    console.warn('Validation:', err.message);
    $('#validation-content').innerHTML = `<div class="placeholder"><span>Validation unavailable</span><small>${esc(err.message)}</small></div>`;
  } finally {
    btn.classList.remove('running');
    btn.disabled = false;
  }
}

function renderValidation(data) {
  const m = data.metrics || {};
  const mccHi = m.mcc_publishable ? ' highlight' : '';
  $('#validation-content').innerHTML = `
    <div class="metrics-display">
      <div class="metric-box${mccHi}"><div class="metric-number">${fmtNum(m.mcc, 3)}</div><div class="metric-label">MCC</div></div>
      <div class="metric-box"><div class="metric-number">${fmtPct(m.accuracy, 1)}</div><div class="metric-label">Accuracy</div></div>
      <div class="metric-box"><div class="metric-number">${fmtPct(m.sensitivity, 1)}</div><div class="metric-label">Sensitivity</div></div>
      <div class="metric-box"><div class="metric-number">${fmtPct(m.specificity, 1)}</div><div class="metric-label">Specificity</div></div>
    </div>
    <div class="confusion-matrix">
      <table>
        <thead><tr><th></th><th>Pred. Pathogenic</th><th>Pred. Benign</th></tr></thead>
        <tbody>
          <tr><th>Actual Pathogenic</th><td class="tp">${m.tp ?? 0} TP</td><td class="fn">${m.fn ?? 0} FN</td></tr>
          <tr><th>Actual Benign</th><td class="fp">${m.fp ?? 0} FP</td><td class="tn">${m.tn ?? 0} TN</td></tr>
        </tbody>
      </table>
    </div>
    <p class="input-hint" style="margin-top:14px">${esc(data.gene)} · ${data.n_variants} ClinVar variants · source: ${esc(data.dataset_source || 'clinvar')}${m.mcc_publishable ? ' · MCC &gt; 0.4 (publishable)' : ''}</p>`;
}

// ─── exports ────────────────────────────────────────────────────────
function exportJSON() {
  if (!state.reportData) { setApiStatus('Run an analysis first', 'err'); return; }
  const blob = new Blob([JSON.stringify(state.reportData, null, 2)], { type: 'application/json' });
  downloadBlob(blob, `${state.gene}_${state.mutation}_report.json`);
}
function exportMarkdown() {
  // Server-rendered markdown (authoritative); falls back to client error.
  const url = `${API}/api/report/${state.gene}/${state.mutation}/markdown`;
  triggerDownload(url);
}
function exportBrief() {
  const url = `${API}/api/report/${state.gene}/${state.mutation}/brief`;
  window.open(url, '_blank');
}
function triggerDownload(url) {
  const a = document.createElement('a');
  a.href = url; a.download = '';
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
}
function downloadBlob(blob, filename) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url; a.download = filename;
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
  URL.revokeObjectURL(url);
}
