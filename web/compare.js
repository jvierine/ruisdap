// Comparison of Millstone analyses of one day, one parameter per row.
// Each row shows the dataset chosen for it. With the row's "Compare" toggle
// on, it also draws a second dataset on the same time and range axes and the
// same colour limits: the second map lies on top and is clipped left of a
// vertical bar, so the left part shows the first dataset and the right part
// the second. Dragging a bar moves the split on every compared row; pan and
// zoom stay shared through the dataset time group.
// Only display happens here: no measurement is changed, and where the two
// analyses use different range gates, both are shown on the finer grid with
// each cell taking the value of the gate it lies in (no interpolation).
import init from './pkg/guisdap_web.js';
import {fieldPanel} from './scientific-plots.js';
import {EISCAT_LIMITS} from './plot-style.js';

const el = id => document.getElementById(id), query = new URLSearchParams(location.search);
const antenna = query.get('antenna') || 'ZENITH';
const json = async u => {const r = await fetch(u, {cache: 'no-store'}); if (!r.ok) throw Error(`${r.status}: ${u}`); return r.json();};
const floats = async u => {const r = await fetch(u); if (!r.ok) throw Error(`${r.status}: ${u}`); return new Float32Array(await r.arrayBuffer());};
// Each map uses two WebGL contexts (field and colour bar) and every row has
// two maps; browsers keep about 16 contexts and silently drop the oldest, so
// at most four rows are shown at once.
const MAX_ROWS = 4;
const ROWS = [
  {name: 'Ne / m⁻³', p: 0, mode: 1, palette: 1, colormap: 'eiscat'},
  {name: 'Te / K', p: 2, mode: 0, palette: 1, colormap: 'eiscat'},
  {name: 'Ti / K', p: 1, mode: 0, palette: 1, colormap: 'eiscat'},
  {name: 'Velocity / m s⁻¹', p: 3, mode: 0, palette: 0, colormap: 'seismic'},
  {name: 'Raw Ne / m⁻³', p: 4, mode: 1, palette: 1, colormap: 'plasma'},
];
// each row's key in the URL
const KEYS = ['ne', 'te', 'ti', 'vi', 'rawne'];

// The eclipse of 8 April 2024 at Millstone Hill, at the ground (isr_analysis
// eclipse_timing.py, astropy; 94 % of the Sun's disk covered at maximum),
// seconds after 00 UTC: drawn as dashed vertical lines on every map.
const ECLIPSE = {
  '2024-04-08': [[18 * 3600 + 15 * 60 + 50, 'first contact'], [19 * 3600 + 29 * 60 + 30, 'maximum, 94 % obscured'],
                 [20 * 3600 + 38 * 60 + 50, 'last contact']],
};
const hhmm = s => `${String(Math.floor(s / 3600)).padStart(2, '0')}:${String(Math.floor(s / 60) % 60).padStart(2, '0')}`;

// One to three sentences per dataset, by analysis id (the info panel).
const ABOUT = [
  [/^madrigal/, "Madrigal analysed fits from Millstone Hill's own processing, as Juha imported them: their original times (mostly 3 min windows, with some 10 s records) and slant-range cells (18 km in mode 115, 4.5 km in mode 97), no new fit. Calibration preliminary (final calibration pending). Velocity positive away from the radar; the sign of isr_analysis's velocity is not independently verified."],
  [/^historical$/, "Products of the older Python pipeline (30 µs gates), imported by Juha and calibrated with the plasma lines. No satellite removal."],
  [/^native-(full-)?20261019$/, "An early test of Juha's Rust pipeline on a sample of pulses only (16 lags, 200–540 km). Incomplete; not for comparisons."],
  [/^native-direct-live/, "Juha's Rust pipeline, reading the raw recordings directly and subtracting every CFAR-detected satellite echo from the voltages. Fitted window by window while it ran; its subtraction is still under Juha's review."],
  [/^native-debris400-.*live/, "The same as the strict Rust run (echoes at least 3 dB above the CFAR threshold, above 400 km), published window by window while it ran; once complete it matches that run."],
  [/^native-debris400/, "Juha's Rust pipeline with a stricter subtraction: only echoes at least 3 dB above the CFAR threshold and above 400 km altitude are subtracted from the voltages."],
  [/^isr-v2rfi-loo/, "Candidate, for Henrik's decision: the final product with variance weights that leave each row out of its own weight (Memo 51), which removes a 17-25 % low bias of the ACF and the bias of the impulse-row removal. The density constant is not yet recalibrated for it, so its Ne is about 10-15 % high."],
  [/^isr-v2rfi-190s/, "Our final product refitted in one window per zenith cycle (about 190 s), the windows of Juha's analyses, for a comparison at matched times (range cells still differ: 9 km against 4.5 km; Memo 51)."],
  [/^isr-nosat-190s/, "Our inversion without satellite handling (variant 3), refitted in one window per zenith cycle (about 190 s), as Juha's analyses."],
  [/^isr-nosat/, "Our variance-weighted lag profile inversion with a background column; satellite echoes are left in (variant 3, Memo 42)."],
  [/^isr-detected/, "Our inversion with each detected satellite echo modelled as an extra unknown, on the pulses where the detector found it (variant 5, Memo 42)."],
  [/^isr-tracked/, "As variant 5, with pass tracking: detections are linked into passes and every pulse of a pass gets a column, also those the detector missed (Memo 45)."],
  [/^isr-v2rfi/, "Our final product: pass tracking with a lower detection threshold (false-alarm probability 1e-6, faint detections kept only when tracking links them into a pass; Memo 48), and the lagged products hit by power-line impulses removed from the inversion (Memo 47; Memo 51 finds that this removal lowers Ne by a few per cent where many rows are removed, through the variance weights)."],
];
const about = id => (ABOUT.find(([re]) => re.test(id)) || [null, ''])[1];

// Ne from the fitted proxy with the finite-Debye ion-power correction, as
// millstone.js computes it for calibrated quicklook analyses.
function calibratedDensity(m, proxy, te, ti) {
  if (!(proxy > 0 && te > 0 && ti > 0)) return NaN;
  const k = 4 * Math.PI * m.radar_frequency_hz / 299792458;
  const a = k * k * 8.8541878188e-12 * 1.380649e-23 * te / (1.602176634e-19 ** 2), ratio = 1 + te / ti;
  const forward = ne => ne * ratio / ((1 + a / ne) * (ratio + a / ne));
  let lo = 0, hi = Math.max(proxy, 1);
  while (forward(hi) < proxy) hi *= 2;
  for (let i = 0; i < 65; i++) {const mid = (lo + hi) / 2; if (forward(mid) < proxy) lo = mid; else hi = mid;}
  return (lo + hi) / 2;
}

async function load(choice) {
  const root = choice.data_root.endsWith('/') ? choice.data_root : choice.data_root + '/';
  const m = await json(root + 'manifest.json');
  if (m.schema === 'ruisdap-millstone-comparison-v1') return loadMadrigal(choice, root, m);
  // Native (Rust) analyses store absolute Ne in their parameters and point
  // their raw Ne and plasma panels at the source export by an absolute URL.
  const native = m.analysis_kind === 'native';
  const url = u => (u.startsWith('/') ? u : root + u);
  const dayStart = Date.parse(m.date + 'T00:00:00Z') / 1000;
  const nt = m.fit_times.length, nr = m.range_m.length;
  const params = await floats(url(m.parameters_url));
  const raw = await floats(url(m.raw_url));
  const gain = native ? m.native_calibration?.gain || null : m.calibration?.gain || null;
  const maps = [];
  for (let p = 0; p < 4; p++) maps.push(params.subarray(p * nt * nr, (p + 1) * nt * nr));
  if (gain && !native) {
    // density map from the proxy, as the main viewer does
    const ne = new Float32Array(nt * nr);
    m.fits.forEach((f, t) => f.parameters.forEach((row, r) => {
      const proxy = f.uncalibrated_density_proxy?.[r] ?? row[0];
      ne[r * nt + t] = calibratedDensity(m, proxy * gain, row[2], row[1]);
    }));
    maps[0] = ne;
  }
  // a native analysis shows its source's raw Ne, not its own: leave it out
  const rawPhysical = !native && m.raw_units === 'm^-3';
  const sorted = a => [...a].sort((x, y) => x - y), median = a => {const b = sorted(a); return b.length ? (b[(b.length - 1) >> 1] + b[b.length >> 1]) / 2 : NaN;};
  // The gates the analysis itself flags, [r * nt + t]: isr_analysis's
  // quality flags (FAILED, BOUND, NOISY; manifest quality_flags), or for a
  // native analysis any fit_status other than "converged". Null when the
  // analysis carries neither (Historical).
  let flagged = null;
  if (m.quality_flags && m.fits.every(f => f.quality_flags)) {
    const bits = m.quality_flags.bits, hide = m.quality_flags.hide.reduce((x, k) => x | bits[k], 0);
    flagged = new Uint8Array(nt * nr);
    m.fits.forEach((f, t) => f.quality_flags.forEach((q, r) => {flagged[r * nt + t] = (q & hide) ? 1 : 0;}));
  } else if (native && m.fits.every(f => f.fit_status)) {
    flagged = new Uint8Array(nt * nr);
    m.fits.forEach((f, t) => f.fit_status.forEach((q, r) => {flagged[r * nt + t] = q === 'converged' ? 0 : 1;}));
  }
  const windowS = median(m.fits.map(f => f.end_unix_s - f.start_unix_s));
  const gateKm = median(m.range_m.slice(1).map((r, i) => (r - m.range_m[i]) / 1000));
  return {
    m, label: choice.label, gain, absolute: !!gain, rawPhysical, windowS, gateKm, fitCount: m.fits.length, flagged, native,
    time: m.fit_times.map(t => t - dayStart), range: m.range_m.map(r => r / 1000),
    intervals: m.fits.map(f => [f.start_unix_s - dayStart, f.end_unix_s - dayStart]),
    maps, raw, rawTime: m.raw_times.map(t => t - dayStart), rawRange: m.raw_range_m.map(r => r / 1000),
  };
}

// Madrigal's analysed fits as imported for comparison: absolute Ne, Ti, Te
// and velocity on Madrigal's own times and slant-range cells, [p][r][t]; no
// raw Ne and no quality flags.
async function loadMadrigal(choice, root, m) {
  const dayStart = Date.parse(m.date + 'T00:00:00Z') / 1000;
  const nt = m.time_unix_s.length, nr = m.range_m.length;
  const params = await floats(root + m.parameters_url);
  const maps = [];
  for (let p = 0; p < 4; p++) maps.push(params.subarray(p * nt * nr, (p + 1) * nt * nr));
  const sorted = a => [...a].sort((x, y) => x - y), median = a => {const b = sorted(a); return b.length ? (b[(b.length - 1) >> 1] + b[b.length >> 1]) / 2 : NaN;};
  return {
    m: {...m, date: m.date}, label: choice.label, gain: null, absolute: true, rawPhysical: false,
    windowS: median(m.intervals.map(([a, b]) => b - a)),
    gateKm: median(m.range_m.slice(1).map((r, i) => (r - m.range_m[i]) / 1000)),
    fitCount: nt, flagged: null, native: false,
    time: m.time_unix_s.map(t => t - dayStart), range: m.range_m.map(r => r / 1000),
    intervals: m.intervals.map(([a, b]) => [a - dayStart, b - dayStart]),
    maps, raw: null, rawTime: [], rawRange: [],
  };
}

// values [r * nt + t] on range grid src, shown on grid dst: each dst cell takes
// the src gate whose cell contains it, NaN outside the src grid.
function regrid(values, nt, src, dst) {
  if (src.length === dst.length && src.every((v, i) => Math.abs(v - dst[i]) < 1e-6)) return values;
  const out = new Float32Array(nt * dst.length).fill(NaN);
  const half = i => (i + 1 < src.length ? src[i + 1] - src[i] : src[i] - src[i - 1]) / 2;
  dst.forEach((y, j) => {
    let best = -1, d = Infinity;
    src.forEach((s, i) => {const e = Math.abs(s - y); if (e < d) {d = e; best = i;}});
    if (best < 0 || d > half(best) + 1e-6) return;
    for (let t = 0; t < nt; t++) out[j * nt + t] = values[best * nt + t];
  });
  return out;
}

const finer = (a, b) => (a.length >= b.length ? a : b);
let split = 0.5;

function applySplit() {
  for (const row of document.querySelectorAll('.swipe-row.compared')) {
    const canvas = row.querySelector('.swipe-after canvas');
    if (!canvas) continue;
    const box = row.getBoundingClientRect(), c = canvas.getBoundingClientRect();
    const x = c.left - box.left + split * c.width;
    row.style.setProperty('--split', x + 'px');
    const handle = row.querySelector('.swipe-handle');
    handle.style.left = x + 'px';
    handle.style.top = (c.top - box.top) + 'px';
    handle.style.height = c.height + 'px';
  }
}

function bindHandle(handle) {
  handle.addEventListener('pointerdown', e => {
    e.preventDefault(); e.stopPropagation(); handle.setPointerCapture(e.pointerId);
    const row = handle.closest('.swipe-row');
    const move = ev => {
      const c = row.querySelector('.swipe-after canvas').getBoundingClientRect();
      split = Math.min(1, Math.max(0, (ev.clientX - c.left) / c.width));
      el('cmp-split').value = Math.round(split * 1000);
      applySplit();
    };
    const up = () => {handle.removeEventListener('pointermove', move); handle.removeEventListener('pointerup', up);};
    handle.addEventListener('pointermove', move); handle.addEventListener('pointerup', up);
  });
}

// Dashed vertical lines at the eclipse's contacts on every map, kept in place
// through pan and zoom by following each map's view on every frame.
let eclipseFrame = null;
function eclipseLines(panels, bounds, date) {
  if (eclipseFrame) cancelAnimationFrame(eclipseFrame);
  const events = ECLIPSE[date] || [];
  el('cmp-eclipse').textContent = events.length ? 'Dashed lines: the eclipse at Millstone Hill (ground): ' +
    events.map(([t, what]) => `${what} ${hhmm(t)}:${String(t % 60).padStart(2, '0')}`).join(', ') + ' UTC.' : '';
  if (!events.length) return;
  const sets = panels.map(p => {
    const box = p.article.querySelector('.ion-raster'), canvas = box.querySelector('canvas');
    const lines = events.map(() => {
      const d = document.createElement('div');
      Object.assign(d.style, {position: 'absolute', width: '0', borderLeft: '1.5px dashed #111', pointerEvents: 'none', zIndex: 3});
      box.append(d); return d;
    });
    return {p, canvas, lines};
  });
  const place = () => {
    for (const {p, canvas, lines} of sets) {
      const v = p.viewport.bounds(), x0 = bounds[0] + (bounds[1] - bounds[0]) * v[0], x1 = bounds[0] + (bounds[1] - bounds[0]) * v[2];
      events.forEach(([t], k) => {
        const f = (t - x0) / (x1 - x0), d = lines[k];
        d.style.display = f < 0 || f > 1 ? 'none' : '';
        Object.assign(d.style, {left: canvas.offsetLeft + f * canvas.offsetWidth + 'px', top: canvas.offsetTop + 'px', height: canvas.offsetHeight + 'px'});
      });
    }
    eclipseFrame = requestAnimationFrame(place);
  };
  place();
}

let hideFlagged = false;

// One row per shown parameter: rows[i] = {a, b, compare} (dataset ids), and
// data maps an id to its loaded analysis. Each render uses fresh element
// ids, so no map reuses a canvas of an earlier render.
let renders = 0;
async function draw(rows, data, choices, bounds, onRowChange) {
  const host = el('cmp-maps');
  host.replaceChildren();
  const gen = ++renders, panels = [];
  const res = d => `${d.label}: ${d.fitCount} fits, median window ${Math.round(d.windowS)} s, range cells ${d.gateKm.toFixed(1)} km`;
  const shown = new Set([...document.querySelectorAll('#cmp-params input:checked')].map(x => Number(x.value)));
  const used = [...new Set([...shown].flatMap(i => (rows[i].compare ? [rows[i].a, rows[i].b] : [rows[i].a])))];
  el('cmp-info-selected').textContent = 'Shown now. ' + used.map(id => res(data.get(id))).join('. ') + '.';
  const flagNote = d => (d.flagged ? (d.native ? 'its fit_status other than converged' : 'its quality flags FAILED, BOUND, NOISY') : 'no flags (shown in full)');
  el('cmp-flag-note').textContent = hideFlagged ? 'Flagged gates hidden, each dataset by its own test: ' +
    used.map(id => `${data.get(id).label}: ${flagNote(data.get(id))}`).join('; ') + '.' : '';
  const last = Math.max(...shown);
  const select = (value, title) => {
    const s = document.createElement('select'); s.title = title;
    for (const c of choices) {const o = document.createElement('option'); o.value = c.id; o.textContent = c.label; s.append(o);}
    s.value = value; return s;
  };
  for (const [i, spec] of ROWS.entries()) {
    if (!shown.has(i)) continue;
    const st = rows[i], A = data.get(st.a), B = st.compare ? data.get(st.b) : null;
    const raw = spec.p === 4;
    // the row's own controls: its dataset, the Compare toggle, the second dataset
    const bar = document.createElement('div');
    bar.className = 'cmp-rowbar';
    const name = document.createElement('span'); name.className = 'cmp-rowname'; name.textContent = spec.name.split(' /')[0];
    const first = select(st.a, B ? 'Left of the bar' : 'The dataset of this row');
    const toggle = document.createElement('button');
    toggle.className = 'cmp-toggle'; toggle.textContent = '⇔ Compare'; toggle.setAttribute('aria-pressed', String(st.compare));
    toggle.title = st.compare ? 'Show one dataset only' : 'Compare with a second dataset, split by a bar';
    bar.append(name, first, toggle);
    let second = null;
    if (st.compare) {second = select(st.b, 'Right of the bar'); bar.append(second);}
    first.onchange = () => onRowChange(i, {a: first.value});
    second && (second.onchange = () => onRowChange(i, {b: second.value}));
    toggle.onclick = () => onRowChange(i, {compare: !st.compare});
    host.append(bar);
    const sides = B ? [['a', A], ['b', B]] : [['a', A]];
    const missing = raw ? sides.filter(([, d]) => !d.rawPhysical).map(([, d]) => d.label) : [];
    if (missing.length) {
      const note = document.createElement('p'); note.className = 'caption cmp-rownote';
      note.textContent = `No raw Ne of its own for ${missing.join(' and ')}.`;
      host.append(note); continue;
    }
    const yA = raw ? A.rawRange : A.range, y = B ? finer(yA, raw ? B.rawRange : B.range) : yA;
    const row = document.createElement('div');
    row.className = 'swipe-row' + (B ? ' compared' : '');
    row.innerHTML = `<div class="swipe-layer swipe-before" id="cmp-${gen}-a-${i}"></div>` +
      (B ? `<div class="swipe-layer swipe-after" id="cmp-${gen}-b-${i}"></div><div class="swipe-handle" title="Drag to move the split"><span></span></div>` : '');
    host.append(row);
    const limits = spec.p === 0 && !sides.every(([, d]) => d.absolute) ? undefined : EISCAT_LIMITS[spec.p];
    for (const [side, d] of sides) {
      const x = raw ? d.rawTime : d.time, nt = x.length;
      let src = raw ? d.raw : d.maps[spec.p];
      if (!raw && hideFlagged && d.flagged) src = src.map((v, k) => (d.flagged[k] ? NaN : v));
      const values = regrid(src, nt, raw ? d.rawRange : d.range, y);
      panels.push(fieldPanel(`cmp-${gen}-${side}-${i}`, `cmp-${gen}-${side}-map-${i}`, spec.name, x, y, values, {
        timeBounds: bounds, intervals: raw ? undefined : d.intervals, mode: spec.mode,
        palette: spec.palette, colormap: spec.colormap, limits, sharedTime: true, sharedTimeLast: i === last,
        ylabel: 'Slant range / km'}));
    }
    if (B) bindHandle(row.querySelector('.swipe-handle'));
  }
  eclipseLines(panels, bounds, el('cmp-day').value);
  // The first map can be drawn before its fonts and labels are ready, and
  // stays blank; draw every map once more after layout.
  await document.fonts?.ready;
  await new Promise(r => requestAnimationFrame(r));
  for (const p of panels) p.redraw();
  applySplit();
}

async function main() {
  await init();
  const config = await json('api/catalog/millstone/datasets');
  const days = config.datasets.filter(d => d.antenna === antenna);
  const date = query.get('date') || days[0]?.date;
  const d = days.find(x => x.date === date);
  if (!d) throw Error(`No ${antenna} dataset for ${date}`);
  for (const x of days) {const o = document.createElement('option'); o.value = o.textContent = x.date; el('cmp-day').append(o);}
  el('cmp-day').value = date;
  // the day's analyses, then its registered comparisons (Madrigal), each id once
  const choices = [{id: 'historical', label: 'Historical', data_root: d.data_root}];
  for (const c of [...(d.analyses || []), ...(d.comparisons || [])]) if (!choices.some(x => x.id === c.id)) choices.push(c);
  for (const id of ['cmp-before', 'cmp-after']) for (const c of choices) {
    const o = document.createElement('option'); o.value = c.id; o.textContent = c.label; el(id).append(o);
  }
  const has = id => choices.some(c => c.id === id);
  for (const c of choices) {
    const dt = document.createElement('dt'), dd = document.createElement('dd');
    dt.textContent = c.label; dd.textContent = about(c.id) || '—';
    el('cmp-info-list').append(dt, dd);
  }
  el('cmp-info-btn').onclick = () => {
    const open = el('cmp-info').hidden; el('cmp-info').hidden = !open;
    el('cmp-info-btn').setAttribute('aria-expanded', String(open));
    requestAnimationFrame(applySplit);
  };
  if (query.get('info') === '1') el('cmp-info-btn').onclick();
  // the "All rows" defaults, then each row's own state from the URL:
  // ne=<id> shows one dataset, ne=<id>~<id> compares two (keys as KEYS)
  const defBefore = has(query.get('before')) ? query.get('before') : (choices.find(c => c.id.startsWith('isr-nosat')) || {id: 'historical'}).id;
  const defAfter = has(query.get('after')) ? query.get('after') : (choices.find(c => c.id.startsWith('isr-v2rfi')) || choices.find(c => c.id.startsWith('isr-tracked')) || choices.at(-1)).id;
  const defCompare = query.get('compare') !== 'off';
  el('cmp-before').value = defBefore; el('cmp-after').value = defAfter;
  const rows = ROWS.map((r, i) => {
    const [a, b] = (query.get(KEYS[i]) || '').split('~');
    if (has(a)) return {a, b: has(b) ? b : defAfter, compare: has(b)};
    return {a: defBefore, b: defAfter, compare: defCompare};
  });
  const data = new Map();
  const dayWindow = [Number(query.get('start')) || 11 * 3600, Number(query.get('end')) || 23 * 3600];
  const show = async () => {
    el('cmp-status').textContent = 'Loading…';
    try {
      const shownRows = [...document.querySelectorAll('#cmp-params input:checked')].map(x => rows[Number(x.value)]);
      const need = [...new Set(shownRows.flatMap(r => (r.compare ? [r.a, r.b] : [r.a])))].filter(id => !data.has(id));
      const got = await Promise.all(need.map(id => load(choices.find(c => c.id === id))));
      need.forEach((id, k) => data.set(id, got[k]));
      await draw(rows, data, choices, dayWindow, (i, change) => {Object.assign(rows[i], change); show();});
      el('cmp-status').textContent = '';
      const u = new URL(location.href);
      for (const k of ['before', 'after', 'compare', ...KEYS]) u.searchParams.delete(k);
      u.searchParams.set('date', date);
      rows.forEach((r, i) => u.searchParams.set(KEYS[i], r.compare ? r.a + '~' + r.b : r.a));
      if (hideFlagged) u.searchParams.set('flags', 'hide'); else u.searchParams.delete('flags');
      history.replaceState(null, '', u);
    } catch (e) {el('cmp-status').textContent = 'Unable to compare: ' + e.message;}
  };
  ROWS.forEach((r, i) => {
    const l = document.createElement('label'), c = document.createElement('input');
    c.type = 'checkbox'; c.value = i; c.checked = i < MAX_ROWS;
    c.onchange = () => {
      const on = document.querySelectorAll('#cmp-params input:checked');
      if (on.length > MAX_ROWS) {c.checked = false; el('cmp-status').textContent = `At most ${MAX_ROWS} parameters at once (the browser's graphics contexts run out).`; return;}
      show();
    };
    l.append(c, ' ' + r.name.split(' /')[0]); el('cmp-params').append(l);
  });
  // "All rows": set every row's first or second dataset at once
  el('cmp-before').onchange = () => {rows.forEach(r => {r.a = el('cmp-before').value;}); show();};
  el('cmp-after').onchange = () => {rows.forEach(r => {r.b = el('cmp-after').value;}); show();};
  el('cmp-swap').onclick = () => {
    const b = el('cmp-before').value; el('cmp-before').value = el('cmp-after').value; el('cmp-after').value = b;
    rows.forEach(r => {[r.a, r.b] = [r.b, r.a];}); show();
  };
  el('cmp-hide-flagged').checked = query.get('flags') === 'hide';
  hideFlagged = el('cmp-hide-flagged').checked;
  el('cmp-hide-flagged').onchange = () => {hideFlagged = el('cmp-hide-flagged').checked; show();};
  el('cmp-day').onchange = () => {const u = new URL(location.href); u.searchParams.set('date', el('cmp-day').value); location.assign(u);};
  el('cmp-split').oninput = () => {split = el('cmp-split').value / 1000; applySplit();};
  window.addEventListener('resize', () => requestAnimationFrame(applySplit));
  await show();
}
main().catch(e => {el('cmp-status').textContent = 'Unable to load: ' + e.message;});
