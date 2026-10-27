// Before/after comparison of two Millstone analyses of one day.
// Each parameter map is drawn twice, from the "before" and the "after"
// analysis, on the same time and range axes and the same colour limits. The
// "after" map lies on top and is clipped left of a vertical bar, so the left
// part shows "before" and the right part "after". Dragging the bar moves the
// split on every map; pan and zoom stay shared through the dataset time group.
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

// One to three sentences per dataset, by analysis id (the info panel).
const ABOUT = [
  [/^historical$/, "Products of the older Python pipeline (30 µs gates), imported by Juha and calibrated with the plasma lines. No satellite removal."],
  [/^native-(full-)?20261019$/, "An early test of Juha's Rust pipeline on a sample of pulses only (16 lags, 200–540 km). Incomplete; not for comparisons."],
  [/^native-direct-live/, "Juha's Rust pipeline, reading the raw recordings directly and subtracting every CFAR-detected satellite echo from the voltages. Fitted window by window while it ran; its subtraction is still under Juha's review."],
  [/^native-debris400-.*live/, "The same as the strict Rust run (echoes at least 3 dB above the CFAR threshold, above 400 km), published window by window while it ran; once complete it matches that run."],
  [/^native-debris400/, "Juha's Rust pipeline with a stricter subtraction: only echoes at least 3 dB above the CFAR threshold and above 400 km altitude are subtracted from the voltages."],
  [/^isr-nosat/, "Our variance-weighted lag profile inversion with a background column; satellite echoes are left in (variant 3, Memo 42)."],
  [/^isr-detected/, "Our inversion with each detected satellite echo modelled as an extra unknown, on the pulses where the detector found it (variant 5, Memo 42)."],
  [/^isr-tracked/, "As variant 5, with pass tracking: detections are linked into passes and every pulse of a pass gets a column, also those the detector missed (Memo 45)."],
  [/^isr-v2rfi/, "Our final product: pass tracking with a lower detection threshold (false-alarm probability 1e-6, faint detections kept only when tracking links them into a pass; Memo 48), and the lagged products hit by power-line impulses removed from the inversion (Memo 47)."],
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
  const windowS = median(m.fits.map(f => f.end_unix_s - f.start_unix_s));
  const gateKm = median(m.range_m.slice(1).map((r, i) => (r - m.range_m[i]) / 1000));
  return {
    m, label: choice.label, gain, rawPhysical, windowS, gateKm, fitCount: m.fits.length,
    time: m.fit_times.map(t => t - dayStart), range: m.range_m.map(r => r / 1000),
    intervals: m.fits.map(f => [f.start_unix_s - dayStart, f.end_unix_s - dayStart]),
    maps, raw, rawTime: m.raw_times.map(t => t - dayStart), rawRange: m.raw_range_m.map(r => r / 1000),
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
  for (const row of document.querySelectorAll('.swipe-row')) {
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

async function draw(before, after, bounds) {
  const host = el('cmp-maps');
  host.replaceChildren();
  const panels = [];
  const res = d => `${d.label}: ${d.fitCount} fits, median window ${Math.round(d.windowS)} s, range cells ${d.gateKm.toFixed(1)} km`;
  el('cmp-info-selected').textContent = `Selected now. Before: ${res(before)}. After: ${res(after)}.`;
  el('cmp-left').textContent = 'Before (left): ' + before.label;
  el('cmp-right').textContent = 'After (right): ' + after.label;
  const shown = new Set([...document.querySelectorAll('#cmp-params input:checked')].map(x => Number(x.value)));
  const last = Math.max(...shown);
  for (const [i, spec] of ROWS.entries()) {
    if (!shown.has(i)) continue;
    const raw = spec.p === 4;
    if (raw && !(before.rawPhysical && after.rawPhysical)) continue;
    const yA = raw ? before.rawRange : before.range, yB = raw ? after.rawRange : after.range;
    const y = finer(yA, yB);
    const row = document.createElement('div');
    row.className = 'swipe-row';
    row.innerHTML = `<div class="swipe-layer swipe-before" id="cmp-a-${i}"></div><div class="swipe-layer swipe-after" id="cmp-b-${i}"></div><div class="swipe-handle" title="Drag to move the split"><span></span></div>`;
    host.append(row);
    const limits = spec.p === 0 && !(before.gain && after.gain) ? undefined : EISCAT_LIMITS[spec.p];
    for (const [side, d] of [['a', before], ['b', after]]) {
      const x = raw ? d.rawTime : d.time, nt = x.length;
      const values = regrid(raw ? d.raw : d.maps[spec.p], nt, raw ? d.rawRange : d.range, y);
      panels.push(fieldPanel(`cmp-${side}-${i}`, `cmp-${side}-map-${i}`, spec.name, x, y, values, {
        timeBounds: bounds, intervals: raw ? undefined : d.intervals, mode: spec.mode,
        palette: spec.palette, colormap: spec.colormap, limits, sharedTime: true, sharedTimeLast: i === last,
        ylabel: 'Slant range / km'}));
    }
    bindHandle(row.querySelector('.swipe-handle'));
  }
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
  const choices = [{id: 'historical', label: 'Historical', data_root: d.data_root}, ...(d.analyses || [])];
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
  el('cmp-before').value = has(query.get('before')) ? query.get('before') : has('isr-nosat-v1') ? 'isr-nosat-v1' : 'historical';
  el('cmp-after').value = has(query.get('after')) ? query.get('after') : has('isr-tracked-v1') ? 'isr-tracked-v1' : choices.at(-1).id;
  const dayWindow = [Number(query.get('start')) || 11 * 3600, Number(query.get('end')) || 23 * 3600];
  const show = async () => {
    el('cmp-status').textContent = 'Loading…';
    const pick = id => choices.find(c => c.id === el(id).value);
    try {
      const [b, a] = await Promise.all([load(pick('cmp-before')), load(pick('cmp-after'))]);
      await draw(b, a, dayWindow);
      el('cmp-status').textContent = '';
      const u = new URL(location.href);
      u.searchParams.set('date', date); u.searchParams.set('before', el('cmp-before').value); u.searchParams.set('after', el('cmp-after').value);
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
  el('cmp-before').onchange = el('cmp-after').onchange = show;
  el('cmp-day').onchange = () => {const u = new URL(location.href); u.searchParams.set('date', el('cmp-day').value); location.assign(u);};
  el('cmp-swap').onclick = () => {const b = el('cmp-before').value; el('cmp-before').value = el('cmp-after').value; el('cmp-after').value = b; show();};
  el('cmp-split').oninput = () => {split = el('cmp-split').value / 1000; applySplit();};
  window.addEventListener('resize', () => requestAnimationFrame(applySplit));
  await show();
}
main().catch(e => {el('cmp-status').textContent = 'Unable to load: ' + e.message;});
