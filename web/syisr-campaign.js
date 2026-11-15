import init from "./pkg/guisdap_web.js";
import { fieldPanel, setViewportBounds, curve } from "./scientific-plots.js";
const $ = (id) => document.getElementById(id),
  root = "data/syisr/full-v2/",
  origin = 1700006400;
const parameters = [
  ["Electron density / m⁻³ (uncalibrated)", [1e10, 1e12], 1, "plasma"],
  ["Ion temperature / K", [200, 2500], 0, "eiscat"],
  ["Electron temperature / K", [300, 4000], 0, "eiscat"],
  ["Spectral velocity / m s⁻¹", [-300, 300], 0, "seismic"],
  ["H⁺ fraction", [0, 1], 0, "plasma"],
  ["He⁺ fraction", [0, 1], 0, "plasma"],
  ["Signal power / ADU²", [1e2, 1e6], 1, "plasma"],
];
let manifest,
  beam,
  data,
  fields = [],
  generation = 0;
async function json(url) {
  const response = await fetch(url, { cache: "no-store" });
  if (!response.ok) throw Error(`${response.status}: ${url}`);
  return response.json();
}
function synchronize(source) {
  const bounds = Array.from(source.viewport.bounds());
  fields.forEach((f) => {
    if (f !== source) {
      setViewportBounds(f.viewport, bounds);
      f.redraw();
    }
  });
}
function draw() {
  const times = beam.time_unix_s.map((t) => t - origin),
    nt = times.length,
    nr = beam.altitude_km.length,
    cells = nt * nr;
  fields = [];
  for (let k = 0; k < parameters.length; k++) {
    const [title, limits, mode, colormap] = parameters[k],
      values = data.slice(k * cells, (k + 1) * cells);
    if ($("reject").checked && k !== 6)
      for (let j = 0; j < cells; j++)
        if (data[7 * cells + j] !== 0) values[j] = NaN;
    let f;
    f = fieldPanel(
      "maps",
      `campaign-${k}`,
      title,
      times,
      beam.altitude_km,
      values,
      {
        limits,
        mode,
        colormap,
        sharedTime: true,
        sharedTimeLast: k === parameters.length - 1,
        timeGap: 0,
        intervals: times.map((t) => [t - 900, t + 900]),
        ylabel: "Altitude / km",
        colorbarTitle: title,
        pick,
        onView: () => synchronize(f),
      },
    );
    fields.push(f);
  }
}
async function selectBeam() {
  const token = ++generation;
  beam = manifest.beams[+$("beam").value];
  const r = await fetch(root + beam.file);
  if (!r.ok) throw Error("Beam tile unavailable");
  const bytes = await r.arrayBuffer();
  if (token !== generation) return;
  const [np, nr, nt] = beam.shape;
  if (bytes.byteLength !== np * nr * nt * 4)
    throw Error("Invalid tile dimensions");
  data = new Float32Array(bytes);
  draw();
  $("status").textContent =
    `${nt} half-hour profiles · azimuth ${beam.azimuth_deg.toFixed(2)}° · elevation ${beam.elevation_deg.toFixed(2)}° · all 41 measured lags · absolute Ne unavailable`;
  $("selection").textContent =
    "Select a pixel to inspect measured lag products and the forward prediction.";
}
async function pick(time, height) {
  const b = beam,
    token = generation,
    nt = b.time_unix_s.length,
    nr = b.altitude_km.length,
    cells = nt * nr;
  let ti = 0,
    ri = 0;
  b.time_unix_s.forEach((t, i) => {
    if (
      Math.abs(t - origin - time) < Math.abs(b.time_unix_s[ti] - origin - time)
    )
      ti = i;
  });
  b.altitude_km.forEach((h, i) => {
    if (Math.abs(h - height) < Math.abs(b.altitude_km[ri] - height)) ri = i;
  });
  const at = (p) => data[p * cells + ri * nt + ti],
    fmt = (v) => (Number.isFinite(v) ? v.toFixed(2) : "unavailable");
  fields.forEach((f) =>
    f.setCursor(b.time_unix_s[ti] - origin, b.altitude_km[ri]),
  );
  $("selection").textContent =
    `${new Date((b.time_unix_s[ti] - 900) * 1000).toISOString().slice(0, 16)}–${new Date((b.time_unix_s[ti] + 900) * 1000).toISOString().slice(11, 16)} UTC · ${b.altitude_km[ri].toFixed(1)} km altitude · Ti ${fmt(at(1))} K · Te ${fmt(at(2))} K · spectral velocity ${fmt(at(3))} m/s · H⁺ ${fmt(at(4))} · He⁺ ${fmt(at(5))} · Ne unavailable · ${at(7) === 0 ? "accepted conditional diagnostic" : at(7) === 1 ? "flagged / prescribed" : "unavailable"}`;
  if (at(7) === 2) {
    $("quality").textContent =
      "No supported fitted parameter at this altitude/time.";
    return;
  }
  const entry = b.diagnostics[ti];
  if (!entry) {
    $("quality").textContent = "No fitted products for this cell.";
    return;
  }
  try {
    const d = await json(root + entry.file);
    if (token !== generation) return;
    const range = b.slant_range_km[ri],
      sample = ((2 * range) / 299792.458 + 240e-6 - 767e-6) / 1e-5;
    const unique = [...new Set(d.products.map((p) => p.sample))],
      gate = unique.reduce(
        (a, k) => (Math.abs(k - sample) < Math.abs(a - sample) ? k : a),
        unique[0],
      );
    const selected = d.products
        .map((p, i) => ({ ...p, predicted: d.predicted[i] }))
        .filter((p) => p.sample === gate),
      all = selected
        .flatMap((p) => [
          ...p.value.map((v, k) => v - p.sigma[k]),
          ...p.value.map((v, k) => v + p.sigma[k]),
          ...p.predicted,
        ])
        .filter((v) => Number.isFinite(v) && Math.abs(v) < 1e20);
    const lo = Math.min(0, ...all),
      hi = Math.max(0, ...all),
      pad = Math.max(1, (hi - lo) * 0.07);
    $("quality").textContent =
      `${d.status} · reduced χ² ${d.reduced_chi_square.toFixed(2)} · receive sample ${gate} · common vertical scale · individual lag-dependent range ambiguities. Error bars: diagonal measurement SEM; response uncertainty is additional.`;
    for (let k = 0; k < 2; k++) {
      const xy = [],
        model = [],
        bars = [];
      for (const p of selected) {
        if (k === 1 && p.lag === 0) continue;
        const x = p.lag * 10;
        xy.push(x, p.value[k]);
        model.push(x, p.predicted[k]);
        bars.push(x, p.value[k] - p.sigma[k], x, p.value[k] + p.sigma[k]);
      }
      curve(k ? "imag" : "real", {
        title: k ? "Imaginary lag products" : "Real lag products",
        bounds: [0, 400, lo - pad, hi + pad],
        xlabel: "Lag / µs",
        ylabel: "Background-subtracted product / ADU²",
        series: [
          { xy: bars, segments: true, color: "#999999" },
          { xy, radius: 2, color: "#444444" },
          { xy: model, line: true, color: "#236b8e" },
        ],
      });
    }
  } catch (e) {
    $("quality").textContent = String(e);
  }
}
try {
  await init();
  manifest = await json(root + "manifest.json");
  manifest.beams.forEach((b, i) =>
    $("beam").add(
      new Option(
        `${i}: az ${b.azimuth_deg.toFixed(1)}°, el ${b.elevation_deg.toFixed(1)}°`,
        i,
      ),
    ),
  );
  $("beam").value = manifest.beams.findIndex((b) => b.elevation_deg > 89);
  $("beam").onchange = () =>
    selectBeam().catch((e) => ($("status").textContent = String(e)));
  $("reject").onchange = draw;
  await selectBeam();
} catch (e) {
  $("status").textContent = String(e);
}
