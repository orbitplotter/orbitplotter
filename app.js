// =============================================================
// app.js
// UI wiring + simulation drivers + plotting for the Orbit Plotter.
// Uses the exact/adaptive propagators defined in physics.js.
// =============================================================

// ---- Plot color palette (mirrors styles.css custom properties) ----
// Kept to the site's actual two-hue identity (blue accent + purple pastel)
// plus neutrals, rather than introducing off-theme reds/greens for markers.
const COLORS = {
  bg: "#191D27",
  border: "#262B37",
  text: "#E9ECF3",
  textMuted: "#8891A3",
  accent: "#5C8DF6",
  accentSoft: "#7FA6FF",
  pastel: "#A99CE0",
};

const PLOTLY_FONT = { family: "Inter, sans-serif", color: COLORS.textMuted, size: 11 };

// ---------- Input helpers ----------
function getVal(name, def = 0) {
  const el = document.querySelector(`[name="${name}"]`);
  const v = el ? parseFloat(el.value) : NaN;
  return Number.isFinite(v) ? v : def;
}

// ---------- Drawing: 2-Body ----------
// markerPoint (optional [x,y,z]): when the elements-mode input specifies
// where on the orbit the spacecraft actually is (a true or mean anomaly),
// this draws a distinct marker there in addition to the orbit itself.
function drawOrbit2Body(points, markerPoint = null) {
  const xs = points.map((p) => p[0]);
  const ys = points.map((p) => p[1]);
  const zs = points.map((p) => p[2]);

  const orbitTrace = {
    x: xs, y: ys, z: zs,
    mode: "lines", type: "scatter3d",
    line: { color: COLORS.accent, width: 4 },
    name: "Orbit",
  };
  const primary = {
    type: "scatter3d", x: [0], y: [0], z: [0],
    mode: "markers", marker: { color: COLORS.text, size: 8 },
    name: "Primary",
  };

  const traces = [orbitTrace, primary];
  if (markerPoint) {
    traces.push({
      type: "scatter3d", x: [markerPoint[0]], y: [markerPoint[1]], z: [markerPoint[2]],
      mode: "markers", marker: { color: COLORS.pastel, size: 7, symbol: "diamond" },
      name: "Spacecraft",
    });
  }

  const maxR = Math.max(
    Math.max(...xs.map((v) => Math.abs(v))),
    Math.max(...ys.map((v) => Math.abs(v))),
    Math.max(...zs.map((v) => Math.abs(v))),
    markerPoint ? Math.max(...markerPoint.map((v) => Math.abs(v))) : 0
  ) || 1;
  const axisRange = [-1.2 * maxR, 1.2 * maxR];

  const layout = {
    scene: {
      xaxis: { title: "x (km)", range: axisRange, color: COLORS.textMuted, gridcolor: COLORS.border, zerolinecolor: COLORS.border },
      yaxis: { title: "y (km)", range: axisRange, color: COLORS.textMuted, gridcolor: COLORS.border, zerolinecolor: COLORS.border },
      zaxis: { title: "z (km)", range: axisRange, color: COLORS.textMuted, gridcolor: COLORS.border, zerolinecolor: COLORS.border },
      aspectmode: "cube",
      bgcolor: COLORS.bg,
    },
    showlegend: true,
    legend: { x: 0.02, y: 0.98, bgcolor: "rgba(0,0,0,0.3)", borderwidth: 0, font: PLOTLY_FONT },
    margin: { t: 10, l: 0, r: 0, b: 0 },
    paper_bgcolor: COLORS.bg,
    plot_bgcolor: COLORS.bg,
    font: PLOTLY_FONT,
  };

  Plotly.newPlot("mainPlot", traces, layout, { displayModeBar: false });
}

// ---------- Drawing: CR3BP ----------
function drawOrbitCR3B(points, mu) {
  const xs = points.map((p) => p[0]);
  const ys = points.map((p) => p[1]);
  const zs = points.map((p) => p[2]);

  const all = [...xs, ...ys, ...zs];
  const minVal = Math.min(...all), maxVal = Math.max(...all);
  const pad = 0.1 * (maxVal - minVal || 1);
  const range = [minVal - pad, maxVal + pad];

  const orbitTrace = {
    x: xs, y: ys, z: zs,
    mode: "lines", type: "scatter3d",
    line: { width: 3, color: COLORS.accent },
    name: "Trajectory",
  };
  const p1 = {
    type: "scatter3d", x: [-mu], y: [0], z: [0],
    mode: "markers", marker: { size: 7, color: COLORS.text },
    name: "Primary 1",
  };
  const p2 = {
    type: "scatter3d", x: [1 - mu], y: [0], z: [0],
    mode: "markers", marker: { size: 5, color: COLORS.pastel },
    name: "Primary 2",
  };

  const layout = {
    scene: {
      xaxis: { title: "x", range, color: COLORS.textMuted, gridcolor: COLORS.border, zerolinecolor: COLORS.border },
      yaxis: { title: "y", range, color: COLORS.textMuted, gridcolor: COLORS.border, zerolinecolor: COLORS.border },
      zaxis: { title: "z", range, color: COLORS.textMuted, gridcolor: COLORS.border, zerolinecolor: COLORS.border },
      aspectmode: "cube",
      bgcolor: COLORS.bg,
    },
    showlegend: true,
    legend: { x: 0.02, y: 0.98, bgcolor: "rgba(0,0,0,0.3)", borderwidth: 0, font: PLOTLY_FONT },
    margin: { t: 10, l: 0, r: 0, b: 0 },
    paper_bgcolor: COLORS.bg,
    plot_bgcolor: COLORS.bg,
    font: PLOTLY_FONT,
  };

  Plotly.newPlot("mainPlot", [orbitTrace, p1, p2], layout, { displayModeBar: false });
}

// ---------- Drawing: Relative Motion ----------
function drawOrbitRelative(r_lvlh) {
  const xs = r_lvlh.map((p) => p[0]); // radial
  const ys = r_lvlh.map((p) => p[1]); // in-track
  const zs = r_lvlh.map((p) => p[2]); // cross-track

  const trace = {
    x: ys, y: xs, z: zs, // LVLH swap for a conventional in-track/radial/cross-track view
    mode: "lines", type: "scatter3d",
    line: { color: COLORS.accent, width: 4 },
    name: "Deputy in LVLH",
  };
  const chief = {
    x: [0], y: [0], z: [0],
    mode: "markers", type: "scatter3d",
    marker: { color: COLORS.text, size: 8 },
    name: "Chief",
  };

  const layout = {
    scene: {
      xaxis: { title: "In-track y (km)", color: COLORS.textMuted, gridcolor: COLORS.border, zerolinecolor: COLORS.border },
      yaxis: { title: "Radial x (km)", color: COLORS.textMuted, gridcolor: COLORS.border, zerolinecolor: COLORS.border },
      zaxis: { title: "Cross-track z (km)", color: COLORS.textMuted, gridcolor: COLORS.border, zerolinecolor: COLORS.border },
      aspectmode: "cube",
      bgcolor: COLORS.bg,
    },
    showlegend: true,
    legend: { x: 0.02, y: 0.98, bgcolor: "rgba(0,0,0,0.3)", borderwidth: 0, font: PLOTLY_FONT },
    margin: { t: 10, l: 0, r: 0, b: 0 },
    paper_bgcolor: COLORS.bg,
    plot_bgcolor: COLORS.bg,
    font: PLOTLY_FONT,
  };

  Plotly.newPlot("mainPlot", [trace, chief], layout, { displayModeBar: false });
}

// ---------- Drawing: conserved-quantity error ----------
function plotConserved(seriesList) {
  if (!seriesList || seriesList.length === 0 || seriesList[0].values.length === 0) {
    Plotly.newPlot("conservedPlot", []);
    return;
  }
  const palette = [COLORS.accent, COLORS.pastel, COLORS.accentSoft];
  const traces = seriesList.map((s, i) => ({
    x: s.values.map((_, idx) => idx),
    y: s.values,
    mode: "lines",
    type: "scatter",
    line: { width: 2, color: palette[i % palette.length] },
    name: s.label,
  }));

  const layout = {
    showlegend: seriesList.length > 1,
    // Explicit height (matching the CSS box exactly) rather than relying
    // solely on Plotly reading the container's size at call-time - that
    // measurement can race with layout/paint timing and previously let
    // the chart render taller than its box, bleeding into the panel below.
    height: 200,
    autosize: false,
    // Vertical stacked legend below the axis (not "h"/horizontal) - the
    // plot itself is narrow (right sidebar column), so a horizontal legend
    // with two fairly long labels could wrap unpredictably; stacking is
    // predictable regardless of width.
    legend: {
      x: 0.5, xanchor: "center",
      y: -0.3, yanchor: "top",
      bgcolor: "rgba(0,0,0,0)", borderwidth: 0,
      font: { ...PLOTLY_FONT, size: 9 },
    },
    paper_bgcolor: COLORS.bg,
    plot_bgcolor: COLORS.bg,
    font: PLOTLY_FONT,
    margin: { t: 12, l: 46, r: 10, b: seriesList.length > 1 ? 52 : 30 },
    xaxis: { title: "sample #", color: COLORS.textMuted, gridcolor: COLORS.border, zerolinecolor: COLORS.border },
    yaxis: { title: "error (%)", color: COLORS.textMuted, gridcolor: COLORS.border, zerolinecolor: COLORS.border, exponentformat: "e" },
  };

  Plotly.newPlot("conservedPlot", traces, layout, { displayModeBar: false });
}

// ---------- Simulators ----------
// ---------- Conserved-quantity note (one sentence, tailored to what was just run) ----------
const CONSERVED_NOTES = {
  twoBody: "Exact closed-form Kepler propagation, so this should stay essentially flat at machine precision (~1e-10%).",
  shapeOnly: "No propagation was performed for a shape-only orbit, so there's nothing to conserve here.",
  cr3bp: "Adaptive numerical integration (no closed-form solution exists), so this stays small but non-zero - typically well under 1e-6%.",
  relativeMotion: "Both orbits use the same exact closed-form Kepler propagation, so this should stay essentially flat at machine precision (~1e-10%).",
};

function setConservedNote(key) {
  const el = document.querySelector(".conserved-note");
  if (el) el.textContent = CONSERVED_NOTES[key] || "";
}

function simulate2Body() {
  const mode = document.querySelector('[name="twoBodyInputMode"]')?.value || "cartesian";
  if (mode === "elements") return simulate2BodyFromElements();
  return simulate2BodyFromCartesian();
}

function simulate2BodyFromCartesian() {
  const MU = getVal("mu", EARTH_MU);
  const tfInput = getVal("tf", 86400);
  const r0 = [getVal("xpos", 7000), getVal("ypos", 0), getVal("zpos", 0)];
  const v0 = [getVal("xvel", 0), getVal("yvel", 7.5), getVal("zvel", 0)];
  const state0 = [...r0, ...v0];

  const { T } = orbitalWindow(state0, MU, tfInput);
  const steps = 1000;
  const dt = T / steps;

  let state = state0.slice();
  const trajectory = [r0];

  const h0 = Math.hypot(...angularMomentum(r0, v0));
  const e0 = specificEnergy(r0, v0, MU);
  const hError = [];
  const energyError = [];

  for (let i = 0; i < steps; i++) {
    state = keplerStep(state, dt, MU);
    trajectory.push(state.slice(0, 3));
    const r = state.slice(0, 3), v = state.slice(3, 6);
    const h = Math.hypot(...angularMomentum(r, v));
    const e = specificEnergy(r, v, MU);
    hError.push(((h - h0) / h0) * 100);
    energyError.push(((e - e0) / e0) * 100);
  }

  drawOrbit2Body(trajectory);
  plotConserved([
    { label: "Angular momentum", values: hError },
    { label: "Specific energy", values: energyError },
  ]);
  setConservedNote("twoBody");
}

// Orbital-elements input path. If no anomaly is given, the elements alone
// only pin down the *shape* of the orbit (its size, shape and orientation),
// not where the spacecraft is along it - so we draw the exact closed conic
// with no time dependence at all (sampled directly from the geometry, not
// propagated). If an anomaly (true or mean) is given, that fixes an actual
// position, so we additionally propagate from there (same exact Kepler
// propagator as the Cartesian path) to get a conserved-quantity diagnostic,
// and mark the spacecraft's position on the shape.
function simulate2BodyFromElements() {
  const MU = getVal("mu", EARTH_MU);
  const a = getVal("oe_a", 7000);
  const e = getVal("oe_e", 0);
  const i = getVal("oe_i", 0);
  const RAAN = getVal("oe_RAAN", 0);
  const argp = getVal("oe_argp", 0);

  const anomalyEl = document.querySelector('input[name="anomalyType"]:checked');
  const anomalyType = anomalyEl ? anomalyEl.value : "none";
  const anomalyValue = getVal("anomalyValue", 0);

  const shape = orbitShapePoints(a, e, i, RAAN, argp, MU, 361);

  if (anomalyType === "none") {
    drawOrbit2Body(shape, null);
    plotConserved([]);
    setConservedNote("shapeOnly");
    return;
  }

  const [r0, v0] =
    anomalyType === "true"
      ? elementsToStateVector(a, e, i, RAAN, argp, anomalyValue, MU)
      : oeToStateVector(a, e, i, RAAN, argp, anomalyValue, MU);
  const state0 = [...r0, ...v0];

  const { T } = orbitalWindow(state0, MU, 86400);
  const steps = 1000;
  const dt = T / steps;
  let state = state0.slice();

  const h0 = Math.hypot(...angularMomentum(r0, v0));
  const e0 = specificEnergy(r0, v0, MU);
  const hError = [];
  const energyError = [];

  for (let k = 0; k < steps; k++) {
    state = keplerStep(state, dt, MU);
    const r = state.slice(0, 3), v = state.slice(3, 6);
    const h = Math.hypot(...angularMomentum(r, v));
    const en = specificEnergy(r, v, MU);
    hError.push(((h - h0) / h0) * 100);
    energyError.push(((en - e0) / e0) * 100);
  }

  drawOrbit2Body(shape, r0);
  plotConserved([
    { label: "Angular momentum", values: hError },
    { label: "Specific energy", values: energyError },
  ]);
  setConservedNote("twoBody");
}

function simulateCR3B() {
  const MU3B = getVal("mu3b", 0.0121505856);
  const r0 = [getVal("xpos", 1.0220215127), getVal("ypos", 0), getVal("zpos", -0.1820967615)];
  const v0 = [getVal("xvel", 0), getVal("yvel", -0.1032563411), getVal("zvel", 0)];
  const state0 = [...r0, ...v0];
  const tf = getVal("tf", 1.51111);

  const traj = integrateAdaptive(CR3B_ode, state0, 0, tf, [MU3B]);
  const trajectory = traj.map((s) => s.slice(0, 3));

  const C0 = jacobiConstant(state0, MU3B);
  const jacobiError = traj.map((s) => ((jacobiConstant(s, MU3B) - C0) / C0) * 100);

  drawOrbitCR3B(trajectory, MU3B);
  plotConserved([{ label: "Jacobi constant error", values: jacobiError }]);
  setConservedNote("cr3bp");
}

function simulateRelativeMotion() {
  const chiefEEl = document.querySelector('input[name="chief_e"]:checked');
  const chief_e = chiefEEl ? parseFloat(chiefEEl.value) : 0;

  const da = getVal("da", 0);
  const de = getVal("de", 0);
  const di = getVal("di", 0);
  const dRAAN = getVal("dRANN", 0);
  const dargp = getVal("dargp", 0);
  const dM = getVal("dM", 0);
  const periods = Math.max(1, Math.round(getVal("orbits", 1)));

  const mu = EARTH_MU;
  const eccentricityMap = { 0.0: 8000, 0.2: 10000, 0.4: 13333.3333333333, 0.6: 20000, 0.8: 40000 };
  const a_chief = eccentricityMap[chief_e] ?? 8000;

  const oe_c = { a: a_chief, e: chief_e, i: 0, RAAN: 0, w: 0, M: 0 };
  const oe_d = {
    a: a_chief + da, e: chief_e + de, i: 0 + di,
    RAAN: 0 + dRAAN, w: 0 + dargp, M: 0 + dM,
  };

  const [r_c0, v_c0] = oeToStateVector(oe_c.a, oe_c.e, oe_c.i, oe_c.RAAN, oe_c.w, oe_c.M, mu);
  const [r_d0, v_d0] = oeToStateVector(oe_d.a, oe_d.e, oe_d.i, oe_d.RAAN, oe_d.w, oe_d.M, mu);

  let xc = [...r_c0, ...v_c0];
  let xd = [...r_d0, ...v_d0];

  const { T } = orbitalWindow(xc, mu, 86400);
  const totalTime = T * periods;
  const samplesPerOrbit = 500;
  const N = samplesPerOrbit * periods;
  const dt = totalTime / N;

  const r_lvlh = [];
  const h0 = Math.hypot(...angularMomentum(r_c0, v_c0));
  const hError = [];

  for (let k = 0; k < N; k++) {
    const rc = xc.slice(0, 3), vc = xc.slice(3, 6);
    const rd = xd.slice(0, 3), vd = xd.slice(3, 6);
    const lv = eciToLvlh(rc, vc, rd, vd);
    r_lvlh.push(lv.r_lvlh);

    const h = Math.hypot(...angularMomentum(rc, vc));
    hError.push(((h - h0) / h0) * 100);

    xc = keplerStep(xc, dt, mu);
    xd = keplerStep(xd, dt, mu);
  }

  drawOrbitRelative(r_lvlh);
  plotConserved([{ label: "Chief angular momentum error", values: hError }]);
  setConservedNote("relativeMotion");
}

function runSimulation() {
  const orbitType = document.getElementById("orbitType")?.value;
  if (orbitType === "2body") simulate2Body();
  else if (orbitType === "cr3bp") simulateCR3B();
  else if (orbitType === "relativemotion") simulateRelativeMotion();
}

// ---------- Example library ----------
const EXAMPLES = {
  "2body": [
    {
      title: "Circular LEO",
      tag: "Elliptical • e ≈ 0",
      desc: "A near-circular 7000 km reference orbit.",
      values: { mu: 398600.4418, xpos: 7000, ypos: 0, zpos: 0, xvel: 0, yvel: 7.5, zvel: 0, tf: 86400 },
    },
    {
      title: "Eccentric Ellipse",
      tag: "Elliptical",
      desc: "A highly inclined, moderately eccentric Earth orbit.",
      values: { mu: 398600.4418, xpos: -1774, ypos: -2787, zpos: 6801, xvel: 6.2, yvel: -1.9, zvel: -0.1, tf: 86400 },
    },
    {
      title: "Hyperbolic Flyby",
      tag: "Hyperbolic",
      desc: "An unbound flyby trajectory with positive specific energy.",
      values: { mu: 398600.4418, xpos: -5700, ypos: -2541, zpos: -14869, xvel: -4.9, yvel: -2.5, zvel: -4.4, tf: 20000 },
    },
    {
      title: "Molniya Shape (Elements)",
      tag: "Orbital elements • shape only",
      desc: "A highly eccentric Molniya-type ellipse defined purely by its elements - no anomaly given, so only the shape is drawn.",
      mode: "elements",
      anomalyType: "none",
      values: { mu: 398600.4418, oe_a: 26600, oe_e: 0.72, oe_i: 63.4, oe_RAAN: 0, oe_argp: 270 },
    },
    {
      title: "Spacecraft at Perigee (Elements)",
      tag: "Orbital elements • true anomaly",
      desc: "Same style of orbit, but with true anomaly ν = 0° - the spacecraft marker sits right at perigee.",
      mode: "elements",
      anomalyType: "true",
      values: { mu: 398600.4418, oe_a: 7000, oe_e: 0.1, oe_i: 28.5, oe_RAAN: 45, oe_argp: 90, anomalyValue: 0 },
    },
    {
      title: "Spacecraft at Given Mean Anomaly (Elements)",
      tag: "Orbital elements • mean anomaly",
      desc: "Mean anomaly M = 120° places the spacecraft marker well past perigee, elsewhere on the same orbit.",
      mode: "elements",
      anomalyType: "mean",
      values: { mu: 398600.4418, oe_a: 7000, oe_e: 0.1, oe_i: 28.5, oe_RAAN: 45, oe_argp: 90, anomalyValue: 120 },
    },
  ],
  cr3bp: [
    {
      title: "Halo Orbit (L2)",
      tag: "Earth–Moon CR3BP",
      desc: "A periodic three-dimensional halo orbit around the L2 point.",
      values: { mu3b: 0.0121505856, xpos: 1.0220215127, ypos: 0, zpos: -0.1820967615, xvel: 0, yvel: -0.1032563411, zvel: 0, tf: 1.51111 },
    },
    {
      title: "Lyapunov Orbit (L1)",
      tag: "Earth–Moon CR3BP",
      desc: "A planar periodic orbit around the L1 point.",
      values: { mu3b: 0.0121505856, xpos: 0.65457084231188, ypos: 0, zpos: 0, xvel: 0, yvel: 0.7413347560791179, zvel: 0, tf: 6.45 },
    },
    {
      title: "2:1 Resonant Orbit",
      tag: "Earth–Moon CR3BP",
      desc: "A planar orbit resonant with the primaries' motion.",
      values: { mu3b: 0.0121505856, xpos: 0.9519486347314083, ypos: 0, zpos: 0, xvel: 0, yvel: -0.952445273435512, zvel: 0, tf: 6.45 },
    },
  ],
  relativemotion: [
    {
      title: "Bounded Safety Ellipse",
      tag: "Circular chief (e=0)",
      desc: "Equal semi-major axis with a small eccentricity/phase offset — a classic bounded relative ellipse.",
      values: { da: 0, de: 0.001, di: 0, dRANN: 0, dargp: 0, dM: 0.05, orbits: 5 },
      radio: "0.0",
    },
    {
      title: "Along-Track Drift",
      tag: "Circular chief (e=0)",
      desc: "A small semi-major-axis offset gives the deputy a different period, so it drifts along-track.",
      values: { da: 5, de: 0, di: 0, dRANN: 0, dargp: 0, dM: 0, orbits: 5 },
      radio: "0.0",
    },
    {
      title: "Cross-Track Weave",
      tag: "Circular chief (e=0)",
      desc: "A small relative inclination produces a bounded cross-track/in-track weave.",
      values: { da: 0, de: 0, di: 0.05, dRANN: 0, dargp: 0, dM: 0, orbits: 5 },
      radio: "0.0",
    },
    {
      title: "Eccentric Chief (e=0.4)",
      tag: "Eccentric chief",
      desc: "A bounded relative orbit about a moderately eccentric chief orbit.",
      values: { da: 0, de: 0.002, di: 0, dRANN: 0, dargp: 0, dM: 0.02, orbits: 5 },
      radio: "0.4",
    },
  ],
};

function renderExamples(orbitType) {
  const list = document.getElementById("examplesList");
  if (!list) return;
  list.innerHTML = "";
  const examples = EXAMPLES[orbitType] || [];
  examples.forEach((ex) => {
    const card = document.createElement("button");
    card.type = "button";
    card.className = "example-card";
    card.innerHTML = `
      <span class="ex-title">${ex.title}</span>
      <span class="ex-desc">${ex.desc}</span>
      <span class="ex-tag">${ex.tag}</span>
    `;
    card.addEventListener("click", () => applyExample(orbitType, ex));
    list.appendChild(card);
  });
}

function applyExample(orbitType, example) {
  if (orbitType === "2body") {
    const mode = example.mode || "cartesian";
    if (mode !== twoBodyMode) {
      twoBodyMode = mode;
      renderInputsNow("2body");
    }
  }
  if (example.radio) {
    const radio = document.querySelector(`input[name="chief_e"][value="${example.radio}"]`);
    if (radio) radio.checked = true;
  }
  if (example.anomalyType) {
    const radio = document.querySelector(`input[name="anomalyType"][value="${example.anomalyType}"]`);
    if (radio) radio.checked = true;
  }
  Object.entries(example.values).forEach(([name, value]) => {
    const el = document.querySelector(`[name="${name}"]`);
    if (el) el.value = value;
  });
  runSimulation();
}

// ---------- Dynamic input templates (same field names/format as before) ----------
const TWO_BODY_CARTESIAN_FIELDS = `
    <div class="field-label">μ (km³/s²)</div>
    <input type="text" name="mu" placeholder="km³/s²">
    <div class="field-label">Position (km)</div>
    <input type="text" name="xpos" placeholder="x (km)">
    <input type="text" name="ypos" placeholder="y (km)">
    <input type="text" name="zpos" placeholder="z (km)">
    <div class="field-label">Velocity (km/s)</div>
    <input type="text" name="xvel" placeholder="vx (km/s)">
    <input type="text" name="yvel" placeholder="vy (km/s)">
    <input type="text" name="zvel" placeholder="vz (km/s)">
    <div class="field-label">Final Time (s)</div>
    <input type="text" name="tf" placeholder="tf (s)">
`;

const TWO_BODY_ELEMENTS_FIELDS = `
    <div class="field-label">μ (km³/s²)</div>
    <input type="text" name="mu" placeholder="km³/s²">
    <div class="field-label">Semi-major axis a (km)</div>
    <input type="text" name="oe_a" placeholder="a (km)">
    <div class="field-label">Eccentricity e (0 ≤ e &lt; 1)</div>
    <input type="text" name="oe_e" placeholder="e">
    <div class="field-label">Inclination i (deg)</div>
    <input type="text" name="oe_i" placeholder="i (deg)">
    <div class="field-label">RAAN Ω (deg)</div>
    <input type="text" name="oe_RAAN" placeholder="Ω (deg)">
    <div class="field-label">Argument of periapsis ω (deg)</div>
    <input type="text" name="oe_argp" placeholder="ω (deg)">
    <fieldset>
      <legend>Spacecraft Position (optional)</legend>
      <label class="radio-row"><input type="radio" name="anomalyType" value="none" checked> Shape only</label>
      <label class="radio-row"><input type="radio" name="anomalyType" value="true"> True anomaly ν</label>
      <label class="radio-row"><input type="radio" name="anomalyType" value="mean"> Mean anomaly M</label>
      <div class="field-label">Anomaly value (deg)</div>
      <input type="text" name="anomalyValue" placeholder="deg">
    </fieldset>
`;

function twoBodyInputsHTML(mode) {
  return `
    <h3>Inputs</h3>
    <div class="field-label">Input Mode</div>
    <select name="twoBodyInputMode">
      <option value="cartesian"${mode === "cartesian" ? " selected" : ""}>Cartesian State Vector</option>
      <option value="elements"${mode === "elements" ? " selected" : ""}>Orbital Elements</option>
    </select>
    ${mode === "elements" ? TWO_BODY_ELEMENTS_FIELDS : TWO_BODY_CARTESIAN_FIELDS}
  `;
}

const INPUT_GROUPS = {
  cr3bp: `
    <h3>Inputs</h3>
    <div class="field-label">μ (mass ratio)</div>
    <input type="text" name="mu3b" placeholder="mass ratio">
    <div class="field-label">Position (nondimensional)</div>
    <input type="text" name="xpos" placeholder="x">
    <input type="text" name="ypos" placeholder="y">
    <input type="text" name="zpos" placeholder="z">
    <div class="field-label">Velocity (nondimensional)</div>
    <input type="text" name="xvel" placeholder="vx">
    <input type="text" name="yvel" placeholder="vy">
    <input type="text" name="zvel" placeholder="vz">
    <div class="field-label">Final Time (nondimensional)</div>
    <input type="text" name="tf" placeholder="tf">
  `,
  relativemotion: `
    <h3>Inputs</h3>
    <fieldset>
      <legend>Chief Orbit Eccentricity</legend>
      <label class="radio-row"><input type="radio" name="chief_e" value="0.0" checked> 0.0</label>
      <label class="radio-row"><input type="radio" name="chief_e" value="0.2"> 0.2</label>
      <label class="radio-row"><input type="radio" name="chief_e" value="0.4"> 0.4</label>
      <label class="radio-row"><input type="radio" name="chief_e" value="0.6"> 0.6</label>
      <label class="radio-row"><input type="radio" name="chief_e" value="0.8"> 0.8</label>
    </fieldset>
    <div class="field-label">Δ semi-major axis (km)</div>
    <input type="text" name="da" placeholder="Δa (km)">
    <div class="field-label">Δ eccentricity</div>
    <input type="text" name="de" placeholder="Δe">
    <div class="field-label">Δ inclination (deg)</div>
    <input type="text" name="di" placeholder="Δi (deg)">
    <div class="field-label">Δ RAAN (deg)</div>
    <input type="text" name="dRANN" placeholder="ΔΩ (deg)">
    <div class="field-label">Δ argument of perigee (deg)</div>
    <input type="text" name="dargp" placeholder="Δω (deg)">
    <div class="field-label">Δ mean anomaly (deg)</div>
    <input type="text" name="dM" placeholder="ΔM (deg)">
    <div class="field-label">Number of orbits</div>
    <input type="text" name="orbits" placeholder="count">
  `,
};

// Which 2-Body input sub-mode is currently shown ("cartesian" | "elements").
let twoBodyMode = "cartesian";

function renderInputsFor(orbitType) {
  if (orbitType === "2body") return twoBodyInputsHTML(twoBodyMode);
  return INPUT_GROUPS[orbitType];
}

// Re-render the inputs panel immediately (no fade), used when an example
// or the mode toggle needs the DOM updated before values are filled in.
function renderInputsNow(orbitType) {
  const inputsFieldsDiv = document.querySelector(".inputs-fields");
  if (!inputsFieldsDiv) return;
  inputsFieldsDiv.innerHTML = renderInputsFor(orbitType);
}

// ---------- Boot ----------
document.addEventListener("DOMContentLoaded", () => {
  const orbitTypeSelect = document.getElementById("orbitType");
  const inputsFieldsDiv = document.querySelector(".inputs-fields");
  const runBtn = document.getElementById("runBtn");

  function updateInputs() {
    const selected = orbitTypeSelect.value;
    if (selected === "2body") twoBodyMode = "cartesian"; // reset to the simple default each time we arrive at 2-Body
    inputsFieldsDiv.style.opacity = 0;
    setTimeout(() => {
      inputsFieldsDiv.innerHTML = renderInputsFor(selected);
      inputsFieldsDiv.style.opacity = 1;
      renderExamples(selected);
      if (window.Plotly?.Plots?.resize) Plotly.Plots.resize(document.getElementById("mainPlot"));
    }, 150);
  }

  orbitTypeSelect.addEventListener("change", updateInputs);

  // Delegated listener: switching the 2-Body Cartesian/Elements toggle
  // re-renders just the inputs fields, without touching orbit type or examples.
  inputsFieldsDiv.addEventListener("change", (e) => {
    if (e.target && e.target.name === "twoBodyInputMode") {
      twoBodyMode = e.target.value;
      renderInputsNow("2body");
    }
  });

  runBtn.addEventListener("click", (e) => {
    e.preventDefault();
    runSimulation();
  });

  updateInputs();

  // Auto-run a sensible default on load
  setTimeout(() => runSimulation(), 200);

  window.addEventListener("resize", () => {
    if (window.Plotly?.Plots?.resize) {
      Plotly.Plots.resize(document.getElementById("mainPlot"));
      Plotly.Plots.resize(document.getElementById("conservedPlot"));
    }
  });
});
