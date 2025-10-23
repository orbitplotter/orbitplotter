// ======================================================
// COMBINED ORBIT VISUALIZER
// Supports: 2-Body, 3-Body, Relative Motion (placeholder)
// ======================================================

// ====== GLOBAL CONSTANTS ======
const MU_EARTH = 398600.4418; // km^3/s^2
const MU3B = 1.0 / (81.30059 + 1); // Earth-Moon system example
let currentMode = "2body"; // default mode

// ======================================================
// ========== 2-BODY FUNCTIONS ==========================
// ======================================================
function twoBodyODE(state) {
  const r = state.slice(0, 3);
  const v = state.slice(3, 6);
  const rNorm = Math.hypot(...r);
  const a = r.map((ri) => -MU_EARTH * ri / (rNorm ** 3));
  return [...v, ...a];
}

function orbitalPeriod(state) {
  const r = state.slice(0, 3);
  const v = state.slice(3, 6);
  const rNorm = Math.hypot(...r);
  const vNorm = Math.hypot(...v);
  const energy = vNorm * vNorm / 2 - MU_EARTH / rNorm;
  const a = -MU_EARTH / (2 * energy);
  if (a <= 0) return null;
  return 2 * Math.PI * Math.sqrt(a ** 3 / MU_EARTH);
}

// ======================================================
// ========== 3-BODY FUNCTIONS ==========================
// ======================================================
function CR3B_ode(state, mu = MU3B) {
  const r = state.slice(0, 3);
  const v = state.slice(3, 6);

  const r1_vec = [r[0] + mu, r[1], r[2]];
  const r2_vec = [r[0] - (1 - mu), r[1], r[2]];

  const r1 = Math.hypot(...r1_vec);
  const r2 = Math.hypot(...r2_vec);

  const ax =
    2 * v[1] +
    r[0] -
    ((1 - mu) * r1_vec[0]) / r1 ** 3 -
    (mu * r2_vec[0]) / r2 ** 3;
  const ay =
    -2 * v[0] +
    r[1] -
    ((1 - mu) * r1_vec[1]) / r1 ** 3 -
    (mu * r2_vec[1]) / r2 ** 3;
  const az =
    -((1 - mu) * r1_vec[2]) / r1 ** 3 -
    (mu * r2_vec[2]) / r2 ** 3;

  return [...v, ax, ay, az];
}

// ======================================================
// ========== RELATIVE MOTION (placeholder) =============
// ======================================================
function simulateRelativeMotion() {
  alert("Relative Motion mode not implemented yet.");
}

// ======================================================
// ========== RK4 INTEGRATOR (shared) ===================
// ======================================================
function rk4Step(f, state, dt, mu = null) {
  const k1 = mu !== null ? f(state, mu) : f(state);
  const k2 = mu !== null ? f(state.map((s, i) => s + 0.5 * dt * k1[i]), mu) : f(state.map((s, i) => s + 0.5 * dt * k1[i]));
  const k3 = mu !== null ? f(state.map((s, i) => s + 0.5 * dt * k2[i]), mu) : f(state.map((s, i) => s + 0.5 * dt * k2[i]));
  const k4 = mu !== null ? f(state.map((s, i) => s + dt * k3[i]), mu) : f(state.map((s, i) => s + dt * k3[i]));
  return state.map((s, i) => s + (dt / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]));
}

// ======================================================
// ========== DRAW FUNCTIONS ============================
// ======================================================
function drawOrbit(points, options = {}) {
  const { title = "3D Orbit Visualization", bodies = [] } = options;

  const xs = points.map((p) => p[0]);
  const ys = points.map((p) => p[1]);
  const zs = points.map((p) => p[2]);

  const orbitTrace = {
    x: xs,
    y: ys,
    z: zs,
    mode: "lines",
    type: "scatter3d",
    line: { color: "#ff5555", width: 4 },
    name: "Trajectory",
  };

  const traces = [orbitTrace, ...bodies];

  const maxR = Math.max(
    Math.max(...xs.map((v) => Math.abs(v))),
    Math.max(...ys.map((v) => Math.abs(v))),
    Math.max(...zs.map((v) => Math.abs(v)))
  );
  const axisRange = [-1.2 * maxR, 1.2 * maxR];

  const layout = {
    title: { text: title, font: { color: "#f1f1f1", size: 20 } },
    scene: {
      xaxis: { title: "x (km)", range: axisRange, color: "#ccc", backgroundcolor: "#111", gridcolor: "#333" },
      yaxis: { title: "y (km)", range: axisRange, color: "#ccc", backgroundcolor: "#111", gridcolor: "#333" },
      zaxis: { title: "z (km)", range: axisRange, color: "#ccc", backgroundcolor: "#111", gridcolor: "#333" },
      aspectmode: "cube",
      bgcolor: "#111",
    },
    paper_bgcolor: "#0d0d0d",
    plot_bgcolor: "#0d0d0d",
    margin: { l: 0, r: 0, b: 0, t: 40 },
    showlegend: true,
    legend: { x: 0, y: 1, font: { color: "#eee" }, bgcolor: "rgba(0,0,0,0)" },
  };

  Plotly.newPlot("mainPlot", traces, layout, { displayModeBar: false });
}

// ======================================================
// ========== SIMULATION WRAPPERS =======================
// ======================================================
function simulate2Body(r0 = [7000, 0, 0], v0 = [0, 7.5, 0]) {
  const state0 = [...r0, ...v0];
  const T = orbitalPeriod(state0);
  if (T === null) {
    alert("Not an elliptical orbit.");
    return;
  }

  const steps = 1000;
  const dt = T / steps;
  let state = state0;
  const trajectory = [r0];

  for (let i = 0; i < steps; i++) {
    state = rk4Step(twoBodyODE, state, dt);
    trajectory.push(state.slice(0, 3));
  }

  drawOrbit(trajectory, {
    title: "2-Body Orbit",
    bodies: [{
      type: "scatter3d",
      x: [0], y: [0], z: [0],
      mode: "markers",
      marker: { color: "#3b82f6", size: 8 },
      name: "Earth"
    }]
  });
}

function simulate3Body() {
  const r0 = [1.0220215127, 0.0, -0.1820967615];
  const v0 = [0.0, -0.1032563411, 0.0];
  const state0 = [...r0, ...v0];
  const steps = 2000;
  const t0 = 0.0, tf = 1.511111111111111;
  const dt = (tf - t0) / steps;

  let state = state0;
  const trajectory = [r0];

  for (let i = 0; i < steps; i++) {
    state = rk4Step(CR3B_ode, state, dt, MU3B);
    trajectory.push(state.slice(0, 3));
  }

  drawOrbit(trajectory, {
    title: "Circular Restricted 3-Body Trajectory",
    bodies: [
      {
        type: "scatter3d",
        x: [-MU3B], y: [0], z: [0],
        mode: "markers",
        marker: { color: "#3b82f6", size: 8 },
        name: "Primary 1",
      },
      {
        type: "scatter3d",
        x: [1 - MU3B], y: [0], z: [0],
        mode: "markers",
        marker: { color: "#eab308", size: 6 },
        name: "Primary 2",
      },
    ]
  });
}

// ======================================================
// ========== MODE HANDLER ==============================
// ======================================================
function simulateOrbit() {
  if (currentMode === "2body") simulate2Body();
  else if (currentMode === "3body") simulate3Body();
  else if (currentMode === "relative") simulateRelativeMotion();
}

// ======================================================
// ========== EVENT HOOKS ===============================
// ======================================================
document.addEventListener("DOMContentLoaded", () => {
  // Run a default 2-body plot with [7000,0,0,0,7.5,0]
  simulate2Body([7000, 0, 0], [0, 7.5, 0]);
});
