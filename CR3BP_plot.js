// 3_body_plot.js
// Circular Restricted 3-Body Problem (CR3BP) Visualization

// CR3B parameters
const MU3B = 1.0 / (81.30059 + 1); // Earth-Moon system example

// CR3B ODE
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

// RK4 integrator
function rk4Step(f, state, dt, mu = MU3B) {
  const k1 = f(state, mu);
  const k2 = f(state.map((s, i) => s + 0.5 * dt * k1[i]), mu);
  const k3 = f(state.map((s, i) => s + 0.5 * dt * k2[i]), mu);
  const k4 = f(state.map((s, i) => s + dt * k3[i]), mu);
  return state.map(
    (s, i) => s + (dt / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
  );
}

// Draw CR3B trajectory with Plotly 3D
function drawOrbit(points, mu = MU3B) {
  const xs = points.map((p) => p[0]);
  const ys = points.map((p) => p[1]);
  const zs = points.map((p) => p[2]);

  const orbitTrace = {
    x: xs,
    y: ys,
    z: zs,
    mode: "lines",
    type: "scatter3d",
    line: { color: "#ff5555", width: 3 },
    name: "Trajectory",
  };

  // Primary bodies
  const earthTrace = {
    type: "scatter3d",
    x: [-mu],
    y: [0],
    z: [0],
    mode: "markers",
    marker: { color: "#3b82f6", size: 8 },
    name: "Primary 1",
  };

  const moonTrace = {
    type: "scatter3d",
    x: [1 - mu],
    y: [0],
    z: [0],
    mode: "markers",
    marker: { color: "#eab308", size: 6 },
    name: "Primary 2",
  };

  const maxR = Math.max(
    ...xs.map((v) => Math.abs(v)),
    ...ys.map((v) => Math.abs(v)),
    ...zs.map((v) => Math.abs(v))
  );
  const axisRange = [-1.2 * maxR, 1.2 * maxR];

  const layout = {
    title: {
      text: "CR3B Trajectory Visualization",
      font: { color: "#f1f1f1", size: 20 },
    },
    scene: {
      xaxis: {
        title: "x",
        range: axisRange,
        color: "#ccc",
        backgroundcolor: "#111",
        gridcolor: "#333",
        zerolinecolor: "#444",
      },
      yaxis: {
        title: "y",
        range: axisRange,
        color: "#ccc",
        backgroundcolor: "#111",
        gridcolor: "#333",
        zerolinecolor: "#444",
      },
      zaxis: {
        title: "z",
        range: axisRange,
        color: "#ccc",
        backgroundcolor: "#111",
        gridcolor: "#333",
        zerolinecolor: "#444",
      },
      aspectmode: "cube",
      bgcolor: "#111",
    },
    paper_bgcolor: "#0d0d0d",
    plot_bgcolor: "#0d0d0d",
    margin: { l: 0, r: 0, b: 0, t: 40 },
    showlegend: true,
    legend: {
      x: 0,
      y: 1,
      font: { color: "#eee" },
      bgcolor: "rgba(0,0,0,0)",
    },
  };

  Plotly.newPlot("mainPlot", [orbitTrace, earthTrace, moonTrace], layout, {
    displayModeBar: false,
  });
}

// Main simulation
function simulateCR3B() {
  const r0 = [
    parseFloat(document.querySelector('input[name="xpos"]').value) || 1.0220215127,
    parseFloat(document.querySelector('input[name="ypos"]').value) || 0.0,
    parseFloat(document.querySelector('input[name="zpos"]').value) || -0.1820967615,
  ];
  const v0 = [
    parseFloat(document.querySelector('input[name="xvel"]').value) || 0.0,
    parseFloat(document.querySelector('input[name="yvel"]').value) || -0.1032563411,
    parseFloat(document.querySelector('input[name="zvel"]').value) || 0.0,
  ];

  const state0 = [...r0, ...v0];
  const t0 = 0.0;
  const tf = 1.511111111111111; // final time
  const steps = 2000;
  const dt = (tf - t0) / steps;

  let state = state0;
  const trajectory = [r0];

  for (let i = 0; i < steps; i++) {
    state = rk4Step(CR3B_ode, state, dt, MU3B);
    trajectory.push(state.slice(0, 3));
  }

  drawOrbit(trajectory, MU3B);
}

// Hook button
document.addEventListener("DOMContentLoaded", () => {
  const btn = document.querySelector(".inputs button");
  if (btn) {
    btn.addEventListener("click", (e) => {
      e.preventDefault();
      simulateCR3B();
    });
  }

  // Run default orbit
  simulateCR3B();
});
