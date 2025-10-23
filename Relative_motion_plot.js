// two_body_stm_plot.js
// Two-body propagation + State Transition Matrix (STM) integration
// Equivalent to the given Python version using RK4

// Earth gravitational parameter [km^3/s^2]
const MU = 398600.4418;

// Two-body + STM derivative
function twoBodySTM(state) {
  const r = state.slice(0, 3);
  const v = state.slice(3, 6);
  const Phi_flat = state.slice(6);
  const Phi = [];
  for (let i = 0; i < 6; i++) Phi.push(Phi_flat.slice(i * 6, (i + 1) * 6));

  const rNorm = Math.hypot(...r);
  const a = r.map((ri) => -MU * ri / (rNorm ** 3));

  // Gravity gradient (G matrix)
  const G = Array(3)
    .fill(0)
    .map(() => Array(3).fill(0));
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      if (i === j) {
        G[i][j] = -MU * (1 / (rNorm ** 3) - 3 * r[i] * r[j] / (rNorm ** 5));
      } else {
        G[i][j] = 3 * MU * r[i] * r[j] / (rNorm ** 5);
      }
    }
  }

  // A matrix (6x6)
  const A = Array(6)
    .fill(0)
    .map(() => Array(6).fill(0));
  for (let i = 0; i < 3; i++) A[i][i + 3] = 1; // upper right = I
  for (let i = 0; i < 3; i++)
    for (let j = 0; j < 3; j++) A[i + 3][j] = G[i][j]; // lower left = G

  // dPhi = A * Phi
  const dPhi = multiplyMatrices(A, Phi);

  // Pack derivative
  const dYdt = new Array(42).fill(0);
  for (let i = 0; i < 3; i++) dYdt[i] = v[i];
  for (let i = 0; i < 3; i++) dYdt[i + 3] = a[i];
  for (let i = 0; i < 6; i++)
    for (let j = 0; j < 6; j++) dYdt[6 + 6 * i + j] = dPhi[i][j];

  return dYdt;
}

// Matrix multiply helper
function multiplyMatrices(A, B) {
  const rowsA = A.length,
    colsA = A[0].length,
    colsB = B[0].length;
  const result = Array(rowsA)
    .fill(0)
    .map(() => Array(colsB).fill(0));
  for (let i = 0; i < rowsA; i++)
    for (let j = 0; j < colsB; j++)
      for (let k = 0; k < colsA; k++) result[i][j] += A[i][k] * B[k][j];
  return result;
}

// RK4 stepper
function rk4Step(f, state, dt) {
  const k1 = f(state);
  const k2 = f(state.map((s, i) => s + 0.5 * dt * k1[i]));
  const k3 = f(state.map((s, i) => s + 0.5 * dt * k2[i]));
  const k4 = f(state.map((s, i) => s + dt * k3[i]));
  return state.map(
    (s, i) => s + (dt / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
  );
}

// Initialize and simulate one orbit
function simulateTwoBodySTM() {
  // Initial circular orbit
  const r0 = [7000, 0, 0];
  const v0 = [0, Math.sqrt(MU / r0[0]), 0];
  const Phi0 = identityMatrix(6).flat();
  const Y0 = [...r0, ...v0, ...Phi0];

  const a = r0[0];
  const T = 2 * Math.PI * Math.sqrt(a ** 3 / MU);
  const steps = 1000;
  const dt = T / steps;

  let Y = Y0;
  const trajectory = [r0];

  for (let i = 0; i < steps; i++) {
    Y = rk4Step(twoBodySTM, Y, dt);
    trajectory.push(Y.slice(0, 3));
  }

  drawOrbit(trajectory);
}

// 6x6 identity matrix
function identityMatrix(n) {
  return Array(n)
    .fill(0)
    .map((_, i) =>
      Array(n)
        .fill(0)
        .map((_, j) => (i === j ? 1 : 0))
    );
}

// 3D orbit plot
function drawOrbit(points) {
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

  const earthTrace = {
    type: "scatter3d",
    x: [0],
    y: [0],
    z: [0],
    mode: "markers",
    marker: { color: "#3b82f6", size: 8 },
    name: "Earth",
  };

  const maxR = Math.max(...xs.map(Math.abs));
  const layout = {
    title: {
      text: "Two-Body STM Orbit Visualization",
      font: { color: "#f1f1f1", size: 20 },
    },
    scene: {
      xaxis: { title: "x (km)", range: [-maxR, maxR], color: "#ccc" },
      yaxis: { title: "y (km)", range: [-maxR, maxR], color: "#ccc" },
      zaxis: { title: "z (km)", range: [-maxR, maxR], color: "#ccc" },
      bgcolor: "#111",
      aspectmode: "cube",
    },
    paper_bgcolor: "#0d0d0d",
    plot_bgcolor: "#0d0d0d",
  };

  Plotly.newPlot("mainPlot", [orbitTrace, earthTrace], layout, {
    displayModeBar: false,
  });
}

document.addEventListener("DOMContentLoaded", simulateTwoBodySTM);
