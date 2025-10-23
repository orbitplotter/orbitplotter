// =============================================================
// MAIN ORBIT VISUALIZATION SCRIPT
// Supports: 2-Body | CR3BP | Relative Motion (placeholder)
// =============================================================

// ========== COMMON UTILS ==========
function rk4Step(f, state, dt, ...extraArgs) {
  const k1 = f(state, ...extraArgs);
  const k2 = f(state.map((s, i) => s + 0.5 * dt * k1[i]), ...extraArgs);
  const k3 = f(state.map((s, i) => s + 0.5 * dt * k2[i]), ...extraArgs);
  const k4 = f(state.map((s, i) => s + dt * k3[i]), ...extraArgs);
  return state.map((s, i) => s + dt / 6 * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]));
}

// ========== 2-BODY ==========

const MU = 398600.4418; // km^3/s^2

function twoBodyODE(state) {
  const r = state.slice(0, 3);
  const v = state.slice(3, 6);
  const rNorm = Math.hypot(...r);
  const a = r.map((ri) => -MU * ri / (rNorm ** 3));
  return [...v, ...a];
}

function orbitalPeriod(state) {
  const r = state.slice(0, 3);
  const v = state.slice(3, 6);
  const rNorm = Math.hypot(...r);
  const vNorm = Math.hypot(...v);
  const energy = vNorm * vNorm / 2 - MU / rNorm;
  const a = -MU / (2 * energy);
  if (a <= 0) return null;
  return 2 * Math.PI * Math.sqrt(a ** 3 / MU);
}

function drawOrbit2Body(points) {
  const xs = points.map(p => p[0]);
  const ys = points.map(p => p[1]);
  const zs = points.map(p => p[2]);

  const orbitTrace = {
    x: xs, y: ys, z: zs,
    mode: 'lines', type: 'scatter3d',
    line: { color: '#ff5555', width: 4 },
    name: 'Orbit Path'
  };

  const earthTrace = {
    type: 'scatter3d',
    x: [0], y: [0], z: [0],
    mode: 'markers',
    marker: { color: '#3b82f6', size: 8 },
    name: 'Earth'
  };

  const maxR = Math.max(
    Math.max(...xs.map(v => Math.abs(v))),
    Math.max(...ys.map(v => Math.abs(v))),
    Math.max(...zs.map(v => Math.abs(v)))
  );
  const axisRange = [-1.2 * maxR, 1.2 * maxR];

  const layout = {
    title: { text: '3D Two-Body Orbit', font: { color: '#f1f1f1', size: 20 } },
    scene: {
      xaxis: { title: 'x (km)', range: axisRange, color: '#ccc', backgroundcolor: '#111' },
      yaxis: { title: 'y (km)', range: axisRange, color: '#ccc', backgroundcolor: '#111' },
      zaxis: { title: 'z (km)', range: axisRange, color: '#ccc', backgroundcolor: '#111' },
      aspectmode: 'cube', bgcolor: '#111'
    },
    paper_bgcolor: '#0d0d0d', plot_bgcolor: '#0d0d0d', margin: { l:0, r:0, b:0, t:40 }
  };

  Plotly.newPlot('mainPlot', [orbitTrace, earthTrace], layout, { displayModeBar: false });
}

function simulate2Body(defaultState = null) {
  let r0, v0;
  if (defaultState) {
    // Use provided initial state (auto-start)
    r0 = defaultState.slice(0, 3);
    v0 = defaultState.slice(3, 6);
  } else {
    // Get values from form inputs
    const getVal = (name) => parseFloat(document.querySelector(`input[name="${name}"]`)?.value) || 0;
    r0 = [getVal("xpos"), getVal("ypos"), getVal("zpos")];
    v0 = [getVal("xvel"), getVal("yvel"), getVal("zvel")];
  }

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

  drawOrbit2Body(trajectory);
}

// ========== CR3BP ==========

const MU3B = 1.0 / (81.30059 + 1); // Earth–Moon system

function CR3B_ode(state, mu = MU3B) {
  const r = state.slice(0, 3);
  const v = state.slice(3, 6);
  const r1_vec = [r[0] + mu, r[1], r[2]];
  const r2_vec = [r[0] - (1 - mu), r[1], r[2]];
  const r1 = Math.hypot(...r1_vec);
  const r2 = Math.hypot(...r2_vec);

  const ax = 2*v[1] + r[0] - ((1 - mu)*r1_vec[0])/r1**3 - (mu*r2_vec[0])/r2**3;
  const ay = -2*v[0] + r[1] - ((1 - mu)*r1_vec[1])/r1**3 - (mu*r2_vec[1])/r2**3;
  const az = -((1 - mu)*r1_vec[2])/r1**3 - (mu*r2_vec[2])/r2**3;

  return [...v, ax, ay, az];
}

function drawOrbitCR3B(points, mu = MU3B) {
  const xs = points.map(p => p[0]);
  const ys = points.map(p => p[1]);
  const zs = points.map(p => p[2]);

  const orbitTrace = { x: xs, y: ys, z: zs, mode: 'lines', type: 'scatter3d', line: { color: '#ff5555', width: 3 }, name: 'Trajectory' };
  const earthTrace = { type: 'scatter3d', x: [-mu], y: [0], z: [0], mode: 'markers', marker: { color: '#3b82f6', size: 8 }, name: 'Primary 1' };
  const moonTrace  = { type: 'scatter3d', x: [1 - mu], y: [0], z: [0], mode: 'markers', marker: { color: '#eab308', size: 6 }, name: 'Primary 2' };

  const maxR = Math.max(...xs.map(v => Math.abs(v)), ...ys.map(v => Math.abs(v)), ...zs.map(v => Math.abs(v)));
  const axisRange = [-1.2 * maxR, 1.2 * maxR];

  const layout = {
    title: { text: 'CR3BP Trajectory Visualization', font: { color: '#f1f1f1', size: 20 } },
    scene: {
      xaxis: { title: 'x', range: axisRange, color: '#ccc', backgroundcolor: '#111' },
      yaxis: { title: 'y', range: axisRange, color: '#ccc', backgroundcolor: '#111' },
      zaxis: { title: 'z', range: axisRange, color: '#ccc', backgroundcolor: '#111' },
      aspectmode: 'cube', bgcolor: '#111'
    },
    paper_bgcolor: '#0d0d0d', plot_bgcolor: '#0d0d0d', margin: { l:0, r:0, b:0, t:40 }
  };

  Plotly.newPlot('mainPlot', [orbitTrace, earthTrace, moonTrace], layout, { displayModeBar: false });
}

function simulateCR3B() {
  const getVal = (name, def) => parseFloat(document.querySelector(`input[name="${name}"]`)?.value) || def;
  const r0 = [getVal("xpos", 1.0220215127), getVal("ypos", 0.0), getVal("zpos", -0.1820967615)];
  const v0 = [getVal("xvel", 0.0), getVal("yvel", -0.1032563411), getVal("zvel", 0.0)];
  const state0 = [...r0, ...v0];

  const t0 = 0.0, tf = getVal("tf", 1.5);
  const steps = 2000;
  const dt = (tf - t0) / steps;

  let state = state0;
  const trajectory = [r0];
  for (let i = 0; i < steps; i++) {
    state = rk4Step(CR3B_ode, state, dt, MU3B);
    trajectory.push(state.slice(0, 3));
  }

  drawOrbitCR3B(trajectory, MU3B);
}

// ========== RELATIVE MOTION ==========
function simulateRelativeMotion() {
  alert("Relative Motion visualization not implemented yet.");
}

// ========== MAIN CONTROLLER ==========
function runSimulation() {
  const orbitType = document.getElementById("orbitType")?.value;
  if (orbitType === "2body") simulate2Body();
  else if (orbitType === "cr3bp") simulateCR3B();
  else if (orbitType === "relativemotion") simulateRelativeMotion();
  else alert("Unknown orbit type selected.");
}

// ========== EVENT HOOK ==========
document.addEventListener("DOMContentLoaded", () => {
  // Attach "Run" button handler (re-attaches on input reload)
  const attachButtonHandler = () => {
    const btn = document.querySelector(".inputs button");
    if (btn) {
      btn.addEventListener("click", (e) => {
        e.preventDefault();
        runSimulation();
      });
    }
  };

  // Observe dynamic input changes
  const inputsDiv = document.querySelector(".inputs");
  const observer = new MutationObserver(attachButtonHandler);
  observer.observe(inputsDiv, { childList: true });

  attachButtonHandler();

  // ✅ Auto-run the default 2-body orbit [7000,0,0,0,7.5,0]
  const defaultState = [7000, 0, 0, 0, 7.5, 0];
  simulate2Body(defaultState);
});
