// =============================================================
// Main_plot.js — Updated (uses mainPlot & conservedPlot ids)
// =============================================================

let conservedData = [];

// ---------- RK4 ----------
function rk4Step(f, state, dt, ...extraArgs) {
  const k1 = f(state, ...extraArgs);
  const k2 = f(state.map((s, i) => s + 0.5 * dt * k1[i]), ...extraArgs);
  const k3 = f(state.map((s, i) => s + 0.5 * dt * k2[i]), ...extraArgs);
  const k4 = f(state.map((s, i) => s + dt * k3[i]), ...extraArgs);
  return state.map((s, i) => s + dt / 6 * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]));
}

// ---------- 2-BODY ----------
let MU = 398600; // default Earth GM

function twoBodyODE(state) {
  const r = state.slice(0,3);
  const v = state.slice(3,6);
  const rNorm = Math.hypot(...r);
  const a = r.map(ri => -MU * ri / (rNorm**3));
  return [...v, ...a];
}

function orbitalPeriod(state) {
  const r = state.slice(0,3), v = state.slice(3,6);
  const rNorm = Math.hypot(...r), vNorm = Math.hypot(...v);
  const energy = (vNorm*vNorm)/2 - MU/rNorm;
  const a = -MU / (2 * energy);
  if (a <= 0) return null;
  return 2 * Math.PI * Math.sqrt(a**3 / MU);
}

/*function drawOrbit2Body(points) {
  const xs = points.map(p => p[0]);
  const ys = points.map(p => p[1]);
  const zs = points.map(p => p[2]);

  const orbitTrace = {
    x: xs, y: ys, z: zs,
    mode: 'lines', type: 'scatter3d',
    line: { color: '#ff5555', width: 3 }, name: 'Orbit Path'
  };

  const center = { type: 'scatter3d', x:[0], y:[0], z:[0], mode:'markers', marker:{color:'#3b82f6', size:6}, name:'Primary' };

  const maxR = Math.max(
    Math.max(...xs.map(v => Math.abs(v))),
    Math.max(...ys.map(v => Math.abs(v))),
    Math.max(...zs.map(v => Math.abs(v)))
  ) || 1;
  const axisRange = [-1.2 * maxR, 1.2 * maxR];

  const layout = {
    title: { text: '3D Two-Body Orbit', font: { color: '#f1f1f1', size: 18 } },
    scene: {
      xaxis: { title: 'x (km)', range: axisRange, color: '#ccc', backgroundcolor: '#111' },
      yaxis: { title: 'y (km)', range: axisRange, color: '#ccc', backgroundcolor: '#111' },
      zaxis: { title: 'z (km)', range: axisRange, color: '#ccc', backgroundcolor: '#111' },
      aspectmode: 'cube', bgcolor: '#111'
    },
    paper_bgcolor: '#0d0d0d', plot_bgcolor: '#0d0d0d', margin: { t:40, l:0, r:0, b:0 }
  };

  Plotly.newPlot('mainPlot', [orbitTrace, center], layout, { displayModeBar: false });
}*/

function drawOrbit2Body(points) {
  const xs = points.map(p => p[0]);
  const ys = points.map(p => p[1]);
  const zs = points.map(p => p[2]);

  const orbitTrace = {
    x: xs, y: ys, z: zs,
    mode: 'lines', type: 'scatter3d',
    line: { color: '#ff5555', width: 3 },
    name: 'Orbit'
  };

  const center = {
    type: 'scatter3d',
    x:[0], y:[0], z:[0],
    mode:'markers',
    marker:{color:'#3b82f6', size:6},
    name:'Primary'
  };

  const maxR = Math.max(
    Math.max(...xs.map(v => Math.abs(v))),
    Math.max(...ys.map(v => Math.abs(v))),
    Math.max(...zs.map(v => Math.abs(v)))
  ) || 1;
  const axisRange = [-1.2 * maxR, 1.2 * maxR];

  const layout = {
    scene: {
      xaxis: { title: 'x (km)', range: axisRange, color: '#ccc' },
      yaxis: { title: 'y (km)', range: axisRange, color: '#ccc' },
      zaxis: { title: 'z (km)', range: axisRange, color: '#ccc' },
      aspectmode: 'cube',
      bgcolor: '#111',
      annotations: [
        {
          showarrow: false,
          text: "<b>3D Two-Body Orbit</b>",
          x: 0.5, y: 1.08, z: 0,
          xref: "paper", yref: "paper",
          font: { color: "#f1f1f1", size: 17 }
        }
      ]
    },
    showlegend: true,
    legend: {
      x: 0.02, y: 0.98,
      bgcolor: "rgba(0,0,0,0.4)",
      borderwidth: 0
    },
    margin: { t: 10, l: 0, r: 0, b: 0 },
    paper_bgcolor: '#0d0d0d',
    plot_bgcolor: '#0d0d0d'
  };

  Plotly.newPlot('mainPlot', [orbitTrace, center], layout, { displayModeBar: false });
}


// ---------- CR3BP ----------
function CR3B_ode(state, mu) {
  const x = state[0], y = state[1], z = state[2];
  const vx = state[3], vy = state[4], vz = state[5];

  const r1_vec = [x + mu, y, z];
  const r2_vec = [x - (1 - mu), y, z];
  const r1 = Math.hypot(...r1_vec);
  const r2 = Math.hypot(...r2_vec);

  const ax = 2*vy + x - ( (1-mu)*r1_vec[0]/r1**3 ) - ( mu*r2_vec[0]/r2**3 );
  const ay = -2*vx + y - ( (1-mu)*r1_vec[1]/r1**3 ) - ( mu*r2_vec[1]/r2**3 );
  const az = - ( (1-mu)*r1_vec[2]/r1**3 ) - ( mu*r2_vec[2]/r2**3 );

  return [vx, vy, vz, ax, ay, az];
}

/*function drawOrbitCR3B(points, mu) {
  const xs = points.map(p => p[0]);
  const ys = points.map(p => p[1]);
  const zs = points.map(p => p[2]);

  const orbitTrace = { x: xs, y: ys, z: zs, mode: 'lines', type: 'scatter3d', line: { width: 2 }, name: 'Trajectory' };
  const p1 = { type:'scatter3d', x:[-mu], y:[0], z:[0], mode:'markers', marker:{size:4}, name:'Primary1' };
  const p2 = { type:'scatter3d', x:[1-mu], y:[0], z:[0], mode:'markers', marker:{size:4}, name:'Primary2' };

  const layout = {
    title: { text: 'CR3BP Trajectory', font: { color: '#f1f1f1', size: 18 } },
    scene: { aspectmode: 'cube', bgcolor: '#111' },
    paper_bgcolor: '#0d0d0d', plot_bgcolor: '#0d0d0d'
  };

  Plotly.newPlot('mainPlot', [orbitTrace, p1, p2], layout, { displayModeBar: false });
}*/

/*function drawOrbitCR3B(points, mu) {
  const xs = points.map(p => p[0]);
  const ys = points.map(p => p[1]);
  const zs = points.map(p => p[2]);

  const orbitTrace = {
    x: xs, y: ys, z: zs,
    mode: 'lines', type: 'scatter3d',
    line: { width: 2, color: '#4da6ff' },
    name: 'Trajectory'
  };

  const p1 = {
    type:'scatter3d', x:[-mu], y:[0], z:[0],
    mode:'markers', marker:{size:5, color:'#ff7f0e'},
    name:'Primary1'
  };

  const p2 = {
    type:'scatter3d', x:[1-mu], y:[0], z:[0],
    mode:'markers', marker:{size:5, color:'#2ca02c'},
    name:'Primary2'
  };

  const layout = {
    scene: {
      aspectmode: 'cube',
      bgcolor: '#111',
      annotations: [
        {
          showarrow: false,
          text: "<b>CR3BP Trajectory</b>",
          x: 0.5, y: 1.08, z: 0,
          xref: "paper", yref: "paper",
          font: { color: "#f1f1f1", size: 17 }
        }
      ]
    },
    showlegend: true,
    legend: {
      x: 0.02, y: 0.98,
      bgcolor: "rgba(0,0,0,0.4)",
      borderwidth: 0
    },
    margin: { t: 10, l: 0, r: 0, b: 0 },
    paper_bgcolor: '#0d0d0d',
    plot_bgcolor: '#0d0d0d'
  };

  Plotly.newPlot('mainPlot', [orbitTrace, p1, p2], layout, { displayModeBar: false });
}*/

/*function drawOrbitCR3B(points, mu) {
  const xs = points.map(p => p[0]);
  const ys = points.map(p => p[1]);
  const zs = points.map(p => p[2]);

  // ---- Correct Axis Scaling (Use true min/max, not abs) ----
  const xmin = Math.min(...xs), xmax = Math.max(...xs);
  const ymin = Math.min(...ys), ymax = Math.max(...ys);
  const zmin = Math.min(...zs), zmax = Math.max(...zs);
  const pad = 0.1;
  const xrange = [xmin - pad, xmax + pad];
  const yrange = [ymin - pad, ymax + pad];
  const zrange = [zmin - pad, zmax + pad];

  // ---- Traces ----
  const orbitTrace = {
    x: xs, y: ys, z: zs,
    mode: 'lines', type: 'scatter3d',
    line: { width: 2, color: '#4da6ff' },
    name: 'Trajectory'
  };

  const p1 = {
    type:'scatter3d', x:[-mu], y:[0], z:[0],
    mode:'markers', marker:{size:5, color:'#ff7f0e'},
    name:'Primary1'
  };

  const p2 = {
    type:'scatter3d', x:[1-mu], y:[0], z:[0],
    mode:'markers', marker:{size:5, color:'#2ca02c'},
    name:'Primary2'
  };

  // ---- Layout w/ Title + Legend Inside, but keep proper scaling ----
  const layout = {
    scene: {
      xaxis: { title: 'x', range: xrange, color: '#ccc' },
      yaxis: { title: 'y', range: yrange, color: '#ccc' },
      zaxis: { title: 'z', range: zrange, color: '#ccc' },
      aspectmode: 'data', // <--- important! maintains proper proportional shape
      bgcolor: '#111',

      annotations: [
        {
          showarrow: false,
          text: "<b>CR3BP Trajectory</b>",
          x: 0.5, y: 1.07, z: 0,
          xref: "paper", yref: "paper",
          font: { color: "#f1f1f1", size: 17 }
        }
      ]
    },

    showlegend: true,
    legend: {
      x: 0.02, y: 0.98,
      bgcolor: "rgba(0,0,0,0.45)",
      borderwidth: 0,
      font: { color: "#ddd" }
    },

    margin: { t: 10, l: 0, r: 0, b: 0 },
    paper_bgcolor: '#0d0d0d',
    plot_bgcolor: '#0d0d0d'
  };

  Plotly.newPlot('mainPlot', [orbitTrace, p1, p2], layout, { displayModeBar: false });
}*/

function drawOrbitCR3B(points, mu) {
  const xs = points.map(p => p[0]);
  const ys = points.map(p => p[1]);
  const zs = points.map(p => p[2]);

  // ---- Compute cube axis scaling ----
  const all = [...xs, ...ys, ...zs];
  const minVal = Math.min(...all), maxVal = Math.max(...all);
  const pad = 0.1 * (maxVal - minVal || 1);
  const min = minVal - pad;
  const max = maxVal + pad;

  const xrange = [min, max];
  const yrange = [min, max];
  const zrange = [min, max];

  // ---- Traces ----
  const orbitTrace = {
    x: xs, y: ys, z: zs,
    mode: 'lines', type: 'scatter3d',
    line: { width: 2, color: '#4da6ff' },
    name: 'Trajectory'
  };

  const p1 = {
    type:'scatter3d', x:[-mu], y:[0], z:[0],
    mode:'markers', marker:{size:5, color:'#ff7f0e'},
    name:'Primary1'
  };

  const p2 = {
    type:'scatter3d', x:[1-mu], y:[0], z:[0],
    mode:'markers', marker:{size:5, color:'#2ca02c'},
    name:'Primary2'
  };

  // ---- Layout: cube, title inside, legend inside ----
  const layout = {
    scene: {
      xaxis: { title: 'x', range: xrange, color: '#ccc' },
      yaxis: { title: 'y', range: yrange, color: '#ccc' },
      zaxis: { title: 'z', range: zrange, color: '#ccc' },
      aspectmode: 'cube',   // <--- back to true cube
      bgcolor: '#111',

      annotations: [
        {
          showarrow: false,
          text: "<b>CR3BP Trajectory</b>",
          x: 0.5, y: 1.07,
          xref: "paper", yref: "paper",
          font: { color: "#f1f1f1", size: 17 }
        }
      ]
    },

    showlegend: true,
    legend: {
      x: 0.02, y: 0.98,
      bgcolor: "rgba(0,0,0,0.45)",
      borderwidth: 0,
      font: { color: "#ddd" }
    },

    margin: { t: 10, l: 0, r: 0, b: 0 },
    paper_bgcolor: '#0d0d0d',
    plot_bgcolor: '#0d0d0d'
  };

  Plotly.newPlot('mainPlot', [orbitTrace, p1, p2], layout, { displayModeBar: false });
}



// ---------- Conserved plotting helper (calls conserved.js function) ----------
function plotConserved(values, label) {
  // if user created conserved.js with a function of same name, this will call that one.
  // otherwise basic fallback here:
  if (typeof window.plotConserved === 'function' && window.plotConserved !== plotConserved) {
    // call the external conserved.js implementation
    window.plotConserved(values, label);
    return;
  }

  // fallback simple plot
  if (!values || values.length === 0) {
    Plotly.newPlot("conservedPlot", []);
    return;
  }
  Plotly.newPlot("conservedPlot", [{
    x: [...Array(values.length).keys()],
    y: values,
    mode: "lines",
    type: "scatter",
    line: { width: 2 }
  }], {
    title: label,
    paper_bgcolor: "#0d0d0d",
    plot_bgcolor: "#0d0d0d",
    font: { color: "#ddd" },
    margin: { t:30, l:40, r:10, b:30 }
  }, { displayModeBar: false });
}

// ---------- SIMULATORS ----------

function simulate2Body(defaultState = null) {
  // Clear plot & conserved
  conservedData = [];
  Plotly.newPlot("conservedPlot", []);

  // build initial state
  let r0, v0;
  if (defaultState) {
    MU = 398600;
    r0 = defaultState.slice(0,3);
    v0 = defaultState.slice(3,6);
  } else {
    const getVal = (name) => parseFloat(document.querySelector(`input[name="${name}"]`)?.value) || 0;
    MU = getVal("mu") || MU;
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

  // conserved: |h|
  const cross = (a,b) => [
    a[1]*b[2] - a[2]*b[1],
    a[2]*b[0] - a[0]*b[2],
    a[0]*b[1] - a[1]*b[0]
  ];
  const h0 = Math.hypot(...cross(r0, v0));
  conservedData = [];

  for (let i = 0; i < steps; i++) {
    state = rk4Step(twoBodyODE, state, dt);
    trajectory.push(state.slice(0,3));

    const r = state.slice(0,3), v = state.slice(3,6);
    const h = Math.hypot(...cross(r, v));
    conservedData.push(((h - h0)/ h0) * 100);
  }

  drawOrbit2Body(trajectory);
  plotConserved(conservedData, "Angular Momentum Over Time");
}

function simulateCR3B() {
  conservedData = [];
  Plotly.newPlot("conservedPlot", []);

  const getVal = (name, def=0) => parseFloat(document.querySelector(`input[name="${name}"]`)?.value) || def;
  const MU3B = getVal("mu3b", 0.0121505856);
  const r0 = [getVal("xpos", 1.0220215127), getVal("ypos", 0), getVal("zpos", -0.1820967615)];
  const v0 = [getVal("xvel", 0), getVal("yvel", -0.1032563411), getVal("zvel", 0)];
  const state0 = [...r0, ...v0];

  const t0 = 0, tf = getVal("tf", 1.5);
  const steps = 2000;
  const dt = (tf - t0) / steps;

  let state = state0;
  const trajectory = [r0];

  function jacobi(s, mu) {
    const x = s[0], y = s[1], z = s[2];
    const vx = s[3], vy = s[4], vz = s[5];
    const r1 = Math.hypot(x + mu, y, z);
    const r2 = Math.hypot(x - (1 - mu), y, z);
    return x*x + y*y + (2*(1-mu))/r1 + (2*mu)/r2 - (vx*vx + vy*vy + vz*vz);
  }

  const C0 = jacobi(state0, MU3B);
  conservedData = [];

  for (let i = 0; i < steps; i++) {
    state = rk4Step(CR3B_ode, state, dt, MU3B);
    trajectory.push(state.slice(0,3));
    const C = jacobi(state, MU3B);
    conservedData.push(((C - C0) / C0) * 100);
  }

  drawOrbitCR3B(trajectory, MU3B);
  plotConserved(conservedData, "Jacobi Constant Over Time");
}

function simulateRelativeMotion() {
  // placeholder — clear conserved plot
  conservedData = [];
  Plotly.newPlot("conservedPlot", []);
  alert("Relative Motion visualization not implemented yet.");
}

// ---------- Controller + Run button wiring ----------
function runSimulation() {
  const orbitType = document.getElementById("orbitType")?.value;
  if (orbitType === "2body") simulate2Body();
  else if (orbitType === "cr3bp") simulateCR3B();
  else if (orbitType === "relativemotion") simulateRelativeMotion();
  else alert("Unknown orbit type selected.");
}

document.addEventListener("DOMContentLoaded", () => {
  // Attach Run button (it is dynamically recreated on orbit-type change in your HTML)
  const inputsDiv = document.querySelector(".inputs");

  const attachButtonHandler = () => {
    const btn = inputsDiv.querySelector("button");
    if (btn) {
      btn.removeEventListener("click", runSimulation); // safe remove
      btn.addEventListener("click", (e) => { e.preventDefault(); runSimulation(); });
    }
  };

  // Observe dynamic input changes and reattach handler
  const observer = new MutationObserver(() => { attachButtonHandler(); });
  observer.observe(inputsDiv, { childList: true, subtree: true });

  attachButtonHandler();

  // Auto-run initial default orbit (use same default state as before)
  const defaultState = [7000, 0, 0, 0, 7.5, 0];
  simulate2Body(defaultState);

  // Resize handler to keep Plotly responsive
  window.addEventListener('resize', () => {
    if (window.Plotly?.Plots?.resize) {
      Plotly.Plots.resize(document.getElementById('mainPlot'));
      Plotly.Plots.resize(document.getElementById('conservedPlot'));
    }
  });
});
