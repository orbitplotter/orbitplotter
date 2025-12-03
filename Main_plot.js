// =============================================================
// Main_plot.js
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
function twoBodyODE(state, MU) {
  const r = state.slice(0,3);
  const v = state.slice(3,6);
  const rNorm = Math.hypot(...r);
  const a = r.map(ri => -MU * ri / (rNorm**3));
  return [...v, ...a];
}

function orbitalPeriod(state, MU, tf) {
  const r = state.slice(0,3), v = state.slice(3,6);
  const rNorm = Math.hypot(...r), vNorm = Math.hypot(...v);
  const energy = (vNorm*vNorm)/2 - MU/rNorm;
  if (energy >= 0) return tf;
  const a = -MU / (2 * energy);
  return 2 * Math.PI * Math.sqrt(a**3 / MU);
}

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
      bgcolor: '#111'
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
      bgcolor: '#111'
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

// ---------- Relative Motion ----------
function oeToStateVector(a, e, i, RAAN, w, M, mu = 398600) {
    // ---- helpers ----
    const deg2rad = d => d * Math.PI / 180;

    i = deg2rad(i);
    RAAN = deg2rad(RAAN);
    w = deg2rad(w);
    M = deg2rad(M);

    // ---- Solve Kepler's Equation: M = E - e sinE ----
    let E = M;
    for (let k = 0; k < 15; k++) {
        E = E + (M - (E - e * Math.sin(E))) / (1 - e * Math.cos(E));
    }

    // ---- True anomaly ----
    const nu = 2 * Math.atan2(
        Math.sqrt(1 + e) * Math.sin(E / 2),
        Math.sqrt(1 - e) * Math.cos(E / 2)
    );

    // ---- Perifocal radius ----
    const r_pf = [
        a * (Math.cos(E) - e),
        a * Math.sqrt(1 - e * e) * Math.sin(E),
        0
    ];

    // ---- Perifocal velocity ----
    const p = a * (1 - e * e);
    const r = Math.sqrt(r_pf[0] ** 2 + r_pf[1] ** 2);
    const v_pf = [
        -Math.sqrt(mu / p) * Math.sin(nu),
        Math.sqrt(mu / p) * (e + Math.cos(nu)),
        0
    ];

    // ---- Rotation matrix from perifocal → ECI ----
    const cosO = Math.cos(RAAN), sinO = Math.sin(RAAN);
    const cosi = Math.cos(i),    sini = Math.sin(i);
    const cosw = Math.cos(w),    sinw = Math.sin(w);

    const R11 =  cosO*cosw - sinO*sinw*cosi;
    const R12 = -cosO*sinw - sinO*cosw*cosi;
    const R13 =  sinO*sini;

    const R21 =  sinO*cosw + cosO*sinw*cosi;
    const R22 = -sinO*sinw + cosO*cosw*cosi;
    const R23 = -cosO*sini;

    const R31 =  sinw*sini;
    const R32 =  cosw*sini;
    const R33 =  cosi;

    // ---- Transform to ECI ----
    const r_eci = [
        R11 * r_pf[0] + R12 * r_pf[1] + R13 * r_pf[2],
        R21 * r_pf[0] + R22 * r_pf[1] + R23 * r_pf[2],
        R31 * r_pf[0] + R32 * r_pf[1] + R33 * r_pf[2]
    ];

    const v_eci = [
        R11 * v_pf[0] + R12 * v_pf[1] + R13 * v_pf[2],
        R21 * v_pf[0] + R22 * v_pf[1] + R23 * v_pf[2],
        R31 * v_pf[0] + R32 * v_pf[1] + R33 * v_pf[2]
    ];

    return [r_eci, v_eci];
}

function eci_to_lvlh(r_c, v_c, r_d, v_d) {
    // --- Helper vector functions ---
    const norm = v => Math.sqrt(v[0]**2 + v[1]**2 + v[2]**2);

    const subtract = (a, b) => [
        a[0] - b[0],
        a[1] - b[1],
        a[2] - b[2]
    ];

    const divide = (v, s) => [v[0]/s, v[1]/s, v[2]/s];

    const cross = (a, b) => [
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    ];

    const matVec = (M, v) => [
        M[0][0]*v[0] + M[0][1]*v[1] + M[0][2]*v[2],
        M[1][0]*v[0] + M[1][1]*v[1] + M[1][2]*v[2],
        M[2][0]*v[0] + M[2][1]*v[1] + M[2][2]*v[2]
    ];

    // --- LVLH basis vectors ---
    const x_hat = divide(r_c, norm(r_c));

    let z_hat = cross(r_c, v_c);
    z_hat = divide(z_hat, norm(z_hat));

    const y_hat = cross(z_hat, x_hat);

    // Rotation matrix: ECI → LVLH
    const R = [
        x_hat,
        y_hat,
        z_hat
    ];

    // Relative position and velocity in ECI
    const r_rel = subtract(r_d, r_c);
    const v_rel = subtract(v_d, v_c);

    // Convert to LVLH
    const r_lvlh = matVec(R, r_rel);
    const v_lvlh = matVec(R, v_rel);

    return { r_lvlh, v_lvlh };
}

function drawOrbitRelative(r_lvlh) {
    const xs = r_lvlh.map(p => p[0]); // radial
    const ys = r_lvlh.map(p => p[1]); // in-track
    const zs = r_lvlh.map(p => p[2]); // cross-track

    const trace = {
        x: ys,   // LVLH swap y,x,z
        y: xs,
        z: zs,
        mode: 'lines',
        type: 'scatter3d',
        line: { color: '#ff5555', width: 3 },
        name: 'Deputy in LVLH'
    };

    const chief = {
        x: [0], y: [0], z: [0],
        mode: 'markers',
        type: 'scatter3d',
        marker: { color: '#3b82f6', size: 6 },
        name: 'Chief'
    };

    const layout = {
        scene: {
            xaxis: { title: 'In-track y [km]', color: '#ccc' },
            yaxis: { title: 'Radial x [km]', color: '#ccc' },
            zaxis: { title: 'Cross-track z [km]', color: '#ccc' },
            aspectmode: 'cube',
            bgcolor: '#111'
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

    Plotly.newPlot('mainPlot', [trace, chief], layout, { displayModeBar: false });
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

  // helper to read inputs with defaults
  const getVal = (name, def = 0) => {
    const el = document.querySelector(`input[name="${name}"]`);
    const v = el ? parseFloat(el.value) : NaN;
    return Number.isFinite(v) ? v : def;
  };

  let MU, tf, r0, v0;
  if (defaultState) {
    MU = 398600;
    tf = 86400;
    r0 = defaultState.slice(0,3);
    v0 = defaultState.slice(3,6);
  } else {
    MU = getVal("mu", 398600);
    tf = getVal("tf", 86400);
    r0 = [getVal("xpos", 7000), getVal("ypos", 0), getVal("zpos", 0)];
    v0 = [getVal("xvel", 0), getVal("yvel", 7.5), getVal("zvel", 0)];
  }

  const state0 = [...r0, ...v0];
  const T = orbitalPeriod(state0, MU, tf);
  const steps = 1000;
  const dt = (T && T > 0) ? (T / steps) : (tf / steps || 1);
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
    state = rk4Step(twoBodyODE, state, dt, MU);
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

  const t0 = 0, tf = getVal("tf", 1.51111);
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

    const getVal = (name, def = 0) => {
        const el = document.querySelector(`[name="${name}"]`);
        const v = el ? parseFloat(el.value) : NaN;
        return Number.isFinite(v) ? v : def;
    };

    // Chief eccentricity from radio buttons
    const chief_e = parseFloat(
        document.querySelector('input[name="chief_e"]:checked').value
    );

    // Deputy orbit element deltas
    const da   = getVal("da", 0);
    const de   = getVal("de", 0);
    const di   = getVal("di", 0);
    const dRAAN = getVal("dRANN", 0);
    const dargp = getVal("dargp", 0);
    const dM    = getVal("dM", 0);

    const periods = getVal("orbits", 1);

    const mu = 398600;

    // --------------------------------------------------------
    // Chief: use the symbolic a values (from your table)
    // --------------------------------------------------------
    const eccentricityMap = {
        0.0: 8000,
        0.2: 10000,
        0.4: 13333.3333333333,
        0.6: 20000,
        0.8: 40000
    };

    const a_chief = eccentricityMap[chief_e];

    // Chief orbital elements (baseline)
    const oe_c = {
        a: a_chief,
        e: chief_e,
        i: 0,
        RAAN: 0,
        w: 0,
        M: 0
    };

    // Deputy orbital elements (chief + deltas)
    const oe_d = {
        a: a_chief + da,
        e: chief_e + de,
        i: 0 + di,
        RAAN: 0 + dRAAN,
        w: 0 + dargp,
        M: 0 + dM
    };

    // --------------------------------------------------------
    // Convert COEs → state vectors (your exact function)
    // --------------------------------------------------------
    const [r_c0, v_c0] = oeToStateVector(
        oe_c.a, oe_c.e, oe_c.i, oe_c.RAAN, oe_c.w, oe_c.M
    );

    const [r_d0, v_d0] = oeToStateVector(
        oe_d.a, oe_d.e, oe_d.i, oe_d.RAAN, oe_d.w, oe_d.M
    );

    const x0chief = [...r_c0, ...v_c0];
    const x0deputy = [...r_d0, ...v_d0];

    // --------------------------------------------------------
    // Determine the orbital period from chief state
    // --------------------------------------------------------
    const T = orbitalPeriod(x0chief, mu, 86400);
    const totalTime = T * periods;

    // Dense sampling (Python used ~2000/orbit)
    const samplesPerOrbit = 2000;
    const N = samplesPerOrbit * periods;

    const dt = totalTime / N;

    // --------------------------------------------------------
    // Propagate both satellites with RK4
    // --------------------------------------------------------
    let xc = x0chief.slice();
    let xd = x0deputy.slice();

    const chiefStates = [];
    const deputyStates = [];

    for (let k = 0; k < N; k++) {
        chiefStates.push(xc.slice());
        deputyStates.push(xd.slice());

        xc = rk4Step(twoBodyODE, xc, dt, mu);
        xd = rk4Step(twoBodyODE, xd, dt, mu);
    }

    // --------------------------------------------------------
    // Convert each sample to LVLH frame
    // --------------------------------------------------------
    const r_lvlh = [];
    const v_lvlh = [];

    for (let i = 0; i < N; i++) {
        const rc = chiefStates[i].slice(0,3);
        const vc = chiefStates[i].slice(3,6);

        const rd = deputyStates[i].slice(0,3);
        const vd = deputyStates[i].slice(3,6);

        const lv = eci_to_lvlh(rc, vc, rd, vd);
        r_lvlh.push(lv.r_lvlh);
        v_lvlh.push(lv.v_lvlh);
    }

    drawOrbitRelative(r_lvlh);
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
