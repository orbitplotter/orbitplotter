/*// 2_body_plot.js

// Earth gravitational parameter [km^3/s^2]
const MU = 398600.4418;

// Two-body ODE
function twoBodyODE(state) {
  const r = state.slice(0,3);
  const v = state.slice(3,6);
  const rNorm = Math.hypot(...r);
  const a = r.map((ri) => -MU * ri / (rNorm**3));
  return [...v, ...a];
}

// RK4 integrator
function rk4Step(f, state, dt) {
  const k1 = f(state);
  const k2 = f(state.map((s,i) => s + 0.5*dt*k1[i]));
  const k3 = f(state.map((s,i) => s + 0.5*dt*k2[i]));
  const k4 = f(state.map((s,i) => s + dt*k3[i]));
  return state.map((s,i) => s + dt/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]));
}

// Compute orbital period for elliptical orbits
function orbitalPeriod(state) {
  const r = state.slice(0,3);
  const v = state.slice(3,6);
  const rNorm = Math.hypot(...r);
  const vNorm = Math.hypot(...v);
  const energy = vNorm*vNorm/2 - MU/rNorm;
  const a = -MU/(2*energy);
  if(a <= 0) return null; // not elliptical
  return 2 * Math.PI * Math.sqrt(a**3 / MU);
}

// Draw orbit (2D x-y projection)
function drawOrbit(points) {
  const canvas = document.getElementById("orbitCanvas");
  const ctx = canvas.getContext("2d");
  ctx.clearRect(0,0,canvas.width,canvas.height);

  // Find max distance to scale
  const xs = points.map(p => p[0]);
  const ys = points.map(p => p[1]);
  const maxR = Math.max(...xs.map(x => Math.abs(x)), ...ys.map(y => Math.abs(y)), 1);

  const scale = canvas.width/2.2 / maxR;
  const cx = canvas.width/2;
  const cy = canvas.height/2;

  // Draw Earth
  ctx.beginPath();
  ctx.arc(cx, cy, 8, 0, 2*Math.PI);
  ctx.fillStyle = "blue";
  ctx.fill();

  // Draw orbit
  ctx.beginPath();
  ctx.moveTo(cx + xs[0]*scale, cy - ys[0]*scale);
  for(let i=1;i<points.length;i++){
    ctx.lineTo(cx + xs[i]*scale, cy - ys[i]*scale);
  }
  ctx.strokeStyle = "red";
  ctx.lineWidth = 2;
  ctx.stroke();
}

// Main simulation
function simulateOrbit() {
  const r0 = [
    parseFloat(document.getElementById("xpos").value),
    parseFloat(document.getElementById("ypos").value),
    parseFloat(document.getElementById("zpos").value)
  ];
  const v0 = [
    parseFloat(document.getElementById("xvel").value),
    parseFloat(document.getElementById("yvel").value),
    parseFloat(document.getElementById("zvel").value)
  ];
  const state0 = [...r0, ...v0];

  const T = orbitalPeriod(state0);
  if(T === null){
    alert("Not an elliptical orbit. Cannot simulate one orbit.");
    return;
  }

  const steps = 1000;
  const dt = T/steps;
  let state = state0;
  const trajectory = [r0];

  for(let i=0;i<steps;i++){
    state = rk4Step(twoBodyODE, state, dt);
    trajectory.push(state.slice(0,3));
  }

  drawOrbit(trajectory);
}

// Hook button
document.addEventListener("DOMContentLoaded",()=>{
  const btn = document.querySelector(".input-container button");
  btn.addEventListener("click",(e)=>{
    e.preventDefault();
    simulateOrbit();
  });

  // Optional: plot default orbit on load
  simulateOrbit();
});
*/

// 2_body_plot.js

// Earth gravitational parameter [km^3/s^2]
const MU = 398600.4418;

// Two-body ODE
function twoBodyODE(state) {
  const r = state.slice(0, 3);
  const v = state.slice(3, 6);
  const rNorm = Math.hypot(...r);
  const a = r.map((ri) => -MU * ri / (rNorm ** 3));
  return [...v, ...a];
}

// RK4 integrator
function rk4Step(f, state, dt) {
  const k1 = f(state);
  const k2 = f(state.map((s, i) => s + 0.5 * dt * k1[i]));
  const k3 = f(state.map((s, i) => s + 0.5 * dt * k2[i]));
  const k4 = f(state.map((s, i) => s + dt * k3[i]));
  return state.map((s, i) => s + dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]));
}

// Compute orbital period for elliptical orbits
function orbitalPeriod(state) {
  const r = state.slice(0, 3);
  const v = state.slice(3, 6);
  const rNorm = Math.hypot(...r);
  const vNorm = Math.hypot(...v);
  const energy = vNorm * vNorm / 2 - MU / rNorm;
  const a = -MU / (2 * energy);
  if (a <= 0) return null; // not elliptical
  return 2 * Math.PI * Math.sqrt(a ** 3 / MU);
}

// Draw orbit with Plotly 3D (dark mode + auto fit)
function drawOrbit(points) {
  const xs = points.map(p => p[0]);
  const ys = points.map(p => p[1]);
  const zs = points.map(p => p[2]);

  const orbitTrace = {
    x: xs,
    y: ys,
    z: zs,
    mode: 'lines',
    type: 'scatter3d',
    line: { color: '#ff5555', width: 4 },
    name: 'Orbit Path'
  };

  // Earth marker
  const earthTrace = {
    type: 'scatter3d',
    x: [0],
    y: [0],
    z: [0],
    mode: 'markers',
    marker: { color: '#3b82f6', size: 8 },
    name: 'Earth'
  };

  // Determine bounds for auto-fit
  const maxR = Math.max(
    Math.max(...xs.map(v => Math.abs(v))),
    Math.max(...ys.map(v => Math.abs(v))),
    Math.max(...zs.map(v => Math.abs(v)))
  );

  const axisRange = [-1.2 * maxR, 1.2 * maxR];

  const layout = {
    title: {
      text: '3D Orbit Visualization',
      font: { color: '#f1f1f1', size: 20 }
    },
    scene: {
      xaxis: {
        title: 'x (km)',
        range: axisRange,
        color: '#cccccc',
        backgroundcolor: '#111',
        gridcolor: '#333',
        zerolinecolor: '#444'
      },
      yaxis: {
        title: 'y (km)',
        range: axisRange,
        color: '#cccccc',
        backgroundcolor: '#111',
        gridcolor: '#333',
        zerolinecolor: '#444'
      },
      zaxis: {
        title: 'z (km)',
        range: axisRange,
        color: '#cccccc',
        backgroundcolor: '#111',
        gridcolor: '#333',
        zerolinecolor: '#444'
      },
      aspectmode: 'cube',
      bgcolor: '#111'
    },
    paper_bgcolor: '#0d0d0d',
    plot_bgcolor: '#0d0d0d',
    margin: { l: 0, r: 0, b: 0, t: 40 },
    showlegend: true,
    legend: {
      x: 0,
      y: 1,
      font: { color: '#eee' },
      bgcolor: 'rgba(0,0,0,0)'
    }
  };

  Plotly.newPlot('mainPlot', [orbitTrace, earthTrace], layout, { displayModeBar: false });
}

// Main simulation
function simulateOrbit() {
  const getVal = (name) =>
    parseFloat(document.querySelector(`input[name="${name}"]`).value) || 0;

  const r0 = [getVal("xpos"), getVal("ypos"), getVal("zpos")];
  const v0 = [getVal("xvel"), getVal("yvel"), getVal("zvel")];
  const state0 = [...r0, ...v0];

  const T = orbitalPeriod(state0);
  if (T === null) {
    alert("Not an elliptical orbit. Cannot simulate one orbit.");
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

  drawOrbit(trajectory);
}

// Hook button
document.addEventListener("DOMContentLoaded", () => {
  const btn = document.querySelector('.inputs button');
  if (btn) {
    btn.addEventListener("click", (e) => {
      e.preventDefault();
      simulateOrbit();
    });
  }

  // Optional: default initial orbit
  simulateOrbit();
});
