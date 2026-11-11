// =============================================================
// physics.js
// Core orbit-propagation math for the Orbit Plotter.
//
// Two-body & relative motion are propagated with the *exact*
// closed-form universal-variable Kepler solution (Curtis, Orbital
// Mechanics for Engineering Students, Alg. 3.4 / Vallado's
// universal Kepler algorithm) instead of a fixed-step numerical
// integrator. This works for elliptic, parabolic and hyperbolic
// orbits alike and conserves energy & angular momentum to machine
// precision - there is no truncation error to accumulate.
//
// The CR3BP has no closed-form solution, so it is propagated with
// an adaptive-step Dormand-Prince RK45 integrator (embedded 4th/5th
// order pair with local error control), which is far more accurate
// per function-evaluation than the fixed-step RK4 it replaces.
// =============================================================

var EARTH_MU = 398600.4418; // km^3/s^2 (standard gravitational parameter, Earth)

// ---------- Stumpff functions ----------
function stumpffC(z) {
  if (z > 1e-6) {
    const sz = Math.sqrt(z);
    return (1 - Math.cos(sz)) / z;
  } else if (z < -1e-6) {
    const sz = Math.sqrt(-z);
    return (Math.cosh(sz) - 1) / (-z);
  }
  return 0.5 - z / 24 + (z * z) / 720;
}

function stumpffS(z) {
  if (z > 1e-6) {
    const sz = Math.sqrt(z);
    return (sz - Math.sin(sz)) / (sz * sz * sz);
  } else if (z < -1e-6) {
    const sz = Math.sqrt(-z);
    return (Math.sinh(sz) - sz) / (sz * sz * sz);
  }
  return 1 / 6 - z / 120 + (z * z) / 5040;
}

// ---------- Exact two-body propagation (universal variables) ----------
// state = [rx, ry, rz, vx, vy, vz]; advances by dt seconds exactly
// (elliptic, parabolic or hyperbolic - alpha = 1/a determines which).
function keplerStep(state, dt, mu) {
  if (dt === 0) return state.slice();

  const r0v = state.slice(0, 3);
  const v0v = state.slice(3, 6);
  const r0 = Math.hypot(...r0v);
  const v0 = Math.hypot(...v0v);
  const vr0 = (r0v[0] * v0v[0] + r0v[1] * v0v[1] + r0v[2] * v0v[2]) / r0;
  const sqrtMu = Math.sqrt(mu);
  const alpha = 2 / r0 - (v0 * v0) / mu; // = 1/a

  // Initial guess for the universal anomaly chi
  let chi;
  if (alpha > 1e-8) {
    chi = sqrtMu * dt * alpha;
  } else if (alpha < -1e-8) {
    const a = 1 / alpha; // negative for hyperbola
    const sign = dt >= 0 ? 1 : -1;
    const num = -2 * mu * alpha * dt;
    const den =
      (r0v[0] * v0v[0] + r0v[1] * v0v[1] + r0v[2] * v0v[2]) +
      sign * Math.sqrt(-mu * a) * (1 - r0 * alpha);
    chi = sign * Math.sqrt(-a) * Math.log(Math.abs(num / den));
  } else {
    // near-parabolic fallback (rare edge case)
    chi = (sqrtMu * dt) / r0;
  }

  // Newton-Raphson solve for chi
  let ratio = 1, iter = 0;
  let C = 0.5, S = 1 / 6, z = 0, r = r0;
  while (Math.abs(ratio) > 1e-10 && iter < 100) {
    z = chi * chi * alpha;
    C = stumpffC(z);
    S = stumpffS(z);
    r = chi * chi * C + ((vr0 * r0) / sqrtMu) * chi * (1 - z * S) + r0 * (1 - z * C);
    const F =
      sqrtMu * dt -
      chi ** 3 * S -
      ((vr0 * r0) / sqrtMu) * chi ** 2 * C -
      r0 * chi * (1 - z * S);
    ratio = F / r;
    chi += ratio;
    iter++;
  }

  const f = 1 - (chi * chi / r0) * C;
  const g = dt - (chi ** 3 / sqrtMu) * S;

  const rVec = [
    f * r0v[0] + g * v0v[0],
    f * r0v[1] + g * v0v[1],
    f * r0v[2] + g * v0v[2],
  ];
  const rNew = Math.hypot(...rVec);

  const fdot = (sqrtMu / (rNew * r0)) * (z * S - 1) * chi;
  const gdot = 1 - (chi * chi / rNew) * C;

  const vVec = [
    fdot * r0v[0] + gdot * v0v[0],
    fdot * r0v[1] + gdot * v0v[1],
    fdot * r0v[2] + gdot * v0v[2],
  ];

  return [...rVec, ...vVec];
}

// Specific orbital energy -> "period" used to size the plotting window.
// Bound (elliptic) orbits: returns the true orbital period.
// Unbound (parabolic/hyperbolic) orbits: returns tf unchanged, since
// there is no period - we just show tf seconds of the flyby.
function orbitalWindow(state, mu, tf) {
  const r = state.slice(0, 3), v = state.slice(3, 6);
  const rNorm = Math.hypot(...r), vNorm = Math.hypot(...v);
  const energy = (vNorm * vNorm) / 2 - mu / rNorm;
  if (energy >= 0) return { T: tf, bound: false };
  const a = -mu / (2 * energy);
  return { T: 2 * Math.PI * Math.sqrt(a ** 3 / mu), bound: true };
}

function angularMomentum(r, v) {
  return [
    r[1] * v[2] - r[2] * v[1],
    r[2] * v[0] - r[0] * v[2],
    r[0] * v[1] - r[1] * v[0],
  ];
}

function specificEnergy(r, v, mu) {
  const rNorm = Math.hypot(...r), vNorm = Math.hypot(...v);
  return (vNorm * vNorm) / 2 - mu / rNorm;
}

// ---------- CR3BP equations of motion (rotating synodic frame) ----------
function CR3B_ode(state, mu) {
  const x = state[0], y = state[1], z = state[2];
  const vx = state[3], vy = state[4], vz = state[5];

  const r1_vec = [x + mu, y, z];
  const r2_vec = [x - (1 - mu), y, z];
  const r1 = Math.hypot(...r1_vec);
  const r2 = Math.hypot(...r2_vec);

  const ax = 2 * vy + x - ((1 - mu) * r1_vec[0]) / r1 ** 3 - (mu * r2_vec[0]) / r2 ** 3;
  const ay = -2 * vx + y - ((1 - mu) * r1_vec[1]) / r1 ** 3 - (mu * r2_vec[1]) / r2 ** 3;
  const az = -((1 - mu) * r1_vec[2]) / r1 ** 3 - (mu * r2_vec[2]) / r2 ** 3;

  return [vx, vy, vz, ax, ay, az];
}

function jacobiConstant(s, mu) {
  const x = s[0], y = s[1], z = s[2];
  const vx = s[3], vy = s[4], vz = s[5];
  const r1 = Math.hypot(x + mu, y, z);
  const r2 = Math.hypot(x - (1 - mu), y, z);
  return x * x + y * y + (2 * (1 - mu)) / r1 + (2 * mu) / r2 - (vx * vx + vy * vy + vz * vz);
}

// ---------- Dormand-Prince adaptive RK45 ----------
const DP = {
  a: [
    [],
    [1 / 5],
    [3 / 40, 9 / 40],
    [44 / 45, -56 / 15, 32 / 9],
    [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729],
    [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656],
    [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84],
  ],
  b5: [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84, 0],
  b4: [5179 / 57600, 0, 7571 / 16695, 393 / 640, -92097 / 339200, 187 / 2100, 1 / 40],
};

function dp45Step(f, state, dt, extra, tol) {
  const n = state.length;
  const k = [];
  for (let s = 0; s < 7; s++) {
    const ys = state.slice();
    for (let j = 0; j < s; j++) {
      const aij = DP.a[s][j];
      if (aij !== 0) for (let m = 0; m < n; m++) ys[m] += dt * aij * k[j][m];
    }
    k.push(f(ys, ...extra));
  }
  const y5 = state.map((s, i) => {
    let sum = 0;
    for (let j = 0; j < 7; j++) sum += DP.b5[j] * k[j][i];
    return s + dt * sum;
  });
  const y4 = state.map((s, i) => {
    let sum = 0;
    for (let j = 0; j < 7; j++) sum += DP.b4[j] * k[j][i];
    return s + dt * sum;
  });
  let errNorm = 0;
  for (let i = 0; i < n; i++) {
    const sc = tol.atol + tol.rtol * Math.max(Math.abs(state[i]), Math.abs(y5[i]));
    errNorm += ((y5[i] - y4[i]) / sc) ** 2;
  }
  errNorm = Math.sqrt(errNorm / n);
  return { y5, errNorm };
}

// Adaptive integration from t0 to tf with local error control.
function integrateAdaptive(f, state0, t0, tf, extra, tol = { rtol: 1e-11, atol: 1e-12 }, maxSteps = 20000) {
  let state = state0.slice();
  let t = t0;
  let dt = (tf - t0) / 1000 || 1e-6;
  const traj = [state.slice()];
  let steps = 0;
  while (t < tf && steps < maxSteps) {
    if (t + dt > tf) dt = tf - t;
    const { y5, errNorm } = dp45Step(f, state, dt, extra, tol);
    if (errNorm <= 1 || dt < 1e-10) {
      t += dt;
      state = y5;
      traj.push(state.slice());
      steps++;
      if (errNorm > 0) dt *= Math.min(5, Math.max(0.2, 0.9 * errNorm ** -0.2));
    } else {
      dt *= Math.max(0.1, 0.9 * errNorm ** -0.25);
    }
  }
  return traj;
}

// ---------- Classical orbital elements -> state vector ----------
// Shared core: given true anomaly nu (radians) and elements in radians,
// build the perifocal-frame state and rotate it into ECI.
function elementsToStateVectorRad(a, e, iRad, RAANRad, wRad, nu, mu) {
  const p = a * (1 - e * e);
  const r = p / (1 + e * Math.cos(nu));
  const r_pf = [r * Math.cos(nu), r * Math.sin(nu), 0];
  const v_pf = [-Math.sqrt(mu / p) * Math.sin(nu), Math.sqrt(mu / p) * (e + Math.cos(nu)), 0];

  const cosO = Math.cos(RAANRad), sinO = Math.sin(RAANRad);
  const cosi = Math.cos(iRad), sini = Math.sin(iRad);
  const cosw = Math.cos(wRad), sinw = Math.sin(wRad);

  const R11 = cosO * cosw - sinO * sinw * cosi;
  const R12 = -cosO * sinw - sinO * cosw * cosi;
  const R13 = sinO * sini;
  const R21 = sinO * cosw + cosO * sinw * cosi;
  const R22 = -sinO * sinw + cosO * cosw * cosi;
  const R23 = -cosO * sini;
  const R31 = sinw * sini;
  const R32 = cosw * sini;
  const R33 = cosi;

  const r_eci = [
    R11 * r_pf[0] + R12 * r_pf[1] + R13 * r_pf[2],
    R21 * r_pf[0] + R22 * r_pf[1] + R23 * r_pf[2],
    R31 * r_pf[0] + R32 * r_pf[1] + R33 * r_pf[2],
  ];
  const v_eci = [
    R11 * v_pf[0] + R12 * v_pf[1] + R13 * v_pf[2],
    R21 * v_pf[0] + R22 * v_pf[1] + R23 * v_pf[2],
    R31 * v_pf[0] + R32 * v_pf[1] + R33 * v_pf[2],
  ];

  return [r_eci, v_eci];
}

// Elements + MEAN anomaly (deg) -> ECI state vector. Solves Kepler's
// equation for the eccentric anomaly, then the true anomaly.
// (Unchanged signature/behavior - relative-motion propagation depends on this.)
function oeToStateVector(a, e, i, RAAN, w, M, mu = EARTH_MU) {
  const deg2rad = (d) => (d * Math.PI) / 180;
  const iRad = deg2rad(i), RAANRad = deg2rad(RAAN), wRad = deg2rad(w);
  const Mrad = deg2rad(M);

  let E = Mrad;
  for (let k = 0; k < 30; k++) {
    E = E + (Mrad - (E - e * Math.sin(E))) / (1 - e * Math.cos(E));
  }
  const nu = 2 * Math.atan2(Math.sqrt(1 + e) * Math.sin(E / 2), Math.sqrt(1 - e) * Math.cos(E / 2));

  return elementsToStateVectorRad(a, e, iRad, RAANRad, wRad, nu, mu);
}

// Elements + TRUE anomaly (deg) -> ECI state vector. Direct - no Kepler
// equation to solve, since true anomaly already fixes the position.
function elementsToStateVector(a, e, i, RAAN, w, nuDeg, mu = EARTH_MU) {
  const deg2rad = (d) => (d * Math.PI) / 180;
  return elementsToStateVectorRad(a, e, deg2rad(i), deg2rad(RAAN), deg2rad(w), deg2rad(nuDeg), mu);
}

// Pure-geometry orbit shape (elliptical, 0 <= e < 1): samples the closed
// conic evenly in true anomaly. No propagation, no time - the exact shape
// of the orbit regardless of where the spacecraft actually is on it.
function orbitShapePoints(a, e, i, RAAN, w, mu = EARTH_MU, numPoints = 361) {
  const eClamped = Math.min(Math.max(e, 0), 0.999);
  const points = [];
  for (let k = 0; k < numPoints; k++) {
    const nuDeg = (360 * k) / (numPoints - 1);
    const [r_eci] = elementsToStateVector(a, eClamped, i, RAAN, w, nuDeg, mu);
    points.push(r_eci);
  }
  return points;
}

// ---------- ECI -> LVLH (Hill/RSW-style) frame transform ----------
function eciToLvlh(r_c, v_c, r_d, v_d) {
  const norm = (v) => Math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2);
  const subtract = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const divide = (v, s) => [v[0] / s, v[1] / s, v[2] / s];
  const cross = (a, b) => [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0],
  ];
  const matVec = (M, v) => [
    M[0][0] * v[0] + M[0][1] * v[1] + M[0][2] * v[2],
    M[1][0] * v[0] + M[1][1] * v[1] + M[1][2] * v[2],
    M[2][0] * v[0] + M[2][1] * v[1] + M[2][2] * v[2],
  ];

  const x_hat = divide(r_c, norm(r_c));
  let z_hat = cross(r_c, v_c);
  z_hat = divide(z_hat, norm(z_hat));
  const y_hat = cross(z_hat, x_hat);

  const R = [x_hat, y_hat, z_hat];
  const r_rel = subtract(r_d, r_c);
  const v_rel = subtract(v_d, v_c);

  return { r_lvlh: matVec(R, r_rel), v_lvlh: matVec(R, v_rel) };
}
