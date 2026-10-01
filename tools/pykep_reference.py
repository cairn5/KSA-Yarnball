"""
Reference results from pykep 2 for the engine's tests (engine/tests/data/pykep_reference.json).

Covers Kepler propagation (pykep.propagate_lagrangian), unpowered flybys (pykep.fb_prop) and the
MGA-1DSM model's burns (pykep.trajopt.mga_1dsm) on random inputs, for the Cassini sequence.
Times in the file are pykep's mjd2000 (days from 2000-01-01 00:00); the engine counts from 12:00.

Run in the conda environment that has pykep 2 (on Windows, with conda's Library\\bin on PATH):
    python tools/pykep_reference.py
"""

import json
import math
import os

import numpy as np
import pykep as pk

rng = np.random.default_rng(2026)
MU = pk.MU_SUN

# --- Kepler propagation: elliptical and hyperbolic, both directions in time ---------------------
kepler = []
for _ in range(200):
    r = rng.normal(0, 1, 3)
    r *= rng.uniform(0.3, 30) * pk.AU / np.linalg.norm(r)
    v_circ = math.sqrt(MU / np.linalg.norm(r))
    v = rng.normal(0, 1, 3)
    v *= rng.uniform(0.2, 1.8) * v_circ / np.linalg.norm(v)
    dt = rng.uniform(-2, 2) * 365.25 * pk.DAY2SEC
    r1, v1 = pk.propagate_lagrangian(tuple(r), tuple(v), dt, MU)
    kepler.append({"r": list(r / 1e3), "v": list(v / 1e3), "dt": dt, "r1": [c / 1e3 for c in r1], "v1": [c / 1e3 for c in v1]})

# --- Unpowered flybys ------------------------------------------------------------------------
fb = []
venus = pk.planet.jpl_lp("venus")
for _ in range(200):
    v_pla = rng.normal(0, 30_000, 3)
    v_in = v_pla + rng.normal(0, 8_000, 3)
    rp = rng.uniform(1.05, 50) * venus.radius
    beta = rng.uniform(-2 * math.pi, 2 * math.pi)
    v_out = pk.fb_prop(tuple(v_in), tuple(v_pla), rp, beta, venus.mu_self)
    fb.append({"v_in": list(v_in / 1e3), "v_pla": list(v_pla / 1e3), "rp": rp / 1e3, "beta": beta, "v_out": [c / 1e3 for c in v_out]})

# --- MGA-1DSM burns on the Cassini sequence --------------------------------------------------
names = ["earth", "venus", "venus", "earth", "jupiter", "saturn"]
seq = [pk.planet.jpl_lp(n) for n in names]
legs = [[100, 400], [100, 500], [30, 300], [400, 1600], [800, 2200]]
udp = pk.trajopt.mga_1dsm(seq=seq, t0=[-3000, 1000], tof=legs, vinf=[2.5, 12.0], add_vinf_dep=False,
                          add_vinf_arr=True, tof_encoding="direct", eta_lb=0.01, eta_ub=0.9, rp_ub=200)
lb, ub = map(np.array, udp.get_bounds())
mga = []
while len(mga) < 300:
    x = lb + rng.uniform(0, 1, len(lb)) * (ub - lb)
    try:
        dv = udp._compute_dvs(list(x))[0]
    except (RuntimeError, ValueError):
        continue
    if all(np.isfinite(dv)):
        mga.append({"x": list(x), "dv": [d / 1e3 for d in dv]})

out = {
    "pykep": pk.__version__,
    "sequence": names,
    "kepler": kepler,
    "flyby": fb,
    "mga_1dsm": mga,
}
path = os.path.join(os.path.dirname(__file__), "..", "engine", "tests", "data", "pykep_reference.json")
os.makedirs(os.path.dirname(path), exist_ok=True)
with open(path, "w") as f:
    json.dump(out, f)
print(f"wrote {path}: {len(kepler)} propagations, {len(fb)} flybys, {len(mga)} chromosomes")
