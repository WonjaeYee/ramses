"""Why MC with llf/hll makes new extrema: a numpy replica of RAMSES's 1-D MUSCL-Hancock step for a scalar.

    acc_replica.py [results.json]

The scalar on a uniform flow (rho = 1, velocity v, c_s = 1): MC slopes; Hancock predictor q_face = q_i +- dq/2
- (dt/2dx) v dq; the scalar flux of llf (F = v(qL+qR)/2 - S(qR-qL)/2, S = v + c_s), hll (S_L = min(v-c_s, 0),
S_R = max(v+c_s, 0)) or an upwind solver (F = v qL); dt = C dx/(v + c_s).  Prints the top-hat's min and max
after one period with and without the predictor and, given acc_run.py's results.json, RAMSES's min.  The replica
reproduces RAMSES; without the predictor the extrema vanish: the predictor moves the face states with v while
llf/hll dissipate at v + c_s, and MC's steep slopes turn the mismatch into over/undershoots.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import acc_cases as A  # noqa: E402


def mc(q):
    dl = q - np.roll(q, 1)
    dr = np.roll(q, -1) - q
    dc = 0.5 * (dl + dr)
    lim = np.minimum(2 * abs(dl), 2 * abs(dr))
    lim[dl * dr <= 0] = 0
    return np.sign(dc) * np.minimum(lim, abs(dc))


def run(N, mach, C, flux, predictor=True):
    h, v, c = 1.0 / N, mach, 1.0
    q = A.exact_on_cells(((np.arange(N) + 0.5) * h)[:, None], h, (1, 0), 1)["tophat"]
    dt, t = C * h / (v + c), 0.0
    while t < 1 / v - 1e-14:                 # RAMSES also takes whole steps past tout
        dq = mc(q)
        pr = 0.5 * dt / h * v * dq if predictor else 0 * dq
        qL = q + 0.5 * dq - pr               # left state at i+1/2
        qR = np.roll(q - 0.5 * dq - pr, -1)  # right state at i+1/2
        if flux == "llf":
            S = v + c
            F = 0.5 * v * (qL + qR) - 0.5 * S * (qR - qL)
        elif flux == "hll":
            SL, SR = min(v - c, 0.0), max(v + c, 0.0)
            F = (SR * v * qL - SL * v * qR + SL * SR * (qR - qL)) / (SR - SL)
        else:
            F = v * qL
        q = q - dt / h * (F - np.roll(F, 1))
        t += dt
    return q.min(), q.max()


def main():
    R = {}
    if len(sys.argv) > 1:
        R = {r["case"]: r for r in json.load(open(sys.argv[1]))["results"] if r["status"] == "OK"}
    print("top-hat, MC, Mach 0.1, N = 64, one period")
    for flux in ("llf", "hll", "hllc"):
        for C in (0.2, 0.8):
            mn, mx = run(64, 0.1, C, flux)
            mn0, mx0 = run(64, 0.1, C, flux, predictor=False)
            r = R.get(f"h1d_s2_{flux}_M0.1_N64" + ("" if C == 0.8 else f"_C{C:g}"))
            ram = f" | RAMSES min {r['tophat']['qmin']:+.3e}" if r else ""
            print(f"  {flux:5s} C={C}: min {mn:+.3e} max-1 {mx - 1:+.3e} | without predictor: min {mn0:+.3e} "
                  f"max-1 {mx0 - 1:+.3e}{ram}")


if __name__ == "__main__":
    main()
