"""Pass/fail of an advection-accuracy run: the scheme's guarantees, and the errors against a reference run.

    acc_compare.py <results.json> [<reference results.json>] [--rtol 1e-4]

Checks (one PASS/FAIL line each, then RESULT: PASS|FAIL; writes compare.json next to results.json):
  complete     every case ran to the end
  mass         box mass of every scalar conserved to 1e-12 (flux form)
  first-order  slope_type 0, 1-D hydro: the sine amplitude equals |G|^n of the first-order scheme with the
               Riemann solver's own scalar flux, to 1e-9 (tests the harness and the solver's scalar flux)
  monotone     no new extrema (> 1e-12) with an upwind scalar flux (hydro hllc/exact/acoustic, MHD hlld/roe)
               and a TVD slope (0, 1, 2, 3, 7, 8): MUSCL-Hancock with an upwind flux is TVD for nu <= 1
  reference    L1, Linf, new extrema, sine amplitude and Gaussian D_num of every case and scalar within
               rtol (relative) + 1e-10 (absolute) of the reference run; cases missing from either side are listed
"""
import argparse
import json
import os

PROFILES = ("tophat", "gauss", "sine")
UPWIND = {"hydro": ("hllc", "exact", "acoustic"), "mhd": ("hlld", "roe")}
TVD_SLOPES = (0, 1, 2, 3, 7, 8)
METRICS = (("L1", "L1"), ("Linf", "Linf"), ("over", "over"), ("under", "under"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results")
    ap.add_argument("ref", nargs="?")
    ap.add_argument("--rtol", type=float, default=1e-4)
    a = ap.parse_args()
    R = json.load(open(a.results))["results"]
    ok = [r for r in R if r["status"] == "OK"]
    lines, bad = [], False

    def say(passed, text):
        nonlocal bad
        bad |= not passed
        lines.append(("PASS: " if passed else "FAIL: ") + text)

    miss = [r["case"] for r in R if r["status"] != "OK"]
    say(not miss, f"complete: {len(ok)}/{len(R)} cases" + (f"; failed: {', '.join(miss[:6])}" if miss else ""))
    if ok:
        w = max(ok, key=lambda r: max(r[p]["mass"] for p in PROFILES))
        m = max(w[p]["mass"] for p in PROFILES)
        say(m <= 1e-12, f"mass: worst |dM/M| {m:.1e} ({w['case']}, tolerance 1e-12)")
    th = [r for r in ok if "amp_theory" in r["sine"]]
    if th:
        w = max(th, key=lambda r: abs(r["sine"]["amp"] - r["sine"]["amp_theory"]))
        d = abs(w["sine"]["amp"] - w["sine"]["amp_theory"])
        say(d <= 1e-9, f"first-order: {len(th)} runs match |G|^n, worst |amp - theory| {d:.1e} ({w['case']})")
    mono = [r for r in ok if r["riemann"] in UPWIND[r["solver"]] and r["slope"] in TVD_SLOPES]
    if mono:
        ext = lambda r: max(max(r[p]["over"], r[p]["under"]) for p in PROFILES)   # noqa: E731
        w = max(mono, key=ext)
        say(ext(w) <= 1e-12, f"monotone: {len(mono)} upwind-flux TVD runs, worst new extremum {ext(w):.1e} "
                             f"({w['case']})")
    out = dict(checks=lines)
    if a.ref:
        ref = {r["case"]: r for r in json.load(open(a.ref))["results"] if r["status"] == "OK"}
        worst, n, fails = (0.0, ""), 0, []
        for r in ok:
            x = ref.get(r["case"])
            if x is None:
                continue
            n += 1
            for p in PROFILES:
                pairs = [(k, r[p][k], x[p][k]) for k, _ in METRICS]
                if p == "sine":
                    pairs.append(("amp", r[p]["amp"], x[p]["amp"]))
                if p == "gauss":
                    pairs.append(("D_num", r[p]["D_num"], x[p]["D_num"]))
                for k, v, v0 in pairs:
                    dv = abs(v - v0)
                    rel = dv / max(abs(v0), 1e-300)
                    if dv > a.rtol * abs(v0) + 1e-10:
                        fails.append(f"{r['case']} {p} {k} {v:.4e} (ref {v0:.4e})")
                    if dv > 1e-10 and rel > worst[0]:
                        worst = (rel, f"{r['case']} {p} {k}")
        only_new = sorted({r["case"] for r in ok} - set(ref))
        only_ref = sorted(set(ref) - {r["case"] for r in ok})
        txt = (f"reference: {n} cases, {len(fails)} metrics beyond rtol {a.rtol:g}; largest relative change "
               f"{worst[0]:.1e}" + (f" ({worst[1]})" if worst[1] else ""))
        if fails:
            txt += "; e.g. " + "; ".join(fails[:4])
        if only_new or only_ref:
            txt += f"; not compared: {len(only_new)} new, {len(only_ref)} only in the reference"
        say(not fails, txt)
        out.update(ref_fails=fails, ref_worst=worst, only_new=only_new, only_ref=only_ref)
    lines.append("RESULT: " + ("FAIL" if bad else "PASS"))
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(a.results)), "compare.json"), "w"), indent=1)
    print("\n".join(lines))
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
