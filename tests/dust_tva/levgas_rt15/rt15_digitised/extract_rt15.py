"""Extract the levitating-atmosphere curves of Rosdahl & Teyssier 2015 (RT15, MNRAS 449, 4380) sec. 3.7
straight from the vector figures in their arXiv source (arXiv:1411.6440, figs/davis_*.pdf).

The figures are vector PDFs, so this is not pixel digitising: `mutool draw -F trace` gives every stroked
path in page coordinates, and `mutool draw -F stext` every tick label with its position. Each panel is
calibrated by snapping its tick labels to the major tick marks and fitting the axis (linear, or log10
for fig. 9); the fit residual is printed and is at the level of the PDF coordinate rounding.

Curves are identified by RT15's colour and dash (captions of figs. 9 and 11-13 and the legends):
  dark green  .333 .420 .184   M1 with trapping (the main run)
  light green .561 .736 .561   M1 without trapping
  red         .697 .133 .133   FLD (Athena, Davis et al. 2014)       [edd_contr: streaming F_y term]
  blue        .255 .412 .883   VET (Athena, Davis et al. 2014)       [edd_contr: trapped dE_t/dy term]
  grey        .439 .502 .564   edd_contr: total f_E,V
Two-point paths are the legend samples (dropped) or the vertical rise at t = 0 clipped by the frame
(kept, flagged clipped=1). Points on the frame edges are clipped values, also flagged.

Writes one CSV per curve, columns t_over_tstar,value,clipped, named <figure>_<panel>_<curve>.csv.
usage: python3 extract_rt15.py          (reads the PDFs in ./pdf, writes ./csv)
"""
import os
import re
import subprocess
import xml.etree.ElementTree as ET

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
COL = {".333008 .419922 .184082": "trap", ".560547 .736328 .560547": "notrap",
       ".697266 .133057 .133057": "FLD", ".254883 .412109 .882812": "VET", ".438965 .501953 .564453": "total"}
DASH = {None: "solid", "56.693 56.693": "dashed", "11.3386 28.3465": "dotted",
        "56.693 28.3465 14.1733 28.3465": "dashdot"}
# per figure: panel names top to bottom, which panels are log10, and per-figure renaming of the curves
FIGS = {
    "davis_eddington": dict(panels=["fEV", "tauV", "tauF_over_tauV"]),       # RT15 fig. 11
    "davis_tau": dict(panels=["tauc"], log=["tauc"],                          # RT15 fig. 9
                      style={"solid": "max", "dashed": "massweighted"}),
    "davis_sigmas": dict(panels=["vy", "sigma"],                              # RT15 fig. 12
                         style={"solid": "sigma", "dashed": "sigmax", "dotted": "sigmay"},
                         style_panels=["sigma"]),
    "davis_edd_contr": dict(panels=["fEV"],                                   # RT15 fig. 13
                            rename={"FLD": "stream", "VET": "trapped", "total": "total"}),
}


def mutool(fmt, pdf):
    return subprocess.run(["mutool", "draw", "-F", fmt, "-o", "-", pdf], capture_output=True, text=True,
                          check=True).stdout


def strokes(pdf):
    out = []
    for el in ET.fromstring(mutool("trace", pdf)).iter("stroke_path"):
        a, b, c, d, e, f = (float(v) for v in el.get("transform").split())
        pts = []
        for p in el:
            k = ("x3", "y3") if p.tag == "curveto" else ("x", "y")
            if p.get(k[0]) is not None:
                x, y = float(p.get(k[0])), float(p.get(k[1]))
                pts.append((a * x + c * y + e, b * x + d * y + f))
        if len(pts) >= 2:
            out.append(dict(color=el.get("color"), dash=el.get("dash"), pts=np.array(pts)))
    return out


def labels(pdf):
    """Numeric tick labels: (value, x_centre, y_centre, bbox)."""
    out = []
    for line in ET.fromstring(mutool("stext", pdf)).iter("line"):
        s = "".join(ch.get("c") for ch in line.iter("char")).replace("−", "-")
        if re.fullmatch(r"-?\d+(\.\d+)?", s):
            x0, y0, x1, y1 = (float(v) for v in line.get("bbox").split())
            out.append((float(s), 0.5 * (x0 + x1), 0.5 * (y0 + y1), (x0, y0, x1, y1)))
    return out


def calibrate(pairs, log):
    """Fit value = A*pix + B (or log10); returns the function and the max residual in axis units."""
    p, v = np.array(pairs).T
    vv = np.log10(v) if log else v
    A, B = np.polyfit(p, vv, 1)
    res = np.max(np.abs(A * p + B - vv))
    return (lambda q: 10 ** (A * q + B)) if log else (lambda q: A * q + B), res


def extract(name, spec):
    pdf = os.path.join(HERE, "pdf", name + ".pdf")
    S = strokes(pdf)
    black = [s for s in S if s["color"] == "0 0 0"]
    seg = [(s["pts"][0], s["pts"][-1]) for s in black if len(s["pts"]) == 2]
    L = lambda a, b: np.hypot(*(b - a))
    # frames: the long vertical black lines are the left/right edges; one (top, bottom) per panel
    vert = sorted({(round(min(a[1], b[1]), 2), round(max(a[1], b[1]), 2)) for a, b in seg
                   if L(a, b) > 50 and abs(a[0] - b[0]) < 0.01})
    xl = min(min(a[0], b[0]) for a, b in seg if L(a, b) > 50)
    xr = max(max(a[0], b[0]) for a, b in seg if L(a, b) > 50)
    assert len(vert) == len(spec["panels"]), (name, vert)
    # major ticks only (the longest short segments): the labels sit ~1 pt off their tick, and IDL nudges
    # the labels at the frame ends inwards, so snapping to any tick would pick up minor ones
    def major(vertical):
        t = [(a, L(a, b)) for a, b in seg if L(a, b) < 10 and abs(a[0 if vertical else 1] - b[0 if vertical else 1]) < 0.01]
        lmax = max(l for _, l in t)
        return [a for a, l in t if l > 0.75 * lmax]
    ticks_v, ticks_h = major(True), major(False)                                     # x ticks, y ticks
    lab = labels(pdf)
    ybot = vert[-1][1]
    # x axis: labels below the bottom frame, snapped to the x tick (on the bottom edge) nearest in x
    xp = []
    for v, cx, cy, bb in lab:
        if bb[1] > ybot:
            t = min((tk for tk in ticks_v if abs(tk[1] - ybot) < 0.05), key=lambda tk: abs(tk[0] - cx))
            xp.append((t[0], v))
    fx, rx = calibrate(xp, False)
    print(f"{name}: x axis from {len(xp)} labels, max residual {rx:.2e} t_*")
    out = {}
    for (top, bot), pan in zip(vert, spec["panels"]):
        log = pan in spec.get("log", [])
        yp = []
        for v, cx, cy, bb in lab:
            if bb[2] < xl and top - 2 < cy < bot + 2:
                t = min((tk for tk in ticks_h if abs(tk[0] - xl) < 0.05), key=lambda tk: abs(tk[1] - cy))
                if abs(t[1] - cy) < 4.0:
                    yp.append((t[1], v))
        fy, ry = calibrate(yp, log)
        print(f"  panel {pan:15s} y from {len(yp)} labels {sorted(v for _, v in yp)}, max residual {ry:.2e}"
              + (" dex" if log else ""))
        groups = {}
        for s in S:
            if s["color"] not in COL:
                continue
            P = s["pts"]
            ym = np.median(P[:, 1])
            if not (top - 0.05 <= ym <= bot + 0.05):
                continue
            if len(P) == 2 and abs(P[0, 0] - P[1, 0]) > 1.0:
                continue                                           # horizontal 2-point path: legend sample
            run = COL[s["color"]]
            run = spec.get("rename", {}).get(run, run)
            sty = DASH.get(s["dash"], s["dash"])
            key = f"{run}_{spec['style'][sty]}" if pan in spec.get("style_panels", spec["panels"]) and "style" in spec else run
            groups.setdefault(key, []).append(P)
        for key, parts in groups.items():
            P = np.concatenate(sorted(parts, key=lambda q: q[:, 0].min()))
            # the vertical clipped rise at t=0 is a 2-point path at the left edge: keep, ordered by x then y
            t = fx(P[:, 0])
            y = fy(P[:, 1])
            clipped = (np.abs(P[:, 1] - top) < 0.05) | (np.abs(P[:, 1] - bot) < 0.05)
            fn = os.path.join(HERE, "csv", f"{name}_{pan}_{key}.csv")
            np.savetxt(fn, np.c_[t, y, clipped.astype(int)], fmt=["%.5f", "%.6g", "%d"], delimiter=",",
                       header=f"RT15 {name}.pdf, panel {pan}, curve {key}; extracted from the vector paths\n"
                              "t_over_tstar,value,clipped", comments="# ")
            out[(pan, key)] = (t, y, clipped)
            print(f"    {key:22s} {len(t):4d} pts  t [{t.min():6.2f},{t.max():6.2f}]  "
                  f"value [{y.min():.4g},{y.max():.4g}]  clipped {clipped.sum()}")
    return out


if __name__ == "__main__":
    os.makedirs(os.path.join(HERE, "csv"), exist_ok=True)
    for name, spec in FIGS.items():
        extract(name, spec)
