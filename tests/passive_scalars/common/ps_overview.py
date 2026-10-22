"""Overview figure of the test matrix: worst error per (test, species set) for each label.

    python ps_overview.py RUNS OUT.png LABEL [LABEL ...]

One panel per metric (sum(top)/rho, ion sums, ion masses); dots per label on a log axis,
round-off tolerance as a dashed line.
"""
import collections
import glob
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

COLORS = ["#eb6834", "#2a78d6", "#1baf7a", "#eda100"]   # categorical slots, fixed order
METRICS = [("top", "|sum(top)/rho - 1|", 1e-12), ("child_abs", "ion sum error / element peak", 1e-12),
           ("ionmass", "ion mass drift (x rho_e)", 1e-12)]


def worst(runs, label):
    w = collections.defaultdict(dict)
    for f in glob.glob(os.path.join(runs, label, "*", "*", "check.json")):
        r = json.load(open(f))
        c = r["case"]
        key = (c["test"], c["set"])
        for m, _, _ in METRICS:
            if m == "ionmass" and (c.get("chem") or c.get("source")):
                continue      # chemistry or injection change the ion masses: no conservation check
            if m in r["values"]:
                w[key][m] = max(w[key].get(m, 0.0), r["values"][m])
    return w


def main():
    runs, out, labels = sys.argv[1], sys.argv[2], sys.argv[3:]
    data = {lab: worst(runs, lab) for lab in labels}
    keys = sorted(set().union(*[set(d) for d in data.values()]))
    fig, axes = plt.subplots(1, len(METRICS), figsize=(14, 0.32 * len(keys) + 1.6), sharey=True)
    y = np.arange(len(keys))
    for ax, (m, title, tol) in zip(axes, METRICS):
        for j, lab in enumerate(labels):
            v = [max(data[lab].get(k, {}).get(m, np.nan), 1e-17) for k in keys]
            ax.scatter(v, y + (j - (len(labels) - 1) / 2) * 0.22, s=36, color=COLORS[j % 4], label=lab,
                       edgecolor="white", linewidth=1.5, zorder=3)
        ax.axvline(tol, color="#888888", ls="--", lw=1, zorder=1)
        ax.set_xscale("log"); ax.set_xlim(1e-17, 1e3)
        ax.set_title(title, fontsize=10)
        ax.grid(axis="x", color="#e6e6e6", lw=0.6); ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[0].set_yticks(y); axes[0].set_yticklabels([f"{t} / {s}" for t, s in keys], fontsize=8)
    axes[0].invert_yaxis()
    axes[-1].legend(loc="lower right", fontsize=8, frameon=False)
    fig.suptitle("Worst error per (test / species set); dashed: round-off tolerance", fontsize=10)
    fig.tight_layout()
    fig.savefig(out, dpi=110)


if __name__ == "__main__":
    main()
