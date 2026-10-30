"""House plotting style for the CALIMA test figures (mirrors pyCALIMA's).

The theme comes from ``pycalima.plotting_style.use_calima_style`` (seaborn
``white``, serif fonts, LaTeX text when a usable installation exists) when
pyCALIMA is importable -- installed, or from its source tree -- and from an
equivalent local definition otherwise.  On top of it go the numeric conventions
of pyCALIMA's plotting code: 300 dpi, inward ticks on all four sides, 16 pt axis
labels, 14 pt ticks and legends, frameless legends, linewidth 2.

Usage, from a test's plot script:

    sys.path.append(str(Path(__file__).resolve().parents[3] / "utils/py/calima"))
    from calima_style import use_style, tidy, save, tex, FIGSIZE
    use_style()
    fig, ax = plt.subplots(figsize=FIGSIZE["single"])
    tidy(ax)
    save(fig, "myfigure")          # writes myfigure.png and .pdf, then closes it

Figures are never displayed: ``use_style`` forces the Agg backend.
"""
from __future__ import annotations

import sys
from pathlib import Path

PYCALIMA_SRC = Path("/Users/currodri/Documents/GitHub/pyCALIMA/src")

# figure sizes matching pyCALIMA's common choices
FIGSIZE = {"single": (7, 5), "square": (6, 5), "wide": (8, 5),
           "two": (8, 5), "row3": (12, 4.2), "row4": (16, 4.6)}

DPI = 300

# figures drawn at their printed size for the paper (openjournal: \columnwidth = 245.3 pt = 3.39 in,
# \textwidth = 513.1 pt, 10 pt Computer Modern body), in the paper's own font.  Use as
# ``with matplotlib.rc_context(PAPER_RC): fig = plt.figure(figsize=(COLUMN_WIDTH, h))``.
COLUMN_WIDTH = 3.4
TEXT_WIDTH = 7.1          # figure* spanning both columns
PAPER_RC = {"font.family": "serif", "font.serif": ["Computer Modern Roman"],
            "axes.labelsize": 10, "axes.titlesize": 10, "xtick.labelsize": 9, "ytick.labelsize": 9,
            "legend.fontsize": 8, "legend.handlelength": 1.6, "legend.labelspacing": 0.25,
            "axes.linewidth": 0.8, "lines.linewidth": 1.0, "lines.markersize": 3,
            "xtick.major.size": 3.5, "ytick.major.size": 3.5, "xtick.minor.size": 2.0, "ytick.minor.size": 2.0,
            "xtick.major.width": 0.7, "ytick.major.width": 0.7, "xtick.minor.width": 0.5, "ytick.minor.width": 0.5,
            "xtick.major.pad": 2.5, "ytick.major.pad": 2.5, "axes.labelpad": 2.5}

# pyCALIMA reaches for these; keep sequential data on one of them
CMAP_SEQ, CMAP_DIV = "viridis", "coolwarm"

# categorical cycle: Okabe & Ito (2008) / Wong (2011, Nat. Methods 8, 441), colour-blind safe,
# without its yellow (too faint on white).  "C0".."C5" in scripts resolve to these, in this order.
PALETTE = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9"]


def _theme(usetex: bool | None = None) -> None:
    """Apply pyCALIMA's theme, or an equivalent fallback."""
    try:
        from pycalima.plotting_style import use_calima_style
    except ModuleNotFoundError:
        if PYCALIMA_SRC.is_dir() and str(PYCALIMA_SRC) not in sys.path:
            sys.path.insert(0, str(PYCALIMA_SRC))
        try:
            from pycalima.plotting_style import use_calima_style
        except ModuleNotFoundError:
            use_calima_style = None
    if use_calima_style is not None:
        use_calima_style(usetex=usetex, force=True)
        return
    import matplotlib.pyplot as plt
    try:
        import seaborn as sns
        sns.set_theme(style="white")
    except ModuleNotFoundError:
        pass
    if usetex is None:
        import shutil
        usetex = all(shutil.which(exe) for exe in ("latex", "dvipng"))
    plt.rcParams.update({"text.usetex": bool(usetex),
                         "font.family": "serif",
                         "font.serif": ["DejaVu Serif", "Times New Roman",
                                        "Times", "serif"]})


def use_style(usetex: bool | None = None) -> None:
    """Set up matplotlib for CALIMA figures (Agg backend: nothing is displayed)."""
    import matplotlib as mpl
    mpl.use("Agg", force=True)
    _theme(usetex)
    mpl.rcParams.update({
        "figure.dpi": 120, "savefig.dpi": DPI,
        "figure.facecolor": "w", "savefig.facecolor": "w",
        "axes.labelsize": 16, "axes.titlesize": 16,
        "xtick.labelsize": 14, "ytick.labelsize": 14,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
        "xtick.bottom": True, "ytick.left": True,     # seaborn "white" turns these two off
        "xtick.minor.visible": True, "ytick.minor.visible": True,
        "legend.fontsize": 14, "legend.frameon": False,
        "lines.linewidth": 2.0, "axes.linewidth": 1.4,
        "axes.edgecolor": "k",
        "xtick.major.size": 6, "ytick.major.size": 6,
        "xtick.minor.size": 3.5, "ytick.minor.size": 3.5,
        "xtick.major.width": 1.2, "ytick.major.width": 1.2,
        "xtick.minor.width": 0.9, "ytick.minor.width": 0.9,
        "xtick.color": "k", "ytick.color": "k",
        "axes.prop_cycle": mpl.cycler(color=PALETTE),
        "image.cmap": CMAP_SEQ,
    })


def tidy(ax):
    """Inward major and minor ticks on all four spines, matching pyCALIMA's tick_params call."""
    ax.minorticks_on()
    ax.tick_params(which="both", axis="both", direction="in",
                   top=True, right=True, bottom=True, left=True)
    return ax


def panel_label(ax, label, loc=(0.03, 0.95), **kw):
    """Bold panel tag, e.g. ``panel_label(ax, "a")`` -> **(a)** in the upper-left corner."""
    return ax.text(*loc, rf"\textbf{{({label})}}" if _usetex() else f"({label})",
                   transform=ax.transAxes, ha="left", va="top", **{"fontsize": mpl_rc("axes.labelsize"), **kw})


def _usetex():
    return bool(mpl_rc("text.usetex"))


def mpl_rc(key):
    import matplotlib as mpl
    return mpl.rcParams[key]


def zoom_inset(ax, bounds, xlim, ylim, **kw):
    """An inset of *ax* at *bounds* (axes fraction: x0, y0, w, h) showing *xlim* x *ylim*.

    Plot into the returned axes with the same calls as into *ax*; the zoomed region is outlined
    on *ax* and joined to the inset. Tick labels are smaller (11 pt) to keep the inset light.
    """
    ins = ax.inset_axes(bounds, **kw)
    ins.set_xlim(*xlim)
    ins.set_ylim(*ylim)
    tidy(ins)
    ins.tick_params(which="both", labelsize=11)
    ax.indicate_inset_zoom(ins, edgecolor="0.3", lw=1.0, alpha=0.9)
    return ins


_TEX_SPECIAL = {"&": r"\&", "%": r"\%", "#": r"\#", "_": r"\_"}


def tex(text: str) -> str:
    """Escape LaTeX-special characters in figure *text* (outside math mode)."""
    return "".join(_TEX_SPECIAL.get(ch, ch) for ch in text)


def save(fig, stem: str, formats=("png", "pdf"), close: bool = True) -> None:
    """Save at 300 dpi as png and pdf, then close the figure.

    A trailing ``.png``/``.pdf`` on *stem* is dropped, so ``save(fig, "a.png")``
    writes ``a.png`` and ``a.pdf``.
    """
    for ext in (".png", ".pdf"):
        if stem.endswith(ext):
            stem = stem[: -len(ext)]
    for ext in formats:
        fig.savefig(f"{stem}.{ext}", format=ext, dpi=DPI, bbox_inches="tight")
    if close:
        import matplotlib.pyplot as plt
        plt.close(fig)
