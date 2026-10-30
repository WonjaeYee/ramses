# RT15 levitating-atmosphere curves, extracted from the vector figures

Rosdahl & Teyssier (2015, MNRAS 449, 4380; arXiv:1411.6440), sec. 3.7. The four figures in `pdf/` are
the vector PDFs from the arXiv source tarball (`figs/davis_*.pdf`). They are RT15's figures and are not
committed: copy them into `pdf/` from the tarball to re-run the extraction. `extract_rt15.py` reads their
stroked paths and tick labels with `mutool` and writes the curves to `csv/`. No pixel digitising is
involved, so the values are as exact as the PDF coordinates allow. The axis fits close to
1.5e-3 in the y-axis units, 1e-6 dex on the log axis and 0.012 t_* in time.

| RT15 figure | file | panels | curves |
|---|---|---|---|
| fig. 9  | `davis_tau`       | `tauc` (log)                      | `{trap,notrap}_{max,massweighted}` |
| fig. 11 | `davis_eddington` | `fEV`, `tauV`, `tauF_over_tauV`   | `trap`, `notrap`, `FLD`, `VET` |
| fig. 12 | `davis_sigmas`    | `vy` (mean v_y / c_*), `sigma`     | `vy`: `trap`, `notrap`, `FLD`, `VET`; `sigma`: the same with `_sigma`, `_sigmax`, `_sigmay` |
| fig. 13 | `davis_edd_contr` | `fEV`                             | `total`, `stream` (kappa rho F_y / c), `trapped` ((1/3) dE_t/dy) |

`trap` is RT15's main run (M1 with trapping, dark green). `notrap` is M1 without trapping (light
green), which ran at c_red throughout in a box half as wide and half as tall (512 h_*, 1024^2 cells, the same
resolution). `FLD` and `VET` are Athena, from
Davis et al. (2014). The VET curves end at 155 t_*, and fig. 13 ends at 182.6 t_*.

Each CSV holds `t_over_tstar,value,clipped`. `clipped=1` marks points on a frame edge, where the
curve leaves the plotted range: the t = 0 spikes of `trap` in f_E,V and tau_V, and of `trapped` in
fig. 13. The values there are the frame limit, not the data.

Note: the fig. 11 caption says "with and without trapping (light green and darker green curves
respectively)". That ordering contradicts both its legend and the fig. 9 caption. The legend is
followed here (dark = with trapping).

```bash
python3 extract_rt15.py     # needs mutool (MuPDF)
```
