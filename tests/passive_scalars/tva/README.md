# tva

CALIMA sets only.  Static uniform gas with a smooth composition; the dust and PAH bins drift
with an imposed velocity (`use_w_drift_test`, w=0.3) relative to the mixture (TVA: rho is the
gas+dust density).  The drift moves dust mass between cells while rho does not change, so the
gas must move the opposite way: without that counter-flux on the gas slots (elements, CO,
ions) `top` grows like the divergence of the drift flux.  PASS for the consistent scheme:
sums and masses to round-off.
