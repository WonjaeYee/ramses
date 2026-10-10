# sink

3D, 32^3, periodic: one individual-star sink (INDIVIDUAL_SINK_STARS, MIST main-sequence
winds, `accrete_sink` in pm/sink_particle.f90) of 40 Msun at the centre of uniform, fully
ionised, 0.1-solar gas.  The wind injects mass, momentum and (neutral) metals every step and
blows a bubble, so the test combines the injection source term with 3D transport of the
composition gradients it creates.  No chemistry, so the box mass of the neutral stage of
every element must grow exactly like the element (`neutral`).

    build.sh <tree> full_nodust 3 <exe> INDI_STAR=1
    python run_sink.py <exe> <rundir> [--label cmaoff] [--np 2]

PASS for the consistent scheme: `top`, `child_abs` to round-off and `neutral` < 1e-8.  With
ion slot = x*rho the injection adds the total injected mass to every element's neutral slot
(consistent with sum(ions)=rho, but the ion masses are wrong) and the hydro adds the
face-state errors.
