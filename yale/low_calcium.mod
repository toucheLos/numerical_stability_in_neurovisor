TITLE Original Low-Threshold Calcium Channel (T-type, Minimal)



COMMENT

Implements a T-type calcium channel with:

  gTbar = 4.0e-4 * 1.0e4 => 0.004 (mS/cm2)

  eca = +120 mV

  instantaneous gate s_inf(V), dynamic gate u

  I_T = gTbar * [s_inf(V)]^2 * u * (v - eca)

ENDCOMMENT



NEURON {

    SUFFIX original_caT

    USEION ca READ eca WRITE ica

    RANGE gTbar, ica

    RANGE u

}



UNITS {

    (mV)     = (millivolt)

    (mA)     = (milliamp)

    (mS)     = (millisiemens)

    (S)      = (siemens)

    (pA)     = (picoamp)

    (um)     = (micron)

    (mA/cm2) = (milliamp / square centimeter)

    (mS/cm2) = (millisiemens / square centimeter)

}



PARAMETER {

    gTbar = 0.004 (mS/cm2)

    Vx = 2 (mV)           : small voltage shift

}



STATE {

    u

}



ASSIGNED {

    v (mV)

    eca (mV)

    ica (mA/cm2)

    s_inf

    u_inf

    tau_u (ms)

}



BREAKPOINT {

    SOLVE states METHOD cnexp

    s_inf = 1.0 / (1.0 + exp(-(v + Vx + 57.0)/6.2))

    ica = gTbar * (s_inf*s_inf) * u * (v - eca)

}



DERIVATIVE states {

    rates(v)

    u' = (u_inf - u) / tau_u

}



PROCEDURE rates(v(mV)) {

    u_inf = 1.0 / (1.0 + exp((v + Vx + 81.0)/4.0))

    tau_u = (30.8 + 211.4 + exp((v + Vx + 113.2)/5.0)) / (3.7 * (1.0 + exp((v + Vx + 84.0)/3.2)))

}

