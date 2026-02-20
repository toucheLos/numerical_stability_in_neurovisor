TITLE Original Slow Potassium (M-current, Minimal)

COMMENT
Implements a slow non-inactivating K+ (M-type) current with:
  gMbar = 7.5e-5 * 1.0e4 => 0.00075 (mS/cm2)
  ek = -90 mV
  I_M = gMbar * p * (v - ek)
  p follows first-order kinetics with:
    p_inf = 1 / (1 + exp(-(v + 35)/10))
    tau_p = tmax / (3.3*exp((v + 35)/20) + exp(-(v + 35)/20))
ENDCOMMENT

NEURON {
    SUFFIX original_kM
    USEION k READ ek WRITE ik
    RANGE gMbar, ik
    RANGE p
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
    gMbar = 0.00075 (mS/cm2)
    tmax = 4000 (ms)
}

STATE {
    p
}

ASSIGNED {
    v (mV)
    ek (mV)
    ik (mA/cm2)
    p_inf
    tau_p (ms)
}

COMMENT
INITIAL {
    rates(v)
    p = p_inf
}
ENDCOMMENT

BREAKPOINT {
    SOLVE states METHOD cnexp
    ik = gMbar * p * (v - ek)
}

DERIVATIVE states {
    rates(v)
    p' = (p_inf - p) / tau_p
}

PROCEDURE rates(v(mV)) {
    p_inf = 1.0 / (1.0 + exp(-(v + 35.0)/10.0))
    tau_p = tmax / (3.3*exp((v + 35.0)/20.0) + exp(-(v + 35.0)/20.0))
}
