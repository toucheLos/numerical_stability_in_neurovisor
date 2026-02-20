TITLE Original Potassium Channel (Minimal, with vshift)

COMMENT
Implements a basic Hodgkin-Huxley style K+ channel with adjustable voltage shift.

Original rate equations:
  alpha_n(V) = 1e3 * 0.032 * (15 - V) / (exp((15 - V)/5) - 1)
  beta_n(V)  = 1e3 * 0.5 * exp((10 - V)/40)

Adding vshift -> use vs = v - vshift

A positive vshift shifts activation rightward (requires more depolarization).
ENDCOMMENT

NEURON {
    SUFFIX original_k
    USEION k READ ek WRITE ik
    RANGE gkbar, ik, n, vshift
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
    gkbar = 0.5 (mS/cm2)
    vshift = 0 (mV)   : positive shifts activation to higher voltages
}

STATE {
    n
}

ASSIGNED {
    v   (mV)
    ek  (mV)
    ik  (mA/cm2)
    alpha_n (1/ms)
    beta_n  (1/ms)
}

INITIAL {
    rates(v)
    n = alpha_n / (alpha_n + beta_n)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ik = gkbar * n^4 * (v - ek)
}

DERIVATIVE states {
    rates(v)
    n' = alpha_n * (1 - n) - beta_n * n
}

PROCEDURE rates(v (mV)) {
    LOCAL vs
    vs = v - vshift  : shifted voltage

    alpha_n = (0.032)*(15.0 - vs) / ( exp((15.0 - vs)/5.0) - 1.0 )
    beta_n  = (0.5)*exp((10.0 - vs)/40.0)
}
