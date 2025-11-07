TITLE Original Potassium Channel (Minimal)

COMMENT
Implements a basic Hodgkin-Huxley style K+ channel with:
  gKbar = 5.0e1 S/m^2 => 0.05 (mS/cm2)
  ek    = -90 mV
Gating variable n with exponent 4, using alpha_n / beta_n from your snippet:
 alpha_n(V) = 1e3 * 0.032*(15 - V)/( exp((15-V)/5)-1 )
 beta_n(V)  = 1e3 * 0.5*exp((10 - V)/40)
ENDCOMMENT

NEURON {
    SUFFIX original_k
    USEION k READ ek WRITE ik
    RANGE gkbar, ik
    RANGE n
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
    gkbar = 0.05 (mS/cm2)   : from 5.0e1 S/m^2 => 0.05 S/cm^2
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

BREAKPOINT {
    SOLVE states METHOD cnexp
    ik = gkbar * n^4 * (v - ek)
}

DERIVATIVE states {
    rates(v)
    n' = alpha_n*(1 - n) - beta_n*n
}

PROCEDURE rates(v(mV)) {
    : Direct from your snippet, no expansions
    alpha_n = (0.032)*(15.0 - v)/( exp((15.0 - v)/5.0) - 1.0 ) * (1e3)
    beta_n  = (0.5)*exp((10.0 - v)/40.0) * (1e3)
}
