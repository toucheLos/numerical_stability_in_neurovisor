TITLE Original Calcium Channel (Minimal)

COMMENT
Matches "CalciumChannel" from your snippet:
  gca = 1.0e1 => 0.01 S/cm^2
  eca = +120 mV
  gating variables: q^2 * r

alpha_q(V) = 0.055*(27 - V) / (exp((27 - V)/3.8) - 1)
beta_q(V)  = 0.94*exp(-(75 + V)/17)
alpha_r(V) = 0.000457*exp(-(13 + V)/50)
beta_r(V)  = 0.0065 / (exp(-(15 + V)/28) + 1)
ENDCOMMENT

NEURON {
    SUFFIX original_ca
    USEION ca READ eca WRITE ica
    RANGE gcabar, ica
    RANGE q, r
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
    gcabar = 0.01 (mS/cm2)
}

STATE {
    q
    r
}

ASSIGNED {
    v (mV)
    eca (mV)
    ica (mA/cm2)
    alpha_q (1/ms)
    beta_q  (1/ms)
    alpha_r (1/ms)
    beta_r  (1/ms)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ica = gcabar * (q*q) * r * (v - eca)
}

INITIAL {
    rates(v)
    q = alpha_q / (alpha_q + beta_q)
    r = alpha_r / (alpha_r + beta_r)
}

DERIVATIVE states {
    rates(v)
    q' = alpha_q*(1 - q) - beta_q*q
    r' = alpha_r*(1 - r) - beta_r*r
}

PROCEDURE rates(v (mV)) {
    alpha_q = 0.055*(-27.0 - v)/(exp((-27.0 - v)/3.8) - 1.0)
    beta_q  = 0.94 * exp(-(75.0 + v)/17.0)

    alpha_r = 0.000457 * exp(-(13.0 + v)/50.0)
    beta_r  = 0.0065 / ( exp(-(15.0 + v)/28.0) + 1.0 )
}
