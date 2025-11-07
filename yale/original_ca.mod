TITLE Original Calcium Channel (Minimal)

COMMENT
Matches "CalciumChannel" from your snippet:
  gca = 1.0e1 => 0.01 S/cm^2
  eca = +120 mV
  gating variables: q^2 * r

 alpha_q(V)= 1e3 * 0.055*(27 - v)/( exp((27 - v)/3.8)-1 )
 beta_q(V) = 1e3 * 0.94*exp((-(75 + v))/17)
 alpha_r(V)= 1e3 * 0.000457*exp((-(13+ v))/50)
 beta_r(V) = 1e3 * 0.0065 / [ exp((-(15+ v))/28)+1 ]
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
    gcabar = 0.01 (mS/cm2)  : from 1.0e1 => 0.01 S/cm^2
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

DERIVATIVE states {
    rates(v)
    q' = alpha_q*(1 - q) - beta_q*q
    r' = alpha_r*(1 - r) - beta_r*r
}

PROCEDURE rates(v (mV)) {
    alpha_q = 0.055*(27.0 - v)/(exp((27.0 - v)/3.8)-1.0)*(1e3)
    beta_q  = 0.94*exp((-(75.0 + v))/17.0)*(1e3)

    alpha_r = 0.000457*exp((-(13.0 + v))/50.0)*(1e3)
    beta_r  = 0.0065 /( exp((-(15.0 + v))/28.0)+1.0 )*(1e3)
}
