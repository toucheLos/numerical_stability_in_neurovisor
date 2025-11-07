TITLE Original Sodium Channel (Minimal)

COMMENT
Implements a sodium channel with:
  gnabar = 60.0e1 S/m^2 => 0.60 (mS/cm2)
  ena = +50 mV
Gating variables m (exponent=3) and h (exponent=1).
From your snippet:

 alpha_m(V)= 1e3*(0.32)*(13 - V) / [exp((13-V)/4)-1]
 beta_m(V) = 1e3*(0.28)*(V - 40)/ [exp((V-40)/5)-1]
 alpha_h(V)= 1e3*(0.128)*exp((17 - V)/18)
 beta_h(V) = 1e3*4.0 / [ exp((40 - V)/5)+1 ]
ENDCOMMENT

NEURON {
    SUFFIX original_na
    USEION na READ ena WRITE ina
    RANGE gnabar, ina
    RANGE m, h
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
    gnabar = 0.60 (mS/cm2)   : from 60.0e1 => 0.60 S/cm^2
}

STATE {
    m
    h
}

ASSIGNED {
    v (mV)
    ena (mV)
    ina (mA/cm2)
    alpha_m (1/ms)
    beta_m  (1/ms)
    alpha_h (1/ms)
    beta_h  (1/ms)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ina = gnabar*(m^3)*h*(v - ena)
}

DERIVATIVE states {
    rates(v)
    m' = alpha_m*(1 - m) - beta_m*m
    h' = alpha_h*(1 - h) - beta_h*h
}

PROCEDURE rates(v(mV)) {
    alpha_m = (0.32)*(13.0 - v)/( exp((13.0 - v)/4.0) - 1.0 )*(1e3)
    beta_m  = (0.28)*(v - 40.0)/( exp((v - 40.0)/5.0) - 1.0 )*(1e3)

    alpha_h = (0.128)*exp((17.0 - v)/18.0)*(1e3)
    beta_h  = (4.0)/( exp((40.0 - v)/5.0)+1.0 )*(1e3)
}
