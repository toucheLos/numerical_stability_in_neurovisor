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
    RANGE gnabar, ina, vshift
    RANGE m, h
}

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (mS) = (millisiemens)
}

PARAMETER {
    gnabar = 0.5 (mS/cm2)
    vshift = 0 (mV) : positive shifts activation to more depolarized voltages
}

STATE {
    m h
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
    ina = gnabar*(m*m*m)*h*(v - ena)
}

DERIVATIVE states {
    rates(v)
    m' = alpha_m*(1 - m) - beta_m*m
    h' = alpha_h*(1 - h) - beta_h*h
}

PROCEDURE rates(v (mV)) {
    LOCAL vs
    vs = v - vshift

    : Rates are in 1/ms (NEURON uses ms). NO *1e3 factor here.
    alpha_m = 0.32 * (13.0 - vs) / (exp((13.0 - vs)/4.0) - 1.0)
    beta_m  = 0.28 * (vs - 40.0) / (exp((vs - 40.0)/5.0) - 1.0)

    alpha_h = 0.128 * exp((17.0 - vs)/18.0)
    beta_h  = 4.0 / (exp((40.0 - vs)/5.0) + 1.0)
}
