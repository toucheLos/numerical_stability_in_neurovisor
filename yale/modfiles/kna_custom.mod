COMMENT
  kna_custom.mod
  A Hodgkin-Huxley style channel containing:
    - Potassium gate n with alpha_n, beta_n from PDE code
    - Sodium gates m, h with alpha_m, beta_m, alpha_h, beta_h from PDE code
  No leak current included. Reversal potentials are read from NEURON's ek, ena.
  Units are consistent with NEURON (mV, ms, mA/cm2, S/cm2).
ENDCOMMENT

NEURON {
  SUFFIX kna_custom
  USEION k  READ ek  WRITE ik
  USEION na READ ena WRITE ina
  RANGE gkbar, gk, gna, gnabar
  RANGE n, m, h
}

PARAMETER {
  gkbar  = 0.005 (S/cm2)  : K conductance (convert from PDE  50 S/m2 => 0.005 S/cm2)
  gnabar = 0.06  (S/cm2)  : Na conductance (convert from PDE 600 S/m2 => 0.06 S/cm2)
}

STATE {
  n
  m
  h
}

ASSIGNED {
  v   (mV)
  ek  (mV)
  ena (mV)
  ik  (mA/cm2)
  ina (mA/cm2)
  alpha_n (1/ms)
  beta_n  (1/ms)
  alpha_m (1/ms)
  beta_m  (1/ms)
  alpha_h (1/ms)
  beta_h  (1/ms)
  gk  (S/cm2)
  gna (S/cm2)
}

BREAKPOINT {
  SOLVE states METHOD cnexp
  gk  = gkbar  * n^4
  gna = gnabar * m^3 * h
  ik  = gk  * (v - ek)
  ina = gna * (v - ena)
}

DERIVATIVE states {
  rates(v)
  n' = alpha_n*(1 - n) - beta_n*n
  m' = alpha_m*(1 - m) - beta_m*m
  h' = alpha_h*(1 - h) - beta_h*h
}

PROCEDURE rates(v (mV)) {
  LOCAL vm  : v in mV
  vm = v

  : alpha_n = 1.0e3 * 0.032*(15 - Vin) / ( exp((15 - Vin)/5) - 1 )
  alpha_n = 1(e3)*0.032 * (15 - vm)/( exp((15 - vm)/5) - 1 )
  beta_n  = 1(e3)*0.5   * exp((10 - vm)/40)

  : alpha_m = 1.0e3 * 0.32*(13 - Vin) / ( exp((13 - Vin)/4) - 1 )
  alpha_m = 1(e3)*0.32 * (13 - vm)/( exp((13 - vm)/4) - 1 )
  beta_m  = 1(e3)*0.28 * (vm - 40)/( exp((vm - 40)/5) - 1 )

  : alpha_h = 1.0e3 * 0.128 * exp((17 - Vin)/18)
  alpha_h = 1(e3)*0.128 * exp( (17 - vm)/18 )
  beta_h  = 1(e3)*4.0   / (exp((40 - vm)/5) + 1 )
}