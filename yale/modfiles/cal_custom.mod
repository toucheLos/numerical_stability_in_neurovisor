TITLE Custom Calcium Channel (cal_custom)
COMMENT
    This mechanism implements a voltage-dependent calcium channel with two 
    gating variables (q and r) using Hodgkin–Huxley-style kinetics. The rate 
    functions are adapted from the sparse solver equations.
    
    Equations used:
      alpha_q(v) = 1e3*0.055*(27-v) / ( exp((-27-v)/3.8) - 1 )
      beta_q(v)  = 1e3*0.94*exp((-75-v)/17)
      alpha_r(v) = 1e3*0.000457*exp((-13-v)/50)
      beta_r(v)  = 1e3*0.0065 / ( exp((-15-v)/28) + 1 )
    
    Steady states and time constants:
      q_inf = alpha_q/(alpha_q+beta_q),   tau_q = 1/(alpha_q+beta_q)
      r_inf = alpha_r/(alpha_r+beta_r),   tau_r = 1/(alpha_r+beta_r)
    
    Current:
      ica = gcalbar * q^2 * r * (v - eca)
ENDCOMMENT

NEURON {
    SUFFIX cal_custom
    USEION ca READ eca WRITE ica
    RANGE gcalbar, eca, ica, qinf, tauq, rinf, taur
}

PARAMETER {
    gcalbar = 0.01 (S/cm2)   : maximum conductance for the Ca channel
    eca = 120 (mV)          : reversal potential for Ca
}

STATE {
    q
    r
}

ASSIGNED {
    v (mV)
    ica (mA/cm2)
    qinf
    tauq (ms)
    rinf
    taur (ms)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ica = gcalbar * (q^2) * r * (v - eca)
}

DERIVATIVE states {
    rate(v)
    q' = (qinf - q)/tauq
    r' = (rinf - r)/taur
}

PROCEDURE rate(v (mV)) {
    LOCAL a_q, b_q, a_r, b_r
    a_q = 1e3 * 0.055 * (27 - v) / ( exp((-27 - v)/3.8) - 1 )
    b_q = 1e3 * 0.94 * exp((-75 - v)/17)
    a_r = 1e3 * 0.000457 * exp((-13 - v)/50)
    b_r = 1e3 * 0.0065 / ( exp((-15 - v)/28) + 1 )
    
    qinf = a_q/(a_q + b_q)
    tauq = 1/(a_q + b_q)
    rinf = a_r/(a_r + b_r)
    taur = 1/(a_r + b_r)
}
