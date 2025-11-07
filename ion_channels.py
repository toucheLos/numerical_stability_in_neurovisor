# Minimal Hodgkin–Huxley-style channel equations (no integrator).
# Each channel has:
#   - rate/steady-time functions (alpha/beta or p_inf/tau)
#   - a current helper I(V, gates, params)
# Units:
#   V in volts; rates in 1/s; conductances in S/m^2; currents in A/m^2.

import numpy as np

# Helpers
def mV(V): 
    # Convert volts -> millivolts (vectorized)
    return 1.0e3 * V

# Default parameters
# Conductances
G_DEFAULT = dict(
    gk   = 5.0e1,            # S/m^2
    gna  = 50.0e1,           # S/m^2
    gl   = 1.0,   # S/m^2
    gcaH = 0.0001 * 1.0E4,   # S/m^2 (High-threshold Ca; "calcium_channel" in your C#)
    gM   = 7.5e-5 * 1.0e4,   # S/m^2 (Slow K / M-current; adjust as in your model)
    gT   = 4.0e-4 * 1.0e4,   # S/m^2 (Low-threshold Ca, T-type)
)

# Reversal potentials (V)
E_DEFAULT = dict(
    Ek  = -90.0e-3,
    Ena =  50.0e-3,
    El  = -70.0e-3,
    Eca = 120.0e-3,
)

# Optional voltage shift
V_SHIFTS = dict(vT=0, Vx=2.0)


# Potassium (KDR) — gate n (Pospischil 2008)
def k_rates(V, vT_mV=V_SHIFTS['vT']):
    # Returns alpha_n(V), beta_n(V) [1/s].
    Vm = mV(V)
    vT = float(vT_mV)
    a = 1.0e3 * (-0.032 * (Vm - vT - 15)) / (np.exp(-(Vm - vT - 15)/5) - 1) 
    b = 1.0e3 * 0.5 * np.exp(-(Vm - vT - 10.0)/40.0)
    return a, b

def k_current(V, n, gk=G_DEFAULT['gk'], Ek=E_DEFAULT['Ek']):
    # I_K = gk * n^4 * (V - Ek)
    return gk * (n**4) * (V - Ek)


# Sodium (Na)
def na_rates(V, vT_mV=V_SHIFTS['vT']):
    # Returns (alpha_m, beta_m, alpha_h, beta_h) [1/s].
    Vm = mV(V)
    vT = float(vT_mV)

    am = 1.0e3 * (-0.32 * (Vm - vT - 13.0)) / (np.exp(-(Vm - vT - 13.0) / 4.0) - 1)
    bm = 1.0e3 * (0.28 * (Vm - vT - 40.0)) / (np.exp((Vm - vT - 40.0) / 5.0) - 1)

    ah = 1.0e3 * 0.128 * np.exp(-(Vm - vT - 17.0) / 18.0)
    bh = 1.0e3 * 4.0 / (1.0 + np.exp(-(Vm - vT - 40.0)/5.0))

    return am, bm, ah, bh
 
def na_current(V, m, h, gna=G_DEFAULT['gna'], Ena=E_DEFAULT['Ena']):
    # I_Na = gna * m^3 * h * (V - Ena).
    return gna * (m**3) * h * (V - Ena)


# Leak
def leak_current(V, gl=G_DEFAULT['gl'], El=E_DEFAULT['El']):
    # I_Leak = gl * (V - El)
    return gl * (V - El)


# High-threshold Calcium (q, r)
def ca_high_rates(V):
    # Returns alpha_q, beta_q, alpha_r, beta_r [1/s].
    Vm = mV(V)

    # alpha_q = 1e3 * 0.055 * (27 - V[mV]) / (exp((-27 - V)/3.8) - 1)
    xq = (-27.0 - Vm)/3.8
    aq = 1.0e3 * 0.055 * np.divide((27.0 - Vm), np.expm1(xq),
                                   out=np.ones_like(Vm), where=(np.expm1(xq)!=0))

    bq = 1.0e3 * 0.94 * np.exp((-75.0 - Vm) / 17.0)

    ar = 1.0e3 * 0.000457 * np.exp((-13.0 - Vm) / 50.0)
    br = 1.0e3 * 0.0065   / (np.exp((-15.0 - Vm) / 28.0) + 1.0)

    return aq, bq, ar, br

def ca_high_current(V, q, r, gca=G_DEFAULT['gcaH'], Eca=E_DEFAULT['Eca']):
    # I_Ca(H) = gca * q^2 * r * (V - Eca)
    return gca * (q**2) * r * (V - Eca)


# Slow Potassium (M-current) — gate p
def m_p_inf(V):
    # p_inf(V) (unitless)
    Vm = mV(V)
    return 1.0 / (1.0 + np.exp(-(Vm + 35.0)/10.0))

def m_tau_p(V, tmax=4):
    # tau_p(V) in seconds
    Vm = mV(V)
    return tmax / (3.3 * np.exp((Vm + 35.0)/20.0) + np.exp(-(Vm + 35.0)/20.0))

def m_alpha_beta(V, tmax=4):
    # Return alpha_p, beta_p from p_inf/tau.
    pinf = m_p_inf(V)
    tau  = m_tau_p(V, tmax)
    a = np.divide(pinf, tau, out=np.zeros_like(pinf), where=(tau!=0))
    b = np.divide(1.0 - pinf, tau, out=np.zeros_like(pinf), where=(tau!=0))
    return a, b

def m_current(V, p, gM=G_DEFAULT['gM'], Ek=E_DEFAULT['Ek']):
    # I_M = gM * p * (V - Ek)
    return gM * p * (V - Ek)


# Low-threshold Calcium (T-type): s is instantaneous; u is dynamic
def t_s_inf(V, Vx_mV=V_SHIFTS['Vx']):
    # s_inf(V) (instantaneous)
    Vm = mV(V)
    Vx = float(Vx_mV)
    return 1.0 / (1.0 + np.exp(-(Vm + Vx + 57.0)/6.2))

def t_u_inf(V, Vx_mV=V_SHIFTS['Vx']):
    Vm = mV(V); Vx=float(Vx_mV)
    return 1.0 / (1.0 + np.exp((Vm + Vx + 81.0)/4.0))

def t_tau_u(V, Vx_mV=V_SHIFTS['Vx']):
    Vm = mV(V); Vx=float(Vx_mV)
    num = (30.8 + 211.4 + np.exp((Vm + Vx + 113.2)/5.0))
    den = 3.7 * (1.0 + np.exp((Vm + Vx + 84.0)/3.2))
    # out = np.full_like(Vm, np.inf)
    out = num / den
    return out

def t_alpha_beta_u(V, Vx_mV=V_SHIFTS['Vx']):
    ui  = t_u_inf(V, Vx_mV)
    tau = t_tau_u(V, Vx_mV)
    a = np.divide(ui, tau, out=np.zeros_like(ui), where=(tau!=0))
    b = np.divide(1.0 - ui, tau, out=np.zeros_like(ui), where=(tau!=0))
    return a, b

def t_current(V, u, gT=G_DEFAULT['gT'], Eca=E_DEFAULT['Eca'], Vx_mV=V_SHIFTS['Vx']):
    """I_T = gT * [s_inf(V)]^2 * u * (V - Eca).  (s is instantaneous)"""
    s = t_s_inf(V, Vx_mV)
    return gT * (s**2) * u * (V - Eca)


# Example usage
if __name__ == "__main__":
    # Test at a few voltages
    V = np.array([-0.070, -0.050, -0.030])  # V

    # K
    an, bn = k_rates(V)
    print("K: alpha_n =", an, "beta_n =", bn)

    # Na
    am, bm, ah, bh = na_rates(V)
    print("Na: alpha_m =", am, "beta_m =", bm, "alpha_h =", ah, "beta_h =", bh)

    # High-threshold Ca
    aq, bq, ar, br = ca_high_rates(V)
    print("Ca-H: alpha_q =", aq, "beta_q =", bq, "alpha_r =", ar, "beta_r =", br)

    # Slow K (M)
    pinf = m_p_inf(V); taup = m_tau_p(V)
    ap, bp = m_alpha_beta(V)
    print("M: p_inf =", pinf, "tau_p =", taup, "alpha_p =", ap, "beta_p =", bp)

    # Low-T Ca (T-type)
    s_inf = t_s_inf(V); ui = t_u_inf(V); tauu = t_tau_u(V)
    au, bu = t_alpha_beta_u(V)
    print("T: s_inf =", s_inf, "u_inf =", ui, "tau_u =", tauu, "alpha_u =", au, "beta_u =", bu)
