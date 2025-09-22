# sbdf2_solver.py
# Minimal single-compartment IMEX–SBDF2 driver that calls channel functions from ion_channels.py
# Mirrors NeuroVISOR structure in tiny form, with explicit references to the original function names.

import numpy as np
import ion_channels as ch

# Channel toggles
CHANNELS = {
    "K":   True,   # KDR (n)
    "Na":  True,   # Sodium (m,h)
    "Leak": True,  # Leak
    "CaH": False,  # High-threshold Ca (q,r)
    "M":   False,  # Slow K (p)
    "T":   True,   # Low-threshold Ca (u, s instantaneous)
}

# Constants
C_m = 1.0e-2     # F/m^2 (1 µF/cm²)
Iapp = 0.0       # A/m^2
E = ch.E_DEFAULT
G = ch.G_DEFAULT

# NeuroVISOR: stateexplicitSBDF2(...)
# S^{n+1} = 4/3 S^n - 1/3 S^{n-1} + dt * [ 4/3 f^n - 2/3 f^{n-1} ],  f = a(V)(1-S) - b(V)S
def f_gate(S, a, b):  # RHS of a single gate ODE
    return a * (1.0 - S) - b * S

def update_gate_sbdf2(Sn, Snm1, Vn, Vnm1, alpha_fn, beta_fn):
    a_n, b_n   = alpha_fn(Vn),   beta_fn(Vn)
    a_m1, b_m1 = alpha_fn(Vnm1), beta_fn(Vnm1)
    fn   = f_gate(Sn,   a_n,  b_n)
    fm1  = f_gate(Snm1, a_m1, b_m1)
    Snp1 = (4.0/3.0)*Sn - (1.0/3.0)*Snm1 + dt*((4.0/3.0)*fn - (2.0/3.0)*fm1)
    return np.clip(Snp1, 0.0, 1.0)

# Special case: M-current uses (p_inf, tau_p)
def update_p_sbdf2(pn, pnm1, Vn, Vnm1):
    fn   = (ch.m_p_inf(Vn)   - pn)   / ch.m_tau_p(Vn)
    fm1  = (ch.m_p_inf(Vnm1) - pnm1) / ch.m_tau_p(Vnm1)
    pnp1 = (4.0/3.0)*pn - (1.0/3.0)*pnm1 + dt*((4.0/3.0)*fn - (2.0/3.0)*fm1)
    return np.clip(pnp1, 0.0, 1.0)

# NeuroVISOR: reactF(...)
# Builds ionic “reaction” in linearized form: sum g_eff(V,gates)*(V - E) = G*V - GE
# We return (G, GE) at a given time level.
def reaction_terms(V, gates):
    Gsum  = np.array([0.0])
    GEsum = np.array([0.0])
    if CHANNELS["K"]:
        n = gates["n"]; gprod = G["gk"] * n**4
        Gsum += gprod; GEsum += gprod * E["Ek"]
    if CHANNELS["Na"]:
        m, h = gates["m"], gates["h"]; gprod = G["gna"] * (m**3) * h
        Gsum += gprod; GEsum += gprod * E["Ena"]
    if CHANNELS["Leak"]:
        gprod = G["gl"]
        Gsum += gprod; GEsum += gprod * E["El"]
    if CHANNELS["CaH"]:
        q, r = gates["q"], gates["r"]; gprod = G["gcaH"] * (q**2) * r
        Gsum += gprod; GEsum += gprod * E["Eca"]
    if CHANNELS["M"]:
        p = gates["p"]; gprod = G["gM"] * p
        Gsum += gprod; GEsum += gprod * E["Ek"]
    if CHANNELS["T"]:
        # s is instantaneous: s_inf(V)
        s = ch.t_s_inf(V); u = gates["u"]; gprod = G["gT"] * (s**2) * u
        Gsum += gprod; GEsum += gprod * E["Eca"]
    return Gsum, GEsum

# NeuroVISOR: SolveStep(...)
# IMEX–SBDF2 (single compartment, scalar solve):
# (3C/(2dt)) V^{n+1} + G_n V^{n+1} =
#   (4C/(2dt)) V^n - (C/(2dt)) V^{n-1} + Iapp
#   + [ 2*(-G_n V^n + GE_n) - (-G_{n-1} V^{n-1} + GE_{n-1}) ]
def step_sbdf2(Vn, Vnm1, gates_n, gates_nm1, dt):
    # 1) update gates explicitly (NeuroVISOR: stateexplicitSBDF2)
    gates_np1 = {k: gates_n[k].copy() for k in gates_n}

    if CHANNELS["K"]:
        gates_np1["n"] = update_gate_sbdf2(gates_n["n"], gates_nm1["n"], Vn, Vnm1, ch.k_rates, lambda V: ch.k_rates(V)[1])
    if CHANNELS["Na"]:
        # m
        gates_np1["m"] = update_gate_sbdf2(gates_n["m"], gates_nm1["m"], Vn, Vnm1, ch.na_rates, lambda V: ch.na_rates(V)[1])
        # h (alpha_h, beta_h are indices 2,3 from na_rates)
        def a_h(V): return ch.na_rates(V)[2]
        def b_h(V): return ch.na_rates(V)[3]
        gates_np1["h"] = update_gate_sbdf2(gates_n["h"], gates_nm1["h"], Vn, Vnm1, a_h, b_h)
    if CHANNELS["CaH"]:
        aq, bq, ar, br = ch.ca_high_rates(Vn);  aqm1, bqm1, arm1, brm1 = ch.ca_high_rates(Vnm1)
        # q
        def a_q(V): return ch.ca_high_rates(V)[0]
        def b_q(V): return ch.ca_high_rates(V)[1]
        gates_np1["q"] = update_gate_sbdf2(gates_n["q"], gates_nm1["q"], Vn, Vnm1, a_q, b_q)
        # r
        def a_r(V): return ch.ca_high_rates(V)[2]
        def b_r(V): return ch.ca_high_rates(V)[3]
        gates_np1["r"] = update_gate_sbdf2(gates_n["r"], gates_nm1["r"], Vn, Vnm1, a_r, b_r)
    if CHANNELS["M"]:
        gates_np1["p"] = update_p_sbdf2(gates_n["p"], gates_nm1["p"], Vn, Vnm1)
    if CHANNELS["T"]:
        # u is dynamic via alpha/beta from t_alpha_beta_u; s is instantaneous inside reaction_terms
        def a_u(V): return ch.t_alpha_beta_u(V)[0]
        def b_u(V): return ch.t_alpha_beta_u(V)[1]
        gates_np1["u"] = update_gate_sbdf2(gates_n["u"], gates_nm1["u"], Vn, Vnm1, a_u, b_u)

    # 2) voltage update (NeuroVISOR: SolveStep + reactF use)
    Gn,  GEn  = reaction_terms(Vn,   gates_n)
    Gm1, GEm1 = reaction_terms(Vnm1, gates_nm1)

    lhs = (3.0*C_m)/(2.0*dt) + Gn
    rhs = ((4.0*C_m)/(2.0*dt))*Vn - ((C_m)/(2.0*dt))*Vnm1 + Iapp \
          + ( -2.0*Gn*Vn + 2.0*GEn ) - ( -Gm1*Vnm1 + GEm1 )

    Vnp1 = rhs / lhs
    return Vnp1, gates_np1

# ---------------------------
# Warm-up step: SBDF1/IMEX-Euler for startup (NeuroVISOR does Euler-like start)
# (C/dt + G_n) V^{1} = (C/dt) V^{0} + Iapp + GE_n
# ---------------------------
def warmup_step(V0, gates0, dt):
    # Euler for gates (explicit)
    g1 = {k: gates0[k].copy() for k in gates0}
    if CHANNELS["K"]:
        a,b = ch.k_rates(V0); g1["n"] += dt * f_gate(g1["n"], a, b); g1["n"] = np.clip(g1["n"],0,1)
    if CHANNELS["Na"]:
        am,bm,ah,bh = ch.na_rates(V0)
        g1["m"] += dt * f_gate(g1["m"], am, bm); g1["m"] = np.clip(g1["m"],0,1)
        g1["h"] += dt * f_gate(g1["h"], ah, bh); g1["h"] = np.clip(g1["h"],0,1)
    if CHANNELS["CaH"]:
        aq,bq,ar,br = ch.ca_high_rates(V0)
        g1["q"] += dt * f_gate(g1["q"], aq, bq); g1["q"] = np.clip(g1["q"],0,1)
        g1["r"] += dt * f_gate(g1["r"], ar, br); g1["r"] = np.clip(g1["r"],0,1)
    if CHANNELS["M"]:
        g1["p"] += dt * ((ch.m_p_inf(V0)-g1["p"]) / ch.m_tau_p(V0)); g1["p"] = np.clip(g1["p"],0,1)
    if CHANNELS["T"]:
        au,bu = ch.t_alpha_beta_u(V0)
        g1["u"] += dt * f_gate(g1["u"], au, bu); g1["u"] = np.clip(g1["u"],0,1)

    # Voltage IMEX-Euler
    Gn, GEn = reaction_terms(V0, g1)
    V1 = ((C_m/dt)*V0 + Iapp + GEn) / ((C_m/dt) + Gn)
    return V1, g1

# Tiny runner for demo/testing
def run(dt=5e-6, T=0.01, V0=-0.065):
    steps = int(round(T/dt))
    # Voltage arrays as 1-element vectors (keeps broadcast simple)
    Vnm1 = np.array([V0], dtype=float)
    # Gates at t0 (use your C# initials)
    gates0 = {}
    if CHANNELS["K"]:   gates0["n"] = np.array([0.0376969])
    if CHANNELS["Na"]:
        gates0["m"] = np.array([0.0147567])
        gates0["h"] = np.array([0.9959410])
    if CHANNELS["CaH"]:
        gates0["q"] = np.array([0.0]); gates0["r"] = np.array([0.0])
    if CHANNELS["M"]:
        gates0["p"] = np.array([0.0])
    if CHANNELS["T"]:
        gates0["u"] = np.array([0.0])

    # Warm-up (NeuroVISOR start)
    Vn, gates_n = warmup_step(Vnm1, gates0, dt)
    gates_nm1 = {k: gates0[k].copy() for k in gates0}

    traceV = np.empty(steps+1); traceV[0]=Vnm1[0]; traceV[1]=Vn[0]

    for i in range(2, steps+1):
        Vnp1, gates_np1 = step_sbdf2(Vn, Vnm1, gates_n, gates_nm1, dt)
        traceV[i] = Vnp1[0]
        # rotate
        Vnm1[:] = Vn; Vn[:] = Vnp1
        gates_nm1 = gates_n; gates_n = gates_np1

        # very simple instability stop
        if not np.isfinite(traceV[i]) or abs(traceV[i]) > 100:
            print(f"[STOP] Instability at step {i}, t={i*dt:.3e}s, V={traceV[i]}")
            return traceV[:i+1], np.linspace(0, i*dt, i+1)

    return traceV, np.linspace(0, T, steps+1)

# Script entry
if __name__ == "__main__":
    dt = 5.0e-6
    for d in [5.0e-6, 5.0e-7]:
        globals()["dt"] = d  # used in update_gate_sbdf2 closure
        V, t = run(dt=d, T=0.01, V0=-0.065)
        print(f"dt={d:.1e}, steps={len(V)-1}, V_end={V[-1]:.4f} V")
