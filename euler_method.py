import numpy as np
import matplotlib.pyplot as plt
import ion_channels as ch  # your channels-only module

# Section 1: Define variables

# Channel toggles (choose what to integrate)
USE_NA  = True # m, h
USE_K   = True # n
USE_LEAK = False
USE_CaH = False  # q, r (high-threshold Ca) -- off here
USE_T   = False # u (T-type Ca). s_inf is instantaneous and NOT integrated.
USE_M   = False # p (slow K / M-current)

# Membrane Capacitance
C_m = 1.0e-2

def I_stim(t):
    return 50.0 if (5e-3 <= t < 10e-3) else 0.0

# Section 2: Update Gating Variables (Forward Euler at V)

# State vector layout: y = [m, n, h, p, q, u]
# If a channel is OFF, we keep its derivative = 0 (state frozen).

def rhs_gates(gates, t, V):
    # This is the RHS for gating variables at fixed voltage V
    # gates = [m, n, h, p, q, u, r]
    m, n, h, p, q, u, r = gates
    Vv = np.array([V])  # ion_channels expects array input

    # Initialize derivatives
    dm = dn = dh = dp = dq = du = dr = 0.0

    if USE_NA:
        am, bm, ah, bh = ch.na_rates(Vv)
        dm = (am*(1.0 - m) - bm*m)[0]
        dh = (ah*(1.0 - h) - bh*h)[0]

    if USE_K:
        an, bn = ch.k_rates(Vv)
        dn = (an*(1.0 - n) - bn*n)[0]        

    if USE_M:
        # M-current gate: dp/dt = (p_inf - p)/tau_p
        dp = ((ch.m_p_inf(Vv) - p) / ch.m_tau_p(Vv))[0]

    if USE_CaH:
        aq, bq, ar, br = ch.ca_high_rates(Vv)
        dq = (aq*(1.0 - q) - bq*q)[0]
        dr = (ar*(1.0 - r) - br*r)[0]  

    if USE_T:
        au, bu = ch.t_alpha_beta_u(Vv)
        du = (au*(1.0 - u) - bu*u)[0]
        # s is also applied, but is instantaneous

    return np.array([dm, dn, dh, dp, dq, du, dr], dtype=float)

# Section 3: Build ionic membrane using gating variables and updated V with Euler:
    # C_m dV/dt = Im - I_stim(t)  =>  V_{k+1} = V_k + dt*(Im - Istim)/C_m

def ionic_current(V, gates, t):
    Im = 0.0
    if USE_K:
        Im += ch.G_DEFAULT['gk'] * (gates['n']**4) * (V - ch.E_DEFAULT['Ek'])
    if USE_NA:
        Im += ch.G_DEFAULT['gna'] * (gates['m']**3) * gates['h'] * (V - ch.E_DEFAULT['Ena'])
    if USE_LEAK:
        Im += ch.G_DEFAULT['gl'] * (V - ch.E_DEFAULT['El'])
    if USE_CaH:
        Im += ch.G_DEFAULT['gcaH'] * (gates['q']**2) * (gates['r']) * (V - ch.E_DEFAULT('Eca'))
    if USE_T:
        s = ch.t_s_inf(np.array([V]))[0]  # instantaneous activation
        Im += ch.G_DEFAULT['gT'] * (s**2) * gates['u'] * (V - ch.E_DEFAULT['Eca'])
    if USE_M:
        Im += ch.G_DEFAULT['gM'] * gates['p'] * (V - ch.E_DEFAULT['Ek'])
    return Im + I_stim(t)

# Section 3: Repeat

def forward_euler(y0, dt, steps, rhs, V0):
    # y0 is the init for gating variables
    y = np.array(y0, dtype=float)
    Y = np.empty((steps + 1, y.size), dtype=float)
    t = np.linspace(0.0, steps*dt, steps + 1)
    Y[0] = y
    for k in range(steps):
        tk = k * dt
        # pass V0 to rhs as your code intends
        y = y + dt * rhs(y, tk, V0)
        # clamp to [0,1] since these are probabilities
        np.clip(y, 0.0, 1.0, out=y)
        Y[k + 1] = y
    return t, Y

# Section 5: Run the program

if __name__ == "__main__":
    # Simulation settings
    T  = 0.02 # seconds
    dt = 1e-6 # seconds
    steps = int(round(T / dt))
    V0 = -0.065 # intial Voltage (in V)

    # m, n, h, p, q, u, r = gates
    y0 = np.array([
        0.0147567 if USE_NA else 0.0, # m
        0.0376969 if USE_K  else 0.0, # n
        0.9959410 if USE_NA else 0.0, # h
        0.0 if USE_M  else 0.0, # p
        0.0, # q
        0.0, # u
        0.0 #r
    ], dtype=float)

    # Integrate gates at fixed voltage V0
    t, Y = forward_euler(y0, dt, steps, rhs_gates, V0)

    labels = ["m", "n", "h", "p", "q", "u", "r"]
    for i, name in enumerate(labels):
        plt.plot(t*1e3, Y[:, i], label=name)
    plt.xlabel("t (ms)")
    plt.ylabel("gate value")
    plt.legend()
    plt.tight_layout()
    plt.show()