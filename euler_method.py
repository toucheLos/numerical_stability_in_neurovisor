import numpy as np
import matplotlib.pyplot as plt
import ion_channels as ch  # channels-only module

# Section 1: Define variables

# Channel toggles (choose what to integrate)
USE_NA = True # m, h
USE_K = True # n
USE_LEAK = False
USE_CaH = False  # q, r (high-threshold Ca) -- off here
USE_T = True # u (T-type Ca). s_inf is instantaneous and NOT integrated.
USE_M = False # p (slow K / M-current)

# Membrane Capacitance
C_m = 1.0e-2

# Soma geometry
# radius = 5.0 # µm  (from NV Debug.Log)
radius = 16.75 # mini-soma
# length = 0.142857142857143 # µm  (Neuron.TargetEdgeLength)
length = 0.25187969924812 # mini-soma
soma_area = 2.0*np.pi*(radius*1e-6)*(length*1e-6)   # m^2
print("soma_area = " + str(soma_area))
C_tot = C_m * soma_area  # F (total capacitance)

def I_stim(t):
    stimAmplitude = 0.15e-11 / soma_area
    # print("stimAmp = " + str(stimAmplitude))
    # return 0.0
    # return 0.03149 if (50e-3 <= t < 100e-3) else 0.0
    return stimAmplitude if (50e-3 <= t < 150e-3) else 0.0

# import math
# r_um = 33.5       # µm
# L_um = 0.142857142857143       # set to your NV TargetEdgeLength (µm)

# A_sphere = 4*math.pi*(r_um*1e-6)**2
# A_cyl    = 2*math.pi*(r_um*1e-6)*(L_um)
# print("A_sphere =", A_sphere, "m^2")
# print("A_cyl    =", A_cyl,    "m^2")

# Section 2: Update Gating Variables (Forward Euler at V)

# State vector layout: y = [m, n, h, p, q, u]
# If a channel is OFF, we keep its derivative = 0 (state frozen).

def rhs(y, t):
    # This is the RHS for gating variables at fixed voltage V
    # gates = [m, n, h, p, q, u, r]
    V, m, n, h, p, q, u, r = y
    Vv = np.array([V])  # ion_channels expects array input

    # Initialize derivatives
    dm = dn = dh = dp = dq = du = dr = 0.0

    if USE_NA:
        am, bm, ah, bh = ch.na_rates(Vv)
        dm = (am[0]*(1.0 - m) - bm[0]*m)
        dh = (ah[0]*(1.0 - h) - bh[0]*h)

    if USE_K:
        an, bn = ch.k_rates(Vv)
        dn = (an[0]*(1.0 - n) - bn[0]*n)

    if USE_M:
        # M-current gate: dp/dt = (p_inf - p)/tau_p
        dp = ((ch.m_p_inf(Vv) - p) / ch.m_tau_p(Vv))[0]

    if USE_CaH:
        aq, bq, ar, br = ch.ca_high_rates(Vv)
        dq = (aq[0]*(1.0 - q) - bq[0]*q)
        dr = (ar[0]*(1.0 - r) - br[0]*r)

    if USE_T:
        au, bu = ch.t_alpha_beta_u(Vv)
        du = (au[0]*(1.0 - u) - bu[0]*u)
        s = ch.t_s_inf(np.array([V]))[0] # instantaneous


    # Section 3: Compute Ionic Current
    Im = 0.0
    if USE_K:
        Im += ch.k_current(V, n)
    if USE_NA:
        Im += ch.na_current(V, m, h)
    if USE_LEAK:
        Im += ch.leak_current(V)
    if USE_CaH:
        Im += ch.ca_high_current(V, q, r)
    if USE_T:
        Im += ch.G_DEFAULT['gT'] * (s**2) * u * (V - ch.E_DEFAULT['Eca'])
    if USE_M:
        Im += ch.m_current(V, p)

    Im_tot = Im / soma_area

    dV = (I_stim(t) - Im) / C_m
    # dV = (I_stim(t) - Im_tot) / C_tot

    # print(str(dV_tot - dV))
    
    return np.array([dV, dm, dn, dh, dp, dq, du, dr], dtype=float)

# Section 3: Repeat

def forward_euler(y0, dt, steps, rhs):
    # y0 is the init for gating variables
    y = np.array(y0, dtype=float)
    Y = np.empty((steps + 1, y.size), dtype=float)
    t = np.linspace(0.0, steps*dt, steps + 1)

    # initialize first element
    Y[0] = y
    t[0] = 0.0

    for k in range(steps):
        tk = k * dt
        # t[k] = tk
        y = y + dt * rhs(y, tk)
        Y[k + 1] = y

    t[steps] = steps * dt
    return t, Y

# Section 5: Run the program

# Change timesteps to capture convergence of a solution

if __name__ == "__main__":
    # Simulation settings
    T  = 200e-3 # seconds
    dt = 50e-7 # seconds
    steps = int(round(T / dt))
    V0 = 0.0 # intial Voltage (in V)

    # m, n, h, p, q, u, r = gates
    y0 = np.array([
        V0,
        0.0147567 if USE_NA else 0.0, # m
        0.0376969 if USE_K else 0.0, # n
        0.9959410 if USE_NA else 0.0, # h
        0.02 if USE_M  else 0.0, # p
        0.02, # q
        0.012 if USE_T else 0.0, # u
        0.01 #r
    ], dtype=float)

    # Integrate the system
    t, Y = forward_euler(y0, dt, steps, rhs)

    # Export time and voltage to CSV
    export_data = np.column_stack((t * 1e3, Y[:, 0] * 1e3))  # time (ms), voltage (mV)
    np.savetxt("./neuron_recordings/euler_trace.csv", export_data, delimiter=",", header="t(ms),V(mV)", comments='')
    print("Exported: euler_trace.csv")

    # Plot voltage
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    ax1.plot(t*1e3, Y[:, 0]*1e3, 'k', linewidth=2)
    ax1.set_xlabel("t (ms)")
    ax1.set_ylabel("V (mV)")
    ax1.set_title("Membrane Voltage")
    ax1.grid(True, alpha=0.3)
    
    # Plot gating variables
    labels = ["m", "n", "h", "p", "q", "u", "r"]
    for i, name in enumerate(labels):
        if Y[:, i+1].max() > 0.01:  # Only plot if gate is active
            ax2.plot(t*1e3, Y[:, i+1], label=name)
    
    ax2.set_xlabel("t (ms)")
    ax2.set_ylabel("gate value")
    ax2.set_title("Gating Variables")
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()

