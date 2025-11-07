from neuron import h
h.load_file("stdrun.hoc")
import numpy as np
import matplotlib.pyplot as plt
import csv

# 1. Geometry and area

radius_um = 16.75       # µm
length_um = 0.251879699 # µm
radius_m = radius_um * 1e-6
length_m = length_um * 1e-6
soma_area_m2 = 2 * np.pi * radius_m * length_m
area_cm2 = soma_area_m2 * 1e4  # convert to cm²


# 2. Create soma

soma = h.Section(name='soma')
soma.L = length_um
soma.diam = 2 * radius_um
soma.cm = 1.0  # µF/cm²  (same as 1e-2 F/m²)
soma.nseg = 1


# 3. Insert channels (Na, K, leak)

for sec in [soma]:
    sec.insert("original_na")
    sec.insert("original_k")
    # sec.insert("pas")
    for seg in sec:
        seg.original_na.gnabar = 0.005      # S/cm²  (matches Python)
        seg.original_k.gkbar = 0.05     # S/cm²
        # seg.pas.g = 0.0001   # S/cm²
        # seg.pas.e = -70      # mV
        seg.ena = 50
        seg.ek  = -90


# 4. Stimulus (match Python’s total current)

# Python: stimAmplitude = 0.15e-11 / soma_area (A/m²)
# So total current = 0.15e-11 A = 0.015 nA

I_total_nA = 0.15  # same total current
stim_start = 50.0
stim_end = 150.0

iclamp = h.IClamp(soma(0.5))
iclamp.delay = stim_start
iclamp.dur = stim_end - stim_start
iclamp.amp = I_total_nA  # nA, same as Python total


# 5. Simulation control
h.v_init = 0
tstop = 200
h.dt = 0.05
h.tstop = tstop
h.finitialize(h.v_init)
h.run()

# 6. Record vectors

t_vec = h.Vector().record(h._ref_t)
v_vec = h.Vector().record(soma(0.5)._ref_v)


# 7. Run simulation

h.finitialize(h.v_init)
h.run()


# 8. Plot

plt.figure(figsize=(9,6))
plt.plot(t_vec, v_vec, 'k', linewidth=2)
plt.xlabel("Time (ms)")
plt.ylabel("Membrane potential (mV)")
plt.title("NEURON Simulation (Aligned with Python Euler)")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()


# 9. Save trace

with open("./neuron_recordings/neuron_trace.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["t(ms)", "V(mV)"])
    for t, v in zip(t_vec, v_vec):
        writer.writerow([float(t), float(v)])
print("Exported: neuron_recordings/neuron_trace.csv")