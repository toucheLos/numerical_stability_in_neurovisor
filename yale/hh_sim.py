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

# Toggle channels on/off here:
USE_NA = True
USE_K = True
USE_CA = False
USE_KM = False
USE_CAT = False
USE_PAS = True

for sec in [soma]:
    if USE_PAS:
        sec.insert("pas")
        for seg in sec:
            seg.pas.g = 1e-8  # S/cm²
            seg.pas.e = -70   # mV

    if USE_NA:
        sec.insert("original_na")
        for seg in sec:
            seg.original_na.gnabar = 50e-3  # S/cm²
            seg.ena = 50

    if USE_K:
        sec.insert("original_k")
        for seg in sec:
            seg.original_k.gkbar = 5e-3  # S/cm²
            seg.ek = -90

    if USE_CA:
        sec.insert("original_ca")
        for seg in sec:
            seg.original_ca.gcabar = 1e-4  # S/cm²
            seg.eca = 120

    if USE_KM:
        sec.insert("original_kM")
        for seg in sec:
            seg.original_kM.gMbar = 7.5e-5

    if USE_CAT:
        sec.insert("original_caT")
        for seg in sec:
            seg.original_caT.gTbar = 4.0e-2
# 
# V Threshold
# soma(0.5).original_na.vshift = -20.0
# soma(0.5).original_k.vshift = -20.0

# 4. Stimulus (match Python’s total current)

# Python: stimAmplitude = 0.15e-11 / soma_area (A/m²)
# So total current = 0.15e-11 A = 0.015 nA

I_total_nA = 0.0015
stim_start = 50.0
stim_end = 150.0

iclamp = h.IClamp(soma(0.5))
iclamp.delay = stim_start
iclamp.dur = stim_end - stim_start
iclamp.amp = I_total_nA  # nA, same as Python total


# 5. Simulation control
h.v_init = 0
tstop = 500
h.dt = 50e-5
h.tstop = tstop
h.finitialize(h.v_init)

# print initial gating values (after initialization, before h.run())
print("Initial gating values:")
for sec in [soma]:
    for seg in sec:
        print(" Section:", sec.name(), " seg:", seg.x)
        if hasattr(seg, "original_na"):
            print("  original_na m,h:", float(seg.original_na.m), float(seg.original_na.h))
        if hasattr(seg, "original_k"):
            print("  original_k n:", float(seg.original_k.n))
        if hasattr(seg, "original_kM"):
            print("  original_kM p (or other):", float(seg.original_kM.p))  
        if hasattr(seg, "original_caT"):
            print("  original_caT u:", float(seg.original_caT.u))
        if hasattr(seg, "original_ca"):
            print("  original_ca r:", float(seg.original_ca.r)) 
            print("  original_ca q:", float(seg.original_ca.q))


# 6. Record vectors

t_vec = h.Vector().record(h._ref_t)
v_vec = h.Vector().record(soma(0.5)._ref_v)


# 7. Run simulation

h.run()

# 9. Save trace

# 9. Save trace

channels = []
if USE_NA:   channels.append("na")
if USE_K:    channels.append("k")
if USE_PAS:  channels.append("leak")
if USE_CA:   channels.append("highca")
if USE_CAT:  channels.append("lowca")
if USE_KM:   channels.append("slowk")

filename = f"../neuron_recordings/{'_'.join(channels)}.csv"
# filename = "../neuron_recordings/neuron_trace.csv"

with open(filename, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["t(ms)", "V(mV)"])
    for t, v in zip(t_vec, v_vec):
        writer.writerow([float(t), float(v)])

print(f"Exported: {filename}")

# 8. Plot

plt.figure(figsize=(9,6))
plt.plot(t_vec, v_vec, 'k', linewidth=2)
plt.xlabel("Time (ms)")
plt.ylabel("Membrane potential (mV)")
plt.title("NEURON Simulation (Aligned with Python Euler)")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()


