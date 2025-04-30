import numpy as np
import matplotlib.pyplot as plt
import math
import functions as f

# Parameters for fused silica and muon
X0_fused_silica = 12.55  # radiation length in cm
mass_muon = 105.66  # MeV/c^2
z_muon = 1  # muon charge
N = 100  # numbr of steps

# Step sizes from 10 µm to 1 mm (in cm)
dx_values_cm = np.linspace(0, 4, N)

# Particle momentum in GeV/c between 0.5 and 1.5 GeV/c
p_values_GeV = np.linspace(0.5, 1.5, N)
p_values_MeV = p_values_GeV * 1000

# Calculate theta values for muon in fused silica
theta_values_mrad = [[
    1 / math.sqrt(2) *
    f.highland_theta(p, f.beta(mass_muon, p), z_muon, dx, X0_fused_silica)
    for p in p_values_MeV
] for dx in dx_values_cm]

# Convert step sizes to mm for plotting/output
dx_values_mm = dx_values_cm * 10

# Plot
fig = plt.figure(figsize=(6, 4))
plt.imshow(theta_values_mrad,
           aspect='auto',
           cmap='gist_rainbow_r',
           origin='lower',
           extent=[0.5, 1.5, 0, 40])
plt.colorbar(label='$\\theta_{0}$ [mrad]')
plt.xlabel('$p$ [GeV/c]')
plt.ylabel('$x$ [mm]')
plt.title('Angular Straggling in Fused Silica')
plt.grid(True)
plt.tight_layout()
plt.savefig('straggling.png', dpi=600)
plt.savefig('straggling.pdf')

np.random.seed(100)  # for reproducibility

# Initial direction angle (in radians)
theta = 0.0
theta_history = [theta]

p = 500  # momentum in MeV/c
L = 4.0  # length of material in cm
dx = L / N  # step size in cm

# Step through the material
for _ in range(N):
    sigma = f.highland_theta(p, f.beta(mass_muon, p), z_muon, dx,
                             X0_fused_silica)
    dtheta = np.random.normal(0, sigma)  # random angular kick
    theta += dtheta  # accumulate angle
    theta_history.append(theta)

# Convert angle history to mrad
theta_history_mrad = np.array(theta_history)

# Plot
fig = plt.figure(figsize=(7, 4))
plt.plot(np.linspace(0, L, N + 1), theta_history_mrad)
plt.xlabel("Depth [cm]")
plt.ylabel("Angle relative to initial direction [mrad]")
plt.title("Accumulated Angular Straggling (1D projection)")
plt.grid()
plt.tight_layout()
plt.savefig('straggling_step.png', dpi=600)
plt.savefig('straggling_step.pdf')
