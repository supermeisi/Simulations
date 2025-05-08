import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import functions as f

cherenkov = f.cherenkov(1, f.rindex(500))
print(cherenkov)

momentum = np.linspace(0.5, 1.5, 100)

angles = [f.cherenkov(f.beta(f.mmu, p), f.rindex(500)) for p in momentum]

# print(angles)

plt.plot(momentum, angles)
plt.title('Cherenkov Angle for Muons')
plt.xlabel('p [GeV/c]')
plt.ylabel('Cherenkov angle [rad]')
plt.savefig('cherenkov_angle.png', dpi=600)
plt.savefig('cherenkov_angle.pdf')
plt.close()

# Grid definitions
momenta_GeV = np.linspace(0.5, 2.0, 100)  # GeV/c
momenta = momenta_GeV * 1000  # MeV/c
thicknesses = np.linspace(0.1, 4.0, 100)  # cm
step_cm = 0.01

# Initialize results
rms_angle_change = np.zeros((len(thicknesses), len(momenta)))
energy_loss = np.zeros((len(thicknesses), len(momenta)))

n_glass = f.rindex(500)

# Main simulation loop
for i, thickness in enumerate(thicknesses):
  steps = int(thickness / step_cm)
  for j, p0 in enumerate(momenta):
    print(f"Thickness: {thickness} cm, Momentum: {p0} MeV/c")
    E = np.sqrt(p0**2 + f.mmu**2)
    E_initial = E
    p = p0
    beta = p / E
    gamma = E / f.mmu
    angles = []
    initial_angle = f.cherenkov(beta, n_glass)
    for _ in range(steps):
      angles.append(f.cherenkov(beta, n_glass))
      dE = f.bethe_bloch(f.mmu, beta, gamma) * step_cm
      E -= dE
      if E <= f.mmu:
        E = f.mmu
        break
      p = np.sqrt(E**2 - f.mmu**2)
      beta = p / E
      gamma = E / f.mmu
    # RMS of Cherenkov angle change (mrad)
    if len(angles) >= 2:
      angles = np.array(angles)
      rms = np.sqrt(np.mean((angles - initial_angle)**2))
    else:
      rms = 0.0
    rms_angle_change[i, j] = rms * 1e3
    energy_loss[i, j] = E_initial - E  # total energy loss in MeV


# Plotting side-by-side heatmaps
fig, axs = plt.subplots(1, 2, figsize=(14, 6))

extent = [momenta_GeV[0], momenta_GeV[-1], thicknesses[0], thicknesses[-1]]

# Cherenkov RMS
im1 = axs[0].imshow(rms_angle_change, aspect='auto', origin='lower', extent=extent, cmap='gist_rainbow_r', norm=colors.LogNorm())
axs[0].set_title('Cherenkov Angle Smearing')
axs[0].set_xlabel('$p$ [GeV/c]')
axs[0].set_ylabel('$x$ [cm]')
cbar1 = fig.colorbar(im1, ax=axs[0])
cbar1.set_label('$\\sigma$ [mrad]')

# Energy Loss
im2 = axs[1].imshow(energy_loss, aspect='auto', origin='lower', extent=extent, cmap='gist_rainbow_r')
axs[1].set_title('Energy Loss')
axs[1].set_xlabel('$p$ [GeV/c]')
axs[1].set_ylabel('$x$ [cm]')
cbar2 = fig.colorbar(im2, ax=axs[1])
cbar2.set_label('$\\Delta E$ [MeV]')

plt.tight_layout()
plt.savefig('energy_loss.png', dpi=600)
plt.savefig('energy_loss.pdf')
plt.close()