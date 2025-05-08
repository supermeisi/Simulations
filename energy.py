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
momenta_GeV = np.linspace(0.5, 1.5, 100)  # GeV/c
momenta = momenta_GeV * 1000  # MeV/c
thicknesses = np.linspace(0.1, 4.0, 100)  # cm
step_cm = 0.01

# Initialize results
rms_angle_change = np.zeros((len(thicknesses), len(momenta)))
energy_loss = np.zeros((len(thicknesses), len(momenta)))

n_glass = f.rindex(500)

m_mu = f.mmu * 1e3

# Main simulation loop
for i, thickness in enumerate(thicknesses):
  steps = int(thickness / step_cm)
  for j, p0 in enumerate(momenta):
    print(f"Thickness: {thickness} cm, Momentum: {p0} MeV/c")
    E = np.sqrt(p0**2 + m_mu**2)
    E_initial = E
    p = p0
    beta = p / E
    gamma = E / m_mu
    angles = []
    initial_angle = f.cherenkov(beta, n_glass)
    for _ in range(steps):
      angles.append(f.cherenkov(beta, n_glass))
      dE = f.bethe_bloch(f.mmu, beta, gamma) * step_cm
      E -= dE
      if m_mu >= E:
        E = m_mu
        break
      p = np.sqrt(E**2 - m_mu**2)
      beta = p / E
      gamma = E / m_mu
    # RMS of Cherenkov angle change (mrad)
    if len(angles) >= 2:
      angles = np.array(angles)
      rms = np.sqrt(np.mean((angles - initial_angle)**2))
    else:
      rms = 0.0
    rms_angle_change[i, j] = rms * 1e3
    energy_loss[i, j] = E_initial - E  # total energy loss in MeV

extent = [momenta_GeV[0], momenta_GeV[-1], thicknesses[0], thicknesses[-1]]

# Cherenkov RMS
fig = plt.figure(figsize=(5, 4))
im1 = plt.imshow(rms_angle_change,
                 aspect='auto',
                 origin='lower',
                 extent=extent,
                 cmap='gist_rainbow_r',
                 norm=colors.LogNorm())
plt.title('Cherenkov Angle Smearing')
plt.xlabel('$p$ [GeV/c]')
plt.ylabel('$x$ [cm]')
cbar1 = fig.colorbar(im1)
cbar1.set_label('$\\sigma$ [mrad]')

plt.tight_layout()
plt.savefig('cherenkov_energy.png', dpi=600)
plt.savefig('cherenkov_energy.pdf')
plt.close()

# Energy Loss
fig = plt.figure(figsize=(5, 4))
im2 = plt.imshow(energy_loss,
                 aspect='auto',
                 origin='lower',
                 extent=extent,
                 cmap='gist_rainbow_r')
plt.title('Energy Loss')
plt.xlabel('$p$ [GeV/c]')
plt.ylabel('$x$ [cm]')
cbar2 = fig.colorbar(im2)
cbar2.set_label('$\\Delta E$ [MeV]')

plt.tight_layout()
plt.savefig('energy_loss.png', dpi=600)
plt.savefig('energy_loss.pdf')
plt.close()
