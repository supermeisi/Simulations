import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import math

import functions as f

momenta = np.linspace(0.1, 6.0, 100)
masses = np.linspace(0.1, 1.0, 100)
polar_angles = np.linspace(0, 40, 100)

polar_angles_rad = polar_angles * math.pi / 180

rindex = np.zeros((len(momenta), len(masses)))

for i, m in enumerate(masses):
  for j, p in enumerate(momenta):
    theta_p = 0
    beta = f.beta(m, p)
    theta_c = f.cherenkov(beta, f.rindex(500))
    if theta_c == 0 and theta_p == 0:
      continue
    n = 1 / math.sin(theta_p + theta_c)
    if n < 1 or n > 2:
      continue
    print(n, theta_c, theta_p)
    rindex[j, i] = n

rindex[rindex == 0] = np.nan

fig = plt.figure(figsize=(5, 4))
plt.imshow(rindex,
           aspect='auto',
           cmap='gist_rainbow_r',
           interpolation='nearest',
           origin='lower',
           extent=[masses[0], masses[-1], momenta[0], momenta[-1]])
plt.colorbar(label='$n$', format='%.1f')
plt.title('Minimum Refractive Index')
plt.xlabel('$m$ [GeV]')
plt.ylabel('$p$ [GeV/c]')
plt.tight_layout()
plt.savefig('rindex_mass.png', dpi=600)
plt.savefig('rindex_mass.pdf')
plt.close()

for i, theta_p in enumerate(polar_angles_rad):
  for j, p in enumerate(momenta):
    beta = f.beta(f.mpi, p)
    theta_c = f.cherenkov(beta, f.rindex(500))
    if theta_c == 0 and theta_p == 0:
      continue
    n = 1 / math.sin(theta_p + theta_c)
    if n < 1 or n > 2:
      continue
    print(n, theta_c, theta_p)
    rindex[i, j] = n

rindex[rindex == 0] = np.nan

fig = plt.figure(figsize=(5, 4))
plt.imshow(rindex,
           aspect='auto',
           cmap='gist_rainbow_r',
           interpolation='nearest',
           origin='lower',
           extent=[momenta[0], momenta[-1], polar_angles[0], polar_angles[-1]])
plt.colorbar(label='$n$', format='%.1f')
plt.title('Minimum Refractive Index')
plt.xlabel('$p$ [GeV/c]')
plt.ylabel('$\\theta_p$ [deg]')
plt.tight_layout()
plt.savefig('rindex_polar_angle.png', dpi=600)
plt.savefig('rindex_polar_angle.pdf')
plt.close()
