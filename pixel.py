import functions as f
import numpy as np
import matplotlib.pyplot as plt
import math

polar_angles = np.linspace(0, 40, 100)
momenta = np.linspace(0.2, 1.5, 100)

polar_angles_rad = polar_angles * math.pi / 180

aa = np.zeros((len(polar_angles), len(momenta)))

beta_min = f.beta(f.mmu, momenta[0])
n = f.rindex(500)
theta_c_min = f.cherenkov(beta_min, f.rindex(500))

for i, theta_p in enumerate(polar_angles_rad):
  for j, p in enumerate(momenta):
    beta = f.beta(f.mmu, p)
    theta_c = f.cherenkov(beta, n)
    phi = (theta_c + theta_p - theta_c_min)
    aa[i, j] = phi / math.sqrt(12) * 1000
    print(theta_c, theta_p)
    print(aa[i, j])

extent = [momenta[0], momenta[-1], polar_angles[0], polar_angles[-1]]

# Pixel resolution
fig = plt.figure(figsize=(7, 4))
im2 = plt.imshow(aa,
                 aspect='auto',
                 origin='lower',
                 extent=extent,
                 cmap='gist_rainbow_r')
plt.title('Number of Pixels')
plt.xlabel('$p$ [GeV/c]')
plt.ylabel('$\\theta_p$ [deg]')
cbar2 = fig.colorbar(im2)
cbar2.set_label('$N$')

plt.tight_layout()
plt.savefig('pixel.png', dpi=600)
plt.savefig('pixel.pdf')
plt.close()
