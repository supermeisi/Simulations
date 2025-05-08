import numpy as np
import math
import matplotlib.pyplot as plt

import functions as f

# Constants
wavelength = 633e-9  # 633 nm in meters
theta_deg = 15
theta = np.radians(theta_deg)
d = 2e-2  # radiator thickness in meters
n = 1.46 # refractive index of fused silica

# Axes for the 2D grid
roughness = np.linspace(0, 10e-9, 100)  # From 0.1 to 5 nm
path_length = np.linspace(0, 10, 100)# From 0 to 1 m

# Create 2D grids
R, L = np.meshgrid(roughness, path_length)

# Calculate reflectivity per bounce due to roughness
reflectivity = np.exp(- (4 * np.pi * R * np.cos(theta) * n / wavelength)**2)

# Calculate total transmission after many bounces
# Number of bounces = L / d
transmission = reflectivity ** (L / d * math.tan(theta))

# Plot using imshow
plt.figure(figsize=(5, 4))
im = plt.imshow(transmission, extent=[roughness[0]*1e9, roughness[-1]*1e9, path_length[0], path_length[-1]],
                origin='lower', aspect='auto', cmap='rainbow')
cbar = plt.colorbar(im)
cbar.set_label('T')

plt.xlabel('$\sigma$ [nm]')
plt.ylabel('$L$ [m]')
plt.title('Surface Roughness')
plt.tight_layout()
#plt.show()
plt.savefig('roughness.png', dpi=600)
plt.savefig('roughness.pdf')