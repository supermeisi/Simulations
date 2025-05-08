import numpy as np
import matplotlib.pyplot as plt

# Constants
wavelength = 633e-9  # 633 nm in meters
n = 1.46             # Refractive index of fused silica
p = 0.22             # Photoelastic coefficient (unitless)
k_B = 1.380649e-23   # Boltzmann constant (J/K)
T_kelvin = 300       # Temperature in Kelvin
beta_T = 7e-11       # Isothermal compressibility (Pa^-1)

# Rayleigh scattering coefficient (alpha_sc in 1/m)
def compute_alpha_sc(wavelength, n, p, T_kelvin, beta_T):
    term1 = (8 * np.pi**3) / (3 * wavelength**4)
    term2 = n**8 * p**2 * k_B * T_kelvin * beta_T / (n**2 - 1)**2
    return term1 * term2

alpha_sc = compute_alpha_sc(wavelength, n, p, T_kelvin, beta_T)

# Path lengths in meters
L = np.linspace(0, 1000, 1000)  # From 0 to 10 meters
transmission = np.exp(-alpha_sc * L)

# Plot
plt.figure(figsize=(5, 4))
plt.plot(L, transmission, label=f'{wavelength*1e9:.0f} nm')
plt.xlabel('Path Length (m)')
plt.ylabel('Transmission')
plt.title('Photon Transmission in Fused Silica vs. Path Length')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('transmission.png', dpi=600)
plt.savefig('transmission.pdf')