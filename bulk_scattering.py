import numpy as np
import matplotlib.pyplot as plt

# Constants
wavelength1 = 800e-9  # 800 nm in meters
wavelength2 = 500e-9  # 633 nm in meters
wavelength3 = 300e-9  # 300 nm in meters

n = 1.46  # Refractive index of fused silica
p = 0.22  # Photoelastic coefficient (unitless)
k_B = 1.380649e-23  # Boltzmann constant (J/K)
T_kelvin = 300  # Temperature in Kelvin
beta_T = 7e-11  # Isothermal compressibility (Pa^-1)


# Rayleigh scattering coefficient (alpha_sc in 1/m)
def compute_alpha_sc(wavelength, n, p, T_kelvin, beta_T):
    term1 = (8 * np.pi**3) / (3 * wavelength**4)
    term2 = n**8 * p**2 * k_B * T_kelvin * beta_T / (n**2 - 1)**2
    return term1 * term2


alpha_sc1 = compute_alpha_sc(wavelength1, n, p, T_kelvin, beta_T)
alpha_sc2 = compute_alpha_sc(wavelength2, n, p, T_kelvin, beta_T)
alpha_sc3 = compute_alpha_sc(wavelength3, n, p, T_kelvin, beta_T)

# Path lengths in meters
L = np.linspace(0, 10, 1000)  # From 0 to 10 meters
transmission1 = np.exp(-alpha_sc1 * L)
transmission2 = np.exp(-alpha_sc2 * L)
transmission3 = np.exp(-alpha_sc3 * L)

# Plot
plt.figure(figsize=(5, 4))
plt.plot(L, transmission1, label=f'{wavelength1*1e9:.0f} nm')
plt.plot(L, transmission2, label=f'{wavelength2*1e9:.0f} nm')
plt.plot(L, transmission3, label=f'{wavelength3*1e9:.0f} nm')
plt.xlim(L[0], L[-1])
plt.xlabel('Path Length (m)')
plt.ylabel('Transmission')
plt.title('Photon Transmission')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('transmission.png', dpi=600)
plt.savefig('transmission.pdf')
