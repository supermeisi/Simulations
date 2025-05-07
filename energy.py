import numpy as np
import matplotlib.pyplot as plt

import functions as f

cherenkov = f.cherenkov(1, f.rindex(500))
print(cherenkov)

momentum = np.linspace(0.5, 1.5, 100)

angles = [f.cherenkov(f.beta(f.mmu, p), f.rindex(500)) for p in momentum]

print(angles)

plt.plot(momentum, angles)
plt.title('Cherenkov Angle for Muons')
plt.xlabel('p [GeV/c]')
plt.ylabel('Cherenkov angle [rad]')
plt.show()