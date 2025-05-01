import matplotlib.pyplot as plt
import numpy as np

import functions as f

m = [f.mmu, f.mpi, f.mK]
colors = ['b', 'g', 'r']
labels = ['mu', 'pi', 'K']

fig = plt.figure(figsize=(5, 4))

for i in range(3):
  xx = []
  yy1 = []
  yy2 = []

  for p in np.arange(0.1, 4., 0.01):
    cherenkov_angle_min = f.cherenkov(f.beta(m[i], p), f.rindex(300.))
    cherenkov_angle_max = f.cherenkov(f.beta(m[i], p), f.rindex(800.))

    xx.append(p)
    yy1.append(cherenkov_angle_min)
    yy2.append(cherenkov_angle_max)

  plt.fill_between(xx, yy1, yy2, color=colors[i], alpha=0.1, label=labels[i])
  plt.plot(xx, yy1, color=colors[i])
  plt.plot(xx, yy2, color=colors[i])

plt.margins(x=0)
plt.ylabel("$\\theta_c$ [rad]")
plt.xlabel("$p$ [GeV/c]")
plt.ylim(0.7, 0.85)
plt.legend()
plt.grid()
#plt.title('Group Velocity')
plt.tight_layout()
plt.savefig('dispersion.pdf')
plt.savefig('dispersion.png', dpi=600)