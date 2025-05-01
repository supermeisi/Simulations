from functions import *

m = [mmu, mpi, mK]
colors = ['b', 'g', 'r']
labels = ['mu', 'pi', 'K']

fig = plt.figure(figsize=(4, 3))

for i in range(3):
  xx = []
  yy1 = []
  yy2 = []

  for p in np.arange(0.1, 4., 0.01):
    cherenkov_angle_min = cherenkov(beta(m[i], p), rindex(300.))
    cherenkov_angle_max = cherenkov(beta(m[i], p), rindex(800.))

    xx.append(p)
    yy1.append(cherenkov_angle_min)
    yy2.append(cherenkov_angle_max)

  plt.fill_between(xx, yy1, yy2, color=colors[i], alpha=0.1, label=labels[i])
  plt.plot(xx, yy1, color=colors[i])
  plt.plot(xx, yy2, color=colors[i])

plt.margins(x=0)
plt.ylabel("$\\theta_c$ [rad]")
plt.xlabel("$p$ [GeV/c]")
plt.ylim(0.5, 0.9)
plt.legend()
plt.grid()
#plt.title('Group Velocity')
plt.tight_layout()
plt.savefig('dispersion.pdf')
plt.savefig('dispersion.png', dpi=600)
plt.show()
plt.close()