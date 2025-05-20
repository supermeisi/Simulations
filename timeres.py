import functions as f
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors

m = [f.mmu, f.mpi, f.mK]

pmin = 0.5
pmax = 1.5

theta_min = 5 * math.pi / 180
theta_max = 35 * math.pi / 180

steps = 100

arr = np.zeros((steps, steps))

n = 3.  # separation power

for i in range(steps):
  for j in range(steps):
    p = pmin + i * 1.0 / 100 * (pmax - pmin)
    theta = theta_min + j * 1.0 / 100 * (theta_max - theta_min)

    top1 = f.top(m[0], p, f.R, 500, theta) * 1e12
    top2 = f.top(m[1], p, f.R, 500, theta) * 1e12

    res = 1 / n * (top2 - top1)

    arr[i, j] = res

    print(f'Time of Propagation: {top1} ps')
    print(f'Time Resolution: {res} ps')

fig = plt.figure(figsize=(5, 4))
plt.imshow(arr,
           aspect='auto',
           cmap='gist_rainbow',
           interpolation='nearest',
           norm=colors.LogNorm(),
           origin='lower',
           extent=[0.5, 1.5, 0, 35])
levels = np.logspace(-1, 2, 10)
plt.colorbar(label='$\Delta t$ [ps]', ticks=levels, format='%.1f')
plt.title('Time Resolution $\mu/\pi$')
plt.xlabel('$p$ [GeV/c]')
plt.ylabel('$\\theta_p$ [deg]')
plt.tight_layout()
plt.savefig('timeres_mupi.png', dpi=600)
plt.savefig('timeres_mupi.pdf')
plt.close()

pmin = 1.5
pmax = 6

theta_min = 5 * math.pi / 180
theta_max = 35 * math.pi / 180

steps = 100

arr = np.zeros((steps, steps))

for i in range(steps):
  for j in range(steps):
    p = pmin + i * 1.0 / 100 * (pmax - pmin)
    theta = theta_min + j * 1.0 / 100 * (theta_max - theta_min)

    top1 = f.top(m[1], p, f.R, 500, theta) * 1e12
    top2 = f.top(m[2], p, f.R, 500, theta) * 1e12

    res = 1 / n * (top2 - top1)

    arr[i, j] = res

    print(f'Time of Propagation: {top1} ps')
    print(f'Time Resolution: {res} ps')

fig = plt.figure(figsize=(5, 4))
plt.imshow(
    arr,
    aspect='auto',
    cmap='gist_rainbow',
    interpolation='nearest',
    norm=colors.LogNorm(),
    origin='lower',
    extent=[pmin, pmax, theta_min / math.pi * 180, theta_max / math.pi * 180])
levels = np.logspace(-1, 3, 20)
plt.colorbar(label='$\Delta t$ [ps]', ticks=levels, format='%.1f')
plt.title('Time Resolution $\pi/K$')
plt.xlabel('$p$ [GeV/c]')
plt.ylabel('$\\theta_p$ [deg]')
plt.tight_layout()
plt.savefig('timeres_piK.png', dpi=600)
plt.savefig('timeres_piK.pdf')
plt.close()
