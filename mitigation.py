from functions import *

aa = []
for i in range(steps):
  bb = []
  for j in range(steps):
    theta = (theta_max - theta_min) / steps * i + theta_min
    p = (pmax - pmin) / steps * j + pmin
    dt = top(mpi, p, distance(R, d, theta, D / 2), wlen_min, theta) - top(
        mpi, p, distance(R, d, theta, D / 2), wlen_max, theta)
    bb.append(dt * 1e12)
  aa.append(bb)

fig = plt.figure(figsize=(5, 4))

#plt.rcParams.update({'font.size': 30})
plt.imshow(
    aa,
    aspect='auto',
    cmap=cm.gist_rainbow,
    origin='lower',
    extent=[pmin, pmax, theta_min / math.pi * 180, theta_max / math.pi * 180])
plt.colorbar(label='$\sigma_t$ [ps]')
plt.xlabel("$p$ [GeV/c]")
plt.ylabel("$\\theta_p$ [deg]")
plt.title('Dipersion Mitigation')
fig.tight_layout()
plt.savefig('mitigation.pdf')
plt.savefig('mitigation.png')
plt.close()
