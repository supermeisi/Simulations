import matplotlib.pyplot as plt
import scipy as sp

# Local imports should be after third-party imports
from functions import cherenkov, franck_tamm, rand_wlen, mK, rindex, mpi, mmu, beta

aa = []

wlen_min = 200
wlen_max = 900

n_bins = int((wlen_max - wlen_min) / 10)

print(franck_tamm(300, rindex(300), mpi, 1))
integral = sp.integrate.quad(lambda x: franck_tamm(x, rindex(x), mpi, 1),
                             wlen_min, wlen_max)
n_photons = int(integral[0]) * 0.02
print(n_photons)

for _ in range(int(n_photons)):
  aa.append(rand_wlen(mpi, 1, wlen_min, wlen_max))

fig = plt.figure(figsize=(5, 4))
plt.hist(aa, bins=n_bins)
plt.xlabel('$\lambda$ [nm]')
plt.ylabel('Entries')
plt.xlim(wlen_min, wlen_max)
plt.tight_layout()
plt.savefig('wavelengths.pdf')
plt.savefig('wavelengths.png', dpi=600)
plt.close()

wlens = []
bb = []
cc = []
dd = []

for i in range(1000):
  wlen = wlen_min + (wlen_max - wlen_min) / 1000 * i
  wlens.append(wlen)
  bb.append(cherenkov(beta(mmu, 1), rindex(wlen)))
  cc.append(cherenkov(beta(mpi, 1), rindex(wlen)))
  dd.append(cherenkov(beta(mK, 1), rindex(wlen)))

fig = plt.figure(figsize=(5, 4))
plt.plot(wlens, bb, label='$\mu$')
plt.plot(wlens, cc, label='$\pi$')
plt.plot(wlens, dd, label='$K$')
plt.xlabel('$\lambda$ [nm]')
plt.ylabel('$\\theta_c$ [rad]')
plt.legend()
plt.tight_layout()
plt.savefig('cherenkov.pdf')
plt.savefig('cherenkov.png', dpi=600)
