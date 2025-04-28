import math
import matplotlib.pyplot as plt

from functions import cherenkov, beta, rindex

thetap = 30 * math.pi / 180
mpi = 0.13957039
wlen = 500
n = rindex(wlen)

thetac = cherenkov(beta(mpi, 1), n)

print(thetac, thetap, wlen, n)

x = []
y = []

for phi in range(0, 360, 1):
  phi = phi * math.pi / 180
  x.append(phi)
  gamma = math.acos(-math.sin(thetac) * math.cos(phi) * math.sin(thetap) +
                    math.cos(thetap) * math.cos(thetac))
  theta = math.asin(1 / n)
  psi = gamma - theta
  y.append(psi)
  print(phi, gamma, theta, psi)

fig = plt.figure(figsize=(6, 4))
plt.plot(x, y)
plt.xlabel("$\phi$ [rad]")
plt.ylabel("$\psi$ [rad]")
plt.grid()
plt.title('Trapping Angle')
plt.tight_layout()
plt.ylim(-0.5, 0.5)
plt.xlim(0, 2 * math.pi)
plt.axhline(y=0, color='r', linestyle='-')
plt.axvline(x=math.pi, color='g', linestyle='-')
plt.savefig('psi.pdf')
plt.savefig('psi.png', dpi=600)
