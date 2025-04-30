import math
import matplotlib.pyplot as plt

from functions import cherenkov, beta, rindex, mmu, mpi, mK

thetap = 30 * math.pi / 180
masses = [mmu, mpi, mK]
wlen = 500
n = rindex(wlen)
p = 1

xx = []
yy = []

for m in masses:
  thetac = cherenkov(beta(m, p), n)
  
  print(thetac, thetap, wlen, n)

  x = []
  y = []
  
  for phi in range(0, 360, 10):
    phi = phi * math.pi / 180
    x.append(phi)
    gamma = math.acos(-math.sin(thetac) * math.cos(phi) * math.sin(thetap) +
                      math.cos(thetap) * math.cos(thetac))
    theta = math.asin(1 / n)
    psi = gamma - theta
    y.append(psi)
    print(phi, gamma, theta, psi)

  xx.append(x)
  yy.append(y)
  
fig = plt.figure(figsize=(6, 4))
plt.plot(xx[0], yy[0], label='mu')
plt.plot(xx[1], yy[1], label='pi')
plt.plot(xx[2], yy[2], label='K')
plt.xlabel("$\phi$ [rad]")
plt.ylabel("$\psi$ [rad]")
plt.grid()
plt.title('Trapping Angle')
plt.tight_layout()
plt.xlim(0, 2 * math.pi)
plt.axhline(y=0, color='r', linestyle='-')
plt.axvline(x=math.pi, color='g', linestyle='-')
plt.legend()
plt.savefig('psi.pdf')
plt.savefig('psi.png', dpi=600)
