import math
import numpy as np
import matplotlib.pyplot as plt

from functions import cherenkov, beta, rindex

thetap = 30 * math.pi / 180
mpi = 0.13957039
wlen = 500
n = rindex(wlen)

thetac = cherenkov(beta(mpi, 1), n)

print(thetac, thetap, wlen, n)

for phi in range(0, 360, 10):
  phi = phi * math.pi / 180
  gamma = math.acos(-math.sin(thetac) * math.cos(phi) * math.sin(thetap) +
                    math.cos(thetap) * math.cos(thetac))
  theta = math.asin(1 / n)

  print(gamma - theta)
