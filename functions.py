import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
import random as rand
from scipy.integrate import quad
from scipy.misc import derivative
import matplotlib.patches as patches

# parameters
steps = 100
c = 3e8  #Speed of light
R = 1  #Size radiator [m]
d = 1.1  #Distance to interaction point [m]
D = 20e-3  #Radiator thickness [m]
d_fel = 0  #Focusing optics length [m]
theta_max = 35 * math.pi / 180
theta_min = 5 * math.pi / 180

pmin = 0.5
pmax = 1.5

mmu = 0.1056583755
mpi = 0.13957039
mK = 0.493677

wlen_min = 200
wlen_max = 900

C = [0.6961663, 0.0684043, 0.4079426, 0.1162414, 0.8974794,
     9.896161]  #Sellmeier coefficients


##calcuation of refractive index
def rindex(wlen):
  wlen = wlen / 1000
  return np.sqrt(1 + C[0] * wlen**2 / (wlen**2 - C[1]**2) + C[2] * wlen**2 /
                 (wlen**2 - C[3]**2) + C[4] * wlen**2 / (wlen**2 - C[5]**2))


# derivation of n with respect to lambda
def drindex(wlen):
  wlen = wlen / 1000

  return ((2 * C[4] * wlen) / (wlen**2 - C[5]**2) - (2 * C[4] * wlen**3) /
          (wlen**2 - C[5]**2)**2 + (2 * C[2] * wlen) / (wlen**2 - C[3]**2) -
          (2 * C[2] * wlen**3) / (wlen**2 - C[3]**2)**2 + (2 * C[0] * wlen) /
          (wlen**2 - C[1]**2) - (2 * C[0] * wlen**3) /
          (wlen**2 - C[1]**2)**2) / (2 * np.sqrt((C[4] * wlen**2) /
                                                 (wlen**2 - C[5]**2) +
                                                 (C[2] * wlen**2) /
                                                 (wlen**2 - C[3]**2) +
                                                 (C[0] * wlen**2) /
                                                 (wlen**2 - C[1]**2) + 1))


# calculating the group velocity
def group(wlen):
  dn = drindex(wlen)

  return c / (rindex(wlen) - dn * wlen / 1000)


# calculating the beta factor
def beta(m, p):
  return p / math.sqrt(m**2 + p**2)


# cherenkov angle
def cherenkov(beta, n):
  try:
    return math.acos(1 / (n * beta))
  except ValueError:
    return 0


# photon path in radiator
def distance(R, d, theta, D):
  return R - d * math.tan(theta) - D * math.tan(theta) + d_fel


# time of propagation
def top(m, p, s, wlen, theta):
  phi = math.pi / 2 - (theta + cherenkov(beta(m, p), rindex(wlen)))
  return s / group(wlen) / math.cos(phi)


def franck_tamm(wlen, n, m, p):
  alpha = 1 / 137.
  z = 1
  return 2 * math.pi * alpha * z**2 * (
      1 / wlen**2 - 1 / (rindex(wlen)**2 * beta(m, p)**2 * wlen**2)) * 1e9


def photons(wlen_min, wlen_max, m, p):
  n = quad(lambda x: franck_tamm(x, rindex(x), mpi, 1), wlen_min, wlen_max)[0]
  return int(n)


def rand_wlen(m, p, wlen_min, wlen_max):
  while (True):
    x = rand.uniform(wlen_min, wlen_max)
    y = rand.random() * 1000
    if y < franck_tamm(x, rindex(x), m, p):
      return x


def avg_wlen(wlen_min, wlen_max, m, p):
  I = quad(lambda x: x * franck_tamm(x, rindex(x), mpi, 1), wlen_min, wlen_max)
  I0 = photons(wlen_min, wlen_max, mpi, 1)
  return I[0] / I0


def std_wlen(wlen_min, wlen_max, m, p):
  mu = avg_wlen(wlen_min, wlen_max, mpi, 1)
  I = quad(lambda x: (x - mu)**2 * franck_tamm(x, rindex(x), mpi, 1), wlen_min,
           wlen_max)
  I0 = photons(wlen_min, wlen_max, mpi, 1)
  std = math.sqrt(I[0] / I0)
  return std


def std_cherenkov(wlen_min, wlen_max, m, p):
  mu = avg_wlen(wlen_min, wlen_max, mpi, 1)
  std = std_wlen(wlen_min, wlen_max, mpi, 1)
  dd = -derivative(lambda x: cherenkov(beta(mpi, 1), rindex(x)), mu)
  sigma = dd * std
  print(std, mu, dd, sigma * 1000)
  return sigma


def std_top(wlen_min, wlen_max, m, p):
  mu = avg_wlen(wlen_min, wlen_max, mpi, 1)
  I = quad(lambda x: (x - mu)**2 * franck_tamm(x, rindex(x), mpi, 1), wlen_min,
           wlen_max)
  I0 = photons(wlen_min, wlen_max, mpi, 1)
  std = math.sqrt(I[0] / I0)
  theta = 0.5 * (theta_min + theta_max)
  dd = -derivative(lambda x: top(mpi, 1, R, x, theta), mu) * std
  return dd


# Highland formula function
def highland_theta(p, beta, z, dx, X0):
  """
    Parameters:
        p: momentum in MeV/c
        beta: v/c
        z: particle charge (in units of e)
        dx: step size in cm
        X0: radiation length of material in cm
    Returns:
        RMS scattering angle in milliradians
    """
  if dx <= 0 or X0 <= 0:
    return 0.0
  theta0_rad = (13.6 / (beta * p)) * z * np.sqrt(
      dx / X0) * (1 + 0.038 * np.log(dx / X0))
  return theta0_rad * 1000  # convert to mrad
