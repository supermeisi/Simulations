from functions import *

steps1 = 200
steps2 = 200

dD = D/steps1
p = 1
theta = 0.5*(theta_max + theta_min)
aa = []
bb = []
wlen_avg = avg_wlen(wlen_min, wlen_max, mpi, 1)
for i in range(steps1):
  x = dD*i
  for j in range(steps2):
    #wlen = (wlen_max - wlen_min)/steps*j + wlen_min
    wlen = rand_wlen(mpi, p, wlen_min, wlen_max)
    dist1 = distance(R, d, theta, x)
    dist2 = distance(R, d, theta, D/2)
    bb.append((cherenkov(beta(mpi, p), rindex(wlen)) - cherenkov(beta(mpi,p),rindex(wlen_avg)))*1e3)
    t1 = top(mpi, p, dist1, wlen, theta)
    t2 = top(mpi, p, dist2, wlen_avg, theta)
    aa.append((t1-t2)*1e12)

fig = plt.figure(figsize=(5, 4))

plt.hist2d(aa, bb, bins=50, cmap=cm.rainbow)

plt.title('Time of Arrival')
plt.colorbar(label='# Entries')
plt.xlabel("$\Delta t$ [ps]")
plt.ylabel("$\Delta \\theta_C$ [mrad]")
fig.tight_layout()
plt.savefig('time.pdf')
plt.savefig('time.png')
plt.close()