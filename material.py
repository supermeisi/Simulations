from functions import *

steps = 100

aa = []

for i in range(steps):
  bb = []
  for j in range(steps):
    n = 1.0 + 1.0 / steps * i
    p = 0.5 + 1.0 / steps * j

    sigma = 1 / 3 * (cherenkov(beta(mmu, p), n) - cherenkov(beta(mpi, p), n))

    #print(n, p, sigma)

    bb.append(sigma * 1000)
  aa.append(bb)

ticks = [0, 0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50, 100]

fig = plt.figure(figsize=(5, 4))
plt.imshow(aa,
           origin='lower',
           norm=colors.LogNorm(),
           cmap=cm.gist_rainbow,
           extent=[0.5, 1.5, 1, 2],
           aspect='auto')
plt.title('$\mu/\pi$ Separation')
ax = plt.gca()
rect = patches.Rectangle((0.5, 1.04),
                         1.0,
                         0.02,
                         linewidth=1,
                         edgecolor='k',
                         facecolor='none')
rect2 = patches.Rectangle((0.5, 1.46),
                          1.0,
                          0.02,
                          linewidth=1,
                          edgecolor='k',
                          facecolor='none')
rect3 = patches.Rectangle((0.5, 1.32),
                          1.0,
                          0.02,
                          linewidth=1,
                          edgecolor='k',
                          facecolor='none')
ax.add_patch(rect)
ax.add_patch(rect2)
ax.add_patch(rect3)
plt.text(1.2, 1.07, 'Aerogel')
plt.text(1.2, 1.49, 'Fused Silica')
plt.text(1.2, 1.35, 'Water')
plt.colorbar(label='$\sigma$ [mrad]', ticks=ticks, format='%.1f')
plt.xlabel('$p$ [GeV/c]')
plt.ylabel('$n$')
plt.tight_layout()
plt.savefig('res_mupi.png', dpi=1200)
plt.savefig('res_mupi.eps')
plt.savefig('res_mupi.pdf')

aa = []

for i in range(steps):
  bb = []
  for j in range(steps):
    n = 1.0 + 1.0 / steps * i
    p = 1.5 + 4.5 / steps * j

    sigma = 1 / 3 * (cherenkov(beta(mpi, p), n) - cherenkov(beta(mK, p), n))

    #print(n, p, sigma)

    bb.append(sigma * 1000)
  aa.append(bb)

ticks = [0, 0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50, 100]

fig = plt.figure(figsize=(5, 4))
plt.imshow(aa,
           origin='lower',
           norm=colors.LogNorm(),
           cmap=cm.gist_rainbow,
           extent=[1.5, 6.0, 1, 2],
           aspect='auto')
plt.title('$\pi/K$ Separation')
ax = plt.gca()
rect = patches.Rectangle((1.5, 1.04),
                         4.5,
                         0.02,
                         linewidth=1,
                         edgecolor='k',
                         facecolor='none')
rect2 = patches.Rectangle((1.5, 1.46),
                          4.5,
                          0.02,
                          linewidth=1,
                          edgecolor='k',
                          facecolor='none')
rect3 = patches.Rectangle((1.5, 1.32),
                          4.5,
                          0.02,
                          linewidth=1,
                          edgecolor='k',
                          facecolor='none')
ax.add_patch(rect)
ax.add_patch(rect2)
ax.add_patch(rect3)
plt.text(2.0, 1.07, 'Aerogel')
plt.text(2.0, 1.49, 'Fused Silica')
plt.text(2.0, 1.35, 'Water')
plt.colorbar(label='$\sigma$ [mrad]', ticks=ticks, format='%.1f')
plt.xlabel('$p$ [GeV/c]')
plt.ylabel('$n$')
plt.tight_layout()
plt.savefig('res_piK.png', dpi=1200)
plt.savefig('res_piK.eps')
plt.savefig('res_piK.pdf')
