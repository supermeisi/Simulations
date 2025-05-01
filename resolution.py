from functions import *

step_width = 5

steps = int((wlen_max - wlen_min) / step_width)

arr = np.zeros((steps + 1, steps + 1))
arr2 = np.zeros((steps + 1, steps + 1))
arr3 = np.zeros((steps + 1, steps + 1))

theta = 0.5 * (theta_min + theta_max)

i1 = 0
for min_wlen in range(wlen_min, wlen_max + step_width, step_width):
           i2 = i1 + 1
           for max_wlen in range(min_wlen + step_width, wlen_max + step_width,
                                 step_width):
                      arr[i2, i1] = std_cherenkov(min_wlen, max_wlen, mpi,
                                                  1) * 1000
                      arr2[i2, i1] = std_top(min_wlen, max_wlen, mpi, 1) * 1e12
                      arr3[i2, i1] = photons(min_wlen, max_wlen, mpi,
                                             1) * 0.05 * 0.02
                      i2 += 1
           i1 += 1

print('finished')

arr[arr == 0] = np.nan
arr2[arr2 == 0] = np.nan
arr3[arr3 == 0] = np.nan

arr4 = arr / np.sqrt(arr3)
arr5 = arr2 / np.sqrt(arr3)

cmap_reversed = cm.get_cmap('gist_rainbow_r')

fig = plt.figure(figsize=(5, 4))
plt.imshow(arr,
           aspect='auto',
           cmap=cmap_reversed,
           origin='lower',
           extent=[wlen_min, wlen_max, wlen_min, wlen_max])
plt.title('Single Photon Resolution')
plt.colorbar(label='$\sigma_{SPR}$ [mrad]')
plt.xlabel("$\lambda_{min}$ [nm]")
plt.ylabel("$\lambda_{max}$ [nm]")
plt.tight_layout()
plt.savefig("single_photon.png")
plt.savefig("single_photon.pdf")
plt.close()

fig = plt.figure(figsize=(5, 4))
plt.imshow(arr3,
           aspect='auto',
           cmap=cmap_reversed,
           origin='lower',
           extent=[wlen_min, wlen_max, wlen_min, wlen_max])
plt.title('Photon Hits')
plt.colorbar(label='$N$')
plt.xlabel("$\lambda_{min}$ [nm]")
plt.ylabel("$\lambda_{max}$ [nm]")
plt.tight_layout()
plt.savefig("photons.png")
plt.savefig("photons.pdf")
plt.close()

fig = plt.figure(figsize=(5, 4))
plt.imshow(arr4,
           aspect='auto',
           cmap=cmap_reversed,
           origin='lower',
           extent=[wlen_min, wlen_max, wlen_min, wlen_max])
plt.title('Detector Resolution')
plt.colorbar(label='$\sigma_{SPR}/\sqrt{N}$ [mrad]')
plt.xlabel("$\lambda_{min}$ [nm]")
plt.ylabel("$\lambda_{max}$ [nm]")
plt.scatter(200, 500, s=50, facecolors='none', edgecolors='black')
plt.scatter(310, 580, s=50, facecolors='none', edgecolors='black')
plt.scatter(500, 800, s=50, facecolors='none', edgecolors='black')
plt.text(220, 500, 'Blue Cathode')
plt.text(330, 580, 'Green Cathode')
plt.text(520, 800, 'Red Cathode')
plt.tight_layout()
plt.savefig("detector.png")
plt.savefig("detector.pdf")
plt.close()

fig = plt.figure(figsize=(5, 4))
plt.imshow(arr5,
           aspect='auto',
           cmap=cmap_reversed,
           origin='lower',
           extent=[wlen_min, wlen_max, wlen_min, wlen_max])
plt.title('Time Difference')
plt.colorbar(label='$\sigma_{t}/\sqrt{N}$ [ps]')
plt.xlabel("$\lambda_{min}$ [nm]")
plt.ylabel("$\lambda_{max}$ [nm]")
plt.tight_layout()
plt.savefig("timediff.png")
plt.savefig("timediff.pdf")
plt.close()
