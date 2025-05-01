from functions import franck_tamm, rand_wlen, rindex, mpi, wlen_min, wlen_max

import matplotlib.pyplot as plt

aa = []

print(franck_tamm(300, rindex(300), mpi, 1))
for _ in range(1000):
  aa.append(rand_wlen(mpi, 1, wlen_min, wlen_max))

fig = plt.figure(figsize=(5, 4))
plt.hist(aa, bins=25)
plt.xlabel('$\lambda$ [nm]')
plt.ylabel('Entries')
plt.xlim(200, 900)
plt.tight_layout()
plt.savefig('wavelengths.pdf')
plt.savefig('wavelengths.png', dpi=600)
plt.close()
