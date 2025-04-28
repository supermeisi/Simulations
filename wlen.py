from functions import *

aa = []
print(franck_tamm(300, rindex(300), mpi, 1))
for i in range(10000):
  aa.append(rand_wlen(mpi, 1, wlen_min, wlen_max))

fig = plt.figure(figsize=(5, 4))
plt.hist(aa, bins=25)
plt.xlabel('$\lambda$ [nm]')
plt.ylabel('Entries')
fig.tight_layout()
plt.savefig('wavelengths.pdf')
plt.savefig('wavelengths.png')
plt.close()
