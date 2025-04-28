from functions import *

aa = np.arange(wlen_min, wlen_max, 1)
fig = plt.figure(figsize=(5, 4))
plt.plot(aa, group(aa))
plt.ylabel("$v_G$ [m/s]")
plt.xlabel("$\lambda$ [nm]")
plt.title('Group Velocity')
fig.tight_layout()
plt.savefig('group.pdf')
plt.savefig('group.png')
plt.close()