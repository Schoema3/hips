import numpy as np
import matplotlib.pyplot as plt

d00 = np.loadtxt('Data_00000.dat')
d01 = np.loadtxt('Data_00001.dat')
d02 = np.loadtxt('Data_00002.dat')
d03 = np.loadtxt('Data_00003.dat')
d04 = np.loadtxt('Data_00004.dat')
d05 = np.loadtxt('Data_00005.dat')
d06 = np.loadtxt('Data_00006.dat')
d07 = np.loadtxt('Data_00007.dat')
d08 = np.loadtxt('Data_00008.dat')
d09 = np.loadtxt('Data_00009.dat')
d10 = np.loadtxt('Data_00010.dat')
d11 = np.loadtxt('Data_00011.dat')

dd = [d00, d01, d02, d03, d04, d05, d06, d07, d08, d09, d10, d11]

fig, ax = plt.subplots(1,1, figsize=(6,7))     # create a figure and axes

n = len(d00[:,0])
i = np.arange(len(d00[:,0]))
for k in range(len(dd)):
    ax.plot(i,dd[k][:,0]+3000*k, lw=1)
    ax.plot(i,np.full(n,300+3000*k), ':', lw=0.5, color='gray')

ax.set_ylim([0,12*3000])
ax.set_xlabel('parcel index')
ax.set_ylabel('T (K), shifted 3000 K for clarity')
ax.set_yticks([300, 3300, 6300, 9300, 12300, 15300, 18300, 21300, 24300, 27300, 30300, 33300])
ax.set_title('time increases from bottom to top')
plt.tight_layout()

plt.show()

