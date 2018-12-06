import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
plt.switch_backend('agg')

data = np.loadtxt('/mnt/home/npanithanpaisal/gaia/mock.txt')

xbin = np.linspace(-1,5,1000)
ybin = np.linspace(-5,17,1000)

H, e1, e2 = np.histogram2d(data[:,4], data[:,3], bins=(xbin, ybin))
H = H.T
fig = plt.figure(figsize=((10,8)))
ax = fig.add_subplot(1,1,1)
X, Y = np.meshgrid(e1, e2)
im = ax.pcolormesh(X, Y, H, cmap='gist_heat')
plt.gca().invert_yaxis()
ax.set_xlabel('Bp - Rp')
ax.set_ylabel('M_G')
ax.set_title('With Extinction')
fig.colorbar(im, ax =ax)
fig.savefig('/mnt/home/npanithanpaisal/gaia/HR_extinction.pdf')



H, e1, e2 = np.histogram2d(data[:,6], data[:,5], bins=(xbin, ybin))
H = H.T
fig = plt.figure(figsize=((10,8)))
ax = fig.add_subplot(1,1,1)
X, Y = np.meshgrid(e1, e2)
im = ax.pcolormesh(X, Y, Hß, cmap='gist_heat')
plt.gca().invert_yaxis()
ax.set_xlabel('Bp - Rp')
ax.set_ylabel('M_G')
ax.set_title('With Extinction')
fig.colorbar(im, ax =ax)
fig.savefig('/mnt/home/npanithanpaisal/gaia/HR_no_extinction.pdf')
