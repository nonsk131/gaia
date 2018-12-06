import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import matplotlib.colors as colors

def pad_withNan(p):
    ii,jj = np.where(p < 1)
    for i, j in zip(ii, jj):
        p[i,j] = np.nan
    return p

data = np.loadtxt('/mnt/home/npanithanpaisal/gaia/mock.txt')

xbin = np.linspace(-1,5,1000)
ybin = np.linspace(-5,17,1000)

H, e1, e2 = np.histogram2d(data[:,4], data[:,3], bins=(xbin, ybin))
H = H.T
#H = pad_withNan(H)
fig = plt.figure(figsize=((10,8)))
ax = fig.add_subplot(1,1,1)
X, Y = np.meshgrid(e1, e2)
im = ax.pcolormesh(X, Y, H, norm=colors.LogNorm(vmin=H.min(), vmax=H.max()), cmap='gist_heat')
ax.set_ylim(-5,17)
plt.gca().invert_yaxis()
ax.set_xlabel('Bp - Rp')
ax.set_ylabel('M_G')
ax.set_title('With Extinction')
fig.colorbar(im, ax =ax)
fig.savefig('/mnt/home/npanithanpaisal/gaia/HR_extinction_log10.png', dpi=300)



H, e1, e2 = np.histogram2d(data[:,6], data[:,5], bins=(xbin, ybin))
H = H.T
#H = pad_withNan(H)
fig = plt.figure(figsize=((10,8)))
ax = fig.add_subplot(1,1,1)
X, Y = np.meshgrid(e1, e2)
im = ax.pcolormesh(X, Y, H, norm=colors.LogNorm(vmin=H.min(), vmax=H.max()), cmap='gist_heat')
ax.set_ylim(-5,17)
plt.gca().invert_yaxis()
ax.set_xlabel('Bp - Rp')
ax.set_ylabel('M_G')
ax.set_title('Without Extinction')
fig.colorbar(im, ax =ax)
fig.savefig('/mnt/home/npanithanpaisal/gaia/HR_no_extinction_log10.png', dpi=300)
