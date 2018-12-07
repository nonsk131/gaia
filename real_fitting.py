import numpy as np
from astropy.coordinates import SkyCoord,Galactocentric
import astropy.units as u
import astropy.coordinates as coord
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
plt.switch_backend('agg')

#Mg,bp_rp,ra,dec,parallax,phot_g_mean_mag
def pick_mstar(data, minmag=2.5, maxmag=3, gmin=5):
    return data[np.where((data[:,0]>gmin) & (data[:,1] > minmag) & (data[:,1] < maxmag))]

def get_bin_edges():
    r_edges = np.linspace(6.5,10.1,37)*u.kpc
    z_edges = np.linspace(-2.5,2.5,101)*u.kpc
    return r_edges, z_edges

def pad_withNan(p):
    ii,jj = np.where(p < 5)
    for i, j in zip(ii, jj):
        p[i,j] = np.nan
    return p

def get_hist(gcentric):
    rr, zz = get_bin_edges()
    H, e1, e2 = np.histogram2d(gcentric.rho, gcentric.z, bins=(rr, zz))
    H = H.T
    print H
    #H_new = pad_withNan(H)
    #print H_new
    fig = plt.figure(figsize=((12,8)))
    ax = fig.add_subplot(1,1,1)
    X, Y = np.meshgrid(e1, e2)
    im = ax.pcolormesh(X, Y, H, cmap='jet')
    fig.colorbar(im, ax =ax)
    fig.savefig('/mnt/home/npanithanpaisal/gaia/star_dist.png', dpi=300)

    l = (len(rr)-1) * (len(zz)-1)
    r_array = np.zeros(l)
    z_array = np.zeros(l)
    val = np.zeros(l)
    r_center = 0.5*(rr[1:] + rr[:-1])
    z_center = 0.5*(zz[1:] + zz[:-1])
    k = 0
    for i, r in enumerate(r_center):
        for j, z in enumerate(z_center):
            r_array[k] = r
            z_array[k] = z
            val[k] = H
            k += 1

    return r_array, z_array, val


def get_volume(rr, zz):
    height = zz[1] - zz[0]
    v = np.zeros(len(rr)-1)
    i = 0
    for r1, r2 in zip(rr[:-1], rr[1:]):
        v[i] = (1/36.) * np.pi * (r2**2 - r1**2) * height
    return v


def fit_func(X, rho, f, l1, h1, l2, h2):
    R, z = X
    return rho* ((np.exp(-R/l1)* np.exp(-z/h1)) + f * (np.exp(-R/l2)* np.exp(-z/h2)))



data = np.loadtxt('/mnt/home/npanithanpaisal/gaia/real_175cut.txt')

data_m = pick_mstar(data)
print 'there are {} m-stars'.format(len(data_m))
c = coord.ICRS(ra=data[:,2] * u.degree,
            dec=data[:,3] * u.degree,
            distance=(1./data[:,4]) * u.kpc)

gcentric = c.transform_to(coord.Galactocentric)
gcentric.representation = 'cylindrical'
get_hist(gcentric)
