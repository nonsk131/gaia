import numpy as np
from astropy.coordinates import SkyCoord,Galactocentric
import astropy.units as u
import astropy.coordinates as coord
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
plt.switch_backend('agg')

#Mg,bp_rp,ra,dec,parallax,phot_g_mean_mag
def pick_mstar(data, minmag=2.5, maxmag=3, gmin=5)
    return np.where((dat[:,0]>gmin) & (dat[:,1] > minmag) & (dat[:,1] < maxmag))

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
    H_new = pad_withNan(H)

    fig = plt.figure(figsize=((12,8)))
    ax = fig.add_subplot(1,1,1)
    X, Y = np.meshgrid(e1, e2)
    im = ax.pcolormesh(X, Y, H_new, cmap='jet')
    fig.colorbar(im, ax =ax)
    fig.savefig('/mnt/home/npanithanpaisal/gaia/star_dist.png', dpi=300)


data = np.loadtxt('/mnt/home/npanithanpaisal/gaia/real_175cut.txt')

data_m = pick_mstar(data)
print 'there are {} m-stars'.format(len(data_m))
c = coord.ICRS(ra=data[:,2] * u.degree,
            dec=data[:,3] * u.degree,
            distance=(1./data[:,4]) * u.kpc)

gcentric = c.transform_to(coord.Galactocentric)
gcentric.representation = 'cylindrical'
get_hist(gcentric)
