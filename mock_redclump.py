import numpy as np
from astropy.coordinates import SkyCoord,Galactocentric
import astropy.units as u
import astropy.coordinates as coord
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
plt.switch_backend('agg')

def pick_clump(data, minmag=1, maxmag=1.4, gmin=0, gmax=1):
    return data[np.where((data[:,5]>gmin) & (data[:,5]<gmax) & (data[:,6] > minmag) & (data[:,6] < maxmag))]

def dispersion(data):
    a = data[:,3] - np.mean(data[:,3])
    b = data[:,4] - np.mean(data[:,4])
    s = (a**2 + b**2).sum()
    return np.sqrt(s/len(a))

def discrepancies(data):
    a = data[:,3]-data[:,5]
    b = data[:,4]-data[:,6]
    s = (a**2 + b**2).sum()
    return np.sqrt(s/len(a))

# data = np.loadtxt('/mnt/home/npanithanpaisal/gaia/mock.txt')
# data = pick_clump(data)
# c = coord.ICRS(ra=data[:,0] * u.degree,
#                 dec=data[:,1] * u.degree,
#                 distance=data[:,2] * u.kpc)
#
# gcentric = c.transform_to(coord.Galactocentric)
# gcentric.representation = 'cylindrical'
#
# # pick out stars within the theta angle
# n = np.append(np.where((gcentric.phi > 175*u.deg))[0], np.where((gcentric.phi < -175*u.deg))[0])
# gcentric = gcentric[n]
# data = data[n]
#
# np.savetxt('/mnt/home/npanithanpaisal/gaia/mock_175cut.txt', data)

data = np.loadtxt('/mnt/home/npanithanpaisal/gaia/mock_175cut_redclump.txt')
c = coord.ICRS(ra=data[:,0] * u.degree,
                dec=data[:,1] * u.degree,
                distance=data[:,2] * u.kpc)

gcentric = c.transform_to(coord.Galactocentric)
gcentric.representation = 'cylindrical'

r_ensemble = np.array([0.1,0.5,1,1.5,2])*u.kpc
#r_ensemble = np.linspace(0.1, 3, 31)*u.kpc
z_ensemble = np.array([2])*u.kpc
dis_array = np.zeros(len(r_ensemble))
i = 0
j = 0
fig = plt.figure(figsize=((10,8)))
ax = fig.add_subplot(1,1,1)
for r_span in r_ensemble:
    j = 0
    #for z_span in z_ensemble:
    z_span=2*u.kpc
    rmin = 8.3*u.kpc-r_span
    rmax = 8.3*u.kpc+r_span
    n = np.where((gcentric.z < z_span) & (gcentric.z > -z_span) & (gcentric.rho > rmin) & (gcentric.rho < rmax))[0]
    print len(n)
    data_cut = data[n]
    dis = discrepancies(data_cut)
    dis_array[j] = dis
    j += 1
ax.plot(r_ensemble, dis_array, linewidth=3)


#ax.plot([0.1,3.5],[0.05, 0.05], linestyle='dashed', linewidth=2)
ax.set_xlabel('r [kpc]')
ax.set_ylabel('dispersion')
#ax.legend(loc='best')
fig.savefig('/mnt/home/npanithanpaisal/gaia/dispersion2.png', dpi=300)
np.savetxt('/mnt/home/npanithanpaisal/gaia/dispersion.txt', dis_array)
