import numpy as np
from astropy.coordinates import SkyCoord,Galactocentric
import astropy.units as u
import astropy.coordinates as coord

def pick_clump(data, minmag=1, maxmag=1.4, gmin=0, gmax=1):
    return data[np.where((data[:,5]>gmin) & (data[:,5]<gmax) & (data[:,6] > minmag) & (data[:,6] < maxmag))]

def dispersion(data):
    a = data[:,3] - np.mean(data[:,3])
    b = data[:,4] - np.mean(data[:,4])
    s = (a**2 + b**2).sum()
    return np.sqrt(s/len(a))

data = np.loadtxt('/mnt/home/npanithanpaisal/gaia/mock.txt')
data = pick_clump(data)
c = coord.ICRS(ra=data[:,0] * u.degree,
                dec=data[:,1] * u.degree,
                distance=data[:,2] * u.kpc)

gcentric = c.transform_to(coord.Galactocentric)
gcentric.representation = 'cylindrical'

# pick out stars within the theta angle
n = np.append(np.where((gcentric.phi > 175*u.deg))[0], np.where((gcentric.phi < -175*u.deg))[0])
gcentric = gcentric[n]
data = data[n]

np.savetxt('/mnt/home/npanithanpaisal/gaia/mock_175cut.txt', data)
# data = np.loadtxt('/mnt/home/npanithanpaisal/gaia/mock_175cut.txt')
# c = coord.ICRS(ra=data[:,0] * u.degree,
#                 dec=data[:,1] * u.degree,
#                 distance=data[:,2] * u.kpc)
#
# gcentric = c.transform_to(coord.Galactocentric(galcen_v_sun = coord.CartesianDifferential((11.1, -232.24, 7.25)*u.km/u.s)))
# gcentric.representation = 'cylindrical'
#
# r_ensemble = np.linspace(0.5.3.5,13)
# z_ensemble = np.lin
# for r_span in np.linspace(0.5)
#     for z_span in (np.array(arange))
