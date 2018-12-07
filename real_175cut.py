import numpy as np
from astropy.coordinates import SkyCoord,Galactocentric
import astropy.units as u
import astropy.coordinates as coord

#Mg,bp_rp,ra,dec,parallax,phot_g_mean_mag
def read_real_data():
    data1 = np.loadtxt('/mnt/home/npanithanpaisal/gaia/full_4kpc-result.csv', skiprows=1,delimiter=',')
    data2 = np.loadtxt('/mnt/home/npanithanpaisal/gaia/backside_2kpc-result.csv', skiprows=1,delimiter=',')
    data3 = np.loadtxt('/mnt/home/npanithanpaisal/gaia/frontside_100pc-result.csv', skiprows=1,delimiter=',')
    data = np.vstack((data1, data2, data3))
    c = coord.ICRS(ra=data[:,2] * u.degree,
                dec=data[:,3] * u.degree,
                distance=(1./data[:,4]) * u.kpc)

    gcentric = c.transform_to(coord.Galactocentric)
    gcentric.representation = 'cylindrical'

    n = np.append(np.where((gcentric.phi > 175*u.deg))[0], np.where((gcentric.phi < -175*u.deg))[0])
    return data[n]

data = read_real_data()
np.savetxt('/mnt/home/npanithanpaisal/gaia/real_175cut.txt', data)
