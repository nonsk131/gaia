import numpy as np
import h5py

data = np.zeros((1,7))

for i in range(10):
    print('doing index {}'.format(i))
    loc = '/mnt/ceph/users/firesims/ananke/GaiaMocks/m12f/test/lsr-{}-rslice-0.m12f-res7100-md-sliced-gcat-dr2.hdf5'.format(i)
    #loc = '/mnt/ceph/users/firesims/ananke/GaiaMocks/m12f/lsr_0/lsr-0-rslice-{}.m12f-res7100-md-sliced-gcat-dr2.hdf5'.format(i)
    f = h5py.File(loc, 'r')

    # pick high latitude stars
    i_high = np.where((f['parallax'][:]>0.25) & (f['parallax_over_error'][:] > 10) & ((f['b'][:] > 20) | (f['b'][:] < -20)) )[0]

    # pick backside
    i_back = np.where((f['parallax'][:]>0.5) & (f['parallax_over_error'][:] > 10) & (f['b'][:] < 20) & (f['b'][:] > -20) & ((f['l'] >90) | (f['l'] < -90)))[0]

    # prick frontside
    i_front = np.where((f['parallax'][:]>10) & (f['parallax_over_error'][:] > 10) & (f['b'][:] < 20) & (f['b'][:] > -20) & (f['l'] <90) & (f['l'] > -90))[0]

    i_all = np.append(i_high, i_back)
    i_all = np.append(i_all, i_front)

    mock_ra = f['ra'][:][i_all]
    mock_dec = f['dec'][:][i_all]
    mock_parallax = f['parallax'][:][i_all]
    mock_bprp = f['bp_rp'][:][i_all]
    mock_bprp_int = f['bp_rp_int'][:][i_all]
    mock_phot_g_mean_mag = f['phot_g_mean_mag'][:][i_all]
    mock_phot_g_mean_mag_int = f['phot_g_mean_mag_int'][:][i_all]

    mock_dist = 1/mock_parallax # in kpc
    mock_G = mock_phot_g_mean_mag - 10 + 5*np.log10(mock_parallax)
    mock_G_int = mock_phot_g_mean_mag_int - 10 + 5*np.log10(mock_parallax)

    dat_i = np.column_stack((mock_ra, mock_dec, mock_dist, mock_G, mock_bprp, mock_G_int, mock_bprp_int))
    data = np.vstack(data, dat_i)

data = data[1:]
np.savetxt('/mnt/home/npanithanpaisal/gaia/mock.txt', data)
