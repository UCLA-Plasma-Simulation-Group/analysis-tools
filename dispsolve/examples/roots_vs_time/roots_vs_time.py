import numpy as np
import copy
import matplotlib
import matplotlib.pyplot as plt
import glob
from scipy import ndimage

# install PyVisOS at https://github.com/UCLA-Plasma-Simulation-Group/pyVisOS
import osh5io
import osh5def
import osh5vis
import osh5utils

import sys
sys.path.insert(0,'../../')
from dispsolve import D_h5data, get_D

font = {'size':18} # 'family' : 'normal', 'weight' : 'bold'
matplotlib.rc('font', **font)

# analytic df/dv of maxwellian dist
def dfMdv(v ,vth):
    return 1./np.sqrt(2*np.pi) * (-v/vth**3) * np.exp(-0.5*(v/vth)**2)


# bohm gross dispersion relation, returns omega/omega_p as func of k * lambda debeye
def omega_bohm_gross(kld):
    return np.sqrt( 1 + 3*kld**2 )


def get_f( h5data ):
    # normalize the distribution func
    norm = np.trapz( h5data.values, np.linspace(h5data.axes[-1].min, h5data.axes[-1].max, h5data.axes[-1].size) )
    h5data.values = h5data.values / norm
    return h5data


def get_dfdv( h5data, sigma=1 ):
    h5data = get_f(h5data)
    dv = h5data.axes[-1].increment
    h5data.values = ndimage.gaussian_filter1d(h5data.values, sigma=sigma, order=1, mode='wrap') / dv
    return h5data


def roots_vs_time():
    # object to facilitate storing the solution for a bunch of k's
    class Fred():
        def __init__(self, time, Wr_vals, Wi_vals, D_vals, D_min):
            self.time = time
            self.Wr_vals = Wr_vals
            self.Wi_vals = Wi_vals
            self.D_vals = D_vals
            self.D_min = D_min

    solns = np.array([])
    for file_idx in range(int(len(file_list))):
        file_path = file_list[file_idx]
        h5data = osh5io.read_h5(file_path)

        time = h5data.run_attrs['TIME'][0]
        time_units = h5data.run_attrs['TIME UNITS']

        print('time {}'.format(time))

        # get a smoothed derivative of the data, dfdv
        dfdv_data = get_dfdv(h5data, sigma=20)

        D_arr = get_D(Wr,Wi,k0,dfdv_data)

        # get the locations of the poles
        # pole_idx = np.where( np.abs(D_arr) < .01 )
        pole_idx = np.where( np.abs(D_arr) < 1.1*np.min(np.abs(D_arr)) )
        # pole_idx = np.where( np.abs(D_arr) == np.min(np.abs(D_arr)) )
        # pole_idx = np.argpartition(np.abs(D_arr.flatten()),5)

        D_min = np.min(np.abs(D_arr))
        print('min {}'.format(D_min))

        # debug
        D_vals = np.array([])
        for i in range(len(pole_idx[0])):
            p0 = pole_idx[0][i]
            p1 = pole_idx[1][i]
            D_vals = np.append( D_vals, np.abs(D_arr[p0, p1]) )
            print( "pole {}, wr {:.3f}, wi {:.4f}, D {}".format(i, Wr[p1], Wi[p0], D_vals[-1]) )

        Wr_vals = Wr[pole_idx[1]]
        Wi_vals = Wi[pole_idx[0]]

        solns = np.append( solns, Fred(time, Wr_vals, Wi_vals, D_vals, D_min) )

    # save the solution so we don't have to recalculate every time we want to edit the plot
    np.save('roots_vs_time.npy',solns)

    return solns


# -------------
# -------------
# get the data
# file_list = sorted(glob.glob('../some_data/MS/PHA/p1/electrons/*.h5'))[::30]
file_list = sorted(glob.glob('../some_other_data/MS/PHA/p1/electrons/*.h5'))[::10]

# set the parameters, you either with the pro ballers or the amateurs
vth = .063
k0 = 2.75

# Note: here we're basically looking for the resonance of plasma waves with negative vph
# For me, the more natural way to do this would be to set k negative and look for positive
# Wr values, but for some reason D_h5data does not like this...Still not sure why.
# This way seems to work though and is equivalent
Wr = np.linspace(-1,-1.2,128) # for this k, 2.75, BG disp rel predicts wr=1.03
Wi = np.linspace(-0.4,0,128) # also, we probably only need to look for damping roots


# -------------
# plot f, dfdv, and theory dfdv as a sanity check
# -------------
fig, ax = plt.subplots(1,3,figsize=(16,8))#,squeeze=False)

h5data = osh5io.read_h5(file_list[0])
f_data = get_f( copy.deepcopy(h5data) )
osh5vis.osplot1d( f_data, ax=ax[0], ylim=[-.5,7], xlim=[-1,1] )
ax[0].axvline(x=-.38, linestyle='--', color='r', label='rescatter EPW $v_{ph}$')

h5data = osh5io.read_h5(file_list[6])
f_data = get_f( copy.deepcopy(h5data) )
osh5vis.osplot1d( f_data, ax=ax[0], ylim=[-.5,7], xlim=[-1,1] )
ax[0].axvline(x=-.38, linestyle='--', color='r', label='rescatter EPW $v_{ph}$')

# vvv = np.linspace(-8*vth,8*vth,len(dfdv_data.values))
# ax[1].plot( vvv, dfMdv(vvv,vth), '--', label='analytic')
# ax[1].legend()

fig.savefig('f_dfdv.png')


# # -------------
# # generate solutions as a a function of time
# # -------------
# solns = roots_vs_time()
# # altenatively just load the solutions from disk if you've already calculated them
# # and just need to edit the plot
# solns = np.load('roots_vs_time.npy',allow_pickle=True)

# # obtain time units for x axis
# file_path = file_list[0]
# h5data = osh5io.read_h5(file_path)
# time_units = h5data.run_attrs['TIME UNITS']

# # setup plot
# plt.figure(figsize=(20,10))
# ax2 = plt.subplot(1,1,1)
# ax2.set_xlabel('$t\ [{}]$'.format(time_units))
# ax2.set_ylabel(r'$\omega_{im} / \omega_p$')
# ax2.grid()

# # plot only the root that mimized D
# x = np.array([])
# y = np.array([])
# for i in solns:
#     time = i.time
#     Wr_vals = i.Wr_vals
#     Wi_vals = i.Wi_vals
#     D_vals = i.D_vals
#     D_min = i.D_min

#     for k in range(len(Wi_vals)):
#         if D_vals[k]==D_min:
#             x = np.append(x,time)
#             y = np.append(y,Wi_vals[k])

# ax2.scatter(x,y)
# ax2.plot(x,y,'--')

# ax2.legend()

# plt.tight_layout()
# plt.savefig( 'roots_vs_time.png' )
