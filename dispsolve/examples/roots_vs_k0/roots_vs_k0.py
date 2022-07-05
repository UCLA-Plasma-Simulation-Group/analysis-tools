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


# object to facilitate storing the solution for a bunch of k's
class Fred():
    def __init__(self, k0, Wr_vals, Wi_vals, D_vals, D_min):
        self.k0 = k0
        self.Wr_vals = Wr_vals
        self.Wi_vals = Wi_vals
        self.D_vals = D_vals
        self.D_min = D_min


def roots_vs_k0():
    solns = np.array([])
    for k0 in k0_vals:
        print('k0 {}'.format(k0))

        D_arr = get_D(Wr,Wi,k0,dfdv_data)

        # save the locations of the best candidates for the poles (ie solns within some tolerance)
        pole_idx = np.where( np.abs(D_arr) < 1.5*np.min(np.abs(D_arr)) )

        # the closer this is to zero, the more confident you can be that this is actually a solution
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

        solns = np.append( solns, Fred(k0, Wr_vals, Wi_vals, D_vals, D_min) )

    # save the solution so we don't have to recalculate every time we want to edit the plot
    np.save('roots_vs_k0.npy',solns)

    return solns


# -------------
# -------------
# get the data
file_list = sorted(glob.glob('../some_data/MS/PHA/p1/electrons/*.h5'))
file_idx = 0 # at the first timestep, the dist fn is pretty much maxwellian, we'll use that for this example
file_path = file_list[file_idx]
h5data = osh5io.read_h5(file_path)

# normalize f
f_data = get_f( copy.deepcopy(h5data) )
f_data.data_attrs['LONG_NAME'] = 'f_e(p_1)'

# get a smoothed derivative of the data, dfdv
dfdv_data = get_dfdv( copy.deepcopy(h5data), sigma=20)
dfdv_data.data_attrs['LONG_NAME'] = '\\frac{df_e(p_1)}{dp_1}'

# set the parameters, you either with the pro ballers or the amateurs
vth = .063                       # for the data grabbed above, vth was .063/c
k0_vals = np.linspace(2.5,6,10)  # limited by the range of momenta for which f is defined in your data
Wr = np.linspace(1,2,128)        # pick this range large enough s.t. resonant mode lies within it, but ideally no larger since solver is somewhat unstable
Wi = np.linspace(-0.4,0,128)     # we only look for damping roots


# -------------
# plot f, dfdv, and theory dfdv as a sanity check
# -------------
fig, ax = plt.subplots(1,2,figsize=(16,8))#,squeeze=False)

osh5vis.osplot1d( f_data, ax=ax[0], ylim=[-.5,7], xlim=[-1,1] )
ax[0].axvline(x=-.38, linestyle='--', color='r', label='rescatter EPW $v_{ph}$')

osh5vis.osplot1d( dfdv_data, ax=ax[1], ylim=[-65,65], xlim=[-1,1] )
ax[1].axvline(x=-.38, linestyle='--', color='r', label='rescatter EPW $v_{ph}$')
ax[1].legend()

vvv = np.linspace(-8*vth,8*vth,len(dfdv_data.values))
ax[1].plot( vvv, dfMdv(vvv,vth), '--', label='analytic')
ax[1].legend()

fig.savefig('f_dfdv.png'.format(file_idx))


# -------------
# compute the value of the dielectric for a particular k0
# as you'll see from the output, the solver has trouble with certain wr-wi pairs, but still
# it seems to be able to predict the roots correctly
# -------------
D_arr = get_D(Wr,Wi,k0_vals[-1],dfdv_data)

# get locations of poles (ie roots, ie values of wr and wi that make D=0)
pole_idx = np.where( np.abs(D_arr) < 1.5*np.min(np.abs(D_arr)) )

# plot the dielectric
plt.figure(figsize=(30,10))

plt.subplot(1,3,1)
plt.title('D Re')
plt.imshow(D_arr.real, origin='lower', extent=[Wr[0], Wr[-1],Wi[0], Wi[-1]], aspect='auto')
plt.xlabel(r'$\omega_R$')
plt.ylabel(r'$\omega_I$')
plt.colorbar()

plt.subplot(1,3,2)
plt.title('D Im')
plt.imshow(D_arr.imag, origin='lower', extent=[Wr[0], Wr[-1],Wi[0], Wi[-1]], aspect='auto')
plt.xlabel(r'$\omega_R$')
plt.ylabel(r'$\omega_I$')
plt.colorbar()

plt.subplot(1,3,3)
plt.title('|D|')
plt.imshow(np.abs(D_arr), origin='lower', extent=[Wr[0], Wr[-1],Wi[0], Wi[-1]], aspect='auto',cmap = 'hot_r')
if len(pole_idx[0])>0:
    plt.scatter(Wr[pole_idx[1]],Wi[pole_idx[0]], c='k')
plt.colorbar()
plt.xlabel(r'$\omega_R$')
plt.ylabel(r'$\omega_I$')

plt.tight_layout()

plt.savefig('dielectric.png',dpi=200)


# -------------
# generate solutions as a function of k0
# -------------
solns = roots_vs_k0()
# altenatively just load the solutions from disk if you've already calculated them
# and just need to edit the plot
solns = np.load('roots_vs_k0.npy',allow_pickle=True)

plt.figure(figsize=(20,10))
ax1 = plt.subplot(1,2,1)
ax1.set_xlabel('$k \lambda_{De}$')
ax1.set_ylabel(r'$\omega_{re}/\omega_p$')
ax1.grid()
ax2 = plt.subplot(1,2,2)
ax2.set_xlabel('$k \lambda_{De}$')
ax2.set_ylabel(r'$\omega_{im} / \omega_p$')
ax2.grid()

# plot only the root that mimized D
x = np.array([])
y = np.array([])
z = np.array([])
for i in solns:
    k0 = i.k0
    Wr_vals = i.Wr_vals
    Wi_vals = i.Wi_vals
    D_vals = i.D_vals
    D_min = i.D_min

    for k in range(len(Wi_vals)):
        # filter out solns not within some tolerance
        # check out how the plot changes if you don't do this
        if D_vals[k]<1.1*D_min:
            x = np.append(x,k0*vth)
            y = np.append(y,Wi_vals[k])
            z = np.append(z,Wr_vals[k])

ax1.scatter(x,z,label='kinetic dispersion')
ax1.plot(x,omega_bohm_gross(x),'--',color='orange',label='Bohm-Gross dispersion')
ax1.legend()

ax2.scatter(x,y)

plt.tight_layout()
plt.savefig( 'roots_vs_k0.png' )


