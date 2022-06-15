import numpy as np

def D_h5data(wr, wi, k0, h5data):
    '''
    landau's dielectric function (eqn 8.4.8 Krall and Trivelpiece)
    numerically integrate to compute, deforming the contour
    accepts a distribution function defined by an array of values
    h5data is of type H5Data from PyVisOS
    normalizing constants are wp and c:
        * wr and wi are in units of wp
        * k0 is k of the EPW in units of wp/c
        * v is in units of c
    '''
    # collect some values of use
    vphi = wr/k0

    dfdv = h5data.values
    dv = h5data.axes[-1].increment
    v_ = np.arange( h5data.axes[-1].min, h5data.axes[-1].max, dv )

    vphi_idx = np.where(v_ < vphi)[0][-1]

    # points along the path of integration
    v0 = v_[0]
    v1 = v_[vphi_idx]
    v2 = v_[vphi_idx+1]
    v3 = v_[-1]

    if wi <= 0:
        # first part of the real axis
        v = v_[0:vphi_idx]
        num = len(v)
        integral = np.trapz( dfdv[0:num] / (v-(wr+1j*wi)/k0), v )

        # deforming contour under pole
        # we need to define an analytic version so we can eval in the complex plane
        # approximate this little bit as just a line
        # using a cannon to shoot a pigeon, or whatever the phrase is
        x = np.linspace(v1,v2,2)
        dfdv_an = np.poly1d(np.polyfit(x,dfdv[vphi_idx:vphi_idx+2],1))

        v_r, v_i = v1, np.linspace(0,wi/k0 - 1e-3,1024)
        v = v_r + 1j*v_i
        integral += np.trapz( dfdv_an(v) / (v-(wr+1j*wi)/k0), 1j*v_i )

        v_r, v_i = np.linspace(v1, v2, 1024), wi/k0-1e-3
        v = v_r + 1j*v_i
        integral += np.trapz( dfdv_an(v) / (v-(wr+1j*wi)/k0), v_r )

        v_r, v_i = v2, np.linspace(wi/k0-1e-3,0,1024)
        v = v_r + 1j*v_i
        integral += np.trapz( dfdv_an(v) / (v-(wr+1j*wi)/k0), 1j*v_i )

        # rest of the real axis
        v = v_[vphi_idx+1:]
        num = len(v)
        integral += np.trapz( dfdv[-num:] / (v-(wr+1j*wi)/k0), v )

    else:
        # all of real axis
        v = v_
        num = len(v)
        integral = np.trapz( dfdv[0:num] / (v-(wr+1j*wi)/k0), v )

    d = 1.-1./(k0**2) * integral

    return d


def get_D(Wr, Wi, k0, h5data):
    '''
    Compute the value of the dielectric in a window of complex omega space for given k0
    and distribution function
    '''
    D_arr = np.zeros((len(Wi), len(Wr)), dtype = complex)
    for i in range(len(Wi)):
        for j in range(len(Wr)):
            D_arr[i,j] = D_h5data(Wr[j], Wi[i], k0=k0, h5data=h5data)
            
    return D_arr
