import numpy as np
from params import *

def solve_mg3():
    '''Solve Poisson equation '''
    pass

def poisson_solver_fft3(dens):
    '''
    Description: solve 3D Poisson equation with FFT.
    Input: 
       dens: density field.
    Return: 
       Phi: gravitational potential field.
    '''
    g_k = np.concatenate((np.mgrid[:Ng//2,:Ng,:Ng][iz], np.mgrid[:Ng//2,:Ng,:Ng][iz] - Ng//2), axis=0) ** 2 + \
          np.concatenate((np.mgrid[:Ng,:Ng//2,:Ng][iy], np.mgrid[:Ng,:Ng//2,:Ng][iy] - Ng//2), axis=1) ** 2 + \
          np.concatenate((np.mgrid[:Ng,:Ng,:Ng//2][ix], np.mgrid[:Ng,:Ng,:Ng//2][ix] - Ng//2), axis=2) ** 2
    
    Phi_k = - np.fft.fftn(dens) / g_k * G * Ng**2 / np.pi
    Phi_k[0][0][0] = 0
    Phi = np.fft.ifftn(Phi_k).real
    
    return Phi