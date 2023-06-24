import numpy as np
from scipy.fft import fftn, ifftn, fftshift
from params import *

def poisson_solver_fft3(dens):
    '''
    Description: solve 3D Poisson equation with FFT.
    Input: 
       dens: density field.
    Return: 
       Phi: gravitational potential field.
    '''
    G = 1
    L = len(dens)
    g_k = np.concatenate((np.mgrid[:L//2,:L,:L][iz], np.mgrid[:L//2,:L,:L][iz] - L/2), axis=0) ** 2 + \
          np.concatenate((np.mgrid[:L,:L//2,:L][iy], np.mgrid[:L,:L//2,:L][iy] - L/2), axis=1) ** 2 + \
          np.concatenate((np.mgrid[:L,:L,:L//2][ix], np.mgrid[:L,:L,:L//2][ix] - L/2), axis=2) ** 2
    
    Phi_k = - fftn(dens) / g_k * G * L**2 / np.pi
    Phi_k[0][0][0] = 0

    #Phi_k1 = fftshift(Phi_k)

    Phi = (ifftn(Phi_k).real)[:Ng, :Ng, :Ng]
    
    return Phi