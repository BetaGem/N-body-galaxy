import numpy as np
from params import *

def sph_sampler(func1, func2, args1, args2, N):
    '''
    Description: produce random points following the given distribution.
    Input:
        - func1
          function, the radial density profile
        - func2
          function, the radial velocity dispersion profile
        - args1(2)
          arguments of func1(2)
        - N
          number of particles
    Return:
        part: position, velocity and mass of particles
    '''
    r_bin = np.arange(0.01, 20, .01)
    
    # (r, theta, phi)
    prob = func1(r_bin, args1) * r_bin**2
    r_rand = np.random.choice(r_bin, N, p=prob/np.sum(prob))
    a_rand = np.random.rand(N) * np.pi * 2
    b_rand = (np.random.rand(N) - 0.5) * np.pi
    # (z, y, x) of the particles
    pos = np.array( [r_rand * np.sin(b_rand),
                     r_rand * np.cos(a_rand) * np.cos(b_rand),
                     r_rand * np.sin(a_rand) * np.cos(b_rand)] ).T
    
    # (sigma_x, sigma_y, sigma_z)
    sigma = func2(r_rand, args2)
    c_rand = np.random.rand(N) * np.pi * 2
    d_rand = (np.random.rand(N) - 0.5) * np.pi
    v_disp = np.array( [sigma * np.sin(d_rand),
                        sigma * np.cos(c_rand) * np.cos(d_rand),
                        sigma * np.sin(c_rand) * np.cos(d_rand)] ).T
    
    # mass = 1
    part = np.hstack( [pos, v_disp, np.full((N, 1), 1)] )
    return part


def set_spheroid(num, center, pec_vel, 
                 rho_pf, rho_pf_args, disp_pf, disp_pf_args):
    '''
    Description: set up a spheroid.
    Input:
        - num
          int, number of particles in this galaxy
        - center
          array-like, coordinate of the galaxy center
        - pec_vel
          array-like, peculiar velocity of the galaxy
        - rho_pf
          function, density profile of the galaxy
          ("const", "power", "double_power", ... other models to be added)
        - disp_pf
          function, velocity disperison profile
    Returns:
        part: particle info
    '''
    particles = sph_sampler( rho_pf, disp_pf, rho_pf_args, disp_pf_args, num )
    particles[:, iz] += center[0]
    particles[:, iy] += center[1]
    particles[:, ix] += center[2]
    particles[:, ivz] += pec_vel[0]
    particles[:, ivy] += pec_vel[1]
    particles[:, ivx] += pec_vel[2]
    return particles

def set_disk():
    pass