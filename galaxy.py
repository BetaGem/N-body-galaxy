import numpy as np
from params import *

def sph_sampler(N, func1, func2, args1, args2, max_radius, label):
    '''
    Description: produce random points following the given distribution.
    Params:
        - N : int
          number of particles
        - func1
          function, the radial density profile
        - func2
          function, the isotropic radial velocity dispersion profile
        - max_radius : float
          maximum radius of the initial posotions for the particles
        - args1(2)
          arguments of func1(2)
        - label : int
          label of particles
    Return:
        part: position, velocity and mass of particles
    '''
    check_positive( (max_radius,) )
    r_bin = np.arange(0.01, max_radius, 0.01)
    
    # (r, theta, phi)
    prob = func1(r_bin, args1) * r_bin**2
    r_rand = np.random.choice(r_bin, N, p=prob/np.sum(prob))
    a_rand = np.random.rand(N) * np.pi * 2
    b_rand = np.arcsin(2*np.random.rand(N) - 1)
    # (z, y, x) of the particles
    pos = np.array( [r_rand * np.sin(b_rand),
                     r_rand * np.cos(a_rand) * np.cos(b_rand),
                     r_rand * np.sin(a_rand) * np.cos(b_rand)] ).T
    
    # (sigma_x, sigma_y, sigma_z)
    sigma = func2(r_rand, args2)
    c_rand = np.random.rand(N) * np.pi * 2
    d_rand = np.arcsin(2*np.random.rand(N) - 1)
    v_disp = np.array( [sigma * np.sin(d_rand),
                        sigma * np.cos(c_rand) * np.cos(d_rand),
                        sigma * np.sin(c_rand) * np.cos(d_rand)] ).T
    
    # mass = 1, label = label
    part = np.hstack( [pos, v_disp,
                       np.full((N, 1), 1), np.full((N, 1), label)] )
    return part


def set_spheroid(num, center, pec_vel,
                 rho_pf, rho_pf_args, disp_pf, disp_pf_args,
                 max_radius=100, label=0):
    '''
    Description:
        set up a spheroid (stellar / dark matter).
    Params:
        - num : int
            number of particles in this galaxy
        - center : array-like
            coordinate of the galaxy center
        - pec_vel : array-like
            peculiar velocity of the galaxy
        - rho_pf : function
            density profile of the galaxy
            ("const", "power", "double_power", ... other models in profiles.py)
        - disp_pf : function
            velocity disperison profile
        - max_radius : float [pixel]
            maximal radius of particle initial positions
        - label : int
            label of particles
    Returns:
        particles: particle info
    '''
    particles = sph_sampler( num, rho_pf, disp_pf,
                             rho_pf_args, disp_pf_args, max_radius, label )
    particles[:, iz] += center[0]
    particles[:, iy] += center[1]
    particles[:, ix] += center[2]
    particles[:, ivz] += pec_vel[0]
    particles[:, ivy] += pec_vel[1]
    particles[:, ivx] += pec_vel[2]
    
    return particles


def disk_sampler(N, func1, args1, incl, min_radius, max_radius, label=0):
    '''
    Description: produce random points following the given distribution.
    Params:
        - N : int
            number of particles
        - func1 : function
            the radial density profile
        - args1
            arguments of func1
        - incl : float [deg]
            inclination angle of the disk (around +X direction)
        - min_radius / max_radius : float [pixel]
            minimal / maximal radius of patricle initial positions
        - label : int
            label of particles
    Return:
        part: position, velocity and mass of particles
    '''
    check_positive( (min_radius, max_radius) )
    r_bin = np.arange(min_radius, max_radius, 1)
    
    # (r, theta, phi)
    prob = func1(r_bin, args1) * r_bin
    r_rand = np.random.choice(r_bin, N, p=prob/np.sum(prob))
    # a_rand = np.random.rand(N) * np.pi * 2
    a_rand = np.linspace(0, 1, N, endpoint=False) * np.pi * 2
    # (z, y, x) of the particles
    pos = np.array( [r_rand * np.cos(a_rand) * np.sin(incl / 180 * np.pi),
                     r_rand * np.cos(a_rand) * np.cos(incl / 180 * np.pi),
                     r_rand * np.sin(a_rand)] ).T
    
    # rotation velocity, mass = 1
    v_0 = 17#np.sqrt( G * np.sum(r_rand[:, np.newaxis] > r_rand, axis=1) * 1 / r_rand ) / 1.5
    v_rot = np.array( [v_0 * np.sin(a_rand) * np.sin(incl / 180 * np.pi),
                       v_0 * np.sin(a_rand) * np.cos(incl / 180 * np.pi),
                       -v_0 * np.cos(a_rand)] ).T
    
    # mass = 1, label = label
    part = np.hstack( [pos, v_rot,
                       np.full((N, 1), 1), np.full((N, 1), label)] )
    return part


def set_disk(num, center, pec_vel,
             r_pf, r_pf_args, incl=0,
             min_radius=0.01, max_radius=100, label=0):
    '''
    Description:
        set up a stellar disk.
    Params:
        - num : int
            number of particles in this galaxy
        - center : array-like
            coordinate of the galaxy center
        - pec_vel : array-like
            peculiar velocity of the galaxy
        - r_pf : function
            radial density profile of the disk
            ("const", "power", "double_power", ... other models in profiles.py)
        - incl : float
            inclination angle of the disk (around +X direction)
        - min_radius / max_radius : float
            minimal / maximal radius of patricle initial positions in pixel
        - label : int
            label of particles
    Returns:
        particles: particle info
    '''
    particles = disk_sampler( num, r_pf, r_pf_args, incl,
                              min_radius, max_radius, label )
    particles[:, iz] += center[0]
    particles[:, iy] += center[1]
    particles[:, ix] += center[2]
    particles[:, ivz] += pec_vel[0]
    particles[:, ivy] += pec_vel[1]
    particles[:, ivx] += pec_vel[2]
    
    return particles


def check_positive(min=0, **params):
    for p in params:
        if p < min:
            raise ValueError("Negative values found. Please check the parameters.")