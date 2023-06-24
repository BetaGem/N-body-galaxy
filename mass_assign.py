# reference: arXiv 1412.5187

import numpy as np
from params import *

def NGP_2D(particle):
    '''
    Description:
        Nearest grid point (NGP) assignment of a particle to a 2D mesh.
    Input: 
        - part
          column 0-1: positions of particles
          column 2-3: velocities of particles
          column 4: masses of particles
    Return: 
        density metrix
    '''
    pass

def NGP_3D(particle):
    '''
    Description:
        Nearest grid point (NGP) assignment of a particle to a 3D mesh.
    Input: 
        - part
          column 0-2: positions of particles
          column 3-5: velocities of particles
          column 6: masses of particles
    Return: 
        density metrix
    '''
    pass

def CIC_2D(particle):
    '''
    Description:
        Clouds-in-cell (CIC) assignment of a particle to a 2D mesh.
    Input: same as NGP_2D
    Return: same as NGP_2D
    '''
    dens = np.zeros((Ng, Ng))
    for i in particle:
        p, q = int(np.floor(i[0] - 1/2)), int(np.floor(i[1] - 1/2))
        ps, qs = i[0] - 1/2 - p, i[1] - 1/2 - q

        dens[q % Ng][p % Ng]             += (1 - ps) * (1 - qs) * i[4]
        dens[q % Ng][(p + 1) % Ng]       += ps * (1 - qs) * i[4]
        dens[(q + 1) % Ng][p % Ng]       += (1 - ps) * qs * i[4]
        dens[(q + 1) % Ng][(p + 1) % Ng] += ps * qs * i[4]
        
    return dens

def CIC_3D(part):
    '''
    Description:
        Clouds-in-cell (CIC) assignment of a particle to a 3D mesh.
    Input: same as NGP_3D
    Return: same as NGP_3D
    '''
    dens = np.zeros((Ng, Ng, Ng))
    for i in part:
        
        m = i[im]
        q, p, r = int(np.floor(i[iz] - 1/2)), int(np.floor(i[iy] - 1/2)), int(np.floor(i[ix] - 1/2))
        qs, ps, rs = i[iz] - 1/2 - q, i[iy] - 1/2 - p, i[ix] - 1/2 - r

        dens[q % Ng][p % Ng][r % Ng]             += (1 - ps) * (1 - qs) * (1 - rs) * m
        dens[q % Ng][(p + 1) % Ng][r % Ng]       += ps * (1 - qs) * (1 - rs) * m
        dens[(q + 1) % Ng][p % Ng][r % Ng]       += (1 - ps) * qs * (1 - rs) * m
        dens[(q + 1) % Ng][(p + 1) % Ng][r % Ng] += ps * qs * (1 - rs) * m
        
        dens[q % Ng][p % Ng][(r + 1) % Ng]             += (1 - ps) * (1 - qs) * rs * m
        dens[q % Ng][(p + 1) % Ng][(r + 1) % Ng]       += ps * (1 - qs) * rs * m
        dens[(q + 1) % Ng][p % Ng][(r + 1) % Ng]       += (1 - ps) * qs * rs * m
        dens[(q + 1) % Ng][(p + 1) % Ng][(r + 1) % Ng] += ps * qs * rs * m
    
    return dens

def zero_padding(dens):
    
    dens_pad = np.zeros((2*Ng, 2*Ng, 2*Ng))
    dens_pad[:Ng, :Ng, :Ng] = dens
    
    return dens_pad