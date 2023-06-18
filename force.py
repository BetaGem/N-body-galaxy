import numpy as np
from params import *

def accel(Phi):
    '''
    Description:
      Calculate acceleration field using potential.
    Input: 
        - Phi
          Gravitational potential field given by the poisson solver.
    Return: acc
      Acceleration field.
    '''
    acc = np.zeros((3, Ng, Ng, Ng))
    
    # main part
    acc[iz][1: Ng - 1]  = (Phi[:Ng-2] - Phi[2:Ng]) / 2
    acc[iy][:,1:Ng - 1] = (Phi[:,:Ng-2] - Phi[:,2:Ng]) / 2
    acc[ix][:,:,1:Ng-1] = (Phi[:,:,:Ng-2] - Phi[:,:,2:Ng]) / 2
    
    # boundaries
    acc[iz][0]     = (Phi[-1] - Phi[1]) / 2
    acc[iy][:,0]   = (Phi[:,-1] - Phi[:,1]) / 2
    acc[ix][:,:,0] = (Phi[:,:,-1] - Phi[:,:,1]) / 2
    
    acc[iz][Ng - 1]   = (Phi[Ng-2] - Phi[-1]) / 2
    acc[iy][:,Ng - 1] = (Phi[:,Ng-2] - Phi[:,-1]) / 2
    acc[ix][:,:,Ng-1] = (Phi[:,:,Ng-2] - Phi[:,:,-1]) / 2
    
    return acc

def Force(pos, acc):
    '''
    Description:
      Calculate specific force for a particle.
    Input: 
        - pos: 
          position of the particle
        - acc:
          acceleration field
    Return:
      F: force
    '''
    
    F = np.array([0., 0., 0.])
    q, p, r = int(np.floor(pos[iz] - 1/2)), int(np.floor(pos[iy] - 1/2)), int(np.floor(pos[ix] - 1/2))
    qs, ps, rs = pos[iz] - 1/2 - q, pos[iy] - 1/2 - p, pos[ix] - 1/2 - r
    
    F += acc[:, q % Ng, p % Ng, (r + 1) % Ng]             * (1 - ps) * (1 - qs) * rs
    F += acc[:, q % Ng, (p + 1) % Ng, (r + 1) % Ng]       * ps * (1 - qs) * rs
    F += acc[:, (q + 1) % Ng, p % Ng, (r + 1) % Ng]       * (1 - ps) * qs * rs
    F += acc[:, (q + 1) % Ng, (p + 1) % Ng, (r + 1) % Ng] * ps * qs * rs
    
    F += acc[:, q % Ng, p % Ng, r % Ng]             * (1 - ps) * (1 - qs) * (1 - rs)
    F += acc[:, q % Ng, (p + 1) % Ng, r % Ng]       * ps * (1 - qs) * (1 - rs)
    F += acc[:, (q + 1) % Ng, p % Ng, r % Ng]       * (1 - ps) * qs * (1 - rs)
    F += acc[:, (q + 1) % Ng, (p + 1) % Ng, r % Ng] * ps * qs * (1 - rs)
        
    return F