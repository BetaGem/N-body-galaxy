import numpy as np
import matplotlib.pyplot as plt

import mass_assign as ms
import plot_me as pm
import galaxy
import Poisson
import force  
import profiles as pf
from params import *

np.seterr(all="ignore")

if __name__ == "__main__":
    
    # initializing particles
    stellar1 = galaxy.set_spheroid(N, [0, 0.6*Ng, 0.6*Ng], [0, -5, 0], 
                                  pf.Gaussian, (1, 6), pf.const, (2,))
    stellar2 = galaxy.set_spheroid(N, [0, 0.5*Ng, 0.5*Ng], [0, 0, 0], 
                                  pf.Gaussian, (1, 3), pf.const, (3,))
    part = np.vstack( (stellar1, stellar2) )
    
    # density field
    dens = ms.CIC_3D(part)
    
    # gravitational potential
    Phi = Poisson.poisson_solver_fft3(dens)
    
    # acceleration field
    a = force.accel(Phi)
    
    for t in np.arange(0, t_tot, dt):
        
        # leapfrog
        for i in range(len(part)):
            v_mid = part[i][ivz:ivx+1] + force.Force(part[i][iz:ix+1], a) * dt / 2 / part[i][im]
            part[i][iz:ix+1] = part[i][iz:ix+1] + v_mid * dt
            part[i][ivz:ivx+1] = v_mid + force.Force(part[i][iz:ix+1], a) * dt / 2 / part[i][im]
        
        # update fields
        dens = ms.CIC_3D(part)
        Phi = Poisson.poisson_solver_fft3(dens)
        a = force.accel(Phi)
        
        # plot
        # img = ms.CIC_3D(part[:N], Ng)
        pm.plot_all(part, dens, Phi, a)