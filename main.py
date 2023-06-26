import numpy as np
import matplotlib.pyplot as plt

import mass_assign as ms
import visualize as vi
import galaxy
import Poisson
import force  
import profiles as pf
from params import *

np.seterr(all="ignore")

if __name__ == "__main__":
    
    # initializing particles
    stellar1 = galaxy.set_spheroid(N, [320, 330, 330], [0, -7.57, 0], 
                                   pf.power, (1, 10, 2), pf.const, (115,))
    stellar2 = galaxy.set_spheroid(N//3, [160, 150, 150], [0, 22.73, 0], 
                                   pf.power, (1, 5, 2), pf.const, (66,))
    part = np.vstack( (stellar1, stellar2) )

    # density field
    dens = ms.CIC_3D(part)
    
    # gravitational potential 
    Phi = Poisson.poisson_solver_fft3(dens)
    
    # acceleration field
    a = force.accel(Phi)

    # time evolution
    for n, t in enumerate(np.arange(0, t_tot, dt)):
        
        # leapfrog
        part[:, iz:ix+1] += part[:, ivz:ivx+1] * dt / 2
        # update fields
        dens = ms.CIC_3D(part)
        Phi = Poisson.poisson_solver_fft3(dens)
        a = force.accel(Phi)
        for i in range(len(part)):
            part[i][ivz:ivx+1] += force.Force(part[i][iz:ix+1], a) * dt / part[i][im]
        part[:, iz:ix+1] += part[:, ivz:ivx+1] * dt / 2
        
        # plot
        # vi.plot_all(part, dens, Phi, a)
        vi.save_data(part, "./data/"+str(n)+".npy")
