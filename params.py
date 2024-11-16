# similation setup
# -------------------
N  =  4096  # number of particles
Ng =  512   # number of grids, with size h = 1
dt =  0.1   # time step
t_tot = 10  # total time

# constants
# -------------------
G = 1e2     # gravitational const.
seed = 42   # random number seed

# visualization parameters
# -------------------
plot_on_the_fly = False   # plot result on-the-fly (slower)
save = False  # if plot_on_the_fly, save the data or not

# indices
# -------------------
iz, iy, ix, ivz, ivy, ivx, im, il = range(8)