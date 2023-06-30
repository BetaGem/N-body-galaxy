# N-body-galaxy
Particle-mesh code for N-body simulations on galaxy structure and evolution.

```main.py```: you should initialize your particles in this file

```params.py```: set the simulation parameters here

```force.py```: calculate acceleration and force

```Poisson.py```: solve the Poisson equation numerically

```mass_assign.py```: get density field from particle distribution

```profiles.py```: custom analytical radial profiles

## Run the code

Run ```main.py``` first to get particle data in each timestep. Then run ```visualize.py``` to see the results.
