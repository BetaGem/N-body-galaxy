# N-body-galaxy
Particle-mesh code for N-body simulations on galaxy structure and evolution.

```main.py```: initialize the particles

```params.py```: set the simulation parameters

```force.py```: calculate acceleration and force

```Poisson.py```: Solve the Poisson equation numerically

```mass_assign.py```: get density field from particle distribution

```profiles.py```: custom analytical radial profiles

Run ```main.py``` first to get particle data in each timestep. Then run ```visualize.py``` to see the results.
