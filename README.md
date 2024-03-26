# N-body-galaxy
A simple particle-mesh code for N-body simulations on galaxy structure and evolution.

```main.py```: **please initialize the particles in this file**

```params.py```: **set the basic simulation parameters here**

```force.py```: calculate acceleration and force

```Poisson.py```: solve the Poisson equation numerically

```mass_assign.py```: get density field from particle distribution

```profiles.py```: custom analytical radial profiles of galaxies

## Run the code

1. Create a new folder ```data``` in your work path.
   ```
   mkdir ./data
   ```
3. Run ```main.py``` to get particle data in each timestep.
   ```
   python main.py
   ```
5. Run ```visualize.py``` to visualize the results.
   ```
   python visualize.py
   ```
