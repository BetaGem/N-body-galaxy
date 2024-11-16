# N-body-galaxy
A simple particle-mesh code for N-body simulations in astrophysics, focusing on galaxy structures and evolution.

## Files:

`main.py`: the main function. **[initialize your particles in this file]**

`params.py`: **[set the basic simulation parameters in this file]**

`force.py`: calculate acceleration and force

`galaxy.py`: galaxy generator

`Poisson.py`: Poisson solver based on FFT

`mass_assign.py`: calculate density field from particle distribution

`profiles.py`: custom analytical radial profiles of galaxies

`major_merger.mp4`: a simulation of elliptical galaxy merger

## Run the code

1. Create a new folder `data` in your work path.  
   ```
   mkdir ./data
   ```
3. Run `main.py` to get particle data in each timestep.  
   ```
   python main.py
   ```
5. Run `visualize.py` to visualize the results.  
   ```
   python visualize.py
   ```

## Reference
*High performance computing and numerical modelling* by Volker Springle: https://arxiv.org/abs/1412.5187
