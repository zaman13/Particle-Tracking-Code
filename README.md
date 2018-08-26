# Particle-Tracking-Code
A Matlab code for tracking colloidal fluorescent nanoparticles. The code automatically compensates for any net drift motion of the nanoparticle and isolates the Brownian behavior. It analyzes the 2D position statistics and fits a Gaussian distribution. 

The current version is written for tracking a single particle. 

# Use
Run the particle_tracker_v_.m file. Make sure that the location of the source file is correctly specified. Also, the initial diameter guess (in pixels) should be close to the size of the particle that is to be tracked. 

## Reference
Based on the IDL particle tracking software. Tracking functions are used from the implementation found in http://site.physics.georgetown.edu/matlab/
