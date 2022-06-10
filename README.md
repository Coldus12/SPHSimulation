# Vltava

Vltava is a work-in-progress SPH-based fluid simulation program. The goal of this project is to create a semi-realtime fluid simulation which harnesses the power of the GPU. To accomplish this, I'm using the Smoothed-Particle-Hydrodynamics computational method.


# Pictures (WIP)

![Pic1]([http://url/to/img.png](https://github.com/Coldus12/SPHSimulation/blob/master/2744particle32container.png))
![Pic2]([http://url/to/img.png](https://github.com/Coldus12/SPHSimulation/blob/master/2744particle32container-2.png))

# To-do list

 - Optimization
	 - neighbourhood search
	 - shader code optimization
- Proper surface-fluid interaction
- Implement other types of SPH
- Implement a soft-body system which works with the project

# Tutorials 

Before I started this project I haven't had any prior experience with either the Vulkan API or the SPH method. So to understand these I followed the tutorials listed below:

- Vulkan tutorials
	- https://vulkan-tutorial.com/
	- https://vkguide.dev/
- SPH tutorials
	- [Eurographics Tutorial 2019](https://interactivecomputergraphics.github.io/SPH-Tutorial/)
	- [Philip Mocz's SPH tutorial](https://philip-mocz.medium.com/create-your-own-smoothed-particle-hydrodynamics-simulation-with-python-76e1cec505f1)
