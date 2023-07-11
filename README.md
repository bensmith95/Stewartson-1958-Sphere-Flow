# Stewartson-1958-Sphere-Flow
MATLAB code that solves the [Stewartson (1958)](https://doi.org/10.1007/978-3-642-45885-9_6) system (2D cartesian Navier Stokes in SF-Vorticity form) for the equatorial flow around a rotating sphere.

[Stewartson (1958)](https://doi.org/10.1007/978-3-642-45885-9_6) derived the equations that he thought governed the flow around a rotating sphere near the equator which turns out to be effectively the same as planar Navier Stokes equations (upto a couple of constants). Although the numerical simulations of the full problem illicits different features, suggesting his model needs some tweaking, the problem is still quite interesting due to how much the resutling PDEs resemble the steady 2D Navier Stokes equations in cartesian coordinates, which are not trivial to solve!

It is hoped that the provided code allows a basis to obtain fast and accurate solutions of steady flows for these types of geometries. 

The PDEs are in Streamfunction-Vorticity form, which decreases the number of equations by one, but also makes boundary conditions difficult to implement (so some care and attention to detail is needed). These are then discretised via finite differences producing a large system of non-linear equations. As a full scale Newton method requires a huge amount of storage, to store and solve the matrix systems, iterative methods are needed. This code utilises a geometric multigrid method to solve large systems of non-linear equations in a reasonably quick manner by iterating on smaller/coarser grids and correcting the solutions on larger/finer grids. The full details are in Smith (2023).

To run the code in its current format simply download all the files into the same directory and then run _Stewartson(**Re**)_ where _**Re**_ is the square root of the Reynolds number (this is the constants I mentioned earlier...). The code will also determine the V component of the velocity and scaled pressure. <ins>Importantly</ins>, make sure that _BL.mat_ is saved in the directory. This is the inlet boundary condition from the [boundary layer region](https://github.com/bensmith95/Rotating-Sphere-Boundary-Layer), otherwise it will not work! After the code runs, all the data will be saved to a file called _Stew_Re=**Re**.mat_. To view the results simply run _figs.m_.

<p align="center">
  <img height="300" src="https://user-images.githubusercontent.com/29705711/176724646-3d343d2e-65ab-404b-a78f-4adf48da0ad0.png">
</p>
