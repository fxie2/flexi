% Flexi Project Documentation
% Numerics Research Group,
  Institute for Aerodynamics and Gas Dynamics
  University of Stuttgart, Germany
% July 22, 2016

\hypertarget{introduction}{}

# Introduction

 [**FLEXI**](http://flexi-project.org) is a high-order discontinuous Galerkin simulation code for the solution of the time-dependent compressible Navier-Stokes equations on unstructured hexahedral elements in three space dimensions. The code was specifically designed for very high order accurate simulations on massively parallel systems. It is licensed under GPLv3, written in Fortran and parallelized with OpenMPI. Implemented features are
 
 * Arbitrary order nodal polynomial tensor product basis using Gauss or Gauss Lobatto collocation points
 * Various Riemann solvers for inter-element coupling
 * Large eddy simulation capabilities through different de-aliasing strategies and subgrid scale models
 * Matching high order curved mesh generation from external mesh formats (CGNS, GMSH) or simple analytic blocks via the open source preprocessor [HOPR](http://hopr-project.org)
 * Nonconforming interfaces based on the mortar approach
 * Non-reflecting boundary conditions and damping zones for direct aeroacoustic computations
 * Automatic domain decomposition for parallel simulations based on a space filling curve
 * High order explicit time Runge-Kutta time integration
 * I/O using the HDF5 library optimized for massively parallel jobs
 
## How to use the user guide

This user guide is organized to both guide the first steps as well as provide a complete overview of the simulation code's features from a user and a developer point of view.

* Chapter \ref{chap:gettingstarted} contains step by step instructions from obtaining the source code up to running a first simulation and visualizing the simulation results. In addition, it provides an overview of the whole simulation framework and the currently implemented features.

* Chapter \ref{chap:codeoptions} lists all compiler and runtime options.

* Chapter \ref{chap:workflow} is meant as a complete user guide with a detailed description how to use and apply the features of **FLEXI** from a user's point of view. This includes setting up solver settings, initial and boundary conditions, the mesh interface, parallel execution and the currently available post processing capabilities.

* Simulation tutorials are contained in Chapter \ref{chap:tutorials}.

* The regression check utility used for automated testing is described in Chapter \ref{chap:regressioncheck}.

* The unit test system used to test key routines with CTest is described in Chapter \ref{chap:unittest}.
