## Lid-driven Cavity

This tutorial describes numerical solution of an  incompressible flow in a three-dimensional square domain. The domain contains of a cavity with isothermal walls (and periodic boundary conditions in the $z-$direction) and prescribed flow state at the open upper face. The following figure shows the used mesh and velocity flow field.

![](tutorials/02_cavity/cavity_mesh.png)   ![](tutorials/02_cavity/cavity_result.png) 

Figure: Mesh and flow field solution with glyph vector plot of the velocity field. View in $x$-$y$-plane.\label{fig:cavity_mesh_and_result}

Copy the ``cavity`` tutorial folder \label{missing:aliases_tutorial_cavity1}

        cp -r $FLEXI_TUTORIALS/cavity .


### Mesh Generation with HOPR

The mesh files used by **FLEXI** are created by supplying an input file *parameter_hopr.ini* with the appropriate information.

    ./hopr parameter_hopr.ini

This creates the mesh file *cavity_mesh.h5* in HDF5 format.

### Flow Simulation with FLEXI

The simulation setup is defined in *parameter_flexi.ini*. The initial condition is selected via the second **RefState=(/1.225,0.,0.,0.,101325./)** which represents the solution vector $(\rho, u, v, w, p)^T$. The first **RefState=(/1.225,1.,0.,0.,101325./)** represents the state for the boundary condition.  

**IniRefState = 2** : the initial condition uses **RefState 2** for the initial flow field solution.

**IniExactFunc = 1** : the used exact function routine uses **RefState 1**, e.g., for the calculation of the $L_2$ error norm.

Constant flow properties like the gas constant are given in table \ref{tab:freestream_flow_prop} 
and define the gas behavior in combination with the ideal gas law $p=\rho R T$.

Table: Numerical settings \label{tab:freestream_flow_prop}

| Property                        | Variable      | Value       |
| ------------------------------- |:-------------:| -----------:|
| dynamic viscosity $\mu$         | mu0           |  0.1        |
| ideal gas constant $R$          | R             |  1          |
| isentropic coefficient $\kappa$ | kappa         |  1.4        |

### Numerical settings

The DG solution on the mesh is represented by piecewise polynomials and the polynomial degree in this tutorial is chosen as $N=3$.

The default settings for these properties are displayed in table \ref{tab:freestream_num_set}. 

Table: Numerical settings \label{tab:freestream_num_set}

| Variable        | Description                            | Value         |
| --------------- |:--------------------------------------:|:-------------:|
| N               | polynomial degree                      | 3             |
| MeshFile        |                                        | cavity_mesh.h5|
| tend            |                                        | 0.2           |
| Analyze_dt      |                                        | 0.2           |
| nWriteData      |                                        | 1             |
| CFLscale        |                                        | 0.99          |
| DFLscale        |                                        | 0.4           |


The command

~~~~~~~
flexi parameter_flexi.ini > std.out
~~~~~~~

runs the code and dumps all output into the file *std.out*. 

