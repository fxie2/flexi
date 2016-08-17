\hypertarget{regressioncheck}{}

# Regression Check \label{chap:regressioncheck}

The regression check (or *Reggie*) is intended to provide continuous code revisal in order to locate bugs.

What checks can be performed?

| Type                            | Modes                 | Description                                         |
| ------------------------------- |:---------------------:| ---------------------------------------------------:|
| Unit-tests                       | run-time compilation  |                                                     |                          
| L2 and Linf norms               | 1, 3                  | Error norms calculated via *ExactFunc*              |
| Record points                   | 1, 2, 3, 4            | time signal of properties recorded over time        |
| Performance                     | 1                     | time measurement of wall time                       |
| No-result                       | -                     | only checks for successful execution                |
| h5diff                          | 4                     | diff multiple HDF5 files                            |
| analyze-tools                   | 1, 2, 3, 4            | high-level analysis (entropy, FFT, time averaging) |
| line plot                       | 2, 3, 4               | plot properties over lines                          |

Table: Regression check types.

Parameters are ...

| Property                        | Value         | Description                                |
| ------------------------------- |:-------------:| ------------------------------------------:|
| Comparison mode                 | 1             | use constant pre-defined value             |
|                                 | 2             | use pre-defined function, e.g., $f(x)=x^2$ |
|                                 | 3             | use supplied data file                     |
|                                 | 4             | HDF5 supplied data file for H5-diff        |

Table: Regression check parameters.

## Using Reggies

The following files are needed for each regression check folder directory

| File                            | Value         | Description                                |
| ------------------------------- |:-------------:| ------------------------------------------:|
| configuration.cmake             | required      | compilation flags: EQNSYS, MPI, etc.       |
| parameter_reggie.ini            | required      | Comparison type, general settings          |
| mesh.h5                         | required      | mesh file in HDF5 format                   |
| parameter_flexi.ini             | required      | IC, BC and numerical settings              |
| parameter_hopr.ini              | optional      | for ANSA, ICEM formats etc.                |

Table: Required files for regression check.

### configuration.cmake

List of builds: MPI F/T, Gauss/GL, EQNSYS (5), BR 1/2, parabolic T/F, viscosity (3)

### parameter_reggie.ini

Comparison mode: const., function, data file, HDF5



~~~~~~~Bash
./regressioncheck parameter_reggie.ini
~~~~~~~

### Running a regression check

The program *regressioncheck* is run

        ./regressioncheck [parameter]

can be done in *run* or *build* mode

        ./regressioncheck run

simply uses the supplied *bin/flexi* binary files and checks alls examples, e.g., */examples/freestream*

        ./regressioncheck build

creates new *build* directories with the supplied compile flags, e.g., *build_GNU_single_parabolic*




## Creating New Reggies

