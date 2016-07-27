# A Simplified 2-D Cloud Resolving Model (CRM)

## Introduction

Hello, world!

The name of the system I want to develope eventually is Super-Parametertization System (SPS).
It is quite big.
And by now, this package is actually a simplified 2-D cloud resolving model (CRM), which is part of my Master of Science (M.S.) work in Chinese Academy of Meteorological Sciences (CAMS) under the guidance of Dr. Guoqiang Xu.

This simplified 2-D CRM consists of a fully compressible dynamic core and mainly a WSM6 microphysics package transplanted from WRF.
The dynamic core is designed on a horizontally stagged Arakawa-C grid with Charney-Phillips vertical variable configuration.
The time integration scheme is Runge-Kutta 3-order, and the spacial discretization scheme is 5-order.
Physical packages are WSM6 microphysics with a simple K-theory subgrid mixing scheme and a simple binary (0-1) cloud fraction scheme.

Some 2-D ideal experiments are conducted in this version, such as density current, convective thermal bubble, internal gravity wave, and an ideal thunderstorm experiment.

The model is written in `Fortran 90/95`, and it is capable of parallel running with the support of `OpenMP`.
If you are interested in developing this model with me, just send me an email.
And suggestions are much appreciated.

Let's hack nature~

## Ideal Experiments

Below are some results of ideal experiments conducted by this model:

### Thermal Bubble
![](http://spsystem.sinaapp.com/results/tb.gif)
### Density Current
![](http://spsystem.sinaapp.com/results/dc.gif)
### Internal Gravity Wave
![](http://spsystem.sinaapp.com/results/igw.gif)

## The modules
+ sp_dynamic.f90                          : the main program
+ sp_module_model.f90                     : settings for the running
+ sp_module_gridvar.f90                   : variables to calculate
+ sp_module_constant.f90                  : some physical constants
+ sp_module_boundary.f90                  : boundary conditions
+ sp_module_initiate.f90                  : initial conditions
+ sp_module_interpolate.f90               : interpolation schemes
+ sp_module_timeschemes.f90               : time integration schemes
+ sp_module_integrate.f90                 : time integration
+ sp_module_advection.f90                 : advection schemes
+ sp_module_tendency.f90                  : tendency calculation
+ sp_module_wsm6.f90                      : WSM6 microphysics scheme
+ sp_module_subgrid.f90                   : subgrid mixing scheme
+ sp_module_cldfra.f90                    : cloud fraction scheme
+ sp_module_output.f90                    : output the results
+ sp_module_debug.f90                     : some debugging tools

## How to use

The basic settings of the model is in `sp_module_model.f90`.

You can make the whole program by `make`.
Note that the default compiler is `ifort`.
And run it with `./sp_dynamic`.
After the running is done, you can plot the theta' field with the NCL scripts in the directory of `script`.

## Copyright

Version: 0.2

Author: Feng Zhu

Email: zhuf.atmos@gmail.com

Date: 2014-06-12

License: MIT

