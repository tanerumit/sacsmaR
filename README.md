
Sacramento Soil Moisture Accounting (SAC-SMA) Model
===================================================

***The model is currently under development***

### Installation

You can install the development version from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("tanerumit/sacsmaR")
```

### R Version

The original model code is written in Fortran and is publicly accessible. The R version of the model is translated from MATLAB code developed by [Sungwook Wi](https://github.com/sungwookwi)

------------------------------------------------------------------------

### Description

    The SAC-SMA is a continuous soil moisture accounting model with spatially 
    lumped parameters that simulates runoff within a basin. The model divides the 
    basin into lower and upper zones at different depths, and defines the distribution 
    of moisture, i.e., tension water components (driven by evapotranspiration and 
    diffusion) and free water components (driven by gravitational forces) in each of 
    these two zones via a set of parameters. The model uses precipitation and temperature 
    variables, along with parameters on soil moisture states and the basin’s relative 
    permeability to estimate the amount of water that enters, is stored in, and leaves 
    the basin. Thus, the model estimates several key hydrologic processes including 
    evapotranspiration, percolation, interflow, and different forms of runoff from a basin. 

    Further information is available at: [NOAA - National Weather Service](http://www.nws.noaa.gov/oh/hrl/nwsrfs/users_manual/part2/_pdf/23sacsma.pdf)

![alt text](http://www.appsolutelydigital.com/ModelPrimer/images/image79.jpeg "SAC-SMA model")

------------------------------------------------------------------------

### Included functions

The package contains for main function, which are described below:

#### Potential Evaporation module based on Hamon equation:

``` r
hamon(par, tavg, lat, jday)
```

#### Sacsma module for run-off generation

``` r
sacsma(par, states.ini = c(0, 0, 5, 5, 5, 0), prcp, pet, lat, elev, verbose = FALSE)
```

#### Lohmann channel routing

``` r
lohmann(par, flength, uh.day = 96, ke = 12)
```

#### Snow17 module for snowmelt calculation

``` r
snow17(par, states.ini = c(0, 0, 0, 0), prcp, tavg, elev, jday, verbose = FALSE)
```

#### HydroSim: wrapper for Sac-sma, Lahmann, and Hamon modules

``` r
hydroSim(par.hamon, par.snow17, par.sacsma, par.lohmann, tavg.grid, prcp.grid, lat.grid, 
  elev.grid, area.grid, flength.grid, subcat.grid = NULL, jday)
```
