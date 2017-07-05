

# Sacramento Soil Moisture Accounting Model in R (SAC-SMA)

[![](http://www.appsolutelydigital.com/ModelPrimer/images/image79.jpeg)

# _DO NOT USE - Currently under development_

## Description

The SAC-SMA is a continuous soil moisture accounting model with spatially lumped 
parameters that simulates runoff within a basin. The model divides the basin into 
lower and upper zones at different depths, and defines the distribution of moisture, 
i.e., tension water components (driven by evapotranspiration and diffusion) and 
free water components (driven by gravitational forces) in each of these two zones 
via a set of parameters. The model uses precipitation and temperature variables, 
along with parameters on soil moisture states and the basin’s relative permeability 
to estimate the amount of water that enters, is stored in, and leaves the basin. 
Thus, the model estimates several key hydrologic processes including evapotranspiration, 
percolation, interflow, and different forms of runoff from a basin. 

The SAC-SMA model is used for various applications that are mainly streamflow or 
runoff centric, for example, river forecasting, water supply forecasting, basin 
hydrologic hazard estimates, and basin climate change assessments. The model is 
ideal for large drainage basins and uses multiple years of records for calibration. 
The SAC-SMA model is a key model used by the U.S. National Weather Service River 
Forecast System (NWSRFS) to issue river forecasts across the country. 

## R Version

The original model code is written in Fortran and is publicly accessible. The R
version of the model is translated from the Matlab code developed by Dr.Sungwook Wi.
