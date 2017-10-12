
#### SACRAMENTO MODEL ALGORITHM

# Notes 
# Each subfunction needs to be a stand-alone function, and then we should be
# able to call them from the main function.


# Input requirements
# >> petHamon2 requires Tavg (time-series), basin_lat (double), coeff (double)
# >> snow17 requires pars, prcp (time-series), temp, elev, state_input, Time


# Standardization
# use lower case argument names 
# use points when coercing argument names


# sacSim <- function {

# CALCULATE PET 
# petHamon2(coeff, basin_lat, Tavg)
# 
# # RUN SNOW MODEL FOR EACH HRU...... 
# snow17(pars, prcp, tavg, elev, statesInput, timeMat) // output = meltNRain
# 
# # RUN SAC-SMA MODEL FOR EACH HRU....
# sac_sma(Prcp, Tavg, Basin_Lat, Basin_Elev, Par, ...)
# 
# ROUT LOCAL FLOWS FOR ALL SUBCATCHMENTS
# rout_lohamann(pars, flowlen, UH_DAY, KE)
# 
# APPLY CONVOLUTION FUNCTION FOR ALL SUBCATCHMENTS
# convol(...)
# 
# 
# return(flow)
# 
# 
# }

f <- function() {
  x <- 1
  y <- 2
  c(x, y)
}


f(3,3)

rm(f)