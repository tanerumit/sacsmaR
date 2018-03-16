

library(lubridate)
library(dplyr)
library(magrittr)

################ SINGLE HRU SIMULATION

#Source all scripts
file_sources = paste0("./R/", list.files(path = "./R/", pattern = "*.R"))
sapply(file_sources, source, .GlobalEnv)

# HRU info file (lat, lon, area, elev, flowlength, id)
grid_info <- readr::read_rds("./data/hru_info_test.Rds")
hru_par   <- readr::read_rds("./data/hru_par_test.Rds")
hru_prcp  <- readr::read_rds("./data/hru_prcp_test.Rds") 
hru_tavg  <- readr::read_rds("./data/hru_tavg_test.Rds") 

# Date-time preferences
str_date <- as.Date("1970/01/1")
end_date <- as.Date("1980/12/31")
seq_date <- seq.Date(str_date, end_date, by = "day") 
sim_date <- seq.Date(str_date, end_date, by = "day") 
sim_per  <- length(sim_date) 

#Parameters
hru_lat     <- unlist(grid_info[,1], use.names = FALSE)
hru_lon     <- unlist(grid_info[,2], use.names = FALSE)
hru_area    <- unlist(grid_info[,3], use.names = FALSE)
hru_elev    <- unlist(grid_info[,4], use.names = FALSE)
hru_flowlen <- unlist(grid_info[,5], use.names = FALSE)

# Set parameters
par_sacsma   <- hru_par[,1:16]
par_petHamon <- hru_par[,17]
par_snow17   <- hru_par[,18:27]
par_routLah  <- hru_par[,28:31]

jday <- lubridate::yday(sim_date)

# #### TEST 1) Check each function individually
pet <- petHamon(coeff = par_petHamon[1], tavg = hru_tavg[[1]],
                lat = hru_lat[[1]], jday = jday)

snowMelt <- snow17(par = par_snow17[1,],
                   prcp = hru_prcp[[1]],
                   tavg = hru_tavg[[1]],
                   elev = hru_elev[[1]],
                   jday = jday)

### SAC MODULE
simflow <- sacSma(par = par_sacsma[1,],
                  prcp = hru_prcp[[1]],
                  pet = pet, lat = hru_lat[1], elev = hru_elev[1])

### ROUTING MODULE
flow <- routeLohmann(par = par_routLah[1,], flowLength = hru_flowlen[1])


################# BASIN-WIDE SIMULATION

hru_par <- hru_par %>% as.matrix()

par_sacsma   = hru_par[,1:16] 
par_petHamon = hru_par[,17]  
par_snow17   = hru_par[,18:27] 
par_routLah  = hru_par[,28:31] 

hru_tavg <- lapply(hru_clim, "[[", "tavg")
hru_prcp <- lapply(hru_clim, "[[", "prcp")

tavg = hru_tavg 
prcp = hru_prcp 
lat  = hru_lat 
elev = hru_elev 
area = hru_area
flowLength = hru_flowlen

flowR <- sacHydroSim(
  par_petHamon = hru_par[,17],
  par_snow17   = hru_par[,18:27],
  par_sacsma   = hru_par[,1:16],
  par_routLah  = hru_par[,28:31],
  hru_tavg = hru_tavg, 
  hru_prcp = hru_prcp, 
  hru_lat  = hru_lat, 
  hru_elev = hru_elev, 
  hru_area = hru_area,
  hru_flowLength = hru_flowlen,
  jday = jday)



 

