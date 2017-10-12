
library(lubridate)
library(dplyr)
library(magrittr)

#Source all scripts
file_sources = paste0("./R/", list.files(path = "./R/", pattern = "*.R"))
sapply(file_sources, source, .GlobalEnv)

# Date-time preferences
str_date <- as.Date("1990/01/1")
end_date <- as.Date("2010/01/1")
seq_date <- seq.Date(str_date, end_date, by = "day") 
sim_date <- seq.Date(str_date, end_date, by = "day") 
sim_per  <- length(sim_date) # Simulation period
sim_dayOfYear     <- lubridate::yday(sim_date)
sim_lastDayOfYear <- yday(as.Date(paste0(sim_doy_mat$year,"/12/31")))

data_folder <- "C:\\Users\\Umit\\Google Drive\\Research\\projects\\stcroix-basin\\input\\"

#Calibration parameters (50 hrus x 31 parameters)
hru_par_dat <- read.table(paste0(data_folder, "stcroix_calib_par_trial1_coord.txt")) %>%
  as_tibble() %>% right_join(grid_info[,c(1,2)], by = c("V1", "V2"))
hru_par <- hru_par_dat %>% select(V3:V33)

write.table(x = hru_par, "stcroix_calib_par_trial1.txt", row.names = F, col.names = F)



# HRU info file (lat, lon, area, elev, flowlength, id)
grid_info <- read.table(paste0(data_folder, "HRUinfo_stcroix_outlet.txt"))
hru_lat     <- grid_info[,1]
hru_lon     <- grid_info[,2]
hru_area    <- grid_info[,3]
hru_elev    <- grid_info[,4]
hru_flowlen <- grid_info[,5]

# HRU climate data
num_hru <- length(hru_flowlen)
hruclim_dir = paste0(data_folder, '\\HRU_clim\\data_')

# Read-in climate data for all HRUs
hru_clim_raw <- list()
for (n in 1:num_hru) {
  lat  <- formatC(grid_info[n,1], format = 'f', flag='0', digits = 4)
  lon  <- formatC(grid_info[n,2], format = 'f', flag='0', digits = 4)
  hru_clim_raw[[n]] <- read.table(paste0(hruclim_dir,lat, "_", lon)) %>% as.matrix()
}

# Subset climate data for the period of interest
hru_clim <- list()
for (n in 1:num_hru) {
  dat <- hru_clim_raw[[n]]
  colnames(dat) <- c("year", "month", "day", "prcp", "tavg")
  grid_date <- as.Date(paste(dat[,1], dat[,2], dat[,3], sep ="/"))
  grid_ind  <- which(str_date == grid_date):which(end_date == grid_date) 
  hru_clim[[n]] <- dat[grid_ind,] %>% as_tibble()
}

# Set parameters
par_sacsma   <- hru_par[,1:16]
par_petHamon <- hru_par[,17]
par_snow17   <- hru_par[,18:27]
par_routLah  <- hru_par[,28:31]

# #### TEST 1) Check each function individually
# pet <- petHamon(coeff = par_petHamon[1], tavg = hru_clim[[1]]$tavg, 
#                 lat = hru_lat[[1]], dayOfYear = sim_dayOfYear)
# 
# snowMelt <- snow17(par = par_snow17[1,], 
#                    prcp = hru_clim[[1]]$prcp, 
#                    tavg = hru_clim[[1]]$tavg, 
#                    elev = hru_elev[[1]],
#                    dayOfYear = sim_dayOfYear,
#                    lastDayOfYear = sim_lastDayOfYear)
# 
# ### SAC MODULE
# simflow <- sacSma(par = par_sacsma[1,], 
#                   prcp = hru_clim[[1]]$prcp, 
#                   pet = pet, lat = hru_lat, elev = hru_elev)
# 
# ### ROUTING MODULE
# flow <- routeLohamann(par = par_routLah[1,], flowLength = hru_flowlen[1])


##### FULL BASIN TEST ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

hru_par <- hru_par %>% as.matrix()

par_sacsma   = hru_par[,1:16] 
par_petHamon = hru_par[,17]  
par_snow17   = hru_par[,18:27] 
par_routLah  = hru_par[,28:31] 
tavg = hru_tavg 
prcp = hru_prcp 
lat  = hru_lat 
elev = hru_elev 
area = hru_area
flowLength = hru_flowlen
dayOfYear = sim_dayOfYear 
lastDayOfYear = sim_lastDayOfYear


par_sacsma[1,]

hru_tavg <- lapply(hru_clim, "[[", "tavg")
hru_prcp <- lapply(hru_clim, "[[", "prcp")

flowR <- sacHydroSim(par_sacsma   = hru_par[,1:16],
                    par_petHamon = hru_par[,17],
                    par_snow17   = hru_par[,18:27],
                    par_routLah  = hru_par[,28:31],
                    tavg = hru_tavg, 
                    prcp = hru_prcp, 
                    lat  = hru_lat, 
                    elev = hru_elev, 
                    area = hru_area,
                    flowLength = hru_flowlen,
                    dayOfYear = sim_dayOfYear, 
                    lastDayOfYear = sim_lastDayOfYear)

### Read-in matlab flow 
matlab_dir <- "C:\\Users\\Umit\\Google Drive\\Research\\fromOthers\\Models\\SACSMA (matlab)"
flowM_dat <- read.table(paste0(matlab_dir, "\\sacsma_stcroix_flow_matlab.txt"))
flowM <- flowM_dat[[1]]

df <- data_frame(date = sim_date, matlab = flowM, R = flowR)

dim(flowM)

#flowM_nosnow <- read.table(paste0(matlab_dir, "\sacsma_stcroix_flow_matlab.txt"))




 

################################################################################
################################################################################

########## PROFILING


l1 <- profvis::profvis({
  out <- sacBasinFlow(sac_par = sac_par, snow_par = snow_par, pet_par = pet_par, 
                      rout_par = rout_par, temp = temp, prcp = prcp, lat = lat, 
                      elev = elev, flowLen = flowLen, dayOfYear = dayOfYear, 
                      lastDayOfYear = lastDayOfYear)})



Ta <- -3.35
a <- microbenchmark(if(Ta < 0) T_snow_new <- Ta else T_snow_new <- 0) 
b <- microbenchmark(T_snow_new <- ifelse(Ta < 0, Ta, 0))

simflow <- c(0,-1,2,3,4,0, -3, 4, -8)
simflow <- ifelse(simflow < 0, 0, simflow)

