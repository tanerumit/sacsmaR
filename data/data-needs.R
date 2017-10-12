
# parameter set #31 parameters for each hru
# initial parameters for each hru

########## Vector data (single value for each HRU)
# flowlength OK
# area  OK
# lat   OK
# elev  OK

########## list data (time-series data for each HRU)
# daily prcp 
# daily tavg

########## matrix data 
# model parameters for each hru OK
# model initial parameters for each hru  OK

########## Other data
# dayOfYear vector OK

#petHamon <- function(coeff, temp, lat, dayOfyear)
#snow17   <- function(par, prcp, temp, elev, statesInput, TIME, verbose = FALSE)
#sacSim   <- function(par, par0, prcp, pet, lat, elev, verbose = FALSE)
#routeLohamann <- function(par, flowLength, UH_DAY = 96, KE  = 12)

############### PREPARE THE DATASET
base_dir <- "C:\\Users\\Umit\\Google Drive\\Research\\Projects\\StCroix-Basin\\"

setwd(base_dir)

# HRU info file for the entire basin (Lat, long, area (%) elev (m), hru_id)
grid_info <- read_csv("./input/coupled/HRUinfo_basin.csv") %>% as.matrix()
grid_key  <- grid_info %>% as_tibble() %>% select(HRU_Lat, HRU_Lon, HRU_id)
num_hru   <- nrow(grid_info)

# Elevation data (vector: 1 x 50 hrus)
hru_elev <- grid_info[,"HRU_Elev(m)"]
hru_area <- grid_info[,"HRU_Area"] / 100
hru_lat  <- grid_info[,"HRU_Lat"]
hru_flowLength <- read_csv("./input/HRUinfo2_subcat.csv") %>% 
  pull(`HRU_FlowLen(m)`)
  
# Calibration parameters
calib_dat <- read.table('./input/coupled/cal_gpcc_trial1.txt', header = T)
calib_par <- calib_dat %>% 
  left_join(grid_key, by = c("HRU_Lat", "HRU_Lon")) %>%
  arrange(HRU_id) %>% 
  select(uztwm:Diff) %>% as.matrix()

# Specify parameters for each model
sac_par  <- calib_par[,1:16]
pet_par  <- calib_par[,17]
snow_par <- calib_par[,18:27]
rout_par <- calib_par[,28:31]

# Sac-sma initial state parameters
sac_iniStates <- c(0,0,5,5,5,0)

# Snow model initial state parameters
snow_iniStates <- c(0,0,0,0)

# Historical climate period
str_date <- as.Date("1990/01/1")
end_date <- as.Date("2010/12/31")
seq_date <- seq.Date(str_date, end_date, by = "day") 
sim_date <- seq.Date(str_date, end_date, by = "day") 
sim_doy  <- lubridate::yday(sim_date)
sim_per  <- length(sim_date) # Simulation period
sim_dayOfYear <- sim_doy

# HRU - climate time-series
# HRU info file for the entire basin (Lat, long, area (%) elev (m), hru_id)
grid_info <- read_csv("./input/HRUinfo_basin.csv") %>% as.matrix()
grid_key  <- grid_info %>% as.tibble() %>% select(HRU_Lat, HRU_Lon, HRU_id)
num_hru   <- nrow(grid_info)

hru_area <- read.table("./input/HRUinfo_basin.txt", header = TRUE) %>% 
  select(HRU_Lat = lat, HRU_Lon = lon, areaPer = area.) %>%
  right_join(grid_key, by = c("HRU_Lat", "HRU_Lon")) %>%
  mutate(areaPer = areaPer / 100) %>% as_tibble()

hru_clim <- list()
for (n in 1:num_hru) {
  lat  <- formatC(grid_info[n,1], format = 'f', flag='0', digits = 4)
  lon  <- formatC(grid_info[n,2], format = 'f', flag='0', digits = 4)
  hru_clim[[n]] <- read.table(
    paste0("./input/HRU_clim/data_",lat, "_", lon)) %>% as.matrix()
  
  # Find starting and ending indices based on the specified dates
  grid_date <- as.Date(paste(hru_clim[[n]][,1], hru_clim[[n]][,2], hru_clim[[n]][,3], sep ="/"))
  grid_ind  <- which(str_date == grid_date):which(end_date == grid_date) 
  
  hru_clim[[n]] <- hru_clim[[n]][grid_ind,] %>% as_tibble()
  colnames(hru_clim[[n]]) <- c("year", "month", "day", "prcp", "tavg")
}

save(list = c('hru_elev', 'hru_area', 'hru_flowLength','hru_lat',
              'pet_par', 'sac_par','rout_par', 'snow_par',
              'calib_par', 'sac_iniStates',
              'snow_iniStates', 'sim_dayOfYear', 'grid_info', 'hru_clim'),
     file = 'sac_testdata.Rdata')




