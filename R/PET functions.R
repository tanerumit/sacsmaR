
### ORIGINAL +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

petHamon <- function(doy, basin_lat, Tavg, Ld, KPEC) {
  
  #Calculate daylight (hours)
  p <- asin(0.39795 * cos(0.2163108 + 2 * atan(0.9671396 * tan(0.00860 * (doy - 186)))))
  d <- 24 - (24 / pi) * acos((sin(0.8333 * pi / 180) + sin(basin_lat * pi / 180) * sin(p)) / 
                               cos(basin_lat * pi / 180) * cos(p))
  
  #Daylight as fraction of 12 hours
  Ld <- d/12
  
  #Calculate PET 
  ESAT <- 6.108 * exp(17.26939 * Tavg/(Tavg + 273.3))
  RHOSAT <- 216.7 * ESAT/(Tavg + 273.3)
  pet <- 0.1651 * Ld * RHOSAT * KPEC
  return(as.numeric(unlist(pet)))

}


### MATLAB BASED EQUATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++

petHamon2 <- function(doy, coeff, basin_lat, Tavg) {
  
  var_theta   <- 0.2163108 + 2*atan(0.9671396 * tan(0.0086 * (doy-186))) 
  var_pi      <- asin(0.39795 * cos(var_theta)) 
  daylighthr  = 24-24/pi * acos((sin(0.8333*pi/180) + sin(basin_lat*pi/180) * sin(var_pi))/(cos(basin_lat*pi/180)*cos(var_pi))) 
  
  esat <- 0.611 * exp(17.27 * Tavg/(237.3 + Tavg))
  pet  <- coeff * 29.8 * daylighthr*(esat/(Tavg + 273.2))
  return(pet)
}
  
  
  





