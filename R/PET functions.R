

dayLightHrs <- function(doy, basin_lat) {
  
  p <- asin(0.39795 * cos(0.2163108 + 2 * atan(0.9671396 * tan(0.00860 * (doy - 186)))))
  d <- 24 - (24 / pi) * acos((sin(0.8333 * pi / 180) + sin(basin_lat * pi / 180) * sin(p)) / 
                               cos(basin_lat * pi / 180) * cos(p))
  return(d)
}

hamonPET <- function (tavg, Ld, KPEC) 
{
  ESAT <- 6.108 * exp(17.26939 * tavg/(tavg + 273.3))
  RHOSAT <- 216.7 * ESAT/(tavg + 273.3)
  PET_daily <- 0.1651 * Ld * RHOSAT * KPEC
  return(PET_daily)
}

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



