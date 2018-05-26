
#' @title Hamon Potential Evapotranspiration Equation
#' @description The Hamon method is also considered as one of the simplest estimates 
#'     of potential Evapotranspiration. 
#' @param par proportionality coefficient (unitless)
#' @param tavg  vector of mean daily temperature (deg C) 
#' @param lat latitude ()
#' @param jday a day number of the year (julian day of the year)
#' @return outputs potential evapotranspiration (mm day-1)
#' @details For details see Haith and Shoemaker (1987) 
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1deg
#'  }
#' }
#' @rdname hamon
#' @export 
hamon <- function(par, tavg, lat, jday) {
  
  var_theta <- 0.2163108 + 2 * atan(0.9671396 * tan(0.0086 * (jday - 186)))
  var_pi <- asin(0.39795 * cos(var_theta))
  daylighthr <- 24 - 24/pi * acos((sin(0.8333 * pi/180) + sin(lat * pi/180) * sin(var_pi))/(cos(lat * 
    pi/180) * cos(var_pi)))
  
  esat <- 0.611 * exp(17.27 * tavg/(237.3 + tavg))
  
  return(par * 29.8 * daylighthr * (esat/(tavg + 273.2)))
  
}