
#' Title
#'
#' @param par 
#' @param tavg 
#' @param lat 
#' @param jdate 
#'
#' @return
#' @export
#'
#' @examples
hamon <- function(par, tavg, lat, jdate) {
  
  var_theta <- 0.2163108 + 2 * atan(0.9671396 * tan(0.0086 * (jdate - 186)))
  var_pi <- asin(0.39795 * cos(var_theta))
  daylighthr <- 24 - 24/pi * acos((sin(0.8333 * pi/180) + sin(lat * pi/180) * sin(var_pi))/(cos(lat * 
    pi/180) * cos(var_pi)))
  
  esat <- 0.611 * exp(17.27 * tavg/(237.3 + tavg))
  
  return(par * 29.8 * daylighthr * (esat/(tavg + 273.2)))
  
}