#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param coeff PARAM_DESCRIPTION
#' @param tavg PARAM_DESCRIPTION
#' @param lat PARAM_DESCRIPTION
#' @param dayOfyear PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname petHamon
#' @export 
petHamon <- function(coeff, tavg, lat, dayOfYear) {
  
  var_theta   <- 0.2163108 + 2*atan(0.9671396 * tan(0.0086 * (dayOfYear-186))) 
  var_pi      <- asin(0.39795 * cos(var_theta)) 
  daylighthr  = 24-24/pi * acos((sin(0.8333*pi/180) + sin(lat*pi/180) * sin(var_pi))/(cos(lat*pi/180)*cos(var_pi))) 
  esat <- 0.611 * exp(17.27 * tavg/(237.3 + tavg))
  
  return(coeff * 29.8 * daylighthr*(esat/(tavg + 273.2)))
}
