
## Routing Model developed by D. Lohmann

route_lohamann <- function(pars, flowlen, KE, UH_DAY) {
  
  # Routing model parameters
  NumRes <- pars[1]       # Grid Unit Hydrograph parameter (number of linear reservoirs)
  K      <- pars[2]       # Grid Unit Hydrograph parameter (reservoir storage constant)
  VELO   <- pars[3]       # wave velocity in the linearized Saint-Venant equation(m/s)
  DIFF   <- pars[4]       # diffusivity in the linearized Saint-Venant equation(m2/s)
  
  DT    <-  3600          # Time step in second for solving Saint-Venant equation. This will affect TMAX
  TMAX  <-  UH_DAY * 24   # Base time of river routing UH in hour because DT is for an hour
  LE		<- 48*50
  
  # Green's function calculation to solve the linearized Saint-Venant
  UH_DAILY <- vector(mode = "numeric", length = UH_DAY)
  if (flowlen == 0) {
    UH_DAILY[1] <- 1 
  } else {
    t <- 0 
    uhm_grid <- vector(mode = "numeric",length = LE)
    for(k in 1:LE) {
      t <- t + DT 
      pot <- ((VELO*t-flowlen)^2)/(4*DIFF*t)
      if(pot <= 69) {
        H <- flowlen/(2*t*sqrt(pi*t*DIFF))*exp(-pot) 
      } else {
        H <- 0 	   
      }
      uhm_grid[k] <- H 
    }
    
    if(sum(uhm_grid) == 0) {
      uhm_grid[1] <-  1.0 
    } else {
      uhm_grid <- uhm_grid/sum(uhm_grid)  
    }
    UHM <- uhm_grid 
      
    # Derive Daily River Impulse Response Function(Green's function)
    FR <- matrix(data = 0, nrow = TMAX, ncol = 2) 
    FR[1:24,1] <- 1/24   

    
    
    
    
        
    #Create the master table for all t * L combinations
    start_time <- Sys.time() 
    loop_tL <- expand.grid(t = 1:TMAX, L = 1:(TMAX+24)) %>%
      as_data_frame() %>%
      filter(t-L > 0) %>%
      group_by(t) %>%
      arrange(t, L) %>%
      mutate(x = cumsum(FR[t-L,1] * UHM[L])) %>%
      slice(n())
 
    FR[,2] <- c(FR[1,2],loop_tL$x)
    Sys.time() - start_time

    UH_DAILY <- sapply(1:UH_DAY, function(t) sum(FR[(24*t-23):(24*t),2]))

    # Nested for loop - Too-slow!
    # Derive Daily River Impulse Response Function(Green's function)
    #FR <- matrix(data = 0, nrow = TMAX, ncol = 2) 
    #FR[1:24,1] <- 1/24   
    #
    #start_time <- Sys.time()
    #for(t in 1:TMAX) {
    #  for(L in 1:(TMAX+24)) {
    #    if((t-L) > 0) {
    #      FR[t,2] <- FR[t,2] + FR[t-L,1] * UHM[L]  
    #    }
    #  }
    #}
    
    #for(t in 1:UH_DAY) {
    #  UH_DAILY[t] <- sum(FR[(24*t-23):(24*t),2]) 
    #}
  
  }
  
  # HRU's Unit Hydrograph represented by Gamma distribution function
  # Matlab code:  @(x) gampdf(x,NumRes,1/K) #(x: values, A = shape parameters, B = scale parameters)
  hruh_fun <- function(x) {
          return(dgamma(x, shape = NumRes, scale = 1/K))
  }

  UH_HRU <- vector(mode = "numeric", length = 12)
  for(i in 1:KE) {
    #integration in matlab: integral(fun,xmin,xmax)
    #integration in R: integrate(f, lower, upper, ...)
    UH_HRU[i] <- integrate(hruh_fun, 24*(i-1), 24*i)$value 
  }
  
  # River UH
  UH_R <- vector(mode = "numeric", length = KE+UH_DAY-1)
  for(k in 1:KE) {
    for(u in 1:UH_DAY) {
      UH_R[k+u-1] <- UH_R[k+u-1] + UH_HRU[k] * UH_DAILY[u] 
    }
  }
  
  UH_R <- UH_R/sum(UH_R) 
  
  return(UH_R)
}