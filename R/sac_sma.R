
# SAC SMA - hydrology simulation model

### Inputs
sac_sma <- function(S_Date, E_Date, Prcp, Tavg, Basin_Lat, Basin_Elev, Par, 
                    IniState, flag_snowmodule = 0) {

  # PARAMETERIZATION -----------------------------------------------------------
  
  # Capacity Thresholds
  uztwm  <-  Par[1]    # Upper zone tension water capacity [mm]
  uzfwm  <-  Par[2]    # Upper zone free water capacity [mm]
  lztwm  <-  Par[3]    # Lower zone tension water capacity [mm]
  lzfpm  <-  Par[4]    # Lower zone primary free water capacity [mm]
  lzfsm  <-  Par[5]    # Lower zone supplementary free water capacity [mm]
  
  # Recession Parameters
  uzk    <-  Par[6]    # Upper zone free water lateral depletion rate [1/day]
  lzpk   <-  Par[7]    # Lower zone primary free water depletion rate [1/day]
  lzsk   <-  Par[8]    # Lower zone supplementary free water depletion rate [1/day]
  
  # Percolation
  zperc  <-  Par[9]    # Percolation demand scale parameter [-]
  rexp   <-  Par[10]   # Percolation demand shape parameter [-]: exponent of the percolation equation
  pfree  <-  Par[11]   # Percolating water split parameter (decimal fraction): Fraction of water percolating from upper zone directly to lower zone free water storage
  
  # Impervious area
  pctim  <-  Par[12]   # Impervious fraction of the watershed area [decimal fraction)
  adimp  <-  Par[13]   # Additional impervious areas (decimal fraction)
  
  # Others
  riva   <-  Par[14]   # Riparian vegetation area (decimal fraction)
  side   <-  Par[15]   # The ratio of deep recharge to channel base flow [-]
  rserv  <-  Par[16]   # Fraction of lower zone free water not transferrable to lower zone tension water (decimal fraction)
  
  # Potential Evapotranspiration Computation
  coeff  <-  Par[17]   # Proportionality Coefficient of Hamon Potential Evapotranspiration
  
  # Snow17
  SCF    <-  Par[18]   # Multiplying factor which adjusts Precipitation (accounts for gage snow catch deficiencies)
  PXTEMP <-  Par[19]   # Temperature that separates rain from snow, deg C
  MFMAX  <-  Par[20]   # Maximum melt factor during non-rain periods - assumed to occur on June 21
  MFMIN  <-  Par[21]   # Minimum melt factor during non-rain periods - assumed to occur on Dec 21
  UADJ   <-  Par[22]   # Average wind function during rain-on-snow periods
  MBASE  <-  Par[23]   # Base temperature for snowmelt computations during non-rain periods, deg C
  TIPM   <-  Par[24]   # Antecedent temperature index parameter (0.01 to 1.0)
  PLWHC  <-  Par[25]   # Percent liquid water holding capacity (maximum value allowed is 0.4)
  NMF    <-  Par[26]   # Maximum negative melt factor
  DAYGM  <-  Par[27]   # A constant daily rate of melt at the soil-snow interface
  snowpar <- c(SCF, PXTEMP, MFMAX, MFMIN, UADJ, MBASE, TIPM, PLWHC, NMF, DAYGM)
  
  # Initial Storage States (SAC-SMA)
  uztwc <- IniState[1]     # Upper zone tension water storage
  uzfwc <- IniState[2]     # Upper zone free water storage
  lztwc <- IniState[3]     # Lower zone tension water storage
  lzfsc <- IniState[4]     # Lower zone supplementary free water storage
  lzfpc <- IniState[5]     # Upper zone primary free water storage
  adimc <- IniState[6]     # Additional impervious area storage
  
  # SNOW17
  W_i     <- IniState[7]   # Accumulated water equivalent of the ice portion of the snow cover (mm)
  ATI     <- IniState[8]   # Antecedent Temperature Index, deg C
  W_q     <- IniState[9]   # Liquid water held by the snow (mm)
  Deficit <- IniState[10]  # Heat Deficit, also known as NEGHS, Negative Heat Storage
  snow_state <- c(W_i, ATI, W_q, Deficit) 

  # PREPARE RESERVOIR ARRAYS ---------------------------------------------------
  
  # RESERVOIR STATE ARRAY INITIALIZATION
  empty_vec <- vector(mode = "numeric", length = length(Prcp))
  # Upper zone states
  uztwc_tot <- empty_vec  # State of Upper zone tension water storage [mm]
  uzfwc_tot <- empty_vec   # State of Upper zone free water storage [mm]
  
  # Lower zone states
  lztwc_tot <- empty_vec   # State of Lower zone tension water storage [mm]
  lzfsc_tot <- empty_vec   # State of Lower zone free water supplementary storage [mm]
  lzfpc_tot <- empty_vec   # State of Lower zone free water primary storage [mm]
  
  # Additional impervious zone states
  adimc_tot <- empty_vec   # State of additional impervious area storages [mm]
  
  # MODEL OUPUT ARRAY INITIALIZATION
  simflow  <- empty_vec  # Simulated Streamflow
  tet_tot  <- empty_vec  # Simulated Actual Evapotranspiration
  base_tot <- empty_vec  # Simulated Base Flow
  surf_tot <- empty_vec  # Simulated Surface&Subsurface water flow
  SWE_tot  <- empty_vec  # Simulated Snow Water Equivalent (SWE)
  
  #CALCULATE PET (HAMON EQUATION) ----------------------------------------------
  #Matlab equivalent
  pet <- petHamon2(doy = yday(seq.Date(S_Date, E_Date, by = "day")), 
                   coeff = coeff, basin_lat = as.numeric(Basin_Lat), Tavg = Tavg)    
  #pet <- petHamon(doy = yday(seq.Date(S_Date, E_Date, by = "day")), 
  #                basin_lat = as.numeric(Basin_Lat), Tavg = Tavg, KPEC = coeff)
  
  #PERFORM HYDROLOGY SIMULATION USING THE SAC-SMA ------------------------------   
  thres_zero  <- 0.00001      # Threshold to be considered as zero
  parea       <- 1 - adimp - pctim
  
  
  #### FIND JULIAN DAY??? 
  #jdate_mat = nan(e_datenum-s_datenum+1,1);
  #datemat   =  datevec(s_datenum : e_datenum);
  #juliandd  =  [datemat(:,1:3) jdate_mat];  % Snow17 input [yr mon day juliandate]
  
  
  for (i in 1:length(pet)) {
    
    # SNOW COMPONENT (IGNORE FOR NOW!!!!!!!!) ***********************
    # Precipitation adjusted by snow process (SNOW17)
    if(flag_snowmodule == 1) {
      
    #  snow_outputs <- SNOW17(pars = snowpar, prcp = Prcp[i], temp = Tavg[i], 
    #    elev = Basin_Elev, states_input = snow_state, TIME = juliandd(i,:))
      
      SWE <- snow_outputs$SWE
      meltNrain  <- snow_outputs$meltNrain
      snow_state <- snow_outputs$snow_state
      
      SWE_tot[i] <- SWE
      pr <- meltNrain
      
    } else {
      pr = Prcp[i]
    }
    
    ### COMPUTE ET LOSS FOR A TIME INTERVAL ++++++++++++++++++++++++++++++++++++
    
    # The amount of maximum ET given by input potential evapotranspiration
    edmnd = pet[i]  
    
    ## Compute for different compnents...
    # ET(1), ET from Upper zone tension water storage
    et1 <- edmnd * uztwc/uztwm
    red <- edmnd - et1  # residual ET demand
    uztwc <- uztwc - et1
    
    # ET(2), ET from upper zone free water storage
    et2 <- 0
    
    # in case et1 > uztws, no water in the upper tension water storage
    if (uztwc <= 0) {
      et1 <- et1 + uztwc #et1 = uztwc
      uztwc <- 0
      red <- edmnd - et1
    
      # when upper zone free water content is less than residual ET
      if (uzfwc < red) {
        
        # all content at upper zone free water zone will be gone as ET
        et2 <- uzfwc
        uzfwc <- 0
        red <- red - et2
        if (uztwc < thres_zero) {uztwc <- 0}
        if (uzfwc < thres_zero) {uzfwc <- 0}  
      
      # when upper zone free water content is more than residual ET
      } else {
        et2 <- red  # all residual ET will be gone as ET
        uzfwc <- uzfwc - et2
        red <- 0 
      }
      
    # in case et1 <= uztws, all maximum et (et1) are consumed at uztwc, 
    # so no et from uzfwc (et2=0)
    } else {
      
      # There's possibility that upper zone free water ratio exceeds 
      #upper zone tension water ratio. If so, free water is transferred to 
      #tension water storage

      if((uztwc / uztwm) < (uzfwc / uzfwm)) {
        uzrat = (uztwc + uzfwc) / (uztwm + uzfwm)
        uztwc = uztwm * uzrat
        uzfwc = uzfwm * uzrat
      }
      
      if(uztwc < thres_zero) {uztwc = 0}
      if(uzfwc < thres_zero) {uzfwc = 0}
      
    }
      
    # ET(3), ET from Lower zone tension water storage when residual ET > 0
    et3 <- red * lztwc / (uztwm + lztwm) #residual ET is always bigger than ET(3)
    lztwc <- lztwc - et3 
    
    # if lztwc is less than zero, et3 cannot exceed lztws
    if(lztwc < 0) {
      et3   <- et3 + lztwc  # et3 = lztwc  
      lztwc <- 0
    }
    #????? not updated red is used later for ET(5) calculation
    
    
    # Water resupply from Lower free water storages to Lower tension water storage
    saved  <- rserv * (lzfpm + lzfsm)
    ratlzt <- lztwc / lztwm
    ratlz  <- (lztwc + lzfpc + lzfsc - saved) / (lztwm + lzfpm + lzfsm - saved)
    
    # water is first taken from supplementary water storage for resupply
    if (ratlzt < ratlz) {
      
      del <- (ratlz - ratlzt) * lztwm
      lztwc <- lztwc + del  # Transfer water from lzfss to lztws
      lzfsc <- lzfsc - del
      
      # if tranfer exceeds lzfsc then remainder comes from lzfps
      if(lzfsc < 0) {   
        lzfpc <- lzfpc + lzfsc
        lzfsc <- 0
      }
    }

    if(lztwc < thres_zero) {lztwc <- 0}
  
    # ET(5), ET from additional impervious (ADIMP) area
    # ????? no idea where this come from, I think there's a possibility that et5 can be negative values
    et5   <- et1 + (red + et2) * (adimc - et1 - uztwc) / (uztwm + lztwm)  
    adimc <- adimc - et5
    if(adimc < 0) { 
        #et5 cannot exceed adims
        et5 <- et5 + adimc # et5 = adimc
        adimc <- 0
    }
    et5 <- et5 * adimp
    
    # Time interval available moisture in excess of uztw requirements
    twx <- pr + uztwc - uztwm
    
    # all moisture held in uztw- no excess
    if(twx < 0) {
      uztwc <- uztwc + pr
      twx <- 0
    # moisture available in excess of uztw storage
    } else {
      uztwc = uztwm 
    }   

    
    # for now twx is excess rainfall after filling the uztwc
    adimc <- adimc + pr - twx
    
    # Compute Impervious Area Runoff
    roimp <- pr * pctim 
    
    # Initialize time interval sums
    sbf   <- 0  # Sum of total baseflow(from primary and supplemental storages)
    ssur  <- 0  # Sum of surface runoff
    sif   <- 0  # Sum of interflow
    sperc <- 0  # Time interval summation of percolation
    sdro  <- 0  # Sum of direct runoff from the additional impervious area
    
    # Determine computational time increments for the basic time interval
    ninc <- floor(1.0 + 0.2*(uzfwc+twx))  # Number of time increments that interval is divided into for further soil-moisture accountng
    dinc <- 1.0 / ninc                    # Length of each increment in days
    pinc <- twx / ninc                    # Amount of available moisture for each increment
    
    # Compute free water depletion fractions for the time increment 
    #(basic depletions are for one day)
    duz   <- 1 - (1 - uzk)^dinc
    dlzp  <- 1 - (1 - lzpk)^dinc
    dlzs  <- 1 - (1 - lzsk)^dinc
    
    # Start incremental for-loop for the time interval
    for(n in 1:ninc) {
      
      adsur <- 0 # Amount of surface runoff. This will be updated.
      
      # Compute direct runoff from adimp area
      ratio <- (adimc - uztwc) / lztwm
      if(ratio < 0) {ratio <- 0} 
      
      # Amount of direct runoff from the additional impervious area     
      addro <- pinc*(ratio^2)
      
      # Compute baseflow and keep track of time interval sum
      # Baseflow from free water primary storage
      bf_p <- lzfpc * dlzp 
      lzfpc <- lzfpc - bf_p
      if(lzfpc <= 0.0001) {
        
        bf_p  <- bf_p + lzfpc
        lzfpc <- 0
      }
      
      sbf <- sbf + bf_p
      
      # Baseflow from free water supplemental storage
      bf_s  <- lzfsc * dlzs  
      lzfsc <- lzfsc - bf_s
      if (lzfsc <= 0.0001) {
        bf_s <- bf_s + lzfsc
        lzfsc <- 0
      }
      
      # Total Baseflow from primary and supplemental storages
      sbf <- sbf + bf_s 
      
      # Compute PERCOLATION- if no water available then skip.
      if((pinc + uzfwc) <= 0.01) {
        uzfwc <- uzfwc + pinc
      } else {
        
        # Limiting drainage rate from the combined saturated lower zone storages
        percm <- lzfpm * dlzp + lzfsm * dlzs 
        perc <- percm * uzfwc / uzfwm

        # DEFR is the lower zone moisture deficiency ratio
        defr <- 1.0 - (lztwc + lzfpc + lzfsc)/(lztwm + lzfpm + lzfsm) 
  
        if(defr < 0) {defr <- 0} 
      
        perc <- perc * (1.0 + zperc * (defr^rexp))
      
        # Note. . . percolation occurs from uzfws before pav is added
      
        # Percolation rate exceeds uzfws
        if(perc >= uzfwc) {perc <- uzfwc}
      
        uzfwc <- uzfwc - perc    # Percolation rate is less than uzfws.
      
        # Check to see if percolation exceeds lower zone deficiency.
        check <- lztwc + lzfpc + lzfsc + perc - lztwm - lzfpm - lzfsm
        if(check > 0) {
          perc <- perc - check
          uzfwc <- uzfwc + check
        }    
  
        # SPERC is the time interval summation of PERC
        sperc <- sperc + perc  
      
      
        # Compute interflow and keep track of time interval sum. Note that PINC has not yet been added.
        del <- uzfwc * duz # The amount of interflow
        sif <- sif + del
        uzfwc <- uzfwc - del
        
        # Distribute percolated water into the lower zones. Tension water
        # must be filled first except for the PFREE area. PERCT is
        # percolation to tension water and PERCF is percolation going to
        # free water.
    
        perct <- perc * (1.0 - pfree)  # Percolation going to the tension water storage
        
        if((perct + lztwc) <= lztwm) {
          
          lztwc <- lztwc + perct
          percf <- 0.0 # Pecolation going to th lower zone free water storages
        
        } else {
          
          percf <- lztwc + perct - lztwm
          lztwc <- lztwm  
          
        }

      # Distribute percolation in excess of tension requirements among the free water storages.
      percf <- percf + (perc * pfree)
      if(percf != 0) {
      
        # Relative size of the primary storage as compared with total lower zone free water storages.
        hpl <- lzfpm / (lzfpm + lzfsm) 
        
        # Relative fullness of each storage.
        ratlp <- lzfpc / lzfpm
        ratls <- lzfsc / lzfsm
        
        # The fraction going to primary
        fracp <- hpl * 2.0 * (1.0 - ratlp) / (1.0 - ratlp + 1.0 - ratls) 
        if(fracp > 1.0) {fracp <- 1.0} 
        
        percp <- percf * fracp # Amount of the excess percolation going to primary
        percs <- percf - percp # Amount of the excess percolation going to supplemental
        lzfsc <- lzfsc + percs
        
        
        if(lzfsc > lzfsm) {
          percs <- percs - lzfsc + lzfsm
          lzfsc <- lzfsm
        }

        lzfpc <- lzfpc + percf - percs
        
        # Check to make sure lzfps does not exceed lzfpm
        if(lzfpc >= lzfpm) {
          excess <- lzfpc - lzfpm
          lztwc <- lztwc + excess
          lzfpc <- lzfpm
        }   
      }
        

      # Distribute PINC between uzfws and surface runoff
      if(pinc != 0) {
        
        # check if pinc exceeds uzfwm
        if((pinc + uzfwc) <= uzfwm) {
          
          uzfwc <- uzfwc + pinc  # no surface runoff
        } else {
          sur <- pinc + uzfwc - uzfwm # Surface runoff
          uzfwc <- uzfwm
          
          ssur = ssur + (sur * parea)
          
          # ADSUR is the amount of surface runoff which comes from
          # that portion of adimp which is not currently generating
          # direct runoff. ADDRO/PINC is the fraction of adimp
          # currently generating direct runoff.
          adsur = sur * (1.0 - addro / pinc) 
          ssur = ssur + adsur * adimp
          
        }
      }  
    }
    
      adimc <- adimc + pinc - addro - adsur
      if(adimc > (uztwm + lztwm)) {
        addro = addro + adimc - (uztwm + lztwm)
        adimc = uztwm + lztwm
      }

      # Direct runoff from the additional impervious area
      sdro  = sdro + (addro * adimp) 
      
      if(adimc < thres_zero) {adimc <- 0}
      
    } # END of incremental for loop
      
      
    # Compute sums and adjust runoff amounts by the area over which they are generated.
    
    # EUSED is the ET from PAREA which is 1.0 - adimp - pctim
    eused <- et1 + et2 + et3 
    sif <- sif * parea 
    
    # Separate channel component of baseflow from the non-channel component
    tbf <- sbf * parea   # TBF is the total baseflow
    bfcc <- tbf * (1.0 / (1.0 + side))     # BFCC is baseflow, channel component  
    
    # Ground flow and Surface flow
    base <- bfcc                       # Baseflow and Interflow are considered as Ground inflow to the channel
    surf <- roimp + sdro + ssur + sif  # Surface flow consists of Direct runoff and Surface inflow to the channel
    
    # ET(4)- ET from riparian vegetation.
    et4 <- (edmnd - eused) * riva  # no effect if riva is set to zero
    
    
    # Compute total evapotransporation - TET
    eused <- eused * parea 
    tet <- eused + et4 + et5 
    
    
    # Check that adims >= uztws
    if(adimc < uztwc) {adimc <- uztwc} 

    # Total inflow to channel for a timestep 
    tot_outflow <- surf + base 
    
    tot_outflow <- tot_outflow - et4 
    if(tot_outflow < 0) {tot_outflow <- 0} 

    # OUTPUTS FOR A TIME STEP
    
    # Main output: Watershed out flow
    simflow[i] <- tot_outflow 
    
    # Optional outputs 
    tet_tot[i]   <- tet 
    surf_tot[i]  <- surf 
    base_tot[i]  <- base 
    uztwc_tot[i] <- uztwc 
    uzfwc_tot[i] <- uzfwc 
    lztwc_tot[i] <- lztwc 
    lzfpc_tot[i] <- lzfpc 
    lzfsc_tot[i] <- lzfsc 
    adimc_tot[i] <- adimc
    
  }
  
  return(simflow)

}