
#///////////////////////////////////////////////////////////////////////////////
# hydrology component: generates streamflow from given conditions
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

#Inputs to sac_sma function


sac_sma <- function(S_Date, E_Date, Prcp, Tavg, Basin_Lat, Basin_Elev, Par, IniState, flag_snowmodule) {
    
  # S_date = start date, as scalar date object [yyyy/mm/dd]
  # E_date = end date, as scalar date object [yyyy/mm/dd]
  
  
  
  
  
  # Prcp = time series of precipitation
  # Basin_lat = basin latitude
  # Basin_Elev = basin elevation
  # Par = Parameter set
  # IniState = initial states
  # flag_snowmodule = snow module flag
  
  
  # SAC_SMA PARAMETERS ---------------------------------------------------------
  
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
  
  
  # HAMON ET CALCULATION ---------------------------------------------------------
  
  
  
  
  
  
  
  
  
  
  
  # EXECUTE MODEL FOR EVERY TIME STEP ------------------------------------------
  
  
  # Initial Storage States
  # SAC-SMA
  uztwc <- IniState(1)     # Upper zone tension water storage
  uzfwc <- IniState(2)     # Upper zone free water storage
  lztwc <- IniState(3)     # Lower zone tension water storage
  lzfsc <- IniState(4)     # Lower zone supplementary free water storage
  lzfpc <- IniState(5)     # Upper zone primary free water storage
  adimc <- IniState(6)     # Additional impervious area storage
  # SNOW17
  W_i     <- IniState(7)   # Accumulated water equivalent of the ice portion of the snow cover (mm)
  ATI     <- IniState(8)   # Antecedent Temperature Index, deg C
  W_q     <- IniState(9)   # Liquid water held by the snow (mm)
  Deficit <- IniState(10)  # Heat Deficit, also known as NEGHS, Negative Heat Storage
  
  snow_state <- c(W_i, ATI, W_q, Deficit) 
  
  # RESERVOIR STATE ARRAY INITIALIZATION
  # Upper zone states
  uztwc_tot <- NaN(size(Prcp))   # State of Upper zone tension water storage [mm]
  uzfwc_tot <- NaN(size(Prcp))   # State of Upper zone free water storage [mm]
  # Lower zone states
  lztwc_tot <- NaN(size(Prcp))   # State of Lower zone tension water storage [mm]
  lzfsc_tot <- NaN(size(Prcp))   # State of Lower zone free water supplementary storage [mm]
  lzfpc_tot <- NaN(size(Prcp))   # State of Lower zone free water primary storage [mm]
  # Additional impervious zone states
  adimc_tot <- NaN(size(Prcp))   # State of additional impervious area storages [mm]
  
  # MODEL OUPUT ARRAY INITIALIZATION
  simflow  <- NaN(size(Prcp))  # Simulated Streamflow
  tet_tot  <- NaN(size(Prcp))  # Simulated Actual Evapotranspiration
  base_tot <- NaN(size(Prcp))  # Simulated Base Flow
  surf_tot <- NaN(size(Prcp))  # Simulated Surface&Subsurface water flow
  SWE_tot  <- NaN(size(Prcp))  # Simulated Snow Water Equivalent (SWE)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}


















