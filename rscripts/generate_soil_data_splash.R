library(rsplash)

# GYF
sand = 65.25
clay = 21.50
silt = 13.25
om = 2.37
bulk_density = 1.04 # gm cc-1
depth = 2.50
lat = 5.283333
elev = 0

soil_props_gyf = soil_hydro(sand=sand, clay = clay, OM=om, bd=bulk_density, fgravel = 0)
soil_props_gyf$depth = depth
soil_props_gyf$site = "gyf"

# TNF
sand = 37.27
clay = 60.09
silt = 2.64
om = 2.54
bulk_density = 1.125 # 1.04 # gm cc-1 # FIXME: Actually should be 1.125, 1.04 used in DroughtMIP submission
depth = 16.1
lat = 2.85
elev = 0

soil_props_tnf = soil_hydro(sand=sand, clay = clay, OM=om, bd=bulk_density, fgravel = 0)
soil_props_tnf$depth = depth
soil_props_tnf$site = "tnf"

library(tidyverse)

soil_info_array = 
  rbind(soil_props_gyf |> as.data.frame(), 
        soil_props_tnf |> as.data.frame()) %>% 
  mutate(
    SAT = SAT*depth*1000,
    WP = WP*depth*1000,
    FC = FC*depth*1000,
    RES = RES*depth*1000,
    Wmax = theta_c*depth*1000,
    KWm = AWC*depth*1000,
    lambda = 1/B,
    bub_press = bubbling_p,
    Au = 0,
    Ai = 250*250,
    cellin = 3,
    cellout = 3,
    AI = 0
  ) %>% 
  select(
    SAT, WP, FC, Ksat, lambda, depth, bub_press, RES, Au, Ai, cellin, cellout, AI
  )
    
soil_info_array

