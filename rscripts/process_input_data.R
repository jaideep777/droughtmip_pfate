library(tidyverse)
library(reshape2)
library(rphydro)

setwd("~/codes/Drought_MIP/drivers_data")

dir.create("data_gen_figures", showWarnings = F)

# for (site in c("GYF", "TNF")){
# for (scenario in c("DY", "HB")){
    
# site = "TNF"
# scenario = "ND16s1"
site = "GYF"
scenario = "ND13s3"

figures_prefix = paste0("data_gen_figures/", site, "_", scenario)

message("- reading hourly data")
hhdf <- read.delim(paste0("PF-CWM_",site,"_",scenario,"hourly.txt"))

ystart = min(hhdf$year)
yend = max(hhdf$year)

# https://glossary.ametsoc.org/wiki/Gas_constant
# https://en.wikipedia.org/wiki/Gas_constant 

Rd = 287.052874 # Gas constant for dry air [J kg-1 K-1]
Rv = 461.52 # Gas constant for water vapour [J kg-1 K-1]

# REF for SH calculation:
# https://earthscience.stackexchange.com/questions/5076/how-to-calculate-specific-humidity-with-relative-humidity-temperature-and-pres
hhdf <- hhdf |> 
  mutate(esat = purrr::map2_dbl(.x=temp.C., .y=pres.Pa., ~calc_esat(.x,.y))) |>
  mutate(eair = esat*rh.../100) |>
  mutate(vpd = esat-eair) |>
  mutate(sh = eair*Rd/Rv/(pres.Pa. - eair)) |> 
  mutate(time = lubridate::as_datetime(paste0(year, "-", month, "-", day, " ", hour-1, ":30:00"))) |>
  mutate(date = lubridate::as_date(time)) |>
  select(-hour, -day, -month, -year)
                       
message("- downsampling - 24 hr means")
ddf_24hr_mean <-
    hhdf |>
      group_by(date) |>
      summarize_all(.funs = mean)


# Check aggregation
cairo_pdf(filename = paste0(figures_prefix, "_fig1_ddf_24hr_mean.pdf"))
p1 = ddf_24hr_mean %>%
  select(-time) %>%
  melt("date") %>%
  ggplot(aes(y=value, x=date)) +
  geom_line(col="aquamarine4") +
  theme_classic() +
  theme(strip.background = element_rect(color = "white", size = 1))+
  facet_wrap(~variable, scales = "free")
p1 %>% print()
dev.off()


# Aggregate around daily maximum ppfd for acclimating model
# ---------------------------------------------------------
test.3day = hhdf %>% filter(date >= as_date(paste0(floor((ystart+yend)/2),"-06-01")) &
                            date <= as_date(paste0(floor((ystart+yend)/2),"-06-03")) )

aggregate_daily_3hr_maxima = function(df){
  # Get the time at which SW_IN is maximum
  maxppfd <- df %>% filter(shortwave.W.m2. == max(shortwave.W.m2.))
  max_t <- maxppfd$time[1]
  
  # Select times that lie in 3-hr interval around max_t
  df_aroundmax <- df %>% filter(time < (max_t + 1.5*3600) &
                                time > (max_t - 1.5*3600) )
  
  # take mean of selected entries
  df_mean <- df_aroundmax |>
    summarize_all(.funs = mean)
  
  df_mean
}

# Test aggregation
# ----------------
test.3day.3hr = test.3day %>% group_by(date) %>% do(aggregate_daily_3hr_maxima(.)) %>% ungroup()

cairo_pdf(filename = paste0(figures_prefix, "_fig2_3hr_maxima_sample.pdf"))
p2 = test.3day %>% 
  select(-date) %>%
  melt("time") %>%
  mutate(type="hourly") %>%
  rbind(test.3day.3hr %>%
          select(-date) %>%
          melt("time") %>%
          mutate(type="daily")
  ) %>%
  ggplot(aes(y=value, x=as.POSIXct(time))) +
  geom_line(data = . %>% filter(type == "hourly"), col="aquamarine4") +
  geom_point(data = . %>% filter(type == "daily")) +
  theme_classic() +
  theme(strip.background = element_rect(color = "white", size = 1))+
  facet_wrap(~variable, scales = "free")
p2 %>% print()
dev.off()



## --------------------------------------------------------------------
# Apply 3hr maxima aggregation to all data
# ----------------------------------------
message("- downsampling - daily 3-hr means around max ppfd")
ddf_3hr_maxima <- hhdf |>
  group_by(date) |>
  do(aggregate_daily_3hr_maxima(.)) |>
  ungroup()


## --------------------------------------------------------------------
aggregate_daily_daylength = function(df){
  # Get the time at which SW_IN > 0
  pos_ppfd <- df %>% filter(shortwave.W.m2. > 10)

  tmax <- max(pos_ppfd$time)
  tmin <- min(pos_ppfd$time)
  
  # Select times that lie in 3-hr interval around max_t
  df_aroundmax <- df %>% filter(time <= tmax &
                                  time >= tmin )
  
  # take mean of selected entries
  df_mean <- df_aroundmax |>
    summarize_all(.funs = mean) |>
    mutate(daylength = difftime(tmax, tmin, units="hours") |> as.numeric())
  
  df_mean
}

# Test aggregation
# ----------------
test.3day.daylen = test.3day %>% group_by(date) %>% do(aggregate_daily_daylength(.)) %>% ungroup()

cairo_pdf(filename = paste0(figures_prefix, "_fig3_daytime_sample.pdf"))
p3 = test.3day %>% 
  select(-date) %>%
  mutate(daylength = NA) %>%
  melt("time") %>%
  mutate(type="hourly") %>%
  rbind(test.3day.daylen %>%
          select(-date) %>%
          melt("time") %>%
          mutate(type="daily")
  ) %>%
  ggplot(aes(y=value, x=as.POSIXct(time))) +
  geom_line(data = . %>% filter(type == "hourly"), col="aquamarine4") +
  geom_point(data = . %>% filter(type == "daily")) +
  theme_classic() +
  theme(strip.background = element_rect(color = "white", size = 1))+
  facet_wrap(~variable, scales = "free")
p3 %>% print()
dev.off()

#'
## --------------------------------------------------------------------
# Apply daytime mean aggregation to all data
# ------------------------------------------
message("- downsampling FLUXNET format - daytime means")
ddf_daytime_mean <- hhdf |>
  group_by(date) |>
  do(aggregate_daily_daylength(.)) |>
  ungroup()


# Check daylength seasonality
cairo_pdf(filename = paste0(figures_prefix, "_fig4_daylenth_seasonality.pdf"))
p4 = ddf_daytime_mean %>%
  ggplot(aes(y=daylength, x=date)) +
  geom_line()+
  theme_classic()
p4 %>% print()
dev.off()


## --------------------------------------------------------------------
# Calculate daily tmax and tmin from hh data
# ------------------------------------------
message("- downsampling - daily tmax & tmin")
tmaxmin <-
  hhdf |>
  group_by(date) |>
  summarize(
    tmax = max(temp.C.),
    tmin = min(temp.C.)
  )


## -----------------------------------------
# Write processed data
# ------------------------------------------

ddf <- read.delim(paste0("PF-CWM_",site,"_",scenario,"daily.txt"))
ddf = ddf |>
  mutate(date = lubridate::as_date(paste0(year, "-", month, "-", day))) |>
  select(-year, -month, -day)
  
ddf_inst_24hr = 
  ddf_24hr_mean %>% 
  rename(sw24hrmean = shortwave.W.m2.) %>% 
  left_join(ddf) %>% 
  left_join(tmaxmin)

cairo_pdf(filename = paste0(figures_prefix, "_fig5_compare_computed_provided_sw24hr.pdf"))
p3 <- ddf_inst_24hr %>%
  filter(date >= as.Date("2000-1-1") & date < as.Date("2000-12-31")) |>
  ggplot()+
  geom_line(aes(y=sw24hrmean, x=date))+
  geom_line(aes(y=shortwave.W.m2., x=date), col="cyan3")
p3 %>% print()
dev.off()

ddf_inst_daytime = 
  ddf_daytime_mean %>% 
  rename(swdaytimemean = shortwave.W.m2.) %>% 
  left_join(ddf) %>% 
  left_join(tmaxmin)


## inst forcing
ddf_inst_24hr %>% 
  write.csv(file = paste0("PF-CWM_",site,"_",scenario,"_inst_forcing_24hr_means.csv"))

ddf_inst_daytime %>% 
  write.csv(file = paste0("PF-CWM_",site,"_",scenario,"_inst_forcing_daytime_means.csv"))

## Acclim forcing
ddf_3hr_maxima %>% 
  write.csv(file = paste0("PF-CWM_",site,"_",scenario,"_acclim_forcing_3hr_maxima.csv"))

##--------------------------------------------
# A final plot of acclim and inst forcing
# --------------------------------------------

png(filename = paste0(figures_prefix, "_fig6_pfate_drivers.png"), height=5*300, width=7*300, res=300)
p4 = ddf_inst_24hr %>%
  select(-time, -sw24hrmean) %>% 
  melt("date") %>%
  mutate(type="24-hr mean") %>%
  rbind(ddf_3hr_maxima %>%
          select(-time) %>% 
          melt("date") %>%
          mutate(type="3-hr maxima")
  ) %>%
  rbind(ddf_daytime_mean %>%
          select(-time) %>% 
          melt("date") %>%
          mutate(type="daytime mean")
  ) %>%
  group_by(variable, type) %>% 
  mutate(value_smooth = loess(value~as.numeric(date), span=1000/length(value)) %>% fitted()) %>% 
  ggplot(aes(x=date)) +
  geom_line(aes(y=value, group=type, col=type), alpha=0.3) +
  geom_line(aes(y=value_smooth, group=type, col=type), alpha=1) +
  theme_classic() +
  theme(strip.background = element_rect(color = "white", size = 1))+
  facet_wrap(~variable, scales = "free")
p4 %>% print()
dev.off()

# }
# }

### write csv file in old format 
library(tidyverse)
library(reshape2)
library(rphydro)

setwd("~/codes/Drought_MIP/drivers_data")

site = "TNF"
scenario = "ND16s1"

## inst forcing
ddf_inst_24hr = read.csv(file = paste0("PF-CWM_",site,"_",scenario,"_inst_forcing_24hr_means.csv"))
ddf_inst_daytime = read.csv(file = paste0("PF-CWM_",site,"_",scenario,"_inst_forcing_daytime_means.csv"))
ddf_3hr_maxima = read.csv(file = paste0("PF-CWM_",site,"_",scenario,"_acclim_forcing_3hr_maxima.csv"))

forcing_inst = 
  ddf_inst_24hr %>% 
  mutate(par = shortwave.W.m2.*2.04,
         swp = 0.05,
         year = lubridate::year(date),
         month = lubridate::month(date),
         doy = as.Date(date) - as.Date(paste0(year, "-", 1, "-", 1)) + 1,
         decimal_year = year + (doy-1)/365.2425,
         vpd = vpd/100 # convert o hPa
         ) %>%
  select(date, year, month, decimal_year, temp.C., vpd, par, swp, precip.mm.) %>% 
  left_join(
    ddf_3hr_maxima %>% 
      select(date, shortwave.W.m2.) %>% 
      mutate(par_max = shortwave.W.m2.*2.04) %>% 
      select(-shortwave.W.m2.)
  ) %>% 
  select(
    year, month, decimal_year, temp.C.,
    vpd, par, par_max, swp, precip.mm.
  ) 

forcing_inst %>% 
  write.csv(file = paste0("PF-CWM_",site,"_",scenario,"_Plant-FATE_oldformat_24hr_par_max.csv"), row.names = F)

forcing_acclim = 
  ddf_3hr_maxima %>% 
  mutate(par = shortwave.W.m2.*2.04,
         swp = 0.05,
         year = lubridate::year(date),
         month = lubridate::month(date),
         doy = as.Date(date) - as.Date(paste0(year, "-", 1, "-", 1)) + 1,
         decimal_year = year + (doy-1)/365.2425,
         vpd = vpd/100, # convert o hPa
         par_max = shortwave.W.m2.*2.04,
         precip.mm. = 0
  ) %>%
  select(
    year, month, decimal_year, temp.C.,
    vpd, par, par_max, swp, precip.mm.
  ) 

forcing_acclim %>% 
  write.csv(file = paste0("PF-CWM_",site,"_",scenario,"_Plant-FATE_oldformat_3hrmax_par_max.csv"), row.names = F)

p5 = forcing_inst %>%
  select(-year, -month) %>% 
  melt("decimal_year") %>%
  mutate(type="24-hr mean") %>%
  rbind(forcing_acclim %>%
          select(-year, -month) %>% 
          melt("decimal_year") %>%
          mutate(type="3-hr maxima")
  ) %>%
  group_by(variable, type) %>% 
  mutate(value_smooth = loess(value~as.numeric(decimal_year), span=1000/length(value)) %>% fitted()) %>% 
  ggplot(aes(x=decimal_year)) +
  geom_line(aes(y=value, group=type, col=type), alpha=0.3) +
  geom_line(aes(y=value_smooth, group=type, col=type), alpha=1) +
  theme_classic() +
  theme(strip.background = element_rect(color = "white", size = 1))+
  facet_wrap(~variable, scales = "free")
p5 %>% print()


#### Analyse data #### 

setwd("~/codes/Drought_MIP/drivers_data")

site = rep(c("GYF", "TNF"), each=3)
scenario = c("HB", "DY", "ND13s3", "HB", "DY", "ND16s1")

par(mfrow = c(2,3), mar=c(4,4,1,1))

df = tibble(site = site, scenario = scenario) %>% 
  mutate(file = paste0("PF-CWM_",site,"_",scenario,"_Plant-FATE_oldformat_24hr_par_max.csv")) %>% 
  group_by(site, scenario) %>% 
  reframe(read.csv(file))

df %>% group_by(site, scenario) %>% 
  summarize(mean_par = mean(par),
            mean_parmax = mean(par_max),
            mean_precip = mean(precip.mm.))

# output of above shows that mean par is higher in ND scenario, 
# which explains the higher GPP
# site  scenario mean_par mean_parmax
# <chr> <chr>       <dbl>       <dbl>
# 1 GYF   DY           418.       1385.
# 2 GYF   HB           401.       1338.
# 3 GYF   ND13s3       424.       1403.
# 4 TNF   DY           431.       1403.
# 5 TNF   HB           420.       1378.
# 6 TNF   ND16s1       432.       1410.

df %>% ggplot() + 
  geom_hex(aes(x=precip.mm., y=par)) + 
  scale_fill_viridis_c(trans="log") +
  facet_wrap(facets = c("site","scenario"))
