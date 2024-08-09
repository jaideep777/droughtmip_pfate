library(tidyverse)
rm(list=ls())


input_dir  = "~/codes/Drought_MIP/input_data/"
output_dir = "~/codes/Drought_MIP/pfate_output/"

read_pfate_outputs = function(input_dir, output_dir, expt_dir){
  wd_back = getwd()
  setwd(paste0(output_dir,"/",expt_dir))

  l = list(
    # seeds1 = read.delim("seeds.csv", header=F, col.names = paste0("V", 1:(n_species+2)))
    # Zp = read.csv("z_star.csv", header=F, col.names = paste0("V", 1:50)),
    # BA1 = read.csv("basal_area.csv", header=F, col.names = paste0("V", 1:(n_species+2)))
    # co = read.csv("canopy_openness.csv", header=F, col.names = paste0("V", 1:50)),
    # lai_v = read.csv("lai_profile.csv", header=F, col.names = paste0("V", 1:27)),
    traits = read.csv("traits.csv"),
    dat_d = readr::read_csv("D_PFATE.csv"),
    # dat$YEAR = decimal_date(as_date(dat$YEAR, format = "%Y-%m-%d %H:%M:%S GMT (doy = %j)"))
    dat2 = read.csv("Y_PFATE.csv"),
    dat3 = read.csv("Y_mean_PFATE.csv"),
    # dist = readr::read_csv("size_distributions.csv", col_names = F),
    x = exp(seq(log(0.01), log(10), length.out=100))
  )

  # l$dist = l$dist[,-ncol(l$dist)]
  # names(l$dist)[1:2] = c("YEAR", "SPP")
  # names(l$Zp)[1] = c("YEAR")
  # names(l$co)[1] = c("YEAR")
  # names(l$lai_v)[1] = c("YEAR")

  l$dat = l$dat_d %>%
    mutate(YEAR = as.integer(YEAR)) %>%
    group_by(YEAR) %>%
    summarize_all(mean)

  l$traits = l$traits %>% filter(!grepl("probe", .$SPP))

  n_species = l$dat2 %>% filter(!grepl("probe", .$PID)) %>% pull(PID) %>% unique() %>% length()
  n_year = length(unique(l$dat2$YEAR))

  setwd(wd_back)

  l
}

cat_outputs = function(list1, list2){
  keys <- unique(c(names(list1), names(list2)))
  l <- lapply(setNames(keys, keys), function(x) {rbind(list1[[x]], list2[[x]])})
  l$x = list1$x
  l
}

subsample = function(l, interval = 10){
  keys <- unique(c(names(l)))
  keys <- keys[keys != "x"]
  keys <- keys[keys != "dat_d"]

  years = unique(l$dat2$YEAR)

  lsub = lapply(setNames(keys, keys), function(x) {l[[x]] = l[[x]] %>% filter(as.integer(YEAR) %% interval == 0)})
  lsub$x = l$x
  lsub$dat_d = l$dat_d
  lsub
}

slice_time = function(l, ymin, ymax){
  keys <- unique(c(names(l)))
  keys <- keys[keys != "x"]

  years = unique(l$dat2$YEAR)

  lsub = lapply(setNames(keys, keys), function(x) {l[[x]] = l[[x]] %>% filter(YEAR < ymax & YEAR > ymin)})
  lsub$x = l$x
  lsub
}


dirs = list(
  sitename = c("GYF", "TNF"),
  scenario = c("SPIN", "HB", "DY", "ND"),
  div = c("High", "Low")) |>
  cross_df() |>
  left_join(data.frame(sitename = c("GYF", "TNF"),
                       nd_suffix = c("13s3", "16s1")
  )) |>
  mutate(scenario = ifelse(scenario=="ND",
                           no=scenario,
                           yes=paste0(scenario, nd_suffix))) |>
  select(-nd_suffix) |>
  mutate(dirs = paste0(
    sitename, "_",
    case_match(scenario,
               "SPIN"~"HB",
               .default = scenario
    ), "_",
    case_when(
      scenario == "SPIN" & div == "High" ~ "20ky",
      scenario == "SPIN" & div == "Low" ~ "1ky",
      substr(scenario, 1,2) %in% c("ND","DY") & div == "High" ~ "10ky",
      substr(scenario, 1,2) %in% c("ND","DY") & div == "Low" ~ "500y",
      scenario == "SPIN" & div == "Low" ~ "1ky",
      .default = "100y"), "_",
    ifelse(div=="High", yes="evol", no="ld"), "_",
    case_when(
      scenario == "SPIN" & div == "High" ~ "7.0",
      .default = "7.1"
    ),
    ifelse(scenario=="SPIN", no = "_cont", yes = "")
  )) |>
  distinct() |>
  filter(scenario != "SPIN") |>
  mutate(filename_prefix = paste0(
    "PF-SPL_",
    sitename, "_",
    div, "_",
    scenario
  ))


data = dirs %>%
  mutate(raw_eco = purrr::map(.x = dirs, .f = function(x){
                  read_pfate_outputs(input_dir, output_dir, x)
    })) %>%
  mutate(raw_ind = purrr::map(.x = dirs, .f = function(x){
                  cat("reading ", x, "\n")
                  readr::read_csv(paste0(output_dir, "/", x, "/cohort_props.csv"))
    }))


data_proc = data %>%
  mutate(processed_eco = purrr::map(.x = raw_eco, .f = function(x){
    x$dat_d %>%
      left_join(x$dat3 %>%
                  mutate(IYEAR = as.integer(YEAR)) %>%
                  select(IYEAR, LAI)
                ) %>%
      mutate(ET = TRANS + AESOIL,
             SW1 = SWCV,
             LFLIT = -9999,
             SW2 = -9999,
             SW3 = -9999,
             SW4 = -9999) %>%
      mutate(GPP = GPP*1e3, # convert kgC m-2 day-1 ---> gC m-2 day-1
             NPP = NPP*1e3, # convert kgC m-2 day-1 ---> gC m-2 day-1
             ET = ET # already in mm day-1
             ) %>%
      select(IYEAR, MON, DAY, GPP, NPP, ET, LAI, LFLIT, SW1, SW2, SW3, SW4) %>%
      rename(YEAR = IYEAR)
  })) %>%
  mutate(processed_ind = purrr::map2(.x = raw_ind, .y = raw_eco, .f = function(x,y){
    x %>%
      mutate(YEAR = as.integer(YEAR)) %>%
      mutate(NLIVE= density*1e4,  # convert m-2 ---> ha-1
             DBH = dbh*100   # convert m --> cm
      ) %>%
      rename(SP = speciesID,
             ID = cohortID,
             HT = height,
             TB = total_biomass,
             AGB = agb) %>%
      select(YEAR, SP, ID, NLIVE, DBH, HT, TB, AGB) %>%
      left_join(y$traits %>%
                  mutate(YEAR = as.integer(YEAR)) %>%
                  mutate(SLA = 1/LMA) %>%
                  select(YEAR, SPP, WD, SLA, HMAT, P50X) %>%
                  rename(SP = SPP)
                )
  })) %>%
  select(-starts_with("raw"), -dirs)

setwd("~/codes/Drought_MIP/data_submitted_v3")

## Filter and write short-term outputs
data_proc %>%
  mutate(processed_eco = purrr::map(.x = processed_eco, .f = function(x){
    x %>% filter(YEAR <= 2500)
    }
    )) %>%
  mutate(processed_ind = purrr::map(.x = processed_ind, .f = function(x){
    x %>% filter(YEAR <= 2500)
    }
  )) %>%
  select(filename_prefix, processed_eco, processed_ind) %>%
  pivot_longer(processed_eco:processed_ind, names_sep = "_", names_to = c("blurb", "level")) %>%
  mutate(filename = paste0(filename_prefix, "_", level, ".csv")) %>%
  mutate(file = here::here("data_submitted_v3", "short_term", filename)) %>%
  group_by(file) %>%
  do(a = write.csv(.$value, file = .$file, row.names = F))


## Write long-term outputs
data_proc %>%
  select(filename_prefix, processed_eco, processed_ind) %>%
  pivot_longer(processed_eco:processed_ind, names_sep = "_", names_to = c("blurb", "level")) %>%
  mutate(filename = paste0(filename_prefix, "_", level, ".csv")) %>%
  mutate(file = here::here("data_submitted_v3", "long_term", filename)) %>%
  group_by(file) %>%
  do(a = write_csv(as.data.frame(.$value), file = .$file))


#
# data_proc %>%
#   mutate(npp_avg = purrr::map_dbl(.x = processed_eco, ~mean(.x$NPP*1e-6*1e4*365)),
#          et_avg = purrr::map_dbl(.x = processed_eco, ~mean(.x$ET)),
#          wd_avg = purrr::map_dbl(.x = processed_ind, ~mean(head(.x$WD))),
#          p50_avg = purrr::map_dbl(.x = processed_ind, ~mean(head(.x$P50))))
#
#
# add_band = function(start = year(as.Date("2000-1-1"))-1e7, end=year(as.Date("2000-1-1"))){
#   polygon(x=c(start,end,end,start), y=c(-1e20,-1e20,1e20,1e20), border = NA, col=scales::alpha("grey",0.2))
# }
#
# add_hband = function(ylim, col="grey30", alpha=0.7, xlim=c(-1e20,200020)){
#   polygon(y=c(ylim[1],ylim[2],ylim[2],ylim[1]), x=c(xlim[1],xlim[1],xlim[2],xlim[2]), border = NA, col=scales::alpha(col, alpha))
# }
#
# # Values for calibration
# obs_values = data.frame(
#   site = c("GYF", "TNF"),
#   agb = c(201, 197), # MgC ha-1
#   bgb = c(42.2, 41.4), # MgC ha-1
#   tb = c(243.2, 238.4), #
#   npp = c(11.6, 13.9), # MgC ha-1 yr-1
#   mort = c(1.13, 1.83), # % yr-1
#   lai = c(6.95, 6.0),
#   litterfall = c(2.1, 4.2) # MgC ha-1 yr-1
# )
#
# plot_all = function(hd, ld, sitename){
#   hd = hd |> group_by(YEAR) |> summarize_all(mean)
#   ld = ld |> group_by(YEAR) |> summarize_all(mean)
#
#   div_cols = c("purple", "goldenrod")
#   matplot(y = cbind(hd$NPP/1e3,
#                     ld$NPP/1e3)*1e3*1e-6*1e4*365.2425, # reconvert
#           x = hd$YEAR,
#           type="l", lty=1, col=div_cols,
#           xlab = "Year", ylab="NPP\n(MgC ha-1 yr-1)",
#           ylim=c(10,15))
#   add_hband((obs_values %>% filter(site==sitename) %>% pull(npp))*c(0.95,1.05), col = "aquamarine1", alpha = 0.3)#, col=scales::alpha("red3", 0.2))
#   add_band()
#   mtext(side=2, line=6, text=sitename)
#
#   matplot(y = cbind(hd$ET,
#                     ld$ET),
#           x = hd$YEAR,
#           type="l", lty=1, col=div_cols,
#           xlab = "Year", ylab="ET\n(mm day-1)")
#   add_band()
#
# }
#
# data_wide = data_proc |>
#   select(sitename, scenario, div, processed_eco) |>
#   pivot_wider(names_from = scenario, values_from = processed_eco) |>
#   pivot_longer(DY:ND16s1, names_to="scenario", values_to="data", values_drop_na = T) |>
#   mutate(hb_sc = purrr::map2(.x=HB, .y=data, .f = rbind)) |>
#   select(sitename, scenario, div, hb_sc) |>
#   pivot_wider(names_from = div, values_from = hb_sc)
#
# par(mfrow=c(4,2), mar=c(3,5,2,1), oma=c(1,5,1,1), mgp=c(3,1,0))
# data_wide %>% filter(scenario == "DY") %>% with(plot_all(High[[1]], Low[[1]], sitename[[1]]))
# mtext(side=3, line=1, outer=F, text="Scenario: HB + DY")
# data_wide %>% filter(scenario == "DY") %>% with(plot_all(High[[2]], Low[[2]], sitename[[2]]))
#
# data_wide %>% filter(scenario != "DY") %>% with(plot_all(High[[1]], Low[[1]], sitename[[1]]))
# mtext(side=3, line=1, outer=F, text="Scenario: HB + ND")
# data_wide %>% filter(scenario != "DY") %>% with(plot_all(High[[2]], Low[[2]], sitename[[2]]))
#


