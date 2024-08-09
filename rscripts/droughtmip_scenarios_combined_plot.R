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


extract_plot_vars = function(l){
  bind_rows(
    l$dat %>%
      mutate(ET = TRANS + AESOIL) %>%
      mutate(MORT = MORT/GPP*100, # convert to % of GPP
             NPP = NPP*1e3*1e-6*1e4*365.2425, # convert kg m-2 day-1 ---> MgC ha-1 yr-1
             GPP = GPP*1e3*1e-6*1e4*365.2425, # convert kg m-2 day-1 ---> MgC ha-1 yr-1
      ) %>%
      select(YEAR, NPP, ET, TRANS, AESOIL, SWCV, MORT) %>%
      pivot_longer(-YEAR),

    l$dat3 %>%
      mutate(AGB = CL+CW,
             BGB = CFR+CCR) %>%
      mutate(AGB = AGB *1e3*1e-6*1e4, # convert kgC m-2 ---> to MgC ha-1
             BGB = BGB *1e3*1e-6*1e4, # convert kgC m-2 ---> to MgC ha-1
             BA = BA*1e4  # convert m2 m-2 ---> to m2 ha-1
             ) %>%
      select(YEAR, BA, LAI, AGB, BGB) %>%
      pivot_longer(-YEAR),

    l$traits %>%
      filter(!grepl("probe", SPP)) %>%
      select(YEAR, SPP, WD, HMAT, P50X) %>%
      pivot_longer(-c(YEAR, SPP))
  )
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
  distinct()


data_lt = dirs %>%
  mutate(data = purrr::map(
    .x = dirs,
    .f = ~read_pfate_outputs(input_dir, output_dir, .x) %>%
            # subsample(interval = 24) %>%
            extract_plot_vars()
    )
  )


data_lt = data_lt %>% unnest(data)

# traits for LD runs

data_lt %>%
  filter(scenario == "HB") %>%
  filter(name %in% c("WD", "HMAT", "P50X")) %>%
  select(sitename, YEAR, name, value, SPP) %>%
  group_by(sitename, name, SPP) %>%
  summarize(value = mean(value)) %>%
  ungroup()


# Values for calibration
obs_values = data.frame(
  sitename = c("GYF", "TNF"),
  AGB = c(201, 197), # MgC ha-1
  # BGB = c(42.2, 41.4), # MgC ha-1
  # TB = c(243.2, 238.4), #
  NPP = c(11.6, 13.9), # MgC ha-1 yr-1
  # MORT = c(1.13, 1.83), # % yr-1
  LAI = c(6.95, 6.0)
  # LITTER = c(2.1, 4.2) # MgC ha-1 yr-1
) %>% pivot_longer(-sitename) %>%
  mutate(
    max = value*1.025,
    min = value*0.975
  )


p_st = data_lt %>%
  filter(div == "High") %>%
  filter(name %in% c("NPP", "AGB", "LAI", "MORT", "ET", "WD", "HMAT", "P50X", "SWCV")) %>%
  filter(YEAR > 1500 & YEAR < 2500) %>%
  ggplot(aes(x=YEAR, y=value))+
  geom_line(aes(col=scenario), alpha=0.7)+
  geom_rect(data = obs_values, aes(ymin=min, ymax=max, xmin=-Inf, xmax=2000), inherit.aes=F,  fill="turquoise", alpha=0.5) +
  ggh4x::facet_grid2(sitename~name, scales="free_y", independent="y")+
  annotate("rect", ymin=-Inf, ymax=Inf, xmin=1900, xmax=2000, fill="grey", col=NA, alpha=0.3)+
  annotate("rect", ymin=-Inf, ymax=Inf, xmin=-Inf, xmax=1900, fill="grey90", col=NA, alpha=0.3)+
  scale_x_continuous(n.breaks=2)+
  theme_bw()+
  scale_color_manual(values=c(
    SPIN="grey40",
    HB="grey30",
    DY="goldenrod",
    ND13s3="orange4",
    ND16s1="orange4"
  ))


p_lt = data_lt %>%
  filter(div == "High") %>%
  filter(name %in% c("NPP", "AGB", "LAI", "MORT", "ET", "WD", "HMAT", "P50X", "SWCV")) %>%
  filter(as.integer(YEAR) %in% seq(min(YEAR), max(YEAR), by=23)) %>%
  ggplot(aes(x=YEAR, y=value))+
  geom_line(aes(col=scenario), alpha=0.7)+
  geom_rect(data = obs_values, aes(ymin=min, ymax=max, xmin=-Inf, xmax=2000), inherit.aes=F,  fill="turquoise", alpha=0.5) +
  ggh4x::facet_grid2(sitename~name, scales="free_y", independent="y")+
  annotate("rect", ymin=-Inf, ymax=Inf, xmin=1900, xmax=2000, fill="grey", col=NA, alpha=0.3)+
  annotate("rect", ymin=-Inf, ymax=Inf, xmin=-Inf, xmax=1900, fill="grey90", col=NA, alpha=0.3)+
  scale_x_continuous(n.breaks=2)+
  theme_bw()+
  scale_color_manual(values=c(
    SPIN="grey40",
    HB="grey30",
    DY="goldenrod",
    ND13s3="orange4",
    ND16s1="orange4"
  ))

library(patchwork)
png(here::here("figures/HD_allscenarios_shortterm_vs_longterm.png"), width=1600*3, height=800*3, res=300)
print(
p_st/p_lt + plot_layout(guides="collect")
)
dev.off()

p_nd = data_lt %>%
  filter(grepl("ND", scenario) | scenario == "HB") %>%
  filter(name %in% c("NPP", "AGB", "LAI", "MORT", "ET", "WD", "HMAT", "P50X", "SWCV")) %>%
  filter(YEAR > 1901 & YEAR < 2500) %>%
  ggplot(aes(x=YEAR, y=value))+
  geom_line(aes(col=div), alpha=0.7)+
  geom_rect(data = obs_values, aes(ymin=min, ymax=max, xmin=-Inf, xmax=2000), inherit.aes=F,  fill="turquoise", alpha=0.5) +
  ggh4x::facet_grid2(sitename~name, scales="free_y", independent="y")+
  annotate("rect", ymin=-Inf, ymax=Inf, xmin=1900, xmax=2000, fill="grey", col=NA, alpha=0.3)+
  annotate("rect", ymin=-Inf, ymax=Inf, xmin=-Inf, xmax=1900, fill="grey90", col=NA, alpha=0.3)+
  scale_x_continuous(n.breaks=2)+
  theme_bw()+
  scale_color_manual(values=c(
    High="purple",
    Low="green3"
  )) +
  ggtitle("ND")

p_dy = data_lt %>%
  filter(grepl("DY", scenario) | scenario == "HB") %>%
  filter(name %in% c("NPP", "AGB", "LAI", "MORT", "ET", "WD", "HMAT", "P50X", "SWCV")) %>%
  filter(YEAR > 1901 & YEAR < 2500) %>%
  ggplot(aes(x=YEAR, y=value))+
  geom_line(aes(col=div), alpha=0.7)+
  geom_rect(data = obs_values, aes(ymin=min, ymax=max, xmin=-Inf, xmax=2000), inherit.aes=F,  fill="turquoise", alpha=0.5) +
  ggh4x::facet_grid2(sitename~name, scales="free_y", independent="y")+
  annotate("rect", ymin=-Inf, ymax=Inf, xmin=1900, xmax=2000, fill="grey", col=NA, alpha=0.3)+
  annotate("rect", ymin=-Inf, ymax=Inf, xmin=-Inf, xmax=1900, fill="grey90", col=NA, alpha=0.3)+
  scale_x_continuous(n.breaks=2)+
  theme_bw()+
  scale_color_manual(values=c(
    High="purple",
    Low="green3"
  )) +
  ggtitle("DY")

png(here::here("figures/SHORT_allscenarios_hd_vs_ld.png"), width=1600*3, height=800*3, res=300)
print(
p_nd/p_dy
)
dev.off()
