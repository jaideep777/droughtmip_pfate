library(tidyverse)
rm(list=ls())

### To convert old format txt outputs to new csv format:
# sed 's/\t/,/g' AmzFACE_D_PFATE_ELE_HD.txt > D_PFATE.csv
# sed 's/\t/,/g' AmzFACE_Y_PFATE_ELE_HD.txt > Y_PFATE.csv
# sed 's/\t/,/g' AmzFACE_Y_mean_PFATE_ELE_HD.txt > Y_mean_PFATE.csv
# sed 's/\t/,/g' canopy_openness.txt > canopy_openness.csv
# sed 's/\t/,/g' lai_profile.txt > lai_profile.csv
# sed 's/\t/,/g' size_distributions.txt > size_distributions.csv
# sed 's/\t/,/g' traits_ELE_HD.txt > traits.csv
# sed 's/\t/,/g' z_star.txt > z_star.csv

sitename = "TNF"
# sitename = "GYF"

input_dir = "~/codes/Drought_MIP/input_data/"
output_dir = "~/codes/Drought_MIP/pfate_output/"

# expt_dir = "calib_GYF_HB_2"
# expt_dir = "calib_TNF_HB_2"

# HB - High
# expt_dir = "TNF_HB_20ky_evol_5.1"
# expt_dir = "GYF_HB_20ky_evol_5.1"
expt_dir = paste0(sitename,"_HB_20ky_evol_7.0")
# expt_dir = "GYF_HB_20ky_evol_7.0"

expt_dir_cont = paste0(sitename,"_HB_100y_evol_7.1_cont")

expt_dir_cont = ifelse(sitename == "GYF",
                       yes = paste0("GYF_ND13s3_10ky_evol_7.1_cont"),
                       no  = paste0("TNF_ND16s1_10ky_evol_7.1_cont")
                       )

expt_dir_cont = ifelse(sitename == "GYF",
                       yes = paste0("GYF_HB_10ky_evol_7.1_cont"),
                       no  = paste0("TNF_HB_10ky_evol_7.1_cont_control")
)


# LD
# expt_dir = "GYF_HB_1ky_ld_6.0"
# expt_dir = "TNF_HB_100y_ld_6.1_cont"


# expt_dir = "cont_test_GYF_HB_100y_evol_4.7"

# expt_dir_cont = "TNF_HB_100y_evol_5.1"
# expt_dir_cont = "GYF_HB_100y_evol_5.1"


# expt_dir_cont = "TNF_ND16s1_500y_ld_5.1_cont"
# expt_dir_cont = "GYF_ND13s3_500y_ld_5.1_cont"

# expt_dir_cont = "corrected_soil_TNF_HB_500y_ld_6.1_cont"
# expt_dir_cont = "GYF_HB_100y_ld_6.0_cont"

# expt_dir_cont = "cont_test_GYF_HB_100y_evol_4.7_cont"


plot_to_file = F
plot_trait_space = F


read_pfate_outputs = function(input_dir, output_dir, expt_dir){
  wd_back = getwd()
  setwd(paste0(output_dir,"/",expt_dir))

  l = list(
    # seeds1 = read.delim("seeds.csv", header=F, col.names = paste0("V", 1:(n_species+2)))
    Zp = read.csv("z_star.csv", header=F, col.names = paste0("V", 1:50)),
    # BA1 = read.csv("basal_area.csv", header=F, col.names = paste0("V", 1:(n_species+2)))
    co = read.csv("canopy_openness.csv", header=F, col.names = paste0("V", 1:50)),
    lai_v = read.csv("lai_profile.csv", header=F, col.names = paste0("V", 1:27)),
    traits = read.csv("traits.csv"),
    dat_d = readr::read_csv("D_PFATE.csv"),
    # dat$YEAR = decimal_date(as_date(dat$YEAR, format = "%Y-%m-%d %H:%M:%S GMT (doy = %j)"))
    dat2 = read.csv("Y_PFATE.csv"),
    dat3 = read.csv("Y_mean_PFATE.csv"),
    dist = readr::read_csv("size_distributions.csv", col_names = F),
    x = exp(seq(log(0.01), log(10), length.out=100))
  )

  l$dist = l$dist[,-ncol(l$dist)]
  names(l$dist)[1:2] = c("YEAR", "SPP")
  names(l$Zp)[1] = c("YEAR")
  names(l$co)[1] = c("YEAR")
  names(l$lai_v)[1] = c("YEAR")

  l$dat = l$dat_d %>%
    mutate(YEAR = as.integer(YEAR)) %>%
    group_by(YEAR) %>%
    summarize_all(mean)

  n_species = l$dat2 %>% filter(!grepl("probe", .$PID)) %>% pull(PID) %>% unique() %>% length()
  n_year = length(unique(l$dat2$YEAR))

  setwd(wd_back)

  l
}

cat_outputs = function(list1, list2){
  keys <- unique(c(names(list1), names(list2)))
  keys <- keys[keys != "x"]
  l <- lapply(setNames(keys, keys), function(x) {
      bind_rows(list1[[x]], list2[[x]])
    })
  l$x = list1$x
  l
}

# aggregate_annual = function(l){
#   keys <- unique(c(names(l)))
#   keys <- keys[keys != "x"]
#
#   years = unique(l$dat2$YEAR)
#
#   lagg = lapply(setNames(keys, keys), function(x) {l[[x]] = l[[x]] %>% group_by(as.integer(YEAR)) %>% summarize_all(mean)})
#   lagg$x = l$x
#   lagg
# }

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

l1 = read_pfate_outputs(input_dir, output_dir, expt_dir)
l2 = read_pfate_outputs(input_dir, output_dir, expt_dir_cont)

as.Date("2000-1-1")-6934855.000000

# l1_slice = slice_time(l1, -20000, 17650)
# l2_slice = slice_time(l2, -20000, 26000)

# l1_slice = slice_time(l1, 1500, 2500)
# l2_slice = slice_time(l2, 1500, 2500)

# l1_slice = slice_time(l1, -2000, 10000)
# l2_slice = slice_time(l2, -2000, 10000)

# l1_sub = subsample(l1, interval = 10)
# # l2_sub = subsample(l2)
#
# l1_sub_slice = slice_time(l1_sub, 2500, 3500)

# l_comb = cat_outputs(l1_slice, l2_slice)

# l_full = slice_time(l1, 2500, 4200)
# l_full = slice_time(l_comb, 2500, 4200)

l1_sub = subsample(l1_slice, interval = 1)
l2_sub = subsample(l2_slice, interval = 1)

l = cat_outputs(l1_sub, l2_sub)
# l = l1_sub

setwd(paste0(output_dir,"/",expt_dir))

add_band = function(start = -20000, end=year(as.Date("1901-1-1"))){
  polygon(x=c(start,end,end,start), y=c(-1e20,-1e20,1e20,1e20), border = NA, col=scales::alpha("yellow2",0.2))
}

add_hband = function(ylim, col="grey30", alpha=0.7, xlim=c(-1e20,200020)){
  polygon(y=c(ylim[1],ylim[2],ylim[2],ylim[1]), x=c(xlim[1],xlim[1],xlim[2],xlim[2]), border = NA, col=scales::alpha(col, alpha))
}

# Values for calibration
obs_values = data.frame(
  site = c("GYF", "TNF"),
  agb = c(201, 197), # MgC ha-1
  bgb = c(42.2, 41.4), # MgC ha-1
  tb = c(243.2, 238.4), #
  npp = c(11.6, 13.9), # MgC ha-1 yr-1
  mort = c(1.13, 1.83), # % yr-1
  lai = c(6.95, 6.0),
  litterfall = c(2.1, 4.2) # MgC ha-1 yr-1
)

traits_obs = read.csv(file = paste0("../../Amz_trait_orig.csv"))
# traits_used = read.csv(file = paste0(input_dir, "/Traits_random_HD2.csv"))

# To get avg size distribution, sum over species and average over years
dist_amb = l$dist %>% filter(YEAR > min(YEAR)) %>% filter(YEAR>max(YEAR)-100) %>% pivot_longer(cols=-(YEAR:SPP), names_to="size_class") %>% group_by(YEAR,size_class) %>% summarize(de = sum(value, na.rm=T)) %>% pivot_wider(names_from = size_class, values_from = de) %>% colMeans(na.rm=T)

# dist_amb = dist %>% filter(YEAR == 1100) %>% pivot_longer(cols=-(YEAR:SPP), names_to="size_class") %>% group_by(YEAR,size_class) %>% summarize(de = sum(value, na.rm=T)) %>% pivot_wider(names_from = size_class, values_from = de) %>% colMeans(na.rm=T)
# dist_ele = dist %>% filter(YEAR == 1101) %>% pivot_longer(cols=-(YEAR:SPP), names_to="size_class") %>% group_by(YEAR,size_class) %>% summarize(de = sum(value, na.rm=T)) %>% pivot_wider(names_from = size_class, values_from = de) %>% colMeans(na.rm=T)

n_species = l$dat2 %>% filter(!grepl("probe", .$PID)) %>% pull(PID) %>% unique() %>% length()
n_year = length(unique(l$dat2$YEAR))
col_species = rainbow(n = n_species, start = 0, end = 0.85, alpha = min(10/n_species, 1))
filter_span = 20/n_year

plot_gauss_mix = function(x, means, sds, wts, add=T, col, add_polygon=F, ...){
  y = x*0
  if (length(sds)==1) sds = rep(sds, length(means))
  for (i in 1:length(means)){
    y = y + wts[i]*dnorm(x, mean = means[i], sd=sds[i])
  }
  if (add) points(y~x, col=col, ...)
  else plot(y~x, col=col, ...)
  if (add_polygon){
    polygon(y~x, col=scales::alpha(col, 0.5), border = NA)
  }
  data.frame(x=x, y=y)
}


if (plot_to_file) png("master_plot.png", width=2814*1.5, height = 1472*1.5, res=300)

par(mfcol=c(4,7), mar=c(4.5,6,.5,1), oma=c(1,1,2,1), cex.lab=1.1, cex.axis=1.1, mgp=c(3.2,1,0), las=1)
seeds = l$dat2 %>% filter(!grepl("probe", .$PID)) %>% select(YEAR, PID, SEEDS) %>% spread(value = "SEEDS", key = "PID")
seeds_smooth = seeds %>% pivot_longer(-YEAR) %>% drop_na() %>% group_by(name) %>% mutate(value = loess(value~YEAR, span=filter_span) %>% fitted()) %>% pivot_wider(names_from=name)
seeds_total = rowSums(seeds[,-1,drop=FALSE], na.rm=T)
matplot(seeds$YEAR, cbind(seeds[,-1], seeds_total), lty=1, col=scales::alpha(c(col_species, "black"), 0.7), type="l",
        las=1, xlab="Time (years)", ylab="Species seed\noutput", log="")
# matplot(seeds_smooth$YEAR, seeds_smooth[,-1], lty=1, type="l",
#         # col=col_species,
#         col=scales::alpha(scales::muted(col_species), alpha=min(30/n_species, 1)),
#         las=1, xlab="Time (years)", ylab="Species seed\noutput", log="",
#         add=T
#         )
mtext(line=0.5, side=3, outer = T, text=expt_dir)
add_band()
abline(v=3100, col="grey")

# matplot(seeds1$V1, seeds1[,-1], lty=1, col=rainbow(n = n_species+1, start = 0, end = 0.85), type="l",
#         las=1, xlab="Time (years)", ylab="Species Seed output", log="")
# mtext(line=0.5, side=3, text=expt_dir)

BA = l$dat2 %>% filter(!grepl("probe", .$PID)) %>% select(YEAR, PID, BA) %>% spread(value = "BA", key = "PID")
matplot(BA$YEAR, cbind(BA[,-1], rowSums(BA[,-1,drop=FALSE], na.rm=T))*1e4, lty=1, col=c(col_species, "black"), type="l",
        las=1, xlab="Time (years)", ylab="Basal area", log="")
add_band()

matplot(l$Zp$YEAR, l$Zp[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
        las=1, xlab="Time (years)", ylab="Z*")
# matplot(co$V1, co[,-1], lty=1, col=rainbow(n = 10, start = 0, end = 0.85), type="l",
#         las=1, xlab="Time (years)", ylab="Io")
# matplot(y=1:24, x=t(-lai_v[,3:26]+lai_v[,2:25]), lty=1, col=rainbow(n = n_year, start = 0, end = 0.85, alpha=0.05), type="l",
#         las=1, xlab="Leaf area density", ylab="Height")
add_band()

matplot(y=cbind(l$dat$DPSI), x=l$dat$YEAR, type="l", lty=1, col=c("cyan4"), ylab="Dpsi\n(MPa)", xlab="Time (years)")
matlines(y=cbind(fitted(loess(l$dat$DPSI~l$dat$YEAR, span=filter_span))), x=l$dat$YEAR, type="l", lty=1, col="black", lwd=c(1,0.5))
# add_hband(c(20,50)) #, col=scales::muted("green4"))
add_band()


matplot(y=1:25, x=t(l$lai_v[,2:26]), lty=1, col=rainbow(n = n_year, start = 0, end = 0.85, alpha=0.05), type="l",
        las=1, xlab="Cumulative LAI", ylab="Height")


plot(l$dat3$LAI~l$dat3$YEAR, type="l", col="red3", ylim=c(0,max(l$dat3$LAI,10.5)), xlab="Time (years)", ylab="Total LAI")
# abline(h=c(5.3, 6.2), col=scales::muted("red"))
add_hband((obs_values %>% filter(site==sitename) %>% pull(lai)*c(0.95,1.05)), col = "aquamarine1")#, col=scales::alpha("red3", 0.2))
# abline(h=c(3.5), col=scales::muted("grey100"))
add_band()

agb = cbind(l$dat3$CL+l$dat3$CW)*1e3*1e-6*1e4
bgb = cbind(l$dat3$CFR+l$dat3$CCR)*1e3*1e-6*1e4

matplot(y=cbind(l$dat$GPP, l$dat$NPP)*1e3*1e-6*1e4*365.2425, x=l$dat$YEAR, type="l", lty=1, col=c("green4", "green3", "brown"), ylab="GPP, NPP\n(MgC/ha/yr)", xlab="Time (years)")
matlines(y=cbind(fitted(loess(l$dat$GPP~l$dat$YEAR, span=filter_span)),
                 fitted(loess(l$dat$NPP~l$dat$YEAR, span=filter_span)))*1e3*1e-6*1e4*365.2425, x=l$dat$YEAR, type="l", lty=1, col="black", lwd=c(1,0.5))
# points(y=l$dat$NPP/l$dat$GPP*4, x=l$dat$YEAR, type="l", lty=1, col=c("yellow1"))
# abline(h=c(3,3.5), col="grey")
add_hband((obs_values %>% filter(site==sitename) %>% pull(npp))*c(0.95,1.05), col = "aquamarine1")#, col=scales::alpha("red3", 0.2))
add_band()

tb = l$dat3$CCR+l$dat3$CFR+l$dat3$CL+l$dat3$CW
nye = min(length(l$dat$MORT), length(tb))
mort_rate = l$dat$MORT[1:nye]/tb[1:nye]*100*365.2425
matplot(y=cbind(mort_rate), x=l$dat$YEAR[1:nye], type="l", lty=1, col=c("brown"), ylim=c(0,6), ylab="MORT\n(% yr-1)", xlab="Time (years)")
add_hband((obs_values %>% filter(site==sitename) %>% pull(mort))*c(0.9,1.1), col = "aquamarine1")#, col=scales::alpha("red3", 0.2))
# abline(h=c(1.31), col=scales::muted("green3"))
add_band()

matplot(y=cbind(l$dat$GS), x=l$dat$YEAR, type="l", lty=1, col=c("cyan3"), ylab="Stomatal conductance\n(mol/m2/s)", xlab="Time (years)")
matlines(y=cbind(fitted(loess(l$dat$GS~l$dat$YEAR, span=filter_span))), x=l$dat$YEAR, type="l", lty=1, col="black", lwd=c(1,0.5))
add_hband(c(0.16, 0.16555))#, col=scales::alpha("cyan4", 0.6))
# abline(h=c(0.16), col=scales::muted("cyan3"))
add_band()

matplot(y=agb, x=l$dat3$YEAR, type="l", lty=1, col=c("yellow4"), ylab="AGB\n(MgC/ha)", xlab = "Time (years)")
add_hband((obs_values %>% filter(site==sitename) %>% pull(agb))*c(0.95, 1.05), col = "aquamarine1")#, col=scales::alpha("red3", 0.2))
add_band()

matplot(y=bgb, x=l$dat3$YEAR, type="l", lty=1, col=c("brown"), ylab="BGB\n(MgC/ha)", xlab = "Time (years)")
add_hband((obs_values %>% filter(site==sitename) %>% pull(bgb))*c(0.95, 1.05), col = "aquamarine1")#, col=scales::alpha("red3", 0.2))
add_band()

matplot(y=cbind(l$dat$VCMAX), x=l$dat$YEAR, type="l", lty=1, col=c("green3"), ylab="Vcmax\n(umol/m2/s)", xlab="Time (years)", ylim=c(0,60))
matlines(y=cbind(fitted(loess(l$dat$VCMAX~l$dat$YEAR, span=filter_span))), x=l$dat$YEAR, type="l", lty=1, col="black", lwd=c(1,0.5))
add_hband(c(20,50)) #, col=scales::muted("green4"))
add_band()

plot_size_dist = function(){
  # matplot(y=cbind(as.numeric(colMeans(filter(dist, V1>1100 & V1<2000)[, -c(1,2)], na.rm = T)),
  #                 as.numeric(colMeans(filter(dist, V1>2100 & V1<3000)[, -c(1,2)], na.rm = T))
  #                 )*1e-2*1e4,
  #         x=x, type="l", log="y", lty=1, col=c("black", "brown"),
  #         xlim=c(0.01, 2), ylim=c(1e-4, 200), ylab="Density (stems/cm/ha)", xlab="Diameter (m)")
  matplot(y=cbind(as.numeric(dist_amb[gtools::mixedsort(names(dist_amb))][-1])
  )*1e-2*1e4, # Convert stems m-1 m-2 --> stems cm-1 ha-1
  x=l$x, type="l", log="y", lty=1, col=c("black", "yellow3"),
  xlim=c(0.01, 1.2), ylim=c(1e-4, 1000), ylab="Density\n(stems/cm/ha)", xlab="Diameter (m)", las=0)

  abline(v=1, col=scales::alpha("red", 0.2))

  xobs = c(15,25,35,45,55,65,75,85,95,105)/100
  # Data for Manaus from https://link.springer.com/article/10.1007/s00442-004-1598-z
  yobs=c(350.5221340921042,
         132.41927918860426,
         62.62503296462008,
         29.61724892214378,
         15.095996574802413,
         5.702923697662178,
         2.3219542502889836,
         1.5968055466971947,
         0.7006940913385968,
         0.5597156879584093)/10
  points(yobs~xobs, pch=20, col=scales::alpha("grey30", 0.4), cex=1.7)
}
try(plot_size_dist())

p50dist_xy = traits_obs %>% select(P50..Mpa., Total.BasalArea_2017.cm2.) %>% drop_na %>%
  # with(density(x =Leaf.LMA..g.m2., weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% plot(ylim=c(0,0.02), las=0, main="", xlab="LMA", col=NA, lwd=2)
  with(plot_gauss_mix(x=seq(-5,0, length.out=1000), means =P50..Mpa., wts=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.), sds=2*0.05, add=F, col="grey", add_polygon=T, type="l", lwd=1.5,  las=0, main="", xlab="P50", ylab="density", ylim=c(0, 2)))
try(
  l$dat2 %>% select(YEAR, PID, BA) %>%
    left_join(l$traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
    filter(!grepl("probe", PID)) %>%
    filter(YEAR == min(20000, l$dat2$YEAR[length(l$dat2$YEAR)-1])) %>%
    with(plot_gauss_mix(x=seq(-5,0, length.out=1000), means =P50X, wts=BA/sum(BA), sds=2*0.2, add=T, col="black", type="l", lwd=1.5))
)


traits_obs %>% select(meanWoodDensity..g.cm3., Total.BasalArea_2017.cm2.) %>% drop_na %>%
  #with(density(x =meanWoodDensity..g.cm3.*1000, weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% plot(ylim=c(0,0.005), las=0, main="", xlab="Wood density", col=NA, lwd=2)
  with(plot_gauss_mix(x=seq(200,1200, length.out=1000), means =meanWoodDensity..g.cm3.*1000, wts=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.), sds=800*0.1, add=F, col="grey", add_polygon=T, type="l", lwd=1.5,  las=0, main="", xlab="Wood density", ylab="density", ylim=c(0, 0.005)))
try(
  l$dat2 %>% select(YEAR, PID, BA) %>%
    left_join(l$traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
    filter(!grepl("probe", PID)) %>%
    filter(YEAR == min(20000, l$dat2$YEAR[length(l$dat2$YEAR)-5])) %>%
    # with(density(x =WD, adjust=1, weights=BA/sum(BA))) %>% points(col="black", type="l", lwd=1.5)
    with(plot_gauss_mix(x=seq(200,1200, length.out=1000), means =WD, wts=BA/sum(BA), sds=800*0.1, add=T, col="black", type="l", lwd=1.5))
)

#
traits_obs %>% select(Height_Max.m., Total.BasalArea_2017.cm2.) %>% drop_na %>%
  # with(density(x =Height_Max.m., weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% plot(ylim=c(0,0.18), las=0, main="", xlab="Max. height", col=NA, lwd=2)
  with(plot_gauss_mix(x=seq(0,50, length.out=1000), means =Height_Max.m., wts=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.), sds=25*0.1, add=F, col="grey", add_polygon=T, type="l", lwd=1.5,  las=0, main="", xlab="Max. Height", ylab="density", ylim=c(0, 0.2)))
try(
  l$dat2 %>% select(YEAR, PID, BA) %>%
    left_join(l$traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
    filter(!grepl("probe", PID)) %>%
    filter(YEAR == min(20000, l$dat2$YEAR[length(l$dat2$YEAR)-5])) %>%
    # with(density(x =HMAT, weights=BA/sum(BA))) %>% points(col="black", type="l", lwd=1.5)
    with(plot_gauss_mix(x=seq(0,50, length.out=1000), means =HMAT, wts=BA/sum(BA), sds=25*0.1, add=T, col="black", type="l", lwd=1.5))
)


# traits_obs %>% select(P50..Mpa., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =P50..Mpa., weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% plot(ylim=c(0,0.5), las=0, main="", xlab="Max. height", col=NA, lwd=2)
# traits_obs %>% select(P50..Mpa., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =P50..Mpa., weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% polygon(col=scales::alpha("grey30", 0.2), border=scales::alpha("grey30",0.4))
# l$dat2 %>% select(YEAR, PID, BA) %>%
#   filter(YEAR == min(2000, max(l$dat2$YEAR)-1)) %>%
#   left_join(traits_obs %>% select(-BA), by = c("PID"="Species")) %>%
#   drop_na %>%
#   with(density(x =P50..Mpa., weights=BA/sum(BA))) %>% points(col="black", type="l", lwd=1.5)
# try(
#   l$dat2 %>% select(YEAR, PID, BA) %>%
#     filter(YEAR == 3000) %>%
#     left_join(traits_obs %>% select(-BA), by = c("PID"="Species")) %>%
#     drop_na %>%
#     with(density(x =P50..Mpa., weights=BA/sum(BA))) %>% points(col="yellow3", type="l", lwd=1.5)
# )


#
# traits_obs %>% select(meanWoodDensity..g.cm3., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =meanWoodDensity..g.cm3.*1000, weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% plot(ylim=c(0,0.005), las=0, main="", xlab="Wood density", col=NA, lwd=2)
# traits_obs %>% select(meanWoodDensity..g.cm3., Total.BasalArea_2017.cm2.) %>% drop_na %>% with(density(x =meanWoodDensity..g.cm3.*1000, weights=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.))) %>% polygon(col=scales::alpha("grey30", 0.2), border=scales::alpha("grey30",0.4))
# l$dat2 %>% select(YEAR, PID, BA) %>%
#   filter(YEAR == min(2000, max(l$dat2$YEAR)-1)) %>%
#   left_join(traits_obs %>% select(-BA), by = c("PID"="Species")) %>%
#   left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
#   # drop_na %>%
#   with(density(x =meanWoodDensity..g.cm3.*1000, weights=BA/sum(BA))) %>% points(col="black", type="l", lwd=1.5)
# try(
#   l$dat2 %>% select(YEAR, PID, BA) %>%
#     filter(YEAR == min(2000, max(l$dat2$YEAR)-1)) %>%
#     left_join(traits_obs %>% select(-BA), by = c("PID"="Species")) %>%
#     left_join(traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
#     #drop_na %>%
#     with(density(x =WD, weights=BA/sum(BA))) %>% points(col="yellow3", type="l", lwd=1.5)
# )
#
#

cwm_wd = traits_obs %>% select(P50..Mpa., Total.BasalArea_2017.cm2.) %>% drop_na %>% summarise(cwm=sum(P50..Mpa.*Total.BasalArea_2017.cm2.)/sum(Total.BasalArea_2017.cm2.))
cwm_wd_pred = l$dat2 %>% select(YEAR, PID, BA) %>%
  left_join(l$traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
  drop_na %>% group_by(YEAR) %>%
  summarize(cwm_wd = sum(P50X*BA)/sum(BA))

hmat = l$traits %>% filter(!grepl("probe", .$SPP)) %>% select(YEAR, SPP, HMAT) %>% pivot_wider(names_from = "SPP", values_from = "HMAT")
wd = l$traits %>% filter(!grepl("probe", .$SPP)) %>% select(YEAR, SPP, WD) %>% pivot_wider(names_from = "SPP", values_from = "WD")
p50x = l$traits %>% filter(!grepl("probe", .$SPP)) %>% select(YEAR, SPP, P50X) %>% pivot_wider(names_from = "SPP", values_from = "P50X")
smx = l$traits %>% filter(!grepl("probe", .$SPP)) %>% select(YEAR, SPP, SMX) %>% pivot_wider(names_from = "SPP", values_from = "SMX")

# matplot(y=hmat[,-1], x=wd[,-1], lty=1, type="o", pch=20, cex=0.5, col=col_species, xlab="Wood density", ylab="Max height")
# matplot(y=p50x[,-1], x=smx[,-1], lty=1, type="o", pch=20, cex=0.5, col=col_species, xlab="SMx", ylab="P50x")

matplot(y=hmat[,-1], x=hmat[,1], lty=1, type="l", pch=20, cex=0.5, col=col_species, ylab="Max height", xlab="YEAR")
add_band()
matplot(y=wd[,-1], x=hmat[,1], lty=1, type="l", pch=20, cex=0.5, col=col_species, ylab="Wood density", xlab="YEAR")
add_band()
matplot(y=p50x[,-1], x=hmat[,1], lty=1, type="l", pch=20, cex=0.5, col=col_species, ylab="P50x", xlab="YEAR")
add_band()
matplot(y=smx[,-1], x=hmat[,1], lty=1, type="l", pch=20, cex=0.5, col=col_species, ylab="SMx", xlab="YEAR")
add_band()
# matplot(x=wd[,1], y=wd[,-1], col=col_species, lty=1, type="l", ylab="Wood density", xlab="Year")
# matplot(x=hmat[,1], y=hmat[,-1], col=col_species, lty=1, type="l", ylab="Max. height", xlab="Year")

zz = l$traits %>%
  select(YEAR, SPP, r0_avg) %>%
  filter(!grepl(x = SPP, "probe")) %>%
  pivot_wider(names_from = SPP, values_from = r0_avg) %>%
  as.matrix()
zz_smooth = zz %>% as.data.frame() %>% pivot_longer(-YEAR) %>% drop_na() %>% group_by(name) %>% mutate(value = loess(value~YEAR, span=60/length(value)) %>% fitted()) %>% pivot_wider(names_from=name) %>% as.matrix

# matplot(y=tanh(zz[,-1]*20)/20, x=zz[,1], type="l", lty=1, ylab="r0", col=col_species)
matplot(y=tanh(zz_smooth[,-1]*20)/20, x=zz[,1], type="l", lty=1, ylab="r0", col=(col_species))
abline(h=0, col="black", lwd=0.2)
# abline(v=1000, col="grey")

# cwm_p50 = traits_obs %>% select(P50..Mpa., Total.BasalArea_2017.cm2.) %>% drop_na %>% summarise(cwm=sum(P50..Mpa.*Total.BasalArea_2017.cm2.)/sum(Total.BasalArea_2017.cm2.))
# cwm_p50_pred = dat2 %>% select(YEAR, PID, BA) %>%
#   left_join(traits_obs, by = c("PID"="Species")) %>%
#   drop_na %>% group_by(YEAR) %>%
#   summarize(cwm_p50 = sum(P50..Mpa.*BA)/sum(BA))
#
# cwm_p50_pred %>% with(plot(cwm_p50~YEAR, type="l")) #, ylim=c(200,900)))
# add_hband(c(cwm_wd,cwm_wd+5))


# p50 = traits %>% filter(!grepl("probe", .$SPP)) %>% select(YEAR, SPP, P50X) %>% pivot_wider(names_from = "SPP", values_from = "P50X")
# matplot(x=p50[,1], y=p50[,-1], col=col_species, lty=1, type="l", ylab="P50", xlab="Year")
#
# p50 %>% tail() %>% print()

plot(l$dat$SWP~l$dat$YEAR, type="l", col="blue", ylab="Mean SWP\n(annual, MPa)", xlab="YEAR")
add_band()

l$dat_d %>% filter(YEAR > max(YEAR)-100) %>% with(plot(SWP~YEAR, type="l", col="blue", ylab="SWP\n(MPa)", xlab="YEAR"))
l$dat_d %>% filter(YEAR > max(YEAR)-100) %>% with(points(SWPA~YEAR, type="l", col="cyan3", ylab="SWP\n(MPa)", xlab="YEAR"))
add_band()
l1_sub$dat_d %>% filter(YEAR > max(YEAR)-50) %>% with(plot(GPP~SWP, pch=".", col="green3", ylab="GPP", xlab="SWP"))
add_band()
# l2_sub$dat_d %>% filter(YEAR > max(YEAR)-50) %>% with(points(GPP~SWP, pch=".", col="green4", ylab="GPP", xlab="SWP"))


matplot(y=cbind(l$dat$TRANS, l$dat$AESOIL), x=l$dat$YEAR, type="l", lty=1, col=c("cyan3", "yellow3"), ylab="E,T\n(mm day-1)", xlab="Time (years)")
matlines(y=cbind(fitted(loess(l$dat$TRANS~l$dat$YEAR, span=filter_span)),
                 fitted(loess(l$dat$AESOIL~l$dat$YEAR, span=filter_span))), x=l$dat$YEAR, type="l", lty=1, col="black", lwd=c(1,0.5))
# abline(h=c(0.16), col=scales::muted("cyan3"))
add_band()

matplot(y=cbind(l$dat$SWCV), x=l$dat$YEAR, type="l", lty=1, col=c("cyan3", "yellow3"), ylab="SWC (m3 m-3)", xlab="Time (years)")
matlines(y=cbind(fitted(loess(l$dat$SWCV~l$dat$YEAR, span=filter_span))), x=l$dat$YEAR, type="l", lty=1, col="black", lwd=c(1,0.5))
# abline(h=c(0.16), col=scales::muted("cyan3"))
add_band()


if (plot_to_file) dev.off()


### Cohorts plot

# dat_coh_spin = readr::read_csv(paste0(output_dir,"/",expt_dir, "/cohort_props.csv"))
# dat_coh_spin %>% filter(YEAR > min(YEAR)+2) %>% ggplot(aes(x=YEAR, y=height, col=log(1e-3+density), group=cohortID)) + geom_line() + scale_color_viridis_c()

# dat_coh_cont = read.csv("~/codes/Drought_MIP/pfate_output/cont_test_GYF_HB_100y_evol_4.5_cont/cohort_props.csv")
# dat_coh_cont %>% filter(YEAR > min(YEAR)+2) %>% ggplot(aes(x=YEAR, y=height, col=log(1e-3+density), group=cohortID)) + geom_line() + scale_color_viridis_c()
#
# dat_coh_spin %>%
#   filter(YEAR > min(YEAR)+2) %>%
#   rbind(dat_coh_cont) %>%
#   ggplot(aes(x=YEAR, y=height, col=log(1e-3+density), group=cohortID)) + geom_line() + scale_color_viridis_c()
#
# dat_coh_spin %>% with(plot(mort~))

# dat_coh_spin %>% select(mort, dbh, YEAR) %>% mutate(YEAR=as.integer(YEAR)) %>%
#   left_join(l1$dat %>% select(YEAR, SWCV, SWP)) %>%
#   left_join(l1$traits %>% filter(SPP=="spp1") %>% select(YEAR, P50X) %>% mutate(YEAR=as.integer(YEAR))) %>%
#   ggplot(aes(x=dbh, y=mort, col=SWCV)) +
#   geom_point()
