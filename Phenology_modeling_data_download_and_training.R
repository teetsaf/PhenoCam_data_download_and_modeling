
###### LIBRARIES - check for and install packages #######

list.of.packages <- c("phenocamr","phenor","ggplot2", "data.table","tidyr","ggpmisc","plyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages,require,character.only=TRUE)

###### SET WORKING DRIVE - note: the folder where this script is located will be populated with downloaded data ######

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

###### VIEW SITE NAMES AND REGION OF INTERESTS (ROIs) BY SITE #######

sites <- list_sites()
(sites$site)
rois <- list_rois()
subset(rois,rois$site=='harvard')
subset(rois,rois$site=='hubbardbrook')
subset(rois,rois$site=='hubbardbrooksfws')
subset(rois,rois$site=='hubbardbrooknfws')
subset(rois,rois$site=='howland2')
subset(rois,rois$site=='proctor')
subset(rois,rois$site=='caryinstitute')
subset(rois,rois$site=='laurentides')
subset(rois,rois$site=='queens')
subset(rois,rois$site=='turkeypointdbf')

###### USE PHENOCAMR TO DOWNLOAD AND FORMAT DATA FROM SELECTED SITES AND ROIs ######

### download data
download_phenocam(site = "bartlett",
                  veg_type = "DB",
                  roi_id = "3000",
                  frequency = 3,
                  trim=2021,
                  outlier_detection = FALSE,
                  smooth = FALSE,
                  out_dir = getwd())

download_phenocam(site = "bartlettir",
                  veg_type = "DB",
                  roi_id = "2000",
                  frequency = 3,
                  trim=2021,
                  outlier_detection = FALSE,
                  smooth = FALSE,
                  out_dir = getwd())

download_phenocam(site = "hubbardbrook",
                  veg_type = "DB",
                  roi_id = "1000",
                  frequency = 3,
                  trim=2021,
                  outlier_detection = FALSE,
                  smooth = FALSE,
                  out_dir = getwd())

download_phenocam(site = "hubbardbrook",
                  veg_type = "DB",
                  roi_id = "2000",
                  frequency = 3,
                  trim=2021,
                  outlier_detection = FALSE,
                  smooth = FALSE,
                  out_dir = getwd())

download_phenocam(site = "hubbardbrook",
                  veg_type = "DB",
                  roi_id = "3000",
                  frequency = 3,
                  trim=2021,
                  outlier_detection = FALSE,
                  smooth = FALSE,
                  out_dir = getwd())

download_phenocam(site = "caryinstitute",
                  veg_type = "DB",
                  roi_id = "1000",
                  frequency = 3,
                  outlier_detection = FALSE,
                  smooth = FALSE,
                  out_dir = getwd())

download_phenocam(site = "harvard",
                  veg_type = "DB",
                  roi_id = "1000",
                  frequency = 3,
                  trim=2021,
                  outlier_detection = FALSE,
                  smooth = FALSE,
                  out_dir = getwd())

download_phenocam(site = "howland2",
                  veg_type = "DB",
                  roi_id = "1000",
                  frequency = 3,
                  outlier_detection = FALSE,
                  smooth = FALSE,
                  out_dir = getwd())

download_phenocam(site = "howland2",
                  veg_type = "DB",
                  roi_id = "2000",
                  frequency = 3,
                  outlier_detection = FALSE,
                  smooth = FALSE,
                  out_dir = getwd())

download_phenocam(site = "laurentides",
                  veg_type = "DB",
                  roi_id = "1000",
                  frequency = 3,
                  outlier_detection = FALSE,
                  smooth = FALSE,
                  out_dir = getwd())

download_phenocam(site = "proctor",
                  veg_type = "DB",
                  roi_id = "1000",
                  frequency = 3,
                  trim=2021,
                  outlier_detection = FALSE,
                  smooth = FALSE,
                  out_dir = getwd())

download_phenocam(site = "queens",
                  veg_type = "DB",
                  roi_id = "1000",
                  frequency = 3,
                  outlier_detection = FALSE,
                  smooth = FALSE,
                  out_dir = getwd())

download_phenocam(site = "queens",
                  veg_type = "DB",
                  roi_id = "2000",
                  frequency = 3,
                  outlier_detection = FALSE,
                  smooth = FALSE,
                  out_dir = getwd())

download_phenocam(site = "turkeypointdbf",
                  veg_type = "DB",
                  roi_id = "1000",
                  frequency = 3,
                  outlier_detection = FALSE,
                  smooth = FALSE,
                  out_dir = getwd())

download_phenocam(site = "turkeypointdbf",
                  veg_type = "DB",
                  roi_id = "2000",
                  frequency = 3,
                  outlier_detection = FALSE,
                  smooth = FALSE,
                  out_dir = getwd())

################################################################################

### call in downloaded ROIs
bartlett1=read_phenocam(file.path('bartlett_DB_3000_3day.csv'))
bartlett2=read_phenocam(file.path('bartlettir_DB_2000_3day.csv'))
hubbard1=read_phenocam(file.path('hubbardbrook_DB_1000_3day.csv'))
hubbard2=read_phenocam(file.path('hubbardbrook_DB_2000_3day.csv'))
hubbard3=read_phenocam(file.path('hubbardbrook_DB_3000_3day.csv'))
cary=read_phenocam(file.path('caryinstitute_DB_1000_3day.csv'))
harv=read_phenocam(file.path('harvard_DB_1000_3day.csv'))
howl1=read_phenocam(file.path('howland2_DB_1000_3day.csv'))
howl2=read_phenocam(file.path('howland2_DB_2000_3day.csv'))
laur=read_phenocam(file.path('laurentides_DB_1000_3day.csv'))
proc=read_phenocam(file.path('proctor_DB_1000_3day.csv'))
queens1=read_phenocam(file.path('queens_DB_1000_3day.csv'))
queens2=read_phenocam(file.path('queens_DB_2000_3day.csv'))
turkey1=read_phenocam(file.path('turkeypointdbf_DB_1000_3day.csv'))
turkey2=read_phenocam(file.path('turkeypointdbf_DB_2000_3day.csv'))

### combine into a list of lists
allsite=list(bartlett1,bartlett2,cary,harv,howl1,howl2,laur,proc,queens1,queens2,turkey1,turkey2)

### prep for for loop
site_ids=as.factor(c('bartlett1','bartlett2','cary','harv','howl1','howl2','laur','proc','queens1','queens2','turkey1','turkey2'))

formatted_data=data.frame(gcc_value=0,transition_10=as.Date('1990-01-01'),transition_50=as.Date('1990-01-01'),transition_80=as.Date('1990-01-01'),
                          transition_10_lower_ci=as.Date('1990-01-01'),transition_50_lower_ci=as.Date('1990-01-01'),transition_80_lower_ci=as.Date('1990-01-01'),
                          transition_10_upper_ci=as.Date('1990-01-01'),transition_50_upper_ci=as.Date('1990-01-01'),transition_80_upper_ci=as.Date('1990-01-01'),
                          threshold_10=0,threshold_50=0,threshold_80=0,min_gcc=0,max_gcc=0,V2=0,
                          season='spring',coordinate='bcc',site='fake',lat=44,lon=-71)

### for loop extracts transition dates at designated thresholds (10%, 50%, and 80% the seasonal amplitude) for green and yellow chromatic coordinate
for(h in 1:length(site_ids)){
  thisSite = allsite[[h]]
  thisSiteName = site_ids[h]
  thisSite = expand_phenocam(thisSite)
  thisSite = detect_outliers(thisSite)
  thisSite = smooth_ts(thisSite)
  
  ### extract phenology dates
  thisSite_phen_dates_gcc=phenophases(thisSite, lower_thresh=0.1, middle_thresh=0.5, upper_thresh=0.8, internal = TRUE)
  
  ### subset spring dates and reformat
  thisSite_spring=data.table(thisSite_phen_dates_gcc$rising,thisSite_phen_dates_gcc$rising$gcc_value=='gcc_90')
  thisSite_spring[,2:10 := lapply(.SD,function(x) as.Date(x,origin='1970-01-01')), .SDcols = 2:10]
  thisSite_spring$season='rising'
  thisSite_spring$coordinate='gcc'
  
  ### subset fall dates and reformat
  thisSite_fall=data.table(thisSite_phen_dates_gcc$falling,thisSite_phen_dates_gcc$falling$gcc_value=='gcc_90')
  thisSite_fall[,2:10 := lapply(.SD,function(x) as.Date(x,origin='1970-01-01')), .SDcols = 2:10]
  thisSite_fall$season='falling'
  thisSite_fall$coordinate='gcc'
  
  ### convert gcc to ycc
  thisSite$data$gcc_90=thisSite$data$gcc_90+thisSite$data$rcc_90
  thisSite$data$smooth_gcc_90=thisSite$data$smooth_gcc_90+thisSite$data$smooth_rcc_90
  thisSite$data$gcc_75=thisSite$data$gcc_75+thisSite$data$rcc_75
  thisSite$data$smooth_gcc_75=thisSite$data$smooth_gcc_75+thisSite$data$smooth_rcc_75
  thisSite$data$gcc_50=thisSite$data$gcc_50+thisSite$data$rcc_50
  thisSite$data$smooth_gcc_50=thisSite$data$smooth_gcc_50+thisSite$data$smooth_rcc_50
  
  ### extract phenology dates
  thisSite_phen_dates_ycc=phenophases(thisSite, lower_thresh=0.1, middle_thresh=0.5, upper_thresh=0.8, internal = TRUE)
  
  ### subset spring dates and reformat
  thisSite_spring_ycc=data.table(thisSite_phen_dates_ycc$rising,thisSite_phen_dates_ycc$rising$gcc_value=='gcc_90')
  thisSite_spring_ycc[,2:10 := lapply(.SD,function(x) as.Date(x,origin='1970-01-01')), .SDcols = 2:10]
  thisSite_spring_ycc$season='rising'
  thisSite_spring_ycc$coordinate='ycc'
  
  ### subset fall dates and reformat
  thisSite_fall_ycc=data.table(thisSite_phen_dates_ycc$falling,thisSite_phen_dates_ycc$falling$gcc_value=='gcc_90')
  thisSite_fall_ycc[,2:10 := lapply(.SD,function(x) as.Date(x,origin='1970-01-01')), .SDcols = 2:10]
  thisSite_fall_ycc$season='falling'
  thisSite_fall_ycc$coordinate='ycc'
  
  formatted=rbind(thisSite_spring,thisSite_fall,thisSite_spring_ycc,thisSite_fall_ycc)
  formatted$site=thisSiteName
  formatted$lat=thisSite$lat
  formatted$lon=thisSite$lon
  
  formatted_data=rbind(formatted_data,formatted)
}

head(formatted_data)
write.csv(formatted_data[-1,],'transition_dates.csv')
cc=read.csv('transition_dates.csv')
cc

################################################################################

####### FORMATTING FOR TRAINING MODELS

#### format spring transition data

spring=subset(cc,season=='rising' & gcc_value=='gcc_90' & coordinate=='gcc')
spring_10 <- data.frame(ID=spring$site, id=spring$site, lat=spring$lat, lon=spring$lon, phenophase=spring$season, year=year(spring$transition_10), doy=yday(spring$transition_10))

ggplot(spring_10,aes(x=year,y=doy))+
  geom_line()+
  geom_point()+
  facet_wrap(~ID,ncol=5)

# write data for beginning of spring (10% greenness) to temporary CSV
write.table(
  spring_10,
  file.path(tempdir(), "synthesis_data.csv"),
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  sep = ","
)

# format beginning of spring dataset with phenoR
spring_10 <- pr_fm_csv(
  file = file.path(tempdir(), "synthesis_data.csv"),
  phenophase = "rising",
  internal = TRUE
)

save(spring_10, file="training_data_spring_transition10.RData")

# write data for end of spring (80% greenness) to temporary CSV
write.table(
  spring_80,
  file.path(tempdir(), "synthesis_data.csv"),
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  sep = ","
)

# format end of spring dataset with phenoR
spring_80 <- pr_fm_csv(
  file = file.path(tempdir(), "synthesis_data.csv"),
  phenophase = "rising",
  internal = TRUE
)

save(spring_80, file="training_data_spring_transition80.RData")


#### format fall transition data

fall=subset(cc,season=='falling' & gcc_value=='gcc_90' & coordinate=='ycc' & yday(transition_50)>200)
Y50 <- data.frame(ID=fall$site, id=fall$site, lat=fall$lat, lon=fall$lon, phenophase=fall$season, year=year(fall$transition_50), doy=yday(fall$transition_50))

fall=subset(cc,season=='falling' & gcc_value=='gcc_90' & coordinate=='gcc' & yday(transition_50)>250)
### exclude the erroneous date at Cary Institute in 2009
G50 <- data.frame(ID=fall$site, id=fall$site, lat=fall$lat, lon=fall$lon, phenophase=fall$season, year=year(fall$transition_50), doy=yday(fall$transition_50))

ggplot(G50,aes(x=year,y=doy))+
  geom_line()+
  geom_point()+
  facet_wrap(~ID,ncol=5)

# write fake data to temporary CSV
write.table(
  Y50,
  file.path(tempdir(), "synthesis_data.csv"),
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  sep = ","
)

# format a "fake" HB dataset (internally, not saved to disk)
fall_Y50 <- pr_fm_csv(
  file = file.path(tempdir(), "synthesis_data.csv"),
  phenophase = "falling",
  offset=365,
  internal = TRUE
)
head(fall_Y50)

save(fall_Y50, file="training_data_fall_transition_Ycc50.RData")

# write fake data to temporary CSV
write.table(
  G50,
  file.path(tempdir(), "synthesis_data.csv"),
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  sep = ","
)

# format a "fake" HB dataset (internally, not saved to disk)
fall_G50 <- pr_fm_csv(
  file = file.path(tempdir(), "synthesis_data.csv"),
  phenophase = "falling",
  offset=365,
  internal = TRUE
)
head(fall_G50)

save(fall_G50, file="training_data_fall_transition_Gcc50.RData")


######################################################################

########## Run for loops

### load training data
load("training_data_spring_transition10.RData")
load("training_data_fall_transition_Ycc50.RData")


#### loop that trains spring models with phenoR

sp_models=factor(levels=c('LIN','TT','PTT','M1','SQ','PA','SM1','PM1'))
sp_par=list(0)
d1=data.frame(model='test',iteration=0,aic=0,rmse=0)

#### set number of chains and iterations
#### this spring model loop can take several hours to run, to test, reduce the number of chains and iterations to single digits

#chains=25
#itertations=40000

for(h in levels(sp_models)){
  for(i in 1:1){
    pars = pr_fit(
      data = spring_10,
      model = paste(h),
      method = "GenSA",
      par_ranges = file.path('parameter_ranges.csv'),
      control = list(max.call = 4))
    d1=rbind(d1,data.frame(model=paste(h),iteration=paste(i),aic=(pars$aic$AIC),rmse=(pars$rmse)))
    sp_par=append(sp_par,list(pars$par))
  }
}

d1[-1,]
write.csv(d1[-1,],"spring_model_fits_G10.csv")
sp_par2=sp_par[-1]
save(sp_par2, file="spring_parameters_G10.RData")
sp_par2


#### loop that trains fall models with phenoR

fall_models = factor(levels=c("CDD", "CDDP"))
fall_par=list(0)
df1=data.frame(model='test',iteration=0,aic=0,rmse=0)

for(h in levels(fall_models)){
  for(i in 1:25){
    pars = pr_fit(
      data = fall_Y50,
      model = paste(h),
      method = "GenSA",
      par_ranges = file.path('parameter_ranges.csv'),
      control = list(max.call = 40000))
    df1=rbind(df1,data.frame(model=paste(h),iteration=paste(i),aic=(pars$aic$AIC),rmse=(pars$rmse)))
    fall_par=append(fall_par,list(pars$par))
  }}

df1
write.csv(df1[-1,],"fall_model_fits_Y50.csv")
fall_par2=fall_par[-1]
save(fall_par2, file="fall_parameters_Y50.RData")

load("fall_parameters_Y50.RData")

######################################################################

#### Use parameters to predict transition dates and visualize results

#### Spring

load("spring_parameters_G10.RData")
sp_par2

d1=read.csv("spring_model_fits_G10.csv")
ddt=data.table(d1)
ddt[ddt[,.I[which.min(aic)],model]$V1][order(aic)]
s10=ddt[ddt[,.I[which.min(aic)],model]$V1]
s10_pars=list(LIN=sp_par2[s10[model=='LIN']$iteration],TT=sp_par2[s10[model=='TT']$iteration+25],
              PTT=sp_par2[s10[model=='PTT']$iteration+50],M1=sp_par2[s10[model=='M1']$iteration+75],
              SQ=sp_par2[s10[model=='SQ']$iteration+100],PA=sp_par2[s10[model=='PA']$iteration+125],
              SM1=sp_par2[s10[model=='SM1']$iteration+150],PM1=sp_par2[s10[model=='PM1']$iteration+175])

#### clunky way to predict and populate a dataframe with observed vs. predicted
sp_spring_10 <- data.frame(LIN=pr_predict(spring_10, model = "LIN", par = s10_pars$LIN[[1]]),
                          TT=pr_predict(spring_10, model = "TT", par = s10_pars$TT[[1]]),
                          PTT=pr_predict(spring_10, model = "PTT", par = s10_pars$PTT[[1]]),
                          M1=pr_predict(spring_10, model = "M1", par = s10_pars$M1[[1]]),
                          SQ=pr_predict(spring_10, model = "SQ", par = s10_pars$SQ[[1]]),
                          PA=pr_predict(spring_10, model = "PA", par = s10_pars$PA[[1]]),
                          SM1=pr_predict(spring_10, model = "SM1", par = s10_pars$SM1[[1]]),
                          PM1=pr_predict(spring_10, model = "PM1", par = s10_pars$PM1[[1]]),
                    obs=c(spring_10$bartlett1$transition_dates,spring_10$bartlett2$transition_dates,
                          #spring_10$hubbard1$transition_dates,spring_10$hubbard2$transition_dates,spring_10$hubbard3$transition_dates,
                          spring_10$cary$transition_dates,spring_10$harv$transition_dates,spring_10$howl1$transition_dates,spring_10$howl2$transition_dates,
                          spring_10$laur$transition_dates,spring_10$proc$transition_dates,spring_10$queens1$transition_dates,spring_10$queens2$transition_dates,
                          spring_10$turkey1$transition_dates,spring_10$turkey2$transition_dates),
                    year=c(spring_10$bartlett1$year,spring_10$bartlett2$year,
                           #spring_10$hubbard1$year,spring_10$hubbard2$year,subset(spring_10$hubbard3$year,spring_10$hubbard3$year<2022),
                           subset(spring_10$cary$year,spring_10$cary$year<2022),spring_10$harv$year,spring_10$howl1$year,spring_10$howl2$year,
                           spring_10$laur$year,spring_10$proc$year,spring_10$queens1$year,subset(spring_10$queens2$year,spring_10$queens2$year<2022),
                           spring_10$turkey1$year,subset(spring_10$turkey2$year,spring_10$turkey2$year<2022)),
                    site=c(rep(spring_10$bartlett1$site,length(spring_10$bartlett1$year)),rep(spring_10$bartlett2$site,length(spring_10$bartlett2$year)),
                           #rep(spring_10$hubbard1$site,length(spring_10$hubbard1$year)),rep(spring_10$hubbard2$site,length(spring_10$hubbard2$year)),rep(spring_10$hubbard3$site,length(spring_10$hubbard3$year)-1),
                           rep(spring_10$cary$site,length(spring_10$cary$year)-1),rep(spring_10$harv$site,length(spring_10$harv$year)),
                           rep(spring_10$howl1$site,length(spring_10$howl1$year)),rep(spring_10$howl2$site,length(spring_10$howl2$year)),
                           rep(spring_10$laur$site,length(spring_10$laur$year)),rep(spring_10$proc$site,length(spring_10$proc$year)),
                           rep(spring_10$queens1$site,length(spring_10$queens1$year)),rep(spring_10$queens2$site,length(spring_10$queens2$year)-1),
                           rep(spring_10$turkey1$site,length(spring_10$turkey1$year)),rep(spring_10$turkey2$site,length(spring_10$turkey2$year)-1)))

sp_10=gather(sp_spring_10,model,modeled,-year,-site,-obs)
sp_10$SITE=substr(sp_10$site,1,4)

sp10_dt=data.table(sp_10)
sp10_dt$season='Spring'
sp10_dt[,.(r2=cor(modeled,obs)^2,rmse=sqrt(mean((modeled-obs)^2))),.(model)]

ggplot(sp10_dt,aes(x=obs,y=modeled,color=SITE,shape=SITE))+
  geom_point()+
  geom_abline()+
  scale_shape_manual(values=c(21,22,23,24,25,15,16,17,18,19))+
  facet_wrap(~model)

write.csv(sp10_dt,'spring_mod_vs_obs_G10.csv')

#### Fall

load("fall_parameters_Y50.RData")
fall_par=fall_par2
df1=read.csv("fall_model_fits_Y50.csv")
#df1=read.csv("fall_model_fits_G50.csv")

dfdt=data.table(df1)
dfdt[dfdt[,.I[which.min(aic)],model]$V1]
fY50=dfdt[dfdt[,.I[which.min(aic)],model]$V1]
fY50_pars=list(CDD=fall_par[fY50[model=='CDD']$iteration],CDDP=fall_par[fY50[model=='CDDP']$iteration+25])
fY50_pars
#### clunky way to predict and populate a dataframe with observed vs. predicted
fa_fall_Y50 <- data.frame(CDD=pr_predict(fall_Y50, model = "CDD", par = fY50_pars$CDD[[1]]),
                    CDDP=pr_predict(fall_Y50, model = "CDDP", par = fY50_pars$CDDP[[1]]),
                    obs=c(fall_Y50$bartlett1$transition_dates,fall_Y50$bartlett2$transition_dates,
                          fall_Y50$hubbard1$transition_dates,fall_Y50$hubbard2$transition_dates,fall_Y50$hubbard3$transition_dates,
                          fall_Y50$cary$transition_dates,fall_Y50$harv$transition_dates,fall_Y50$howl1$transition_dates,fall_Y50$howl2$transition_dates,
                          fall_Y50$laur$transition_dates,fall_Y50$proc$transition_dates,fall_Y50$queens1$transition_dates,fall_Y50$queens2$transition_dates,
                          fall_Y50$turkey1$transition_dates,fall_Y50$turkey2$transition_dates),
                    year=c(fall_Y50$bartlett1$year,fall_Y50$bartlett2$year,
                           fall_Y50$hubbard1$year,fall_Y50$hubbard2$year,subset(fall_Y50$hubbard3$year,fall_Y50$hubbard3$year<2022),
                           subset(fall_Y50$cary$year,fall_Y50$cary$year<2022),fall_Y50$harv$year,fall_Y50$howl1$year,fall_Y50$howl2$year,
                           fall_Y50$laur$year,fall_Y50$proc$year,fall_Y50$queens1$year,subset(fall_Y50$queens2$year,fall_Y50$queens2$year<2022),
                           fall_Y50$turkey1$year,subset(fall_Y50$turkey2$year,fall_Y50$turkey2$year<2022)),
                    site=c(rep(fall_Y50$bartlett1$site,length(fall_Y50$bartlett1$year)),rep(fall_Y50$bartlett2$site,length(fall_Y50$bartlett2$year)),
                           rep(fall_Y50$hubbard1$site,length(fall_Y50$hubbard1$year)),rep(fall_Y50$hubbard2$site,length(fall_Y50$hubbard2$year)),rep(fall_Y50$hubbard3$site,length(fall_Y50$hubbard3$year)),
                           rep(fall_Y50$cary$site,length(fall_Y50$cary$year)),rep(fall_Y50$harv$site,length(fall_Y50$harv$year)),
                           rep(fall_Y50$howl1$site,length(fall_Y50$howl1$year)),rep(fall_Y50$howl2$site,length(fall_Y50$howl2$year)),
                           rep(fall_Y50$laur$site,length(fall_Y50$laur$year)),rep(fall_Y50$proc$site,length(fall_Y50$proc$year)),
                           rep(fall_Y50$queens1$site,length(fall_Y50$queens1$year)),rep(fall_Y50$queens2$site,length(fall_Y50$queens2$year)),
                           rep(fall_Y50$turkey1$site,length(fall_Y50$turkey1$year)),rep(fall_Y50$turkey2$site,length(fall_Y50$turkey2$year))))



fa_Y50=gather(fa_fall_Y50,model,modeled,-year,-site,-obs)
fa_Y50$SITE=substr(fa_Y50$site,1,4)

fa_Y50_dt=data.table(fa_Y50)
fa_Y50_dt$season='Fall'
fa_Y50_dt[,.(r2=cor(modeled,obs)^2,rmse=sqrt(mean((modeled-obs)^2))),.(model)]

write.csv(fa_Y50_dt,'fall_mod_vs_obs_Y50.csv')



######## fit models to predict the Fcrit required to end spring and start fall using other parameters fixed (from previous parameterization)
###### this step requires manually changing the parameter file (CSV file) to "fix" all parameters other than Fcrit


### load spring 80
load("training_data_spring_transition80.RData")
spring_80

sp_models2=factor(levels=c('LIN','TT','PTT','M1','SQ','PA','SM1','PM1'))
sp_par80=list(0)
d2=data.frame(model='test',iteration=0,aic=0,rmse=0)

for(h in levels(sp_models2)){
  for(i in 1:1){
    pars = pr_fit(
      data = spring_80,
      model = paste(h),
      method = "GenSA",
      par_ranges = file.path('parameter_ranges_constrained_AFT.csv'),
      control = list(max.call = 40000))
    d2=rbind(d2,data.frame(model=paste(h),iteration=paste(i),aic=(pars$aic$AIC),rmse=(pars$rmse)))
    sp_par80=append(sp_par80,list(pars$par))
  }
}

d2[-1,]
write.csv(d2[-1,],"spring_model_fits_G80.csv")
sp_G80=sp_par80[-1]
save(sp_G80, file="spring_parameters_G80.RData")
sp_G80

ddt80=data.table(d2[-1,])
ddt80[ddt80[,.I[which.min(aic)],model]$V1][order(aic)]
s80=ddt80[ddt80[,.I[which.min(aic)],model]$V1]
s80_pars=list(LIN=sp_G80[1],TT=sp_G80[2],
              PTT=sp_G80[3],M1=sp_G80[4],
              SQ=sp_G80[5],PA=sp_G80[6],
              SM1=sp_G80[7],PM1=sp_G80[8])
s80_pars

#### loop that trains fall models with phenoR

fY50_pars$CDD
fY50_pars$CDDP

load("training_data_fall_transition_Gcc50.RData")
fall_G50

fall_models = factor(levels=c("CDD", "CDDP"))
fall_parG50=list(0)
df2=data.frame(model='test',iteration=0,aic=0,rmse=0)

for(h in levels(fall_models)){
  for(i in 1:1){
    pars = pr_fit(
      data = fall_G50,
      model = paste(h),
      method = "GenSA",
      par_ranges = file.path('parameter_ranges_constrained_AFT.csv'),
      control = list(max.call = 40000))
    df2=rbind(df2,data.frame(model=paste(h),iteration=paste(i),aic=(pars$aic$AIC),rmse=(pars$rmse)))
    fall_parG50=append(fall_parG50,list(pars$par))
  }}

df2
write.csv(df2[-1,],"fall_model_fits_G50_.csv")
fG50_pars=fall_parG50[-1]
save(fG50_pars, file="fall_parameters_G50.RData")
fG50_pars

df2=read.csv("fall_model_fits_G50.csv")
load("fall_parameters_G50.RData")

ddtG50=data.table(df2[-1,])
ddtG50[ddtG50[,.I[which.min(aic)],model]$V1][order(aic)]
fG50_pars=list(CDD=fG50_pars[1],CDDP=fG50_pars[2])
fG50_pars


######################################################################


#parameters

### Start of spring parameters

load("spring_parameters_G10.RData")
sp_par2
d1=read.csv("spring_model_fits_G10.csv")
ddt=data.table(d1)
s10=ddt[ddt[,.I[which.min(aic)],model]$V1]
s10_pars=list(LIN=sp_par2[s10[model=='LIN']$iteration],TT=sp_par2[s10[model=='TT']$iteration+25],
              PTT=sp_par2[s10[model=='PTT']$iteration+50],M1=sp_par2[s10[model=='M1']$iteration+75],
              SQ=sp_par2[s10[model=='SQ']$iteration+100],PA=sp_par2[s10[model=='PA']$iteration+125],
              SM1=sp_par2[s10[model=='SM1']$iteration+150],PM1=sp_par2[s10[model=='PM1']$iteration+175])
s10_pars

### End of spring parameters
load("spring_parameters_G80.RData")
sp_G80
s80_pars=list(LIN=sp_G80[1],TT=sp_G80[2],PTT=sp_G80[3],M1=sp_G80[4],SQ=sp_G80[5],PA=sp_G80[6],SM1=sp_G80[7],PM1=sp_G80[8])
s80_pars

### End of fall parameters
load("fall_parameters_Y50.RData")
fall_par=fall_par2
df1=read.csv("fall_model_fits.csv")
dfdt=data.table(df1)
fY50=dfdt[dfdt[,.I[which.min(aic)],model]$V1]
fY50_pars=list(CDD=fall_par[fY50[model=='CDD']$iteration],CDDP=fall_par[fY50[model=='CDDP']$iteration+25])
fY50_pars

### Start of fall parameters
load("fall_parameters_G50.RData")
fG50_pars=list(CDD=fG50_pars[1],CDDP=fG50_pars[2])
fG50_pars

######################################################################


### predict phenology at Hubbard Brook

year <- seq(1993,2020,1)
lat <- 43.941613
lon <- -71.747016
id <- "HB"
ID <- "HB"
phenophase <- "rising"
doy <- 0
df <- data.frame(ID, id, lat, lon, phenophase, year, doy)
str(df)

# write fake data to temporary CSV
write.table(
  df,
  file.path(tempdir(), "HB.csv"),
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  sep = ","
)

# format a "fake" HB dataset (internally, not saved to disk)
HB_spring <- pr_fm_csv(
  file = file.path(tempdir(), "HB.csv"),
  phenophase = "rising",
  offset = 264,
  internal = TRUE
)

##### predict spring transition dates using G10 parameters
HB_sp <- data.frame(year=1993:2020,
                         LIN=pr_predict(HB_spring, model = "LIN", par = s10_pars$LIN[[1]]),
                         TT=pr_predict(HB_spring, model = "TT", par = s10_pars$TT[[1]]),
                         PTT=pr_predict(HB_spring, model = "PTT", par = s10_pars$PTT[[1]]),
                         M1=pr_predict(HB_spring, model = "M1", par = s10_pars$M1[[1]]),
                         SQ=pr_predict(HB_spring, model = "SQ", par = s10_pars$SQ[[1]]),
                         PM1=pr_predict(HB_spring, model = "PM1", par = s10_pars$PM1[[1]]),
                         PA=pr_predict(HB_spring, model = "PA", par = s10_pars$PA[[1]]),
                         SM1=pr_predict(HB_spring, model = "SM1", par = s10_pars$SM1[[1]]))
    

HB_sp              
HBspring=gather(HB_sp,model,prediction,2:9)

##### predict fall transition dates using Y50 parameters
year <- rev(seq(1993,2020,1))
phenophase_fall <- "falling"
df <- data.frame(ID, id, lat, lon, phenophase_fall, year, doy)

# write fake data to temporary CSV
write.table(
  df,
  file.path(tempdir(), "HB.csv"),
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  sep = ","
)

HB_fall <- pr_fm_csv(
  file = file.path(tempdir(), "HB.csv"),
  phenophase = "falling",
  offset = 365,
  internal = TRUE
)

######################################################################

### predict future phenology at Bartlett Experimental Forest using data downloaded from NA-CORDEX (https://na-cordex.org/)
### downloaded and formatted data available on GitHub


year <- seq(1980,2020,1)
lat <- 44.0646
lon <- -71.2881
id <- "bart0"
ID <- "bart0"
phenophase <- "rising"
doy <- 0
df <- data.frame(ID, id, lat, lon, phenophase, year, doy)
str(df)

# write fake data to temporary CSV
write.table(
  df,
  file.path(tempdir(), "bart0.csv"),
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  sep = ","
)

# format a "fake" HB dataset (internally, not saved to disk)
bart0 <- pr_fm_csv(
  file = file.path(tempdir(), "bart0.csv"),
  phenophase = "rising",
  offset = 264,
  internal = TRUE
)


### predict future 
str(bart0)
bart0$bart0$transition_dates=rep(120,94)
bart0$bart0$year=2007:2100
bart0$bart0$Ti=cbind(bart0$bart0$Ti,bart0$bart0$Ti,bart0$bart0$VPDi[,1:14])
bart0$bart0$Tmini=cbind(bart0$bart0$Tmini,bart0$bart0$Tmini,bart0$bart0$VPDi[,1:14])
bart0$bart0$Tmaxi=cbind(bart0$bart0$Tmaxi,bart0$bart0$Tmaxi,bart0$bart0$VPDi[,1:14])
bart0$bart0$Li=cbind(bart0$bart0$Li,bart0$bart0$Li,bart0$bart0$Li[,1:14])
bart0$bart0$Pi=cbind(bart0$bart0$Pi,bart0$bart0$Pi,bart0$bart0$Pi[,1:14])
bart0$bart0$VPDi=cbind(bart0$bart0$VPDi,bart0$bart0$VPDi,bart0$bart0$VPDi[,1:14])

### read in CODEX
codex_min=read.csv('Tmin_codex_2022_v2.csv')
codex_max=read.csv('Tmax_codex_2022_v2.csv')

LIN_future <- data.frame(year=2007:2100,model='LIN')
M1_future <- data.frame(year=2007:2100,model='M1')
PA_future <- data.frame(year=2007:2100,model='PA')
PM1_future <- data.frame(year=2007:2100,model='PM1')
PTT_future <- data.frame(year=2007:2100,model='PTT')
SM1_future <- data.frame(year=2007:2100,model='SM1')
SQ_future <- data.frame(year=2007:2100,model='SQ')
TT_future <- data.frame(year=2007:2100,model='TT')


for(i in 1:12){
  bart0$bart0$Tmaxi=matrix(codex_max[264:34573,6+i],nrow=365,ncol=94)
  bart0$bart0$Tmini=matrix(codex_min[264:34573,6+i],nrow=365,ncol=94)
  bart0$bart0$Ti=matrix((codex_min[264:34573,6+i]+codex_max[264:34573,6+i])/2,nrow=365,ncol=94)
  LIN_future[2+i]=pr_predict(bart0, model = 'LIN', par = s10_pars$LIN[[1]])
  M1_future[2+i]=pr_predict(bart0, model = 'M1', par = s10_pars$M1[[1]])
  PA_future[2+i]=pr_predict(bart0, model = 'PA', par = s10_pars$PA[[1]])
  PM1_future[2+i]=pr_predict(bart0, model = 'PM1', par = s10_pars$PM1[[1]])
  PTT_future[2+i]=pr_predict(bart0, model = 'PTT', par = s10_pars$PTT[[1]])
  SM1_future[2+i]=pr_predict(bart0, model = 'SM1', par = s10_pars$SM1[[1]])
  SQ_future[2+i]=pr_predict(bart0, model = 'SQ', par = s10_pars$SQ[[1]])
  TT_future[2+i]=pr_predict(bart0, model = 'TT', par = s10_pars$TT[[1]])
  colnames(LIN_future)[2+i]=colnames(codex_min)[6+i]
  colnames(M1_future)[2+i]=colnames(codex_min)[6+i]
  colnames(PA_future)[2+i]=colnames(codex_min)[6+i]
  colnames(PM1_future)[2+i]=colnames(codex_min)[6+i]
  colnames(PTT_future)[2+i]=colnames(codex_min)[6+i]
  colnames(SM1_future)[2+i]=colnames(codex_min)[6+i]
  colnames(SQ_future)[2+i]=colnames(codex_min)[6+i]
  colnames(TT_future)[2+i]=colnames(codex_min)[6+i]
}


for(i in 1:12){
  LIN_future[14+i]=0
  M1_future[14+i]=0
  PA_future[14+i]=0
  PM1_future[14+i]=0
  PTT_future[14+i]=0
  SM1_future[14+i]=0
  SQ_future[14+i]=0
  TT_future[14+i]=0
}


for(i in 1:12){
  bart0$bart0$Tmaxi=matrix(codex_max[264:34573,6+i],nrow=365,ncol=94)
  bart0$bart0$Tmini=matrix(codex_min[264:34573,6+i],nrow=365,ncol=94)
  bart0$bart0$Ti=matrix((codex_min[264:34573,6+i]+codex_max[264:34573,6+i])/2,nrow=365,ncol=94)
  LIN_future[14+i]=pr_predict(bart0, model = 'LIN', par = s80_pars$LIN[[1]])
  M1_future[14+i]=pr_predict(bart0, model = 'M1', par = s80_pars$M1[[1]])
  PA_future[14+i]=pr_predict(bart0, model = 'PA', par = s80_pars$PA[[1]])
  PM1_future[14+i]=pr_predict(bart0, model = 'PM1', par = s80_pars$PM1[[1]])
  PTT_future[14+i]=pr_predict(bart0, model = 'PTT', par = s80_pars$PTT[[1]])
  SM1_future[14+i]=pr_predict(bart0, model = 'SM1', par = s80_pars$SM1[[1]])
  SQ_future[14+i]=pr_predict(bart0, model = 'SQ', par = s80_pars$SQ[[1]])
  TT_future[14+i]=pr_predict(bart0, model = 'TT', par = s80_pars$TT[[1]])
  colnames(LIN_future)[14+i]=colnames(codex_min)[6+i]
  colnames(M1_future)[14+i]=colnames(codex_min)[6+i]
  colnames(PA_future)[14+i]=colnames(codex_min)[6+i]
  colnames(PM1_future)[14+i]=colnames(codex_min)[6+i]
  colnames(PTT_future)[14+i]=colnames(codex_min)[6+i]
  colnames(SM1_future)[14+i]=colnames(codex_min)[6+i]
  colnames(SQ_future)[14+i]=colnames(codex_min)[6+i]
  colnames(TT_future)[14+i]=colnames(codex_min)[6+i]
}


spring_fut=rbind(LIN_future,M1_future,PA_future,PM1_future,PTT_future,SM1_future,SQ_future,TT_future)
head(future)

fut_start=gather(spring_fut[1:14],projection,doy_start,-year,-model)
fut_start$scenario=paste("RCP_",substr(fut_start$projection, nchar(fut_start$projection)-1, nchar(fut_start$projection)),sep="")

fut_end=gather(spring_fut[,c(1,2,15:26)],projection,doy_end,-year,-model)
fut_end$scenario=paste("RCP_",substr(fut_end$projection, nchar(fut_end$projection)-1, nchar(fut_end$projection)),sep="")

fut=data.table(merge(fut_start,fut_end))
fut

sft_summ=fut[,.(mean_start=mean(doy_start),mean_end=mean(doy_end),sd_start=sd(doy_start),sd_end=sd(doy_end)),.(year,model,scenario)]
sft_summ

spring_future=merge(fut,sft_summ)
spring_future


### predict future fall senescence

year <- 2020:1980
lat <- 44.0646
lon <- -71.2881
id <- "bart0"
ID <- "bart0"
phenophase <- "falling"
doy <- 0
df <- data.frame(ID, id, lat, lon, phenophase, year, doy)
str(df)

# write fake data to temporary CSV
write.table(
  df,
  file.path(tempdir(), "bart0_data_future.csv"),
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  sep = ","
)

# format a "fake" bart0 dataset (internally, not saved to disk)
bart0f <- pr_fm_csv(
  file = file.path(tempdir(), "bart0_data_future.csv"),
  phenophase = "falling",
  offset = 365,
  internal = TRUE
)

str(bart0f)
bart0f$bart0$transition_dates=rep(240,94)
bart0f$bart0$year=2100:2007
bart0f$bart0$Ti=cbind(bart0f$bart0$Ti,bart0f$bart0$Ti,bart0f$bart0$VPDi[,1:14])
bart0f$bart0$Tmini=cbind(bart0f$bart0$Tmini,bart0f$bart0$Tmini,bart0f$bart0$VPDi[,1:14])
bart0f$bart0$Tmaxi=cbind(bart0f$bart0$Tmaxi,bart0f$bart0$Tmaxi,bart0f$bart0$VPDi[,1:14])
bart0f$bart0$Li=cbind(bart0f$bart0$Li,bart0f$bart0$Li,bart0f$bart0$Li[,1:14])
bart0f$bart0$Pi=cbind(bart0f$bart0$Pi,bart0f$bart0$Pi,bart0f$bart0$Pi[,1:14])
bart0f$bart0$VPDi=cbind(bart0f$bart0$VPDi,bart0f$bart0$VPDi,bart0f$bart0$VPDi[,1:14])
str(bart0f)
### read in CODEX
codex_min=read.csv('Tmin_codex_2022_v2.csv')
codex_max=read.csv('Tmax_codex_2022_v2.csv')

CDD_future <- data.frame(year=2100:2007,model='CDD')
CDDP_future <- data.frame(year=2100:2007,model='CDDP')

codex_min_rev=codex_min[order(-codex_min$year,codex_min$doy),]
head(codex_min_rev)
codex_max_rev=codex_max[order(-codex_max$year,codex_max$doy),]
head(codex_max_rev)


for(i in 1:12){
  bart0f$bart0$Tmaxi=matrix(codex_max_rev[1:34310,6+i],nrow=365,ncol=94)
  bart0f$bart0$Tmini=matrix(codex_min_rev[1:34310,6+i],nrow=365,ncol=94)
  bart0f$bart0$Ti=matrix((codex_min_rev[1:34310,6+i]+codex_max_rev[1:34310,6+i])/2,nrow=365,ncol=94)
  CDD_future[2+i]=pr_predict(bart0f, model = 'CDD', par = fY50_pars$CDD[[1]])
  CDDP_future[2+i]=pr_predict(bart0f, model = 'CDDP', par = fY50_pars$CDDP[[1]])
  colnames(CDD_future)[2+i]=colnames(codex_min)[6+i]
  colnames(CDDP_future)[2+i]=colnames(codex_min)[6+i]
} 

for(i in 1:12){
  CDD_future[14+i]=0
  CDDP_future[14+i]=0
}

for(i in 1:12){
  bart0f$bart0$Tmaxi=matrix(codex_max_rev[1:34310,6+i],nrow=365,ncol=94)
  bart0f$bart0$Tmini=matrix(codex_min_rev[1:34310,6+i],nrow=365,ncol=94)
  bart0f$bart0$Ti=matrix((codex_min_rev[1:34310,6+i]+codex_max_rev[1:34310,6+i])/2,nrow=365,ncol=94)
  CDD_future[14+i]=pr_predict(bart0f, model = 'CDD', par = fG50_pars$CDD[[1]])
  CDDP_future[14+i]=pr_predict(bart0f, model = 'CDDP', par = fG50_pars$CDDP[[1]])
  colnames(CDD_future)[14+i]=colnames(codex_min)[6+i]
  colnames(CDDP_future)[14+i]=colnames(codex_min)[6+i]
} 


fall_fut=rbind(CDD_future,CDDP_future)

fall_start=gather(fall_fut[1:14],projection,doy_start,-year,-model)
fall_start$scenario=paste("RCP_",substr(fall_start$projection, nchar(fall_start$projection)-1, nchar(fall_start$projection)),sep="")

fall_end=gather(fall_fut[,c(1,2,15:26)],projection,doy_end,-year,-model)
fall_end$scenario=paste("RCP_",substr(fall_end$projection, nchar(fall_end$projection)-1, nchar(fall_end$projection)),sep="")

fut=data.table(merge(fall_start,fall_end))
fut_summ=fut[,.(mean_start=mean(doy_start),mean_end=mean(doy_end),sd_start=sd(doy_start),sd_end=sd(doy_end)),.(year,model,scenario)]
fall_future=merge(fut,fut_summ)




### merge future spring and fall projections

spring_future$season='Spring'
fall_future$season='Fall'

future=rbind(spring_future,fall_future)

write.csv(future,'future_phenology_projections.csv')








