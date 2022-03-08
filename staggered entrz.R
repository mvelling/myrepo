#load in packages
library(survival)
library(ggfortify)
library(mixPHM)
library(smoothHR)
library(survminer)
library(MuMIn)
library(tidyverse)
library(sm)
# after Henriks email clarification:
#The way I thought t_start will always be smaller than t_end

#For example:
  
#  For the hunting part
#First attack, day = 0
#Lynx shot, day = 5 (this will be the new t_start)
#Second affect, day =10 (this will be the t_end) - as before
#This results in short "exposure time" only 5 days (t_end minus t_start)

#For the non-hunting part
#First attack, day = 0 (this will be the t_start) - as before
#Second affect, day =10 (this will be the t_end) - as before

#For the non-hunting part (with a lynx shot after the second attack)
#First attack, day = 0 (this will be the t_start) - as before
#Second affect, day =10 (this will be the t_end) - as before
#Lynx shot, day = 15 (this data will anyway be "non-hunting part") - as before


##### Data handling ######
dfst <- df

for (i in 1:nrow(dfst)) {
  if(!is.na(dfst$time_to_hunt[i])) {
    dfst$tstart[i] <- dfst$time_to_hunt[i]
  }
}

dfst <- dfst %>% select(
  tstart, timespan, hunt, waterdist, roe_dens, first_event, status, year, famgroup_nr
)

for (i in 1:nrow(dfst)) {
  if(dfst$timespan[i] == 366){
    dfst$status[i] <- 0
  }
}
#I need to change the one time where timespan is 0.7 to 1
dfst$timespan[1061] <- 1


dfst$year <- as.numeric(dfst$year)
hist(dfst$year)


#If I want to omit all the nas
#dfst <- na.omit(dfst)
hist(dfst$year)
#now we can check the new model:
dfstweek <- dfst
dfstweek$week <- as.integer(dfst$timespan/7)+1
dfstweek$weekstart <- as.integer(dfst$tstart/7)+0.9
dfstweek <- subset(dfstweek,status == 1)
kaplan_week <- survfit(Surv(weekstart,week,status)~1, dfstweek)
summary(kaplan_week)
autoplot(kaplan_week)
 ##### Kaplan model #####
kaplan_new <- survfit(Surv(tstart, timespan, status)~ 1, dfst)
summary(kaplan_new)
autoplot(kaplan_new)


kaplan_data <- data.frame(n.risk = rep(0,46), n.event = rep(0,46), surv = rep(0,46), time = rep(0,46))
kaplan_data$n.risk <- as.data.frame(kaplan_week$n.risk)
kaplan_data$n.event <- kaplan_week$n.event
kaplan_data$surv <- kaplan_week$surv
kaplan_data$time <- kaplan_week$time
kaplan_data$risk[1] <- 1 - kaplan_data$surv[1]

for (i in 2:nrow(kaplan_data)) {
  kaplan_data$risk[i] <- kaplan_data$surv[i-1] - kaplan_data$surv[i]
}

plot(kaplan_data$time, kaplan_data$risk)


Surv(dfst$tstart)
hist(dfst$tstart[which(dfst$tstart != 0)])
hist(df$time_to_hunt)

dfst <- dfst %>% mutate(
  famgroup_presence = cut(famgroup_nr, breaks = c(- Inf, 0, 2, Inf), labels = c(0,1,2)),# 0 is absence, 1 is 1-2, 2 is more than 2
)

#what we still need to think about. This analysis only considers data from 2009 till 2019, all other data has been disregarded
df$status[which(df$week == 53)] <- 0
df$year <- as.numeric(df$year)
hist(df$year)
yearstatus <- df %>% group_by(year) %>% summarize(sum(status))
plot(yearstatus$year, yearstatus$`sum(status)`, type = "h")
#### Cox Model #####
cox_st <- coxph(Surv(tstart,timespan,status) ~ hunt + famgroup_presence + waterdist + roe_dens, data = dfst, na.action = "na.fail")
summary(cox_st)
cox.zph(cox_st)
dredge(cox_st) #so seems like the best fit model according to AICc is actually the full model, family groups seem to be the most important, then hunt, then roe dens then water

#### ggadjustedcurves ######
cox_st_graph <- ggadjustedcurves(cox_st, data = dfst,conf.int = TRUE, method = "average", variable = "hunt", palette = "hue", ggtheme = theme_grey(), reference = split1, ylim = c(0.7,1))
cox_st_graph + labs(title = "Survival rate for repeated attacks for hunting events") + xlab("Time (days)") + scale_color_manual(values = c("blue","red"),name = "", labels = c("no hunt", "succesfull hunt"))

cox_st_graph_fam <- ggadjustedcurves(cox_st, data = dfst,conf.int = TRUE, method = "average", variable = "famgroup_presence", palette = "hue", ggtheme = theme_grey(), reference = split1, ylim = c(0,1))
cox_st_graph_fam +labs(title = "Survival rate for repeated attacks for family group presences") + xlab("Time (days)")+ scale_color_manual(values = c("green", "blue", "red")
                                                                                                                   ,name = "Family group presence", labels = c("no presence", "1-2 family groups", "3+ family groups"))

testdft <- dfst %>%  select(tstart, timespan, status, hunt, famgroup_presence, waterdist, roe_dens)

ggcompetingrisks(kaplan_new)

##### Hazard ratio Plots ######
library(mixPHM)
library(smoothHR)
library(Greg)
library(simPH)
hr1 <- smoothHR::smoothHR(data = dfst, time = "tstart", time2 = "timespan", status = "status", formula = ~hunt + famgroup_presence + waterdist + roe_dens )
#hr1 <- smoothHR(data = dfst, coxfit =  cox_st)
smoothHR::plot.HR(hr1, predictor = "roe_dens")
plot.HR(hr1, predictor = "waterdist")
smoothHR::predict.HR(hr1, predictor = "waterdist", prob = 0, conf.level = 0.95)
Greg::plotHR(cox_st, term = "waterdist", rug = "density", xlim = c(0,2400), xlab = "Distance to water in meters")
Greg::plotHR(cox_st, term = "roe_dens", rug = "density", xlim = c(0,16), xlab = "Number of roe-deer harvested per 1000 ha")
simGG(cox_st)
ggforest(cox_st)











###################### This is old stuff that didn't work #########################
######################### Making data #####################
dfst <- workdf
dfst <- as.data.frame(dfst)

dfst$status <- 1

dfst$year <- format(dfst$Damage_date_from, "%Y")
dfst$week[is.na(dfst$week)] <- 53
dfst$status[which(dfst$week == 53)] <- 0
str(dfst)
dfst$hunt <- as.factor(dfst$hunt)
dfst$status <- as.factor(dfst$status)
dfst$log1_un <- log(dfst$un_dist+1)
dfst$log1_ge <- log(dfst$ge_dist+1)
dfst$log1_sett <- log(dfst$sett_dist+1)
dfst$status <- as.numeric(dfst$status)
Surv(dfst$week, dfst$status)
#hunting has been calculated wrong, So I need to recategorize it, hunt is marked at 1 if it occured either in the last
#365 days for no repeated and within attack period  between the attacks

### this part get's taken out since I don't want to exclude hunting after lynx
#df$hunt <- if_else(!is.na(df$time_to_hunt), 1, 0 )
#for (i in 1:nrow(df)) {
#  if(!is.na(df$timespan[i]) & !is.na(df$time_to_hunt[i])){
#    if(df$hunt[i] == 1 & df$timespan[i] < df$time_to_hunt[i]) {
#      df$hunt[i] <- 0
#    }
#  }
#}
#for (i in 1:nrow(df)) {
#  if(df$hunt[i] == 0 & !is.na(df$time_to_hunt[i])){
#    df$time_to_hunt[i] <- NA
#  }
#  
#}


#df$time_to_hunt[which(df == 53)] <- 366
dfst$timespan[which(dfst$week == 53)] <- 366
dfst$tstart <- 0

dfst$tstart[which(dfst$timespan < dfst$time_to_hunt)] <- 5

for (i in 1:nrow(dfst)) {
  if(!is.na(dfst$time_to_hunt[i])) {
    if(dfst$timespan[i] < dfst$time_to_hunt[i]){
    dfst$tstart[i] <- dfst$time_to_hunt[i]
  }}
}

dfst <- dfst %>% mutate(
  famgroup_presence = cut(famgroup_nr, breaks = c(- Inf, 0, 2, Inf), labels = c(0,1,2)),# 0 is absence, 1 is 1-2, 2 is more than 2
)


############# new statistical analysis #################

#Kaplan still seems to show a difference
kaplan_new <- survfit(Surv(tstart, timespan, status)~ hunt, df)
summary(kaplan_new)
autoplot(kaplan_new)

cox_st <- coxph(Surv(tstart,timespan,status) ~ hunt + famgroup_presence + waterdist + roe_dens, data = dfst, id = first_event)
summary(cox_st)
cox.zph(cox_st)


split1 <- survSplit(Surv(timespan,status) ~ hunt + famgroup_category + famgroup_nr + waterdist + roe_dens, data = df,
                    cut = c(50), episode = "timegroup", id = "first_event")
                    