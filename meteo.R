##This script analyzes the impact of the weather and sampling date and time 
##(external factors that had nothing to do with effects investigated) on the
##sampling dates on thenumber of pollinators caught.

setwd("C:/Users/sohe1/Documents/Master General Biology/Master_Thesis/R")
 
data <- read.csv("meteo_extended.csv", h = T) #load meteorological observations and  # pollinators per site

#Visualize potential weather influence
par(mfrow = c(2,2))
plot(data$mean_air_temp_C,data$Poll_per_hr, main = "Temperature", 
     ylab = "Nr pollinators per hour", xlab = "Temperature (°C) during visit",
     col = "orange")
plot(data$mean_cloud_cover_daily_av,data$Poll_per_hr, main = "Cloud cover", 
     ylab = "Nr pollinators per hour", xlab = "Daily average cloud cover (%)",
     col = "darkgray")
plot(data$prec_mm.day,data$Poll_per_hr, main = "Precipitation", 
     ylab = "Nr pollinators per hour", xlab = "Sum precipitation over day (mm)",
     col = "darkblue")
plot(data$max_wind_m.s,data$Poll_per_hr, main = "Wind speed", 
     ylab = "Nr pollinators per hour", xlab = "Average max wind speed 
     per 10 min over day (m/s)", col = "lightblue")


#Visualize influence of date and time of day 
par(mfrow = c(1,2))
plot(data$Date_from_start,data$Poll_per_hr, main = "Date", 
     ylab = "Nr pollinators per hour", xlab = "Date",
     col = "darkgreen")
plot(data$Start_time_R,data$Poll_per_hr, main = "Time of day", 
     ylab = "Nr pollinators per hour", xlab = "Start time of sampling",
     col = "lightgreen")
par(mfrow = c(1,1))

#model selection
m0 <- lm(Poll_per_hr ~ 1, data = data)
m1 <- lm(Poll_per_hr ~ mean_air_temp_C*mean_cloud_cover_daily_av*
           prec_mm.day*max_wind_m.s*Date_from_start*Start_time_R,
         data = data )
m2 <- lm(Poll_per_hr ~ mean_air_temp_C*mean_cloud_cover_daily_av*
           prec_mm.day*max_wind_m.s+Date_from_start+Start_time_R,
         data = data)
m3 <- lm(Poll_per_hr ~ mean_air_temp_C+mean_cloud_cover_daily_av+
           prec_mm.day+max_wind_m.s+Date_from_start+Start_time_R,
         data = data)
m4 <- lm(Poll_per_hr ~ prec_mm.day*max_wind_m.s+Date_from_start+Start_time_R,
         data = data)
m5 <- lm(Poll_per_hr ~ prec_mm.day+max_wind_m.s+Date_from_start+Start_time_R,
         data = data)
m6 <- lm(Poll_per_hr ~ prec_mm.day+Date_from_start+Start_time_R,
         data = data)
m7 <- lm(Poll_per_hr ~ Date_from_start+Start_time_R,
         data = data)
m8 <- lm(Poll_per_hr ~ prec_mm.day*max_wind_m.s+Date_from_start,
         data = data)
m9 <- lm(Poll_per_hr ~ prec_mm.day+max_wind_m.s+Date_from_start,
         data = data)
m10 <- lm(Poll_per_hr ~ prec_mm.day*max_wind_m.s,data = data)
m11 <- lm(Poll_per_hr ~ prec_mm.day+max_wind_m.s,data = data)
m12 <- lm(Poll_per_hr ~ prec_mm.day,data = data)
m13 <- lm(Poll_per_hr ~ Date_from_start,data = data)

mlist = list(m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13)
AICTab = AIC(m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13) 
AICTab$logLik = unlist(lapply(mlist, logLik)) 
AICTab = AICTab[order(AICTab$AIC, decreasing=F),]
AICTab$delta = round(AICTab$AIC - min(AICTab$AIC), 2)
lh = exp(-0.5*AICTab$delta)
AICTab$w = round(lh/sum(lh), 2)
AICTab

#model with all interaction ranked highest: 
summary(m1)
#model I was most interested in (ranked 7th): 
summary(m4)
#weather does not seem to have an influence, but time of the day and sampling date does!!

#model with only weather data (ranked very low), so definitely NOT good explanation):
summary(m10)
#model with only date and time of day (ranked second):
summary(m4)


#sampling time as qualitative:
data$time_quali <- rep(NA, length(data$Visit_id))
data[data$End_time_R < 11, ]$time_quali <- "morning"
data[data$End_time_R >= 11 & data$End_time_R < 13, ]$time_quali <- "early_midday"
data[data$End_time_R >= 13 & data$End_time_R < 15, ]$time_quali <- "late_midday"
data[data$End_time_R >= 15 & data$End_time_R < 17, ]$time_quali <- "afternoon"
data[data$End_time_R >= 17, ]$time_quali <- "evening"

#new model with date from start and time as qualitative factor:
m_new <- lm(Poll_per_hr ~ Date_from_start*time_quali,
            data = data)
summary(m_new)
#without the interaction:
m_new2 <- lm(Poll_per_hr ~ Date_from_start+time_quali,
             data = data)
summary(m_new2)

#model selection on new models:
mlist = list(m0, m_new, m_new2)
AICTab = AIC(m0, m_new, m_new2) 
AICTab$logLik = unlist(lapply(mlist, logLik)) 
AICTab = AICTab[order(AICTab$AIC, decreasing=F),]
AICTab$delta = round(AICTab$AIC - min(AICTab$AIC), 2)
lh = exp(-0.5*AICTab$delta)
AICTab$w = round(lh/sum(lh), 2)
AICTab

#model with interaction ranked higher
#most pollinators found in the morning (significantly), followed by the afternoon,
#then the late, then the early midday, then the evening; this seems consistent 
#with literature

###m_all <- lm(Poll_per_hr ~ mean_air_temp_C*mean_cloud_cover_daily_av*
###              prec_mm.day*max_wind_m.s*Date_from_start*time_quali,
###         data = data )
###
###mlist = list(m0, m_new, m_new2, m_all)
####AICTab = AIC(m0, m_new, m_new2, m_all) 
###AICTab$logLik = unlist(lapply(mlist, logLik)) 
###AICTab = AICTab[order(AICTab$AIC, decreasing=F),]
###AICTab$delta = round(AICTab$AIC - min(AICTab$AIC), 2)
###lh = exp(-0.5*AICTab$delta)
###AICTab$w = round(lh/sum(lh), 2)
###AICTab

library(lme4)
m_all <- lmer(Poll_per_hr ~ Date_from_start*time_quali + (1|mean_air_temp_C) +
                 (1|mean_cloud_cover_daily_av) + (1|prec_mm.day) +
                 (1|max_wind_m.s), data = data)

#model selection overall: 
mlist = list(m0, m_new, m_new2, m_all)
AICTab = AIC(m0, m_new, m_new2, m_all) 
AICTab$logLik = unlist(lapply(mlist, logLik)) 
AICTab = AICTab[order(AICTab$AIC, decreasing=F),]
AICTab$delta = round(AICTab$AIC - min(AICTab$AIC), 2)
lh = exp(-0.5*AICTab$delta)
AICTab$w = round(lh/sum(lh), 2)
AICTab
