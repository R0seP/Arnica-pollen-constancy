##This script analyzes the impact of the weather and sampling date and time 
##(external factors that had nothing to do with effects investigated) on the
##sampling dates on the number of pollinators caught. Second try.

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


#sampling time as qualitative:
data$time_quali <- rep(NA, length(data$Visit_id))
data[data$End_time_R < 11, ]$time_quali <- "morning"
data[data$End_time_R >= 11 & data$End_time_R < 13, ]$time_quali <- "early_midday"
data[data$End_time_R >= 13 & data$End_time_R < 15, ]$time_quali <- "late_midday"
data[data$End_time_R >= 15 & data$End_time_R < 17, ]$time_quali <- "afternoon"
data[data$End_time_R >= 17, ]$time_quali <- "evening"

#sample sizes: 
table(data$time_quali, useNA = 'always')


library(lme4)
library(lmerTest)
#model selection 
m0 <- lm(Poll_per_hr ~ 1, data = data)
m1 <- lmer(Poll_per_hr ~ Date_from_start*time_quali + (1|prec_mm.day) + 
             (1|max_wind_m.s), data = data)
m2 <- lmer(Poll_per_hr ~ Date_from_start*time_quali + (1|prec_mm.day), data = data)
m3 <- lm(Poll_per_hr ~ Date_from_start*time_quali, data = data)
m4 <- lm(Poll_per_hr ~ Date_from_start + time_quali, data = data)
m5 <- lm(Poll_per_hr ~ Date_from_start, data = data)
m6 <- lm(Poll_per_hr ~ time_quali, data = data)
m7 <- lmer(Poll_per_hr ~ prec_mm.day*max_wind_m.s+(1|Date_from_start)+
              (1|time_quali), data = data)
m8 <- lmer(Poll_per_hr ~ max_wind_m.s+(1|Date_from_start)+
              (1|time_quali), data = data)
m9 <- lmer(Poll_per_hr ~ prec_mm.day+(1|Date_from_start)+
              (1|time_quali), data = data)
m10 <- lmer(Poll_per_hr ~ Date_from_start*time_quali + (1|prec_mm.day) +
              (1|max_wind_m.s) + (1|mean_air_temp_C) + 
              (1|mean_cloud_cover_daily_av), data = data)

mlist = list(m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10)
AICTab = AIC(m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10) 
AICTab$logLik = unlist(lapply(mlist, logLik)) 
AICTab = AICTab[order(AICTab$AIC, decreasing=F),]
AICTab$delta = round(AICTab$AIC - min(AICTab$AIC), 2)
lh = exp(-0.5*AICTab$delta)
AICTab$w = round(lh/sum(lh), 2)
AICTab

summary(m2)
##most pollinators found in the morning (significantly), followed by the afternoon,
#then the late, then the early midday, then the evening; this seems consistent 
#with literature
##the further away from start date the more pollinators caught, probably because 
#of better catching ability of me but maybe also because of more flowers flowering
##interactions: the later in sampling period the fewer insects caught in morning in 
#comparison to start of sampling, more at all other times; maybe because of there 
#being more pollinators always and reaching the catching limit in morning, or maybe
#because of me getting better at catching infrequent pollinators

coef <- coef(summary(m2))
coef
#characterize effect size
##morning:
###day 1, average -se
mo1minse <- (coef[1,1]-coef[1,2])+(coef[2,1]-coef[2,2])*1+(coef[6,1]-coef[6,2])+(coef[10,1]-coef[10,2])*1
###day 1, av
mo1av <- coef[1,1]+coef[2,1]*1+coef[6,1]+(coef[10,1])*1
###day 1, av +se
mo1pluse <- (coef[1,1]+coef[1,2])+(coef[2,1]+coef[2,2])*1+(coef[6,1]+coef[6,2])+(coef[10,1]+coef[10,2])*1
###day 24, average -sd
mo24minse <- (coef[1,1]-coef[1,2])+(coef[2,1]-coef[2,2])*24+(coef[6,1]-coef[6,2])+(coef[10,1]-coef[10,2])*24
###day 24, av
mo24av <- coef[1,1]+coef[2,1]*24+coef[6,1]+(coef[10,1])*24
###day 24, av +sd
mo24pluse <- (coef[1,1]+coef[1,2])+(coef[2,1]+coef[2,2])*24+(coef[6,1]+coef[6,2])+(coef[10,1]+coef[10,2])*24
##prediction from 100 (66 - 135) pollinators/hr on day 1 to 11 (-58 - 81) pollinators/hr on day 24
##unrealistic! No extrapolation possible, data only for late period (therefore 1 more realistic)


##early mid-day:
###day 1, average -sd
em1minse <- (coef[1,1]-coef[1,2])+(coef[2,1]-coef[2,2])*1+(coef[3,1]-coef[3,2])+(coef[7,1]-coef[7,2])*1
###day 1, av
em1av <- coef[1,1]+coef[2,1]*1+coef[3,1]+(coef[7,1])*1
###day 1, av +sd
em1pluse <- (coef[1,1]+coef[1,2])+(coef[2,1]+coef[2,2])*1+(coef[3,1]+coef[3,2])+(coef[7,1]+coef[7,2])*1
###day 24, average -sd
em24minse <- (coef[1,1]-coef[1,2])+(coef[2,1])*24-coef[2,2]+(coef[3,1]-coef[3,2])+(coef[7,1])*24-coef[7,2]
###day 24, av
em24av <- coef[1,1]+coef[2,1]*24+coef[3,1]+(coef[7,1])*24
###day 24, av +sd
em24pluse <- (coef[1,1]+coef[1,2])+(coef[2,1]+coef[2,2])*24+(coef[3,1]+coef[3,2])+(coef[7,1]+coef[7,2])*24
##prediction: from 9 (4-14) pollinators/hr on day  to 15 (2 - 26) pollinators/hr on day 24
##might be realistic, but sd very high

##late mid-day:
###day 1, average -sd
lm1minse <- (coef[1,1]-coef[1,2])+(coef[2,1]-coef[2,2])*1+(coef[5,1]-coef[5,2])+(coef[9,1]-coef[9,2])*1
###day 1, av
lm1av <- coef[1,1]+coef[2,1]*1+coef[5,1]+(coef[9,1])*1
###day 1, av +sd
lm1pluse <- (coef[1,1]+coef[1,2])+(coef[2,1]+coef[2,2])*1+(coef[5,1]+coef[5,2])+(coef[9,1]+coef[9,2])*1
###day 24, average -sd
lm24minse <- (coef[1,1]-coef[1,2])+(coef[2,1]-coef[2,2])*24+(coef[5,1]-coef[5,2])+(coef[9,1]-coef[9,2])*24
###day 24, av
lm24av <- coef[1,1]+coef[2,1]*24+coef[5,1]+(coef[9,1])*24
###day 24, av +sd
lm24pluse <- (coef[1,1]+coef[1,2])+(coef[2,1]+coef[2,2])*24+(coef[5,1]+coef[5,2])+(coef[9,1]+coef[9,2])*24
##prediction: from 5 (-1 - 11) pollinators/hr on day 1 to 11 (-3 - 25) pollinators/hr on day 24

##afternoon:
###day 1, average -sd
af1minse <- (coef[1,1]-coef[1,2])+(coef[2,1]-coef[2,2])*1
###day 1, av
af1av <- coef[1,1]+(coef[2,1])*1
###day 1, av +sd
af1pluse <- (coef[1,1]+coef[1,2])+(coef[2,1]+coef[2,2])*1
###day 24, average -sd
af24minse <- (coef[1,1]-coef[1,2])+(coef[2,1]-coef[2,2])*24
###day 24, av
af24av <- coef[1,1]+(coef[2,1])*24
###day 24, av +sd
af24pluse <- (coef[1,1]+coef[1,2])+(coef[2,1]+coef[2,2])*24
##prediction: from 10 (8 - 12) pollinators/hr on day 1 to 13 (8 - 18) pollinators/hr on day 24

##evening:
###day 1, average -sd
ev1minse <- (coef[1,1]-coef[1,2])+(coef[2,1]-coef[2,2])*1+(coef[4,1]-coef[4,2])+(coef[8,1]-coef[8,2])*1
###day 1, av
ev1av <- coef[1,1]+coef[2,1]*1+coef[4,1]+(coef[8,1])*1
###day 1, av +sd
ev1pluse <- (coef[1,1]+coef[1,2])+(coef[2,1]+coef[2,2])*1+(coef[4,1]+coef[4,2])+(coef[8,1]+coef[8,2])*1
###day 24, average -sd
ev24minse <- (coef[1,1]-coef[1,2])+(coef[2,1]-coef[2,2])*24+(coef[4,1]-coef[4,2])+(coef[8,1]-coef[8,2])*24
###day 24, av
ev24av <- coef[1,1]+coef[2,1]*24+coef[4,1]+(coef[8,1])*24
###day 24, av +sd
ev24pluse <- (coef[1,1]+coef[1,2])+(coef[2,1]+coef[2,2])*24+(coef[4,1]+coef[4,2])+(coef[8,1]+coef[8,2])*24
##prediction: from -8 (-27 - 11) pollinators/hr on day 1 to 16 (-26 - 59) pollinators/hr on day 24
## unrealistic because of too small sample size and data only for short time span!

mean(c(mo1av, mo24av)) #mean poll caught in morning: 56
mean(c(em1av, em24av)) #mean poll caught in early midday: 12
mean(c(lm1av, lm24av)) #mean poll caught in late midday: 8
mean(c(af1av, af24av)) #mean poll caught in afternoon: 12
mean(c(ev1av, ev24av)) #mean poll caught in evening: 4


library(MuMIn)
r.squaredGLMM(m2)
#r-squared of model: 0.45

#variancMuMIn#variance partitioning:
anova(m2) #get variance partitioning
total_model <- 16.353 + 118.398 + 97.042
16.353/total_model #7.1% of variance explained by model explained by date from start
118.398/total_model #51.1% of variance explained by model explained by time of day
97.042/total_model #41.9% of variance explained by model explained by interaction


#visualize the data and the model
library(ggplot2)

ggplot(data,aes(Date_from_start, Poll_per_hr, shape = time_quali)) +
  geom_point() +
  geom_smooth(method='lm', aes(col = time_quali)) +
  theme_minimal() +
  xlab("Day from start of sampling") + ylab("Pollinators caught per hour") +
  labs(shape = "Time of day", col = "Time of day") +
  theme(text = element_text(size = 20))

        