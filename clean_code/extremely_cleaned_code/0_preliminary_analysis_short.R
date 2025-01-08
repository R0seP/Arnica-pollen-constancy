##This script analyzes the impact of the weather and sampling date and time 
##(external factors) on the on the number of pollinators caught.

library(lme4)
library(lmerTest)
library(MuMIn)
library(ggplot2)

data <- read.csv("meteo_extended.csv", h = T) #load meteorological observations and  # pollinators per site

#Data preparation----
#sampling time as qualitative:
data$time_quali <- rep(NA, length(data$Visit_id))
data[data$End_time_R < 11, ]$time_quali <- "morning"
data[data$End_time_R >= 11 & data$End_time_R < 13, ]$time_quali <- "early_midday"
data[data$End_time_R >= 13 & data$End_time_R < 15, ]$time_quali <- "late_midday"
data[data$End_time_R >= 15 & data$End_time_R < 17, ]$time_quali <- "afternoon"
data[data$End_time_R >= 17, ]$time_quali <- "evening"


#Model----
m2 <- lmer(Poll_per_hr ~ Date_from_start*time_quali + (1|prec_mm.day), data = data)

summary(m2)

r.squaredGLMM(m2)
#r-squared of model: 0.45


#Effect sizes----
coef <- coef(summary(m2))
coef

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



#Variance partitioning----
anova(m2) #get variance partitioning
total_model <- 16.353 + 118.398 + 97.042
16.353/total_model #7.1% of variance explained by model explained by date from start
118.398/total_model #51.1% of variance explained by model explained by time of day
97.042/total_model #41.9% of variance explained by model explained by interaction


#Model visualization----
ggplot(data,aes(Date_from_start, Poll_per_hr, shape = time_quali)) +
  geom_point() +
  geom_smooth(method='lm', aes(col = time_quali)) +
  theme_minimal() +
  xlab("Day from start of sampling") + ylab("Pollinators caught per hour") +
  labs(shape = "Time of day", col = "Time of day") +
  theme(text = element_text(size = 20))

        