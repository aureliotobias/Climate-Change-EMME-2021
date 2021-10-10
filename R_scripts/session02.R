###########################################################################################################################
###########################################################################################################################
#####  Workshop "Time-series regression for short-term effects of environmental exposures in the EMME".               #####
#####  Aurelio Tobías (aurelio.tobias@idaea.csic.es)                                                                  #####
#####                                                                                                                 #####
#####  2nd International Conference on Climate Change in the Eastern Mediterranean and Middle East.                   ##### 
#####  Cyprus/online, 11th October 2021.                                                                              #####  
#####                                                                                                                 #####
#####  session02.R - Quantifying short-term effects of environmental risk factors.                                    #####
#####  LAST UPDATE: 10/10/2021                                                                                        #####
###########################################################################################################################
###########################################################################################################################

# Remove all previous objects.
rm(list = ls())

###########################################################################################################################
### Load London dataset
###########################################################################################################################

# Load dataset.
data <- read.csv('londondataset.csv', sep=",")

# Formating date.
data$date <- as.Date(data$date)

# Generating time variables
data$time <- seq(nrow(data))
data$year  <- as.factor(substr(data$date, 1, 4))
numyears <- length(unique(data$year))

###########################################################################################################################
### Modeling framework for time-series regression
###########################################################################################################################

library(splines)
library(tsModel)

# Model with natural cubit splits with 6 df/year.
spl <- ns(data$time, df=6*numyears)

model1 <- glm(deaths ~ spl , data, family=quasipoisson)
summary(model1)
pred1 <- predict(model1, type="response")

plot(data$date, data$deaths, ylim=c(100,300), pch=19, cex=0.5, col=grey(0.6),
     ylab="Num. of deaths", xlab="Date")
lines(data$date, pred1, lwd=5, col="blue")

###########################################################################################################################
### Defining crossbasis for environmental exposures 
##########################################################################################################################

install.packages("dlnm")
library(dlnm)
source("findmin.R")

# Crossbasis for temperature.
varknots <- equalknots(data$temp, nk=2)
lagknots <- logknots(21, nk=3)
cbtemp   <- crossbasis(data$temp, lag=21, 
                      argvar=list(fun="bs", degree=2, knots=varknots),
                      arglag=list(knots=lagknots))
summary(cbtemp)

# Crossbasis for ozone.
cbo3 <- crossbasis(data$ozone, lag=7, 
                      argvar=list(fun="lin"),
                      arglag=list(fun="ns", knots=c(1,3)))
summary(cbo3)

###########################################################################################################################
### Estimating temperature-mortality association
###########################################################################################################################

# FIGURE 6. Common environmental exposures.
par(mex=0.8,mfrow=c(2,3))

# Undajusted model.
model1a <- glm(deaths ~ cbtemp, data, family=quasipoisson)
min1a   <- findmin(cbtemp, model1a)
pred1a  <- crosspred(cbtemp, model1a, cen=min1a)
plot(pred1a, "overall", ylim=c(0.9,3), xlim=c(0,30),
     main="Temperature unadjusted", xlab="Temperature (ºC)", ylab="Risk of mortality (RR")

# Time-series regression model.
model1b <- glm(deaths ~ spl + cbtemp, data, family=quasipoisson)
min1b  <- findmin(cbtemp, model1b)
pred1b  <- crosspred(cbtemp, model1b, cen=min1b)
plot(pred1b, "overall", ylim=c(0.9,3), xlim=c(0,30), 
     main="Adjusted for time-trend", xlab="Temperature (ºC)", ylab="Risk of mortality (RR")

# Time-series regression model adjusted for ozone.
model1c <- glm(deaths ~ spl + cbtemp + cbo3, data, family=quasipoisson)
min1c   <- findmin(cbtemp, model1c)
pred1c  <- crosspred(cbtemp, model1c, cen=min1c)
plot(pred1c, "overall", ylim=c(0.9,3), xlim=c(0,30),
     main="Adjusted for time-trend and ozone", xlab="Temperature (ºC)", ylab="Risk of mortality (RR")

###########################################################################################################################
### Estimating ozone-mortality association
###########################################################################################################################

# Undajusted model.
model2a <- glm(deaths ~ cbo3, data, family=quasipoisson)
min2a   <- findmin(cbo3, model2a)
pred2a  <- crosspred(cbo3, model2a, cen=min2a)
plot(pred2a, "overall", ylim=c(0.9,1.3), 
     main="Ozone unadjusted", xlab="Ozone (ug/m3)", ylab="Risk of mortality (RR")

# Time-series regression model.
model2b <- glm(deaths ~ spl + cbo3, data, family=quasipoisson)
min2b   <- findmin(cbo3, model2b)
pred2b  <- crosspred(cbo3, model2b, cen=min2b)
plot(pred2b, "overall", ylim=c(0.9,1.3), 
     main="Adjusted for time-trend", xlab="Ozone (ug/m3)", ylab="Risk of mortality (RR")

# Time-series regression model adjusted for ozone.
model2c <- glm(deaths ~ spl + cbtemp + cbo3, data, family=quasipoisson)
min2c   <- findmin(cbo3, model2c)
pred2c  <- crosspred(cbo3, model2c, cen=min2c)
plot(pred2c, "overall", ylim=c(0.9,1.3), 
     main="Adjusted for time-trend and temperature", xlab="Ozone (ug/m3)", ylab="Risk of mortality (RR")

# END OF FIGURE 6.
layout(1)

###########################################################################################################################
### Estimating health effects of temperature
###########################################################################################################################

# FIGURE 7. Effects of temperature.
par(mex=0.8,mfrow=c(3,2))

temp3d <- plot(pred1c, shade=0.01, col=grey(1),
               xlab="Temperature", ylab="Lag", zlab="RR")
lines(trans3d(x= 2, y=0:21, z=pred1c$matRRfit[as.character( 2),], pmat=temp3d), lwd=3, col="blue")
lines(trans3d(x=22, y=0:21, z=pred1c$matRRfit[as.character(22),], pmat=temp3d), lwd=3, col="red")

temp3d <- plot(pred1c, shade=0.01, col=grey(1),
               xlab="Temperature", ylab="Lag", zlab="RR")
lines(trans3d(x=pred1c$predvar, y=2, z=pred1c$matRRfit[,"lag2"], pmat=temp3d), lwd=3, col="blue")
lines(trans3d(x=pred1c$predvar, y=7, z=pred1c$matRRfit[,"lag7"], pmat=temp3d), lwd=3, col="red")

plot(pred1c, var=2, ylim=c(0.9,1.2), 
     main="Temperature = 2ºC (cold)", ylab="RR", col="blue")

plot(pred1c, lag=2, ylim=c(0.9,1.2), 
     main="lag 2", ylab="RR", xlab="Temperature (ºC)", col="blue")

plot(pred1c, var=22, ylim=c(0.9,1.2),  
     main="Temperature = 22ºC (heat)", ylab="RR", col="red")

plot(pred1c, lag=7, ylim=c(0.9,1.2),  
     main="lag 7", ylab="RR", xlab="Temperature (ºC)", col="red")

# END OF FIGURE 7.
layout(1)

###########################################################################################################################
### Estimating health effects of ozone
###########################################################################################################################

# FIGURE 8. Overall ozone effects.
par(mex=0.8,mfrow=c(2,1))

ozone3d <- plot(pred2c, zlim=c(0.98,1.05), shade=0.01, col=grey(1),
                xlab="Ozone", ylab="Lag", zlab="RR")

plot(pred2c, "overall", ylim=c(0.9,1.3),
     xlab="Ozone (ug/m3)", ylab="RR")

# END OF FIGURE 8.
layout(1)


# FIGURE 9. Ozone effects by lag.
par(mex=0.8,mfrow=c(3,2))

plot(pred2c, lag=0, ylim=c(0.9,1.1),  
     main="lag 0", ylab="RR", xlab="Ozone (ug/m3)", col="red")

plot(pred2c, lag=1, ylim=c(0.9,1.1),  
     main="lag 1", ylab="RR", xlab="Ozone (ug/m3)", col="red")

plot(pred2c, lag=2, ylim=c(0.9,1.1),  
     main="lag 2", ylab="RR", xlab="Ozone (ug/m3)", col="red")

plot(pred2c, lag=3, ylim=c(0.9,1.1),  
     main="lag 3", ylab="RR", xlab="Ozone (ug/m3)", col="red")

plot(pred2c, lag=4, ylim=c(0.9,1.1),  
     main="lag 4", ylab="RR", xlab="Ozone (ug/m3)", col="red")

plot(pred2c, lag=5, ylim=c(0.9,1.1),  
     main="lag 5", ylab="RR", xlab="Ozone (ug/m3)", col="red")

plot(pred2c, lag=6, ylim=c(0.9,1.1),  
     main="lag 6", ylab="RR", xlab="Ozone (ug/m3)", col="red")

plot(pred2c, lag=7, ylim=c(0.9,1.1),  
     main="lag 7", ylab="RR", xlab="Ozone (ug/m3)", col="red")

# END OF FIGURE 9.
layout(1)


# FIGURE 10. Ozone effects by 10 ug/m3
plot(pred2c, var=10, type="p", ci="bars", pch=19, 
     ylim=c(0.99,1.01),  
     ylab="RR", xlab="Lag", col="red")

###########################################################################################################################
### Estimating health effects of desert dust
###########################################################################################################################

# Check the materials:
# Pre-conference workshop on modeling desert dust exposure events for epidemiological short-term health effects studies 
# 31st Annual Conference of the International Society for Environmental Epidemiology 2019.
# https://github.com/aureliotobias/dust/


###########################################################################################################################
###########################################################################################################################
#####                                                                                                                 #####
#####                                             END OF SCRIPT                                                       #####
#####                                                                                                                 #####
###########################################################################################################################
