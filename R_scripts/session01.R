###########################################################################################################################
###########################################################################################################################
#####  Workshop "Time-series regression for short-term effects of environmental exposures in the EMME".               #####
#####  Aurelio Tobías (aurelio.tobias@idaea.csic.es)                                                                  #####
#####                                                                                                                 #####
#####  2nd International Conference on Climate Change in the Eastern Mediterranean and Middle East.                   ##### 
#####  Cyprus/online, 11th October 2021.                                                                              #####  
#####                                                                                                                 #####
#####  session01.R - Basic concepts for the time-series design.                                                       #####
#####  LAST UPDATE: 09/10/2021                                                                                        #####
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

###########################################################################################################################
### Time-series data
###########################################################################################################################

# List first 10 rows of time-series data.
head(data)

# FIGURE 1. Time series plots.
par(mex=0.8,mfrow=c(3,1))

# Deaths. 
plot(data$date, data$deaths, 
    type="l", col= "blue",
    main="London, Jan 2002-Dec 2006", 
    ylab="Num. of deaths", xlab="Date")

# Temperature.
plot(data$date, data$temp, 
     type="l", col= "red",
     ylab="Temperature (ºC)", xlab="Date")

# Ozone.
plot(data$date, data$ozone, 
     type="l", col= "springgreen4",
     ylab="Ozone (ug/m3)", xlab="Date")

# END OF FIGURE 1.
layout(1)


###########################################################################################################################
### Modeling framework for time-series regression
###########################################################################################################################

# Time stratified model
#######################

# Generate month and year.
data$month <- as.factor(months(data$date, abbr=TRUE))
data$year  <- as.factor(substr(data$date, 1, 4))

# FIGURE 2. Fit for time stratified models.
par(mex=0.8,mfrow=c(3,1))

# Model with year.
model1a <- glm(deaths ~ factor(year), data, family=quasipoisson)
summary(model1a)
pred1a <- predict(model1a, type="response")

plot(data$date, data$deaths, ylim=c(100,300), pch=19, cex=0.5, col=grey(0.6),
    ylab="Num. of deaths", xlab="Date")
lines(data$date, pred1a, lwd=5, col="blue")

# Model with year and month.
model1b <- glm(deaths ~ factor(year)+month, data, family=quasipoisson)
summary(model1b)
pred1b <- predict(model1b, type="response")

plot(data$date, data$deaths, ylim=c(100,300), pch=19, cex=0.5, col=grey(0.6),
     ylab="Num. of deaths", xlab="Date")
lines(data$date, pred1b, lwd=5, col="blue")

# Model with month nested in each year.
model1c <- glm(deaths ~ year/month, data, family=quasipoisson)
summary(model1c)
pred1c <- predict(model1c, type="response")

plot(data$date, data$deaths, ylim=c(100,300), pch=19, cex=0.5, col=grey(0.6),
     ylab="Num. of deaths", xlab="Date")
lines(data$date, pred1c, lwd=5, col="blue")

# END OF FIGURE 2.
layout(1)


# Periodic functions (Fourier terms)
####################################

install.packages("tsModel")
library(tsModel)
data$time <- seq(nrow(data))

# FIGURE 3. Fit for periodic functions models.
par(mex=0.8,mfrow=c(3,1))

# Model with one sine-cosine pairs.
fourier <- harmonic(data$time, nfreq=1, period=365.25)

model2a <- glm(deaths ~ fourier + time, data, family=quasipoisson)
summary(model2a)
pred2a <- predict(model2a, type="response")

plot(data$date, data$deaths, ylim=c(100,300), pch=19, cex=0.5, col=grey(0.6),
     ylab="Num. of deaths", xlab="Date")
lines(data$date, pred2a, lwd=5, col="blue")

# Model with two sine-cosine pairs.
fourier <- harmonic(data$time, nfreq=2, period=365.25)

model2b <- glm(deaths ~ fourier + time, data, family=quasipoisson)
summary(model2b)
pred2b <- predict(model2b, type="response")

plot(data$date, data$deaths, ylim=c(100,300), pch=19, cex=0.5, col=grey(0.6),
     ylab="Num. of deaths", xlab="Date")
lines(data$date, pred2b, lwd=5, col="blue")

# Model with four sine-cosine pairs.
fourier <- harmonic(data$time, nfreq=4, period=365.25)

model2c <- glm(deaths ~ fourier + time, data, family=quasipoisson)
summary(model2c)
pred2c <- predict(model2c, type="response")

plot(data$date, data$deaths, ylim=c(100,300), pch=19, cex=0.5, col=grey(0.6),
     ylab="Num. of deaths", xlab="Date")
lines(data$date, pred2c, lwd=5, col="blue")

# END OF FIGURE 3.
layout(1)


# Spline functions
##################

library(splines)

# FIGURE 4. Model fit for spline functions models.
par(mex=0.8,mfrow=c(3,1))
numyears <- length(unique(data$year))

# Model with natural cubit splits with 1 df/year.
spl <- ns(data$time, df=1*numyears)

model3a <- glm(deaths ~ spl , data, family=quasipoisson)
summary(model3a)
pred3a <- predict(model3a, type="response")

plot(data$date, data$deaths, ylim=c(100,300), pch=19, cex=0.5, col=grey(0.6),
     ylab="Num. of deaths", xlab="Date")
lines(data$date, pred3a, lwd=5, col="blue")

# Model with natural cubit splits with 6 df/year.
spl <- ns(data$time, df=6*numyears)

model3b <- glm(deaths ~ spl , data, family=quasipoisson)
summary(model3b)
pred3b <- predict(model3b, type="response")

plot(data$date, data$deaths, ylim=c(100,300), pch=19, cex=0.5, col=grey(0.6),
     ylab="Num. of deaths", xlab="Date")
lines(data$date, pred3b, lwd=5, col="blue")

# Model with natural cubit splits with 12 df/year.
spl <- ns(data$time, df=12*numyears)

model3c <- glm(deaths ~ spl , data, family=quasipoisson)
summary(model3c)
pred3c <- predict(model3c, type="response")

plot(data$date, data$deaths, ylim=c(100,300), pch=19, cex=0.5, col=grey(0.6),
     ylab="Num. of deaths", xlab="Date")
lines(data$date, pred3c, lwd=5, col="blue")

# END OF FIGURE 4.
layout(1)


###########################################################################################################################
### Comparing modeling strategies
###########################################################################################################################

# Load function for Quasi-AIC.
source("qAIC.R")

# FIGURE 5. Autocorrelation functions.
par(mex=0.8,mfrow=c(3,1))

# Time stratified model
#######################

qaic1 <- qAIC(model1b, type="dev")
disp1 <- sqrt(summary(model1b)$dispersion)
res1 <- residuals(model1b, type="response")
acf(res1, na.action=na.omit, 
    main=paste("Time stratified model | ", "qAIC=" , round(qaic1,1) , ", overdispersion=" , round(disp1,2)))

# Periodic functions (Fourier terms)
####################################

qaic2 <- qAIC(model2c, type="dev")
disp2 <- sqrt(summary(model2c)$dispersion)
res2 <- residuals(model2c, type="response")
acf(res2, na.action=na.omit, 
    main=paste("Time stratified model | ", "qAIC=" , round(qaic2,1) , ", overdispersion=" , round(disp2,2)))

# Spline functions
##################

qaic3 <- qAIC(model3b, type="dev")
disp3 <- sqrt(summary(model3b)$dispersion)
res3 <- residuals(model3b, type="response")
acf(res3, na.action=na.omit, 
    main=paste("Time stratified model | ", "qAIC=" , round(qaic3,1) , ", overdispersion=" , round(disp3,2)))

# END OF FIGURE 5.
layout(1)

###########################################################################################################################
###########################################################################################################################
#####                                                                                                                 #####
#####                                             END OF SCRIPT                                                       #####
#####                                                                                                                 #####
###########################################################################################################################
###########################################################################################################################
