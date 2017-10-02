# September 10, 2017
# Lab 1, Pollution example
# Here we will learn how to read data and form
# the linear regression model
# Read data from cvs file with ";" instead of ","
pollutiondata <- read.csv2("AirQualityData.csv",header=TRUE)

# Replace missing values (-200) with Not A number (NA)
pollutiondata[pollutiondata==-200]<-NA

# Remove rows with NA
# Note: Don't remove entire row in practice! Never remove valid data samples!
pollutiondatareduced<-pollutiondata[complete.cases(pollutiondata), ]

# Print the first 10 data rows
head(pollutiondatareduced,10)

# Plot Benzene pollutant vs time
#time is a factor - done hourly
par(mfrow=c(1,2))
plot(pollutiondatareduced$Time,pollutiondatareduced$PT08.S2.NMHC.,xlab="Time", ylab="Benzene")
# Notice that we get a box plot in the scatter plot. Why is this ?

# Plot Benzene pollutant vs Temperature

plot(as.numeric(as.character(pollutiondatareduced$Temperature)),pollutiondatareduced$PT08.S2.NMHC.,xlab="Temperature", ylab="Benzene")

# Let us look at Temperature vs Benzene polutant
# Attach reduced data set
attach(pollutiondatareduced)

# Some of the numbers in dataframe for time are loaded as characters, 
# Force them to numeric
Temperature = as.numeric(as.character(Temperature))

# Make scatter plot
par(mfrow=c(1,1))
plot(Temperature,pollutiondatareduced$PT08.S2.NMHC.,xlab="Temperature", ylab="Benzene")

# Perform linear regression between temperature and Benzene
m1 <- lm(PT08.S2.NMHC. ~ Temperature)
summary(m1)

# Plot figure with least square fit
plot(Temperature,PT08.S2.NMHC.,xlab="Temperature", ylab="Benzene")
abline(lsfit(Temperature,PT08.S2.NMHC.))

#t-value for 95% confidence
tval <- qt(1-0.05/2,825)
tval

#95% confidence intervals
round(confint(m1,level=0.95),3)

# Confidence and prection values for 95%
predict(m1,newdata=data.frame(Temperature=c(5,10,15,20,25,30)),interval="confidence",level=0.95)
predict(m1,newdata=data.frame(Temperature=c(5,10,15,20,25,30)),interval="prediction",level=0.95)
