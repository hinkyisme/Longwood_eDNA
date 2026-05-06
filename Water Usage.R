#USGS Data on Water Use in the United States
#Downloaded from  website on March 12, 2024, https://waterdata.usgs.gov/nwis/water_use/
#Data are the total water used within various sectors for 1950-2015, which is the most recent dataset
#Data are collected by the USGS every 5 years, but it can take 3 years before the data are publicly available
#Units are Billions of Gallons of water per day (Bgal/d)

#Dependencies
install.packages("tidyverse")
install.packages("ggpubr"). # only install these if you haven't already.
require(tidyverse)
require(ggpubr)
#Import data file
water <- read_csv("wateruse.csv")

##Review Data 
head(water) # Look at the top
tail(water) # Look at the bottom
View(water) # Look at how R views the data

###############################################################
###############################################################
###US Population and Total Water Withdrawal

#Plot US Population by Year
ggscatter(water, x = "Population", y = "Year", add = "reg.line", ylab = "Population Size (millions)", main = "Population by Year")

#Linear Regression
pop_model <- lm(Population~Year, data= water)
summary(pop_model)

#Place Equation on Plot
ggscatter(water, x = "Population", y = "Year", add = "reg.line", ylab = "Population Size (millions)", main = "Population by Year") + 
  stat_regline_equation(aes(label = after_stat(c(eq.label, adj.rr.label))))

#Plot Total Water Withdrawal vs Human Population size
plot(Population~Year, data=water, ylim=c(100,650), xlim=c(1950, 2020), pch=19, col="black", type="o", ylab="Population Size (millions)")
points(Total_Withdrawal~Year, data=water, ylim=c(0,650), xlim=c(1950,2020), pch=19, col="blue", type="o")
mtext("Total Water Withdrawal (Bgal per day)", side = 4)
legend(1950,630, c("Total Water Withdrawal","US Human Population"), fill=c("black","blue"))

###Calculate Correlation Coefficient
cor(water$Population, water$Total_Withdrawal)

#################################################################################################
#################################################################################################
### Comparing Surface and Groundwater Use

#Plot Total Water Withdrawal in the US from 1950-2010
plot(Total_Withdrawal~Year, data=water, ylim=c(0,650), xlim=c(1950,2020), pch=19, col="blue", type="o", ylab="Total Water Withdrawl (Bgal per day)")

#Add points for surface and groundwater and legend
points(Total_Surface~Year, data=water, type="o", pch=19, col="green")
points(Total_Ground~Year, data=water, type="o", pch=19, col="brown")
legend(1950,630, c("Total","Surface Waters","Groundwater"), fill=c("blue","green","brown"))

##Calculate percent change over time by year
library(tidyverse)
library(dplyr)
water %>%
    mutate(pct_change_total = (Total_Withdrawal - dplyr::lag(Total_Withdrawal))/dplyr::lag(Total_Withdrawal) * 100) %>%
    mutate(pct_change_surface = (Total_Surface - dplyr::lag(Total_Surface))/dplyr::lag(Total_Surface) * 100) %>%
    mutate(pct_change_ground = (Total_Ground - dplyr::lag(Total_Ground))/dplyr::lag(Total_Ground) * 100)


###############################################################################################
###############################################################################################
##Water Use by Sector
#Plot Water Withdrawal for different Sectors
plot(Irrigation~Year, data=water, ylim=c(0,300), xlim=c(1950,2020), pch=19, col="green", type="o", ylab="Water Withdrawal (Bgal per day)")
points(Municipal~Year, data=water, type="o", pch=19, col="blue")
points(Industry~Year, data=water, type="o", pch=19, col="black")
points(Thermoelectric~Year, data=water, type="o", pch=19, col="orange")
points(Livestock~Year, data=water, type="o", pch=19, col="brown")

legend(1950,300, c("Thermoelectric","Irrigation","Industry", "Municipal", "Livestock"), fill=c("orange","green","black","blue","brown"))


##Calculate percent change over time by year
water <- water %>%
  mutate(pct_change_electric = (Thermoelectric - dplyr::lag(Thermoelectric))/dplyr::lag(Thermoelectric) * 100) %>%
  mutate(pct_change_irrigation = (Irrigation - dplyr::lag(Irrigation))/dplyr::lag(Irrigation) * 100) %>%
  mutate(pct_change_industry = (Industry- dplyr::lag(Industry))/dplyr::lag(Industry) * 100) %>%
  mutate(pct_change_municipal = (Municipal- dplyr::lag(Municipal))/dplyr::lag(Municipal) * 100) %>%
  mutate(pct_change_livestock = (Livestock - dplyr::lag(Livestock))/dplyr::lag(Livestock) * 100) 

water

##Add Irrigation and Livestock to make Agriculture Column
water <- water %>%
  mutate(Agriculture = (Irrigation + Livestock))

water

##Add these points to plot
#Plot Water Withdrawal for different Sectors
plot(Irrigation~Year, data=water, ylim=c(0,300), xlim=c(1950,2020), pch=19, col="green", type="o", ylab="Water Withdrawal (Bgal per day)")
points(Municipal~Year, data=water, type="o", pch=19, col="blue")
points(Industry~Year, data=water, type="o", pch=19, col="black")
points(Thermoelectric~Year, data=water, type="o", pch=19, col="orange")
points(Livestock~Year, data=water, type="o", pch=19, col="brown")
points(Agriculture~Year, data=water,type="o", pch=19, col="yellow" )

legend(1950,300, c("Thermoelectric","Irrigation","Industry", "Municipal", "Livestock", "Agriculture"), fill=c("orange","green","black","blue","brown","yellow"))

