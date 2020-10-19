##Load in Required Packages
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(MASS)
library(AER)
library(geosphere)
library(boot)
library(gridExtra)
library(ggfortify)


##Read in Mark Release Recapture Data Set
MRR_DATA = read_csv("./MRR_DATA.csv")

##Tidy Mark Release Recapture Data Set
##Blue mosquitoes removed from analysis
##Names of Columns Simplified and spaces removed
##New Column for Flight_Distance (Release Point to BG Trap)
##New Column for Flight_Direction (Release Point to BG Trap)

MRR_Tidy = MRR_DATA%>%
  filter(Color != "Blue")%>%
  rename(N.Released = `N released`)%>%
  rename(N.Collected = `N collected`)%>%
  rename(Wolbachia = wb)%>%
  rename(BG_Location = `BG location`)%>%
  rename(Human_Density = `People (N)`)%>%
  rename(Premise_Type = `Premise Type`)%>%
  rename(Wind_Intensity = `Wind intensity (m/s)`)%>%
  rename(Wind_Direction = `Wind direction`)%>%
  rename(Av_Max_Temp = `Average MAX Temp`)%>%
  rename(Av_Min_Temp = `Average MIN Temp`)%>%
  rename(Av_Mean_Temp = `Average MEAN temp`)%>%
  rename(Temp_Max_Absolute = `Temp MAX absolut`)%>%
  rename(Temp_Min_Absolute = `Temp MIN absolut`)%>%
  rename(Absolute_Max_RH = `Average MAX UR`)%>%
  rename(Average_Min_RH = `Average MIN UR`)%>%
  rename(Average_Med_RH = `Average MED UR`)%>%
  rename(Rainfall = `Accumulated rain (mm)`)%>%
  rename(Wild = wild)%>%
  mutate(Wild_and_Wolbachia = Wolbachia + Wild)%>%
  rename(Latitude = lat)%>%
  rename(Longitude = long)%>%
  mutate(BG = as.factor(BG))%>%
  mutate(BG_Location = as.factor(BG_Location))%>%
  mutate(Premise_Type = as.factor(Premise_Type))%>%
  mutate(Color = as.factor(Color))%>%
  group_by(Color, BG)%>%
  mutate(Dispersal_Distance = (mean(Dist, na.rm = TRUE)))%>%  
  mutate(Dispersal_Direction = (mean(Direction, na.rm = TRUE)))
 
##Filter data to have only datapoints which have dist, lat, long and direction
Capture_Locations_Dist_Dir = MRR_Tidy%>%
  filter(Dist > 0)


##make Coordinates data frame to hold loop
Coordinates = data.frame(1,2)

##run loop of destPointRhumb (assumes mercator projection - fine as over short distances)
for(i in 1:119){
  Coordinates[i,] = (destPointRhumb(c(Capture_Locations_Dist_Dir$Longitude[i], 	Capture_Locations_Dist_Dir$Latitude[i]),
                                    b =Capture_Locations_Dist_Dir$Direction [i], d=Capture_Locations_Dist_Dir$Dist[i]))
}

##Rename column names in Coordinates dataset
Coordinates = Coordinates%>%
  rename(Release_Longitude = X1)%>%
  rename(Release_Latitude = X2)
  
##Bind the release Coordinates DF with the MRR_Release DF
MRR_Release = bind_cols(Capture_Locations_Dist_Dir, Coordinates)


##Confirm accuracy of calculations through back calculating. For supplementary information
##Create temp vectors for loops
Distance_Calculated = c()
Direction_Calculated = c()

##run loop for distRhumb and bearingRhumb (assumes mercator projection)
##to get distance and directions between BG and release site (allows checking that
##initial calculation work)
for(j in 1:119){
  Distance_Calculated[j] = distRhumb(c(MRR_Release$Longitude[j], MRR_Release$Latitude[j]), 
                                     c(MRR_Release$Release_Longitude[j], MRR_Release$Release_Latitude[j]))
 Direction_Calculated[j] = bearingRhumb(c(MRR_Release$Longitude[j], MRR_Release$Latitude[j]), 
                                     c(MRR_Release$Release_Longitude[j], MRR_Release$Release_Latitude[j]))
}


##convert to dataframe
Distance_Calculated = as.data.frame(Distance_Calculated)
Direction_Calculated = as.data.frame(Direction_Calculated)

##add back calculated distance and direction to release info DF
All_Release_Info = bind_cols(MRR_Release, Distance_Calculated, Direction_Calculated)

##Plot to confirm reliability of Distance back calculation
ggplot(All_Release_Info, aes(x=Dist, y=Distance_Calculated)) +
  geom_point()


##Plot to confirm reliability of Direction back calculation
ggplot(All_Release_Info, aes(x=Direction, y=Direction_Calculated)) +
         geom_point()

##Both plots have straight line confirming transformations worked. 

###Average Central Release point of each release point by Color and Week.
All_Release_Info_A = All_Release_Info%>%
  group_by(Week, Color)%>%
  mutate(Central_Release_Lat = mean(Release_Latitude))%>%
  mutate(Central_Release_Long = mean(Release_Longitude))


##Add calculated Central release points to main dataset; NAs introduced
All_Data = All_Release_Info_A%>%
  full_join(MRR_Tidy)

##Calculate centre of weekly release points
All_Data2 = All_Data%>%
  group_by(Week, Color)%>%
  mutate(Av_Release_Lat = mean(Central_Release_Lat, na.rm = TRUE))%>%##Calculate centre of weekly release points
  mutate(Av_Release_Long = mean(Central_Release_Long, na.rm = TRUE))%>%##Calculate centre of weekly release points
  ungroup()%>%
  group_by(Color)%>%
  mutate(Aver_Release_Lat = mean(Av_Release_Lat, na.rm = TRUE))%>%##Calculate average of study release points
  mutate(Aver_Release_Long = mean(Av_Release_Long, na.rm = TRUE))##Calculate average of study release points

##Remove NAs and rename release lat and long columns
NAN_Av_Release = All_Data2%>%
  filter(is.na(Av_Release_Lat))%>%
  rename(Average_Release_Lat = "Aver_Release_Lat")%>%
  rename(Average_Release_Long = "Aver_Release_Long")%>%
  dplyr::select(-c(Av_Release_Lat, Av_Release_Long))

##Remove NAs and rename release lat and long columns
Not_NAN_Av_Release = All_Data2%>%
  filter(!is.na(Av_Release_Lat))%>%
  rename(Average_Release_Lat = "Av_Release_Lat")%>%
  rename(Average_Release_Long = "Av_Release_Long")%>%
  dplyr::select(-c(Aver_Release_Lat, Aver_Release_Long))
  
##Combine data frames
MRR_Tidy_Final = Not_NAN_Av_Release%>%
  rbind(NAN_Av_Release)
    
##Calculate Dist and Direction to release sites from BGs
NAN_MRR_Tidy= MRR_Tidy_Final%>%
  filter(is.na(Dist))%>%
  rowwise( )%>%
  mutate(Direction = bearingRhumb(c(Average_Release_Long, Average_Release_Lat), c(Longitude, Latitude)))%>%
  mutate(Dist = distRhumb(c(Average_Release_Long, Average_Release_Lat), c(Longitude, Latitude)))%>%
  ungroup()


##Filter non NAs
Not_NAN_MRR_Tidy = MRR_Tidy_Final%>%
  filter(!is.na(Dist))

##Join dataframes to give dataframe with disersal dist and direction
MRR_Dataset = union(Not_NAN_MRR_Tidy, NAN_MRR_Tidy)
  

ggplot(MRR_Dataset, aes(x=Dist, y=N.Collected)) +
  geom_point()

##SUMMARY STATISTICS##
##Total Number of Mosquitoes Released
Total_Released = MRR_Dataset%>%
  filter(BG == 1)%>%
  summarise(sum(N.Released))

print(sum(Total_Released$`sum(N.Released)`))

##Total Mosquitoes Released per Colour
MRR_Released = MRR_Dataset%>%
  filter(BG == 1)%>%
  group_by(Color)%>%
  summarise(sum(N.Released))

print(MRR_Released)

##Total Marked Mosquitoes Released per Week
MRR_Dataset%>%
  filter(BG == 1)%>%
  group_by(Week)%>%
  summarise(sum(N.Released))


##Total Mark Mosquitoes Captured per Colour
MRR_Captured = MRR_Dataset%>%
  group_by(Color)%>%
  summarise(sum(N.Collected))

print(MRR_Captured)
##Total number marked mosquitoes captured
Total_Marked = sum(MRR_Captured$`sum(N.Collected)`)
print(Total_Marked)

##Total number of wild Aedes aegypti caught
MRR_Yellow = MRR_Dataset%>%
  filter(Color == "Yellow")

Total_Wild = sum(MRR_Yellow$Wild)
print(Total_Wild)

##Total number of Aedes aegypti with Wolbachia
Total_Wolbachia = sum(MRR_Yellow$Wolbachia)
print(Total_Wolbachia)
#Total number mosquitoes caught
Total_Caught = sum(Total_Wild + Total_Wolbachia + Total_Marked)
print(Total_Caught)

##percentages of marked, wild, wolbachia
Total_Marked/Total_Caught * 100
Total_Wild/Total_Caught * 100
Total_Wolbachia/Total_Caught * 100
 
MRR_Trap_Info = MRR_Yellow%>%
  ungroup()%>%
  filter(Week == 15)

##Recapture Success per Color 
MRR_Color = MRR_Released%>%
  inner_join(MRR_Captured)%>%
  rename(Collected = "sum(N.Collected)")%>%
  rename(Released = "sum(N.Released)")%>%
  mutate(Catch_Percentage = Collected/Released * 100)

print(MRR_Color)

##Overall Recapture Success
Overall_Perc_Caught = sum(MRR_Color$Collected) / sum(MRR_Color$Released) * 100
print(Overall_Perc_Caught)

##Maximum and Minimum trap distances
Range_Traps = range(MRR_Dataset$Dist)

print(Range_Traps)

##Maximum and Minimum actual flight distances
Range_Dispersal = range(No_Zeroes_MRR$Dist)

print(Range_Dispersal)

##percentage max dispersal for trap range
Range_Dispersal[2]/Range_Traps[2] * 100


View(MRR_Dataset)

##Calculate sampling area in km^2 using max and min lat/long of BG sentinels
Width = distRhumb(c(-43.22196, 0), c(-43.23019, 0)) / 1000
Length = distRhumb(c(0, -22.78027), c(0, -22.78870)) / 1000

Total_Area = Width*Length
print(Total_Area)
30/Total_Area
 

##Data Visualisation
Marked_Histogram = ggplot(No_Zeroes_MRR, aes(x=N.Collected, fill = Color)) +
  geom_histogram(binwidth = 1) +
  scale_fill_manual(values= c("green", "grey", "purple", "orange", "pink", "yellow"))+
  xlab("Number of Marked Mosquitoes Captured") +
  ylab("Frequency")+
  theme_bw()+
  theme(legend.position = "none")+
  ggtitle("A")

##Calculate frequency of no marked mosquitoes being caught regardless of colour
nrow(MRR_Dataset%>%
  filter(N.Collected == 0))

##Use single colour to prevent duplication of count values in histograms
MRR_Dataset_Yellow = MRR_Dataset%>%
  filter(Color == "Yellow")

#histogram of Wild mosquito counts
Wild_Histogram = ggplot(MRR_Dataset_Yellow, aes(x=Wild)) +
  geom_histogram(binwidth = 1) +
  theme_bw() +
  xlab("Number of Wild Mosquitoes Captured") +
  ylab("Frequency")+
  ggtitle("B") 

##histogram of Wolbachia counts
Wolbachia_Histogram = ggplot(MRR_Dataset_Yellow, aes(x=Wolbachia))+
  geom_histogram(binwidth = 1) +
  theme_bw() +
  xlab("Number of Wolbachia Mosquitoes Captured") +
  ylab("Frequency")+
  ggtitle("C") 

##Histograms of Count Frequencies
Wild_Wolbachia_Histograms = grid.arrange(Wild_Histogram, Wolbachia_Histogram)

grid.arrange(Marked_Histogram, Wild_Wolbachia_Histograms, ncol = 2)


##Summary of trap information  
summary(MRR_Trap_Info)

##Filter data to only include non zero marked mosquito captures 
No_Zeroes_MRR = MRR_Dataset%>%
  filter(N.Collected > 0)


##trap locations
Trap_Locations = MRR_Tidy_Final%>%
   dplyr::select(Latitude, Longitude, BG_Location, Premise_Type)

colnames(MRR_Dataset)
pal = c("#008d00","#a0db00","#dbac2d","#5590bd",
        "#ff0000" ,"#8c8b8b", "#ffffff", "#000000")

##Plot of traps and number collected:
Number_Wild_Wolbachia_Collected = ggplot(MRR_Dataset_Yellow, aes(x=Longitude, y=Latitude, size = (Wild + Wolbachia), shape = BG_Location, color = Premise_Type)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~Week)

Number_Wild_Wolbachia_Collected

##create plot of movement from release sites to traps
ggplot(No_Zeroes_MRR, aes()) +
  geom_segment(aes(x=Average_Release_Long, y=Average_Release_Lat, ##Release Site Coordinates
                   xend = Longitude, yend= Latitude, ##BG Coordinates
                   size = N.Collected,
                   color = Color)) +
  scale_colour_manual(values= c("green", "grey", "purple", "orange", "pink", "yellow")) +   ##colours match marked mosquitoes colours
  geom_point(aes(x=Average_Release_Long, y=Average_Release_Lat, colour = Color, size = 10))+
  geom_point(data = MRR_Dataset_Yellow, aes(x = Longitude, y= Latitude, size = (Wild + Wolbachia), shape = BG_Location))+
  theme_classic()+ 
  theme(legend.position = "none")+
  xlab("Longitude")+
  ylab("Latitude") +
  facet_wrap(~Week)


##Boxplot of Trap Counts (excluding zeroes to allow scale to be visible)
ggplot(No_Zeroes_MRR, aes(x = BG, y=N.Collected)) +
  geom_boxplot()+
  theme_bw() +
  ylab("Number of Marked Mosquitoes Collected") +
  xlab("BG Sentinel Trap ID Number")
  

##plot of number of marked mosquites captured and distance
ggplot(No_Zeroes_MRR, aes(x=Dist, y=N.Collected, color = Color)) +
  geom_point() +
  scale_colour_manual(values= c("green", "grey", "purple", "orange", "pink", "yellow"))+
  xlab("Dispersal Distance (m)") +
  ylab("Number Collected")



###GRAVITY MODELS##
##Testing for linearity; between variables potentially likely to display colinearity
##Not colinear

##Testing Humidity and Temperature for co-linearity
Humid_Temp_Data = MRR_Tidy_Final%>%
  filter(BG == 1)%>%
  filter(Color == "Yellow")
Humid_Temp_Lm = lm(formula = Average_Med_RH ~ Av_Mean_Temp, data = Humid_Temp_Data)
summary(Humid_Temp_Lm)

Humid_Temp_Scatter = ggplot(Humid_Temp_Data, aes(y = Average_Med_RH, x = Av_Mean_Temp)) +
  geom_point()+
  geom_abline(slope = -3.5470, intercept = 169.6682) +
  xlab("Mean Temperature , °C") +
  ylab("Average Relative Humidity, %") +
  ggtitle("A")

Humid_Temp_Lm_Residuals = residuals(Humid_Temp_Lm)


##Broadly normally distributed
Humid_Temp_Hist = ggplot(, aes(x=Humid_Temp_Lm_Residuals)) +
  geom_histogram(binwidth = 4) +
  theme_bw() +
  xlab("Residuals") +
  ggtitle("B")
shapiro.test(Humid_Temp_Lm_Residuals)
autoplot(Humid_Temp_Lm)

grid.arrange(Humid_Temp_Scatter, Humid_Temp_Hist)  

##broadly normally distributed
##Highly co-linear; only need to input either temperature or humidity

##Rainfall and Humdity
Rainfall_Humid_LM = lm(formula = Average_Med_RH ~ Rainfall, data = Humid_Temp_Data)
summary(Rainfall_Humid_LM)


##N.Collected is transformed using log(N.Collected + 1)
##as zeroes cannot be logged. No impact of release size on capture size
ggplot(MRR_Dataset, aes(x=N.Released, y=N.Collected))+
  geom_point() +
  geom_smooth(method = "lm")

##Gravity_Model1 Generalised Linear Model: Poisson, No additional Covariates (only distance)
Poisson_GM1 = glm(
  formula = log(N.Collected + 1) ~ log(Dist),
                     data = MRR_Dataset,
                    family = "poisson") 
##residuals
Poisson_GM1_Residuals = residuals.glm(Poisson_GM1, "pearson")
summary(Poisson_GM1)


##Test Normality of residuals of Poisson Gravity Model 1
ggplot(,aes(x=Poisson_GM1_Residuals)) +
  geom_histogram(binwidth = 0.1) +
  xlab("Residuals") +
  theme_bw()

##Shapiro-Wilk test for Normality on Pearson Residuals, Poisson Grav Model 1
shapiro.test(Poisson_GM1_Residuals)

##diagnostic plots
autoplot(Poisson_GM1)

#Test for overdispersion
dispersiontest(Poisson_GM1) ##Suggests data is not overdispersed, but INF AIC of Poisson is a problem

##sequence generation to make prediction on
Distance = seq(1, 611, 1)
Colour = c("Green", "Grey", "Lilac", "Orange", "Pink", "Yellow")

predict.poisson.count =
  predict(Poisson_GM1, newdata=data.frame(Dist = Distance),
                                      type="response", se.fit = TRUE, interval = "prediction",
              dispersion = 0.971834)

##Calculate 95% CIs
upperCI = as.numeric(predict.poisson.count$fit + (1.96 * predict.poisson.count$se.fit))
lowerCI = as.numeric(predict.poisson.count$fit - (1.96 * predict.poisson.count$se.fit))
fit_values = as.numeric(predict.poisson.count$fit)

##combine into dataframe
Poisson_Predict = data.frame(Distance, upperCI, lowerCI, fit_values)


##Plot prediction of Poisson GM with confidence intervals and actual values
Poisson_Predict_Graph = ggplot(Poisson_Predict, aes(x=log(Distance), y=log(fit_values + 1))) +
  geom_line(colour = "green") +
  geom_line(colour = "red", aes(y=log(upperCI +1), x= log(Distance))) +
  geom_line(colour = "red", aes(y=log(lowerCI + 1), x = log(Distance))) +
  geom_point(data = MRR_Dataset,aes(x=log(Dist), y= log(N.Collected+1), colour = Color)) +
  scale_colour_manual(values= c("green", "grey", "purple", "orange", "pink", "yellow")) +
  xlab("log(Distance (m) )")+
  ylab("log(Number Marked Mosquitoes Collected + 1)") +
  theme_bw()+
  ggtitle("Poisson")+
  theme(legend.position = "none")



##Poisson Gravity Model with additional Explanatory Variables
##Plotting of humidity vs temperature shows high levels of colinearity
##Temperature more likely to play a role in dispersal??

Poisson_GM2 = glm(formula = log(N.Collected + 1) ~ 
                    log(Dist) +
                    (log(Direction) +
                    log(Wind_Intensity) * 
                    log(Wind_Direction) +
                    log(Rainfall+1) + 
                    log(Av_Mean_Temp)) +
                    (log(Wild_and_Wolbachia + 1) +
                    log(Human_Density) +
                    Color + BG_Location + BG),
                    data = MRR_Dataset,
                     family = poisson(link = "log"))
##Summary and diagnostic plots of Poisson gravity model containing all 
##potential diagnostics
summary(Poisson_GM2) ##still INF AIC, do not use

Poisson_GM2_Residuals = residuals.glm(Poisson_GM2, "pearson")

##Poisson Grav Model 2 Diagnostics
dispersiontest(Poisson_GM2)
autoplot(Poisson_GM2)

ggplot(, aes(x=Poisson_GM2_Residuals)) +
  geom_histogram(binwidth = 0.1) +
  xlab("Residuals") +
  theme_bw()

shapiro.test(Poisson_GM2_Residuals)

##Stepwise AIC on Poisson Gravity Model containing all covariates.

step(Poisson_GM2)
##Stepwise AIC not possible due to infinite AIC (overdispersed data)


##Gravity Models with Negative Binomial Distributions

##Gravity Model containing only Distance as a covariate
Neg_Bin_GM1 = glm.nb(formula = log(N.Collected + 1) ~ log(Dist), method = "glm.fit", link = "log",
                     data = MRR_Dataset)

##Summaries and diagnostic plots of residuals of negative binomial gravity
##model containing only distance as a predictor. 
summary(Neg_Bin_GM1)  
autoplot(Neg_Bin_GM1)

Neg_Bin_GM1_Residuals = residuals.glm(Neg_Bin_GM1, "pearson")

ggplot(, aes(x=Neg_Bin_GM1_Residuals)) +
  geom_histogram(binwidth = 0.1) +
  xlab("Residuals") +
  theme_bw() 

shapiro.test(Neg_Bin_GM1_Residuals)


predict.neg.bin.count =
  predict(Neg_Bin_GM1, newdata=data.frame(Dist = Distance),
          type="response", se.fit = TRUE, dispersion = 1.1433)

##Calculate 95% CIs
upperCI_NB = as.numeric(predict.neg.bin.count$fit + (1.96 * predict.neg.bin.count$se.fit))
lowerCI_NB = as.numeric(predict.neg.bin.count$fit - (1.96 * predict.neg.bin.count$se.fit))
fit_values_NB = as.numeric(predict.neg.bin.count$fit)

Neg_Bin_Predict = data.frame(Distance, upperCI_NB, lowerCI_NB, fit_values_NB)

##plot graph of predicts (and CIs) with actual data too
Negative_Binomial_Predict_Graph = ggplot(Neg_Bin_Predict, aes(x=log(Distance), y=log(fit_values_NB + 1))) +
  geom_line(colour = "green") +
  geom_line(colour = "red", aes(y=log(upperCI_NB +1), x= log(Distance))) +
  geom_line(colour = "red", aes(y=log(lowerCI_NB + 1), x = log(Distance))) +
  geom_point(data = MRR_Dataset, aes(x=log(Dist), y=log(N.Collected + 1), colour = Color)) +
  scale_colour_manual(values= c("green", "grey", "purple", "orange", "pink", "yellow")) +
  xlab("log(Distance (m))") +
  ylab("log(Number of Marked Mosquitoes Captured + 1)")+
  ggtitle("Negative Binomial")+
  theme_bw() +
  theme(legend.position = "none")


##Negative Binomial GLm with Multiple Covariates:
Neg_Bin_GM2 = glm.nb(formula = log(N.Collected + 1) ~ 
                        log(Dist) +
                        log(Direction) +
                        log(Wind_Intensity) * 
                        log(Wind_Direction) +
                        log(Rainfall + 1) + 
                        log(Wild_and_Wolbachia + 1) +
                        log(Human_Density) +
                        log(Av_Mean_Temp) +
                        Color + BG_Location + BG, method = "glm.fit", link = "log",
            data = MRR_Dataset)

##Summaries and diagnostic plots of residuals for Negative Binomial gravity model
##containing all potential covariates.
summary(Neg_Bin_GM2)

##residuals
Neg_Bin_GM2_Residuals = residuals.glm(Neg_Bin_GM2, "pearson")

autoplot(Neg_Bin_GM2) ##diagnostic plots
ggplot(Neg_Bin_GM2, aes(x=Neg_Bin_GM2_Residuals)) +
  geom_histogram(binwidth = 0.1) +
  xlab("Residuals") +
  theme_bw()
shapiro.test(Neg_Bin_GM2_Residuals) ##on pearson residuals


##Run stepwise AIC on gravity model containing all potential explanatory covariates
StepAIC_Neg_Bin_GM2 = stepAIC(Neg_Bin_GM2)

##Take final stepped output of Step(Neg_Bin_GM2)
Neg_Bin_GM3 = glm.nb(log(N.Collected + 1) ~ 
                       log(Dist) + 
                       log(Direction) + 
                       log(Wind_Intensity) + 
                       log(Wild_and_Wolbachia + 1) + 
                       Color,
                     method = "glm.fit", link = "log",
                     data = MRR_Dataset)


###Confirm stepwise is completed.
stepAIC(Neg_Bin_GM3) ##no additional covariates eliminated
summary(Neg_Bin_GM3)
autoplot(Neg_Bin_GM3) ##diagnostic plots

Neg_Bin_GM3_Residuals = residuals(Neg_Bin_GM3, "pearson")

ggplot(, aes(x=Neg_Bin_GM3_Residuals)) +
  geom_histogram(binwidth = 0.1) +
  xlab("Residuals") +
  theme_bw()

shapiro.test(Neg_Bin_GM3_Residuals)

##Test to see if using flight direction and wind direction as categorical variables improve
##or affect model.

Wind_NE = MRR_Dataset%>%
  filter(Wind_Direction <= 90)%>%
  mutate(Wind_Direction_NESW = "NE")

Wind_SE = MRR_Dataset%>%
  filter(Wind_Direction > 90 & Wind_Direction <= 180)%>%
  mutate(Wind_Direction_NESW = "SE")
  
Wind_SW = MRR_Dataset%>%
  filter(Wind_Direction > 180 & Wind_Direction <= 270)%>%
  mutate(Wind_Direction_NESW = "SW")

Wind_NW = MRR_Dataset%>%
  filter(Wind_Direction > 270 & Wind_Direction <= 360)%>%
  mutate(Wind_Direction_NESW = "NW")


Wind_NESW_Data = rbind(Wind_NE, Wind_SE, Wind_SW, Wind_NW)

##same for dispersal direction
Flight_NE = MRR_Dataset%>%
  filter(Direction <= 90)%>%
  mutate(Flight_Direction_NESW = "NE")

Flight_SE = MRR_Dataset%>%
  filter(Direction > 90 & Direction <= 180)%>%
  mutate(Flight_Direction_NESW = "SE")

Flight_SW = MRR_Dataset%>%
  filter(Direction > 180 & Direction <= 270)%>%
  mutate(Flight_Direction_NESW = "SW")

Flight_NW = MRR_Dataset%>%
  filter(Direction > 270 & Direction <= 360)%>%
  mutate(Flight_Direction_NESW = "NW")


Flight_NESW_Data = rbind(Flight_NE, Flight_SE, Flight_SW, Flight_NW)


MRR_NESW = inner_join(Wind_NESW_Data, Flight_NESW_Data)

##GLM with categorical direction variables
Neg_Bin_NESW = glm.nb(formula = log(N.Collected + 1) ~ 
                       log(Dist) +
                       log(Rainfall + 1) + 
                       log(Wild_and_Wolbachia + 1) +
                       log(Human_Density) +
                       log(Av_Mean_Temp) +
                       Color + BG_Location + 
                        log(Wind_Intensity)*
                        Wind_Direction_NESW +
                        Flight_Direction_NESW, method = "glm.fit", link = "log",
                     data = MRR_NESW)



summary(Neg_Bin_NESW)
Neg_Bin_NESW_Residuals = residuals.glm(Neg_Bin_NESW, "pearson")

ggplot(, aes(x=Neg_Bin_NESW_Residuals)) +
  geom_histogram(binwidth = 0.1) +
  xlab("Residuals")

autoplot(Neg_Bin_NESW)##diagnostic plot

shapiro.test(Neg_Bin_NESW_Residuals)


##Stepwise elimination on glm.nb() with categorical directional variables
Step_NESW_GM = stepAIC(Neg_Bin_NESW)

Stepped_NB_GM_NESW = glm.nb(log(N.Collected + 1) ~ 
                         log(Dist) + 
                         log(Rainfall + 1) + 
                         log(Wild_and_Wolbachia + 1) + 
                         log(Human_Density) + 
                         log(Av_Mean_Temp) +
                         Color + 
                         log(Wind_Intensity) + 
                         Wind_Direction_NESW,  
                         method = "glm.fit", 
                         link = "log",
                         data = MRR_NESW)
                        


summary(Stepped_NB_GM_NESW)

Neg_Bin_NESW_2_Residuals = residuals.glm(Stepped_NB_GM_NESW, "pearson")

#histogram of pearson residuals
ggplot(, aes(Neg_Bin_NESW_2_Residuals)) +
  geom_histogram(binwidth = 0.1) +
  xlab("Residuals") +
  theme_bw()

autoplot(Stepped_NB_GM_NESW) ##diagnostic plot
shapiro.test(Neg_Bin_NESW_2_Residuals)


##Convert N.Collected to Binary and test if makes better model
MRR_Dataset_Collected = MRR_NESW%>%
  filter(N.Collected > 0)%>%
  mutate(Collection = 1)

MRR_Dataset_Not_Collected = MRR_NESW%>%
  filter(N.Collected == 0)%>%
  mutate(Collection = 0)

##rejoin binary subsets
MRR_Dataset_Binary = rbind(MRR_Dataset_Collected, MRR_Dataset_Not_Collected)

##binomial glm with distance only covariate
Binary_Dist_GLM = glm(formula = cbind(log(N.Collected+1), log(N.Released)) ~ log(Dist),
                      family = "binomial", data = MRR_Dataset_Binary)

summary(Binary_Dist_GLM)
Binary_Dist_GLM_Residuals = residuals.glm(Binary_Dist_GLM, "pearson") ##pearson residuals

##histogram of pearson residuals
ggplot(, aes(x=Binary_Dist_GLM_Residuals)) +
         geom_histogram(binwidth = 0.1) +
          xlab("Residuals")

shapiro.test(Binary_Dist_GLM_Residuals)
autoplot(Binary_Dist_GLM) ##diagnostic plot

##predictions of binomial glm with distance as predictor
predict.binom.count =
  predict(Binary_Dist_GLM, newdata=data.frame(Dist = Distance),
          type="response", se.fit = TRUE)

##Calculate 95% CIs
upperCI_Binom = as.numeric(predict.binom.count$fit + (1.96 * predict.binom.count$se.fit))
lowerCI_Binom = as.numeric(predict.binom.count$fit - (1.96 * predict.binom.count$se.fit))
fit_values_Binom = as.numeric(predict.binom.count$fit)

Binom_Predict = data.frame(Distance, upperCI_Binom, lowerCI_Binom, fit_values_Binom)

##graph of binomial predictions with CIs and data points
Binomial_Predict_Graph = ggplot(Binom_Predict, aes(x=log(Distance), y=fit_values_Binom)) +
  geom_line(colour = "green") +
  geom_line(colour = "red", aes(y=upperCI_Binom, x= log(Distance))) +
  geom_line(colour = "red", aes(y=lowerCI_Binom, x = log(Distance))) +
  geom_point(data = MRR_Dataset_Binary, aes(x=log(Dist), y=Collection, colour = Color))+
  scale_colour_manual(values= c("green", "grey", "purple", "orange", "pink", "yellow"))+
  ggtitle("Binomial") +
  theme_bw() +
  xlab("log(Distance (m) )")+
  ylab("Marked Mosquito Presence") +
  theme(legend.position = "none")
  
Binomial_Predict_Graph

##Binomial glm with multiple covariates
Binary_GLM = glm(formula = Collection ~ 
                        log(Dist) +
                        log(Rainfall + 1) + 
                        log(Wild_and_Wolbachia + 1) +
                        log(Human_Density) +
                        log(Av_Mean_Temp) +
                        Color + BG_Location +
                        log(Direction) +
                        log(Wind_Direction)*
                        log(Wind_Intensity), 
                      family = "binomial",
                      data = MRR_Dataset_Binary)

summary(Binary_GLM)
##stepwise elimination
step(Binary_GLM)

##run stepwise result
Stepped_Binary = glm(formula = Collection ~ 
                       log(Dist) + 
                       log(Wild_and_Wolbachia + 1) + 
                       log(Human_Density) + 
                       log(Av_Mean_Temp) + 
                       Color + BG_Location + 
                       log(Wind_Direction) +
                       log(Wind_Intensity) + 
                       log(Wind_Direction):log(Wind_Intensity), 
                     family = "binomial", data = MRR_Dataset_Binary)


summary(Stepped_Binary)

Stepped_Binary_Residuals = residuals.glm(Stepped_Binary, "pearson")
shapiro.test(Stepped_Binary_Residuals)
#histogram of residuals
ggplot(, aes(x=Stepped_Binary_Residuals)) +
  geom_histogram(binwidth = 0.1) +
  xlab("Residuals")
autoplot(Stepped_Binary)


##binomial glm using categorical directions
Binary_GLM_NESW = glm(formula = Collection ~ 
                   log(Dist) +
                   log(Rainfall + 1) + 
                   log(Wind_Intensity)*
                   Wind_Direction_NESW +    
                   log(Wild_and_Wolbachia + 1) +
                   log(Human_Density) +
                   log(Av_Mean_Temp) +
                   Color + BG_Location +
                   Flight_Direction_NESW,
                   family = "binomial",
                 data = MRR_Dataset_Binary)

summary(Binary_GLM_NESW)
##stepwise elimination
step(Binary_GLM_NESW)
##run result from stepwise
Stepped_Binary_NESW = glm(formula = Collection ~ 
                            log(Dist) + 
                            log(Wild_and_Wolbachia + 1) +
                            log(Human_Density) + 
                            Color + 
                            BG_Location + 
                            Flight_Direction_NESW, 
                          family = "binomial", data = MRR_Dataset_Binary)
summary(Stepped_Binary_NESW)

Stepped_Binary_NESW_Residuals = residuals.glm(Stepped_Binary_NESW, "pearson")

ggplot(, aes(x=Stepped_Binary_NESW_Residuals)) +
  geom_histogram(binwidth = 0.1) +
  xlab("Residuals")

autoplot(Stepped_Binary_NESW)

shapiro.test(Stepped_Binary_NESW_Residuals)

##Plot Predictions from GLMs using distance only  as predictor
grid.arrange(Poisson_Predict_Graph, Negative_Binomial_Predict_Graph, Binomial_Predict_Graph, ncol = 3)
             


