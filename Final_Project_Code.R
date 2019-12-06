#Libraries
install.packages("sf")
install.packages("plyr")
install.packages("dplyr")
install.packages("spdep")
install.packages("GISTools")
install.packages("raster")
install.packages("maptools")
install.packages("rgdal")
install.packages("spatstat")
install.packages("sp")
install.packages("tmap")
install.packages("gstat")
install.packages("leaflet")
install.packages("spgwr")
library(sf)
library(plyr)
library(dplyr)
library(spdep)
library(GISTools)
library(raster)
library(maptools)
library(rgdal)
library(spatstat)
library(sp)
library(tmap)
library(gstat)
library(leaflet)
library(spgwr)

#Cleaning data
#Set working directory
dir <- "/Users/mikenewton/desktop/Geog418Lab/FinalProject/Working"
setwd(dir)

#Reading in particulate matter dataset
pm25 <- read.csv("PM25.csv") #Read in PM2.5 data
#Select only columns 1 and 2
pm25 <- pm25[,1:2]
#Change the column names 
colnames(pm25) <- c("POSTALCODE", "PM25")
pm25 <- na.omit(pm25)

#Reading in postal code shapefile
postalcodes <- shapefile("./BC_PostalCodes/BC_Postal_Codes.shp") #Read in related postal code data
postalcodes <- spTransform(postalcodes, CRS("+init=epsg:32610")) #change projection to NAD83 UTM (Zone 10)
crs(postalcodes)
#Reading in dissemination tract and income data
income <- read.csv("Income.csv") #Read in census income data  
colnames(income) <- c("DAUID", "Income") #Select only ID and Income columns
census.tracts <- shapefile("./BC_DA/BC_DA.shp") #Read in dissemination tract shapefile
census.tracts <- spTransform(census.tracts, CRS("+init=epsg:32610")) #change projection to NAD83 UTM (Zone 10)
crs(census.tracts)
income.tracts <- merge(census.tracts,income, by = "DAUID") #Merge income and dissemination data
nrow(income.tracts) #Determine the number of columns in the dataframe
income.tracts <- income.tracts[!is.na(income.tracts$Income),]
nrow(income.tracts)

#Create choropleth map of income
png("./Figures/Median_Income_Chloropleth.png", width = 4, height = 3, units = "in", res = 300) #create a png file called Median_Income.png
tm_shape(income.tracts) + 
  tm_polygons(col = "Income", 
              title = "Median Income", #title of legend
              style = "jenks", #use natural break Jenks classification
              palette = "viridis", n = 6, lwd = 0.2)
#tm_layout((title = "Median Income Census Tracts in Vancouver"), title.position = c("left", "TOP")) #lowercase title posision includes margins, uppercase title position without margins. makes it tighter to the frame. can also use 0-1 for X and Y coordinates
dev.off() #close the png file

#Select postal codes that fall within dissemination tracts)
postalcodes <- intersect(postalcodes,income.tracts)
#plot(postalcodes) #See what the data looks like spatially
head(postalcodes) #See what the data looks like in tabular form

#Join PM2.5 data with postal code data
pm25.spatial <- merge(postalcodes,pm25,by = "POSTALCODE")

#Aggregate the PM2.5 values in each DA in order to have a single value per DA. Here we aggregate based on the max.
pm25.aggregate <- aggregate((as.numeric(pm25.spatial$PM25)/10)~pm25.spatial$DAUID,FUN=max)

#Re-join aggregated data to the income.tracts layer.
colnames(pm25.aggregate) <- c("DAUID", "PM25AGG") #Select only ID and Income columns
income.pm25 <- merge(income.tracts,pm25.aggregate, by = "DAUID") #Merge income and dissemination data

#Re-join aggregated data to the pm25.spatial points layer.
pm25.points.aggregate <- merge(pm25.spatial, pm25.aggregate, by = "DAUID")

#Create a subsample of the datapoints provided in the PM2.5 dataset using the sample n provided on CourseSpaces
set.seed(240)
sampleSize=240
spSample <- pm25.points.aggregate[sample(1:length(pm25.points.aggregate),sampleSize),]

#Create a grid called grd to use in your interpolation
# Create an empty grid where n is the total number of cells
grd <- as.data.frame(spsample(income.tracts, "regular", n=20000)) #n defines how many cells are in the created grid. balance between looking smooth and not taking too much time to compute (not too fine but not too coarse) #income.tracts for full study area, spSample for boundary box of points
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object
proj4string(grd) <- proj4string(spSample)

######## mapping subset points
png("./Figures/SubsetPoints_PM25.png", width = 6.5, height = 6, units = "in", res = 300)
tm_shape(income.tracts) + 
  tm_polygons() +
  tm_style(style = 'white') +
  tm_shape(spSample) +
  tm_dots(col = "PM25AGG", palette = "YlOrBr",
          title="Sampled PM (<2.5µm)\nin (µg/m^3)", size=0.5, shape = 21) + 
  tm_layout(legend.title.size = 0.9, legend.text.size = 0.7, legend.position = c("left", "bottom")) +
  tm_compass(position = c(0.75,0.85)) +
  tm_scale_bar(position = c(0.65,0.75), breaks = c(0,5,10,15,20), text.size = 0.7)
dev.off()

############ Morans I (Global and Local) for Income
CRS <- crs(census.tracts)
van.nb <- poly2nb(income.tracts) #neighbour weights matrix
van.net <- nb2lines(van.nb,coords=coordinates(income.tracts), proj4string = CRS) #convert that neighbour matrix to lines that we can plot

png("./Figures/Queen_Weight_Neighbourhoods.png") #create a png file called "" in current directory
tm_shape(income.tracts) + tm_borders(col='darkgrey') + 
  tm_shape(van.net) + tm_lines(col='blue') + #shows neighbours of the polygons (queen)
  #tm_layout((title = "Queen Weight Neighbourhoods\nin Vancouver"), title.position = c("left", "TOP")) +
  tm_add_legend(type = "line", title = "Weight Neighbourhood Classification", labels = "Queen", col = 'blue')
dev.off()

########################
van.lw <- nb2listw(van.nb, zero.policy = TRUE, style = "W") #weight matrix
print.listw(van.lw, zero.policy = TRUE)
########################
income.tracts$IncLagMeans = lag.listw(van.lw, income.tracts$Income, zero.policy = TRUE) # calculate lag means. how a polygon compares to the mean of its neighbours

png("./Figures/Median_Income_Lagged_Mean.png")
tm_shape(income.tracts) + 
  tm_polygons(col = "IncLagMeans", 
              title = "Median Income\nLagged Means", 
              style = "fisher", 
              palette = "viridis", n = 6, lwd = 0.2)
#tm_layout(title = "Median Income Lagged Means Across Census\nTracts in Victoria Census Municipality Area 2016", title.position = c("left", "TOP"))
dev.off() #close the png file

########################
mi <- moran.test(income.tracts$Income, van.lw, zero.policy = TRUE) #global morans I test. run it on median income and give it the weights matrix
mi #results are what we use to calculate a Z score (Moran I statistic vs expected)

moran.range <- function(lw) { #build your own function
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}
moran.range(van.lw)

mI <- mi$estimate[[1]] #grabs the statistic
eI <- mi$estimate[[2]] #grabs the expected
var <- mi$estimate[[3]] #grabs the variance

z <- (mI - eI)/sqrt(var) #z score formula

########################  

lisa.test <- localmoran(income.tracts$Income, van.lw) #local Moran I #attach to van data

income.tracts$Ii <- lisa.test[,1]
income.tracts$E.Ii<- lisa.test[,2]
income.tracts$Var.Ii<- lisa.test[,3]
income.tracts$Z.Ii<- lisa.test[,4]
income.tracts$P<- lisa.test[,5]
########################
png("./Figures/Median_Income_Local_Morans_I_Z_Score.png", width = 6.5, height = 6, units = "in", res = 300)
tm_shape(income.tracts) + 
  tm_polygons(col = "Z.Ii", 
              title = "Local Moran's I\nZ Score", 
              style = "fisher", 
              palette = "viridis", n = 6, lwd = 0.2) +
  tm_compass(position = c(0.75,0.85)) +
  tm_scale_bar(position = c(0.65,0.75), breaks = c(0,5,10,15,20), text.size = 0.7) +
  tm_layout(legend.title.size = 0.9, legend.text.size = 0.7, legend.position = c("left", "bottom"))
#tm_layout(title = "Local Moran's I for Median Income Across Census\nTracts in Vancouver", title.position = c("left", "TOP"))
dev.off() 


########################

#scatterplot of positive autocorrelation
png("./Figures/Median_Income_Moran_Plot.png") 
moran.plot(income.tracts$Income, van.lw, zero.policy=NULL, spChk=NULL, labels=NULL, xlab="Median Income", 
           ylab="Spatially Lagged Median Income", quiet=NULL, title("Comparing Median Income against Spatially Lagged\nMedian Income to Determine Autocorrelation"))
dev.off()

###### Map P Values
png("./Figures/Median_Income_PVal.png")
tm_shape(income.tracts) + 
  tm_polygons(col = "P", 
              title = "P Value", 
              style = "fixed", breaks = c(0,0.01,0.05,0.1,1),
              palette = "viridis", lwd = 0.2)
#tm_layout(title = "LISA P Values for Median Income Across Census\nTracts in Vancouver", title.position = c("left", "TOP"))
dev.off() 

############ Interpolation for PM2.5
#Universal Kriging 2nd order trend model
#Add X and Y to P
spSample$X <- coordinates(spSample)[,1] #creating coordinates X and Y
spSample$Y <- coordinates(spSample)[,2]

f.2 <- as.formula(PM25AGG ~ X + Y + I(X*X)+I(Y*Y) + I(X*Y)) 

var.smpl <- variogram(f.2, spSample, cloud = FALSE)
dat.fit  <- fit.variogram(var.smpl, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(psill=0.75, model="Sph", range=7500, nugget=0))
#plot(var.smpl, dat.fit)
png("./Figures/PM25_Universal_Kriging_semivariogram_2nd_order.png", width = 4, height = 3, units = "in", res = 300)
plot(var.smpl, dat.fit)
dev.off()
# Perform the krige interpolation
dat.krg <- krige(f.2, spSample, grd, dat.fit)

# Convert kriged surface to a raster object for clipping
r.UK <- raster(dat.krg)
r.m.UK <- mask(r.UK, census.tracts)

# Plot the map
png("./Figures/PM25_Kriging_2nd_order.png", width = 4, height = 2.5, units = "in", res = 300)
tm_shape(r.m.UK) + 
  tm_raster(n=10, palette="YlOrBr",  
            title="Predicted PM (<2.5 µm)\n(in µg/m^3)") +
  tm_shape(spSample) + tm_dots(size=0.07) +
  tm_compass(position = c(0.05,0.55), size = 1, text.size = 0.4) +
  tm_scale_bar(position = c(0.02, 0.01), breaks = c(0,5,10), text.size = 0.3) +
  tm_layout(legend.title.size = 0.4, legend.text.size = 0.3, legend.position = c(0.02,0.08)) +
  tm_shape(income.tracts) +
  tm_borders(col = "grey40", lwd=0.1)
#tm_legend(legend.outside=TRUE)
dev.off()

# Variance and 95% CI 2nd order
r.UKvar   <- raster(dat.krg, layer="var1.var")
r.m.UKvar <- mask(r.UKvar, census.tracts)

png("./Figures/PM25_Kriging_Variance_2nd_order.png", width = 4, height = 2.5, units = "in", res = 300)
tm_shape(r.m.UKvar) + 
  tm_raster(n=7, palette ="Reds",
            title="Variance map \n(in squared µg/m^3)") +tm_shape(spSample) + tm_dots(size=0.07) +
  tm_compass(position = c(0.05,0.55), size = 1, text.size = 0.4) +
  tm_scale_bar(position = c(0.02, 0.01), breaks = c(0,5,10), text.size = 0.3) +
  tm_layout(legend.title.size = 0.4, legend.text.size = 0.3, legend.position = c(0.02,0.08)) +
  tm_shape(income.tracts) +
  tm_borders(col = "grey40", lwd=0.1)
dev.off()

r.UK.CI   <- sqrt(raster(dat.krg, layer="var1.var")) * 1.96
r.m.UK.CI <- mask(r.UK.CI, census.tracts)

png("./Figures/PM25_Kriging_95CI_2nd_order.png", width = 4, height = 2.5, units = "in", res = 300)
tm_shape(r.m.UK.CI) + 
  tm_raster(n=7, palette ="Reds",
            title="95% CI map \n(in µg/m^3)") +tm_shape(spSample) + tm_dots(size=0.07) +
  tm_compass(position = c(0.05,0.55), size = 1, text.size = 0.4) +
  tm_scale_bar(position = c(0.02, 0.01), breaks = c(0,5,10), text.size = 0.3) +
  tm_layout(legend.title.size = 0.4, legend.text.size = 0.3, legend.position = c(0.02,0.08)) +
  tm_shape(income.tracts) +
  tm_borders(col = "grey40", lwd=0.1)
dev.off()

############ Linear Regression

#combining Income and PM2.5 for regression
pm.income.poly <- extract(r.UK, income.tracts, fun=mean, sp=TRUE)
# View(pm.income.poly@data)
colnames(pm.income.poly@data)[31] <- "PM25"

#Plot income and PM2.5 from the pm.income.poly dataset you created
#plot(pm.income.poly$Income~pm.income.poly$PM25, xlab = "PM 2.5", ylab = "Median Income ($)")

#Notice that there are a lot of 0's in this dataset. If you decide to remove them, use the following line:
#pm.income.poly <-  pm.income.poly[pm.income.poly$PM25 != 0, ] #don't want to remove 0's as they can still occur
pm.income.poly <- pm.income.poly[!is.na(pm.income.poly$PM25),] #remove null values
pm.income.poly$PM25[pm.income.poly$PM25 < 0] <- 0 #make negative values, zero
#View(pm.income.poly@data)
#Now plot the data again
plot(pm.income.poly$Income~pm.income.poly$PM25, xlab = "PM 2.5 (µg/m^3)", ylab = "Median Income ($)")

#Perform a linear regression on the two variables. You should decide which one is dependent.
lm.model <- lm(pm.income.poly$Income~pm.income.poly$PM25)
#Add the regression model to the plot you created
png("./Figures/PM25_Income_Plot.png", width = 4, height = 4, units = "in", res = 300)
plot(pm.income.poly$Income~pm.income.poly$PM25, xlab = "PM 2.5 (µg/m^3)", ylab = "Median Income ($)", cex=0.4) +
  abline(lm.model, lwd = 2, col = "red")
dev.off()

#Get the summary of the results. 
summary(lm.model)
#par(mfrow = c(2,2)) #sets up result plots in 2 by 2 pane
#plot(lm.model) 

#You want to determine if the model residuals are spatially clustered. 
#First obtain the residuals from the model
model.resids <- as.data.frame(residuals.lm(lm.model))
#Then add the residuals to your spatialpolygon dataframe
pm.income.poly$residuals <- residuals.lm(lm.model)
#Observe the result to make sure it looks correct
head(pm.income.poly)

#Now, create choropleth map of residuals
# resids <- pm.income.poly$residuals
# shades <- auto.shading(resids, n=6, cols = brewer.pal(6, 'Greens'))

png("./Figures/chloro_residuals.png", width = 4, height = 2.5, units = "in", res = 300)
tm_shape(pm.income.poly) + 
  tm_polygons(col = "residuals", 
              title = "Residuals", 
              style = "fisher",
              palette = "viridis", n = 6, lwd = 0.2)
dev.off()

############ Morans I (Global and Local) for the residuals from the regression
############
van.nb.res <- poly2nb(pm.income.poly) #neighbour weights matrix
van.net.res <- nb2lines(van.nb.res,coords=coordinates(pm.income.poly)) #convert that neighbour matrix to lines that we can plot

png("./Figures/Residuals_Queen_Weight_Neighbourhoods.png") 
tm_shape(pm.income.poly) + tm_borders(col='darkgrey') + 
  tm_shape(van.net.res) + tm_lines(col='blue') + #shows neighbours of the polygons (queen)
  #tm_layout((title = "Queen Weight Neighbourhoods\nin Vancouver"), title.position = c("left", "TOP")) +
  tm_add_legend(type = "line", title = "Weight Neighbourhood Classification", labels = "Queen", col = 'blue')
dev.off()

########################
van.lw.res <- nb2listw(van.nb.res, zero.policy = TRUE, style = "W") #weight matrix
print.listw(van.lw.res, zero.policy = TRUE)

########################
pm.income.poly$ResLagMeans = lag.listw(van.lw.res, pm.income.poly$residuals, zero.policy = TRUE) # calculate lag means. how a polygon compares to the mean of its neighbours

png("./Figures/Residuals_Lagged_Mean.png")
tm_shape(pm.income.poly) + 
  tm_polygons(col = "ResLagMeans", 
              title = "Residual\nLagged Means", 
              style = "fisher", 
              palette = "viridis", n = 6, lwd = 0.2)
dev.off()

########################
mi.res <- moran.test(pm.income.poly$residuals, van.lw.res, zero.policy = TRUE) #global morans I test. run it on median income and give it the weights matrix
mi.res #results are what we use to calculate a Z score (Moran I statistic vs expected)

moran.range.res <- function(lw) { #build your own function
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}
moran.range.res(van.lw.res)


mI.res <- mi.res$estimate[[1]] #grabs the statistic
eI.res <- mi.res$estimate[[2]] #grabs the expected
var.res <- mi.res$estimate[[3]] #grabs the variance

z.res <- (mI.res - eI.res)/sqrt(var.res) #z score formula

########################  
lisa.test <- localmoran(pm.income.poly$residuals, van.lw.res) #local Moran I #attach to van data

pm.income.poly$Ii.res <- lisa.test[,1]
pm.income.poly$E.Ii.res<- lisa.test[,2]
pm.income.poly$Var.Ii.res<- lisa.test[,3]
pm.income.poly$Z.Ii.res<- lisa.test[,4]
pm.income.poly$P.res<- lisa.test[,5]
########################

png("./Figures/Residuals_Local_Morans_I_Z_Score.png", width = 4, height = 2.5, units = "in", res = 300)
tm_shape(pm.income.poly) + 
  tm_polygons(col = "Z.Ii.res", 
              title = "Local Moran's I\nZ Score (residuals)", 
              style = "fisher", 
              palette = "viridis", n = 6, lwd = 0.2) +
  tm_compass(position = c(0.05,0.55), size = 1, text.size = 0.4) +
  tm_scale_bar(position = c(0.02, 0.01), breaks = c(0,5,10), text.size = 0.3) +
  tm_layout(legend.title.size = 0.4, legend.text.size = 0.3, legend.position = c(0.02,0.08))
#tm_layout(title = "Local Moran's I for Residuals Across Census\nTracts in Vancouver", title.position = c("left", "TOP"))
dev.off() 

#scatterplot of positive autocorrelation
png("./Figures/Residuals_Moran_Plot.png") 
moran.plot(pm.income.poly$residuals, van.lw.res, zero.policy=NULL, spChk=NULL, labels=NULL, xlab="Residuals", 
           ylab="Spatially Lagged Residuals", quiet=NULL, title("Comparing Residuals against Spatially Lagged\nResiduals to Determine Autocorrelation"))
dev.off() 

###### Map P Values
png("./Figures/Residuals_PVal.png")
tm_shape(pm.income.poly) + 
  tm_polygons(col = "P.res", 
              title = "P Value", 
              style = "fixed", breaks = c(0,0.01,0.05,0.1,1),
              palette = "viridis", lwd=0.2)
#tm_layout(title = "LISA P Values for Residuals Across Census\nTracts in Vancouver", title.position = c("left", "TOP"))
dev.off() 

############ Geographically Weighted Regression (get values for every dissemination area)

#add the polygon coordinates to the spatialpolygondataframe.
#You can obtain the coordinates using the "coordinates" function from the sp library
pm.income.poly.coords <- sp::coordinates(pm.income.poly)
#Observe the result
head(pm.income.poly.coords)
#Now add the coordinates back to the spatialpolygondataframe
pm.income.poly$X <- pm.income.poly.coords[,1]
pm.income.poly$Y <- pm.income.poly.coords[,2]
head(pm.income.poly)

###Determine the bandwidth for GWR: this will take a while
GWRbandwidth <- gwr.sel(pm.income.poly$Income~pm.income.poly$PM25, 
                        data=pm.income.poly, coords=cbind(pm.income.poly$X,pm.income.poly$Y),adapt=T) 

###Perform GWR on the two variables with the bandwidth determined above
###This will take a looooooong while
gwr.model = gwr(pm.income.poly$Income~pm.income.poly$PM25, 
                data=pm.income.poly, coords=cbind(pm.income.poly$X,pm.income.poly$Y), 
                adapt=GWRbandwidth, hatmatrix=TRUE, se.fit=TRUE) 

#Print the results of the model
gwr.model

#Look at the results in detail
results<-as.data.frame(gwr.model$SDF)
head(results)

#Now for the magic. Let's add our local r-square values to the map
pm.income.poly$localr <- results$localR2

#Create choropleth map of r-square values
png("./Figures/choroGWR_r2.png", width = 4, height = 2.5, units = "in", res = 300)
tm_shape(pm.income.poly) + 
  tm_polygons(col = "localr", 
              title = "R^2 of GWR", 
              style = "fisher", 
              palette = "viridis", n = 6, lwd=0.4)+
  tm_compass(position = c(0.05,0.55), size = 1, text.size = 0.4) +
  tm_scale_bar(position = c(0.02, 0.01), breaks = c(0,5,10), text.size = 0.3) +
  tm_layout(legend.title.size = 0.5, legend.text.size = 0.4, legend.position = c(0.02,0.08))
#tm_layout(title = "R^2 Values for GWR Across Census\nTracts in Vancouver", title.position = c("left", "TOP"))
dev.off() 

#Map the coefficients
pm.income.poly$coeff <- results$pm.income.poly.PM25

#Create choropleth map of the coefficients
png("./Figures/choroGWR_coef.png", width = 4, height = 2.5, units = "in", res = 300)
tm_shape(pm.income.poly) + 
  tm_polygons(col = "coeff", 
              title = "Coefficient of GWR", 
              style = "fisher", 
              palette = "RdYlBu", n = 7, lwd=0.4, midpoint = NA)+
  tm_compass(position = c(0.05,0.55), size = 1, text.size = 0.4) +
  tm_scale_bar(position = c(0.02, 0.01), breaks = c(0,5,10), text.size = 0.3) +
  tm_layout(legend.title.size = 0.5, legend.text.size = 0.4, legend.position = c(0.02,0.08))
#tm_layout(title = "Coefficient Values for GWR Across Census\nTracts in Vancouver", title.position = c("left", "TOP"))
dev.off() 

############ Point Pattern Analysis
#check for and remove duplicated points
zd <- zerodist(spSample)
zd
#remove duplicates if needed
spSample <- remove.duplicates(spSample)

#create an "extent" object which can be used to create the observation window for spatstat
spSample.ext <- as.matrix(extent(spSample)) 

#observation window
window <- as.owin(list(xrange = spSample.ext[1,], yrange = spSample.ext[2,]))
window

#create ppp oject from spatstat
spSample.ppp <- ppp(x = spSample$X, y = spSample$Y, window = window)

##QUADRAT ANALYSIS
##First, determine the number of quadrats
quads <- 10
qcount <- quadratcount(spSample.ppp, nx = quads, ny = quads)

# png(paste("./Figures/quadrat_analysis.png"))
# plot(spSample.ppp, pch = "+", cex = 0.7, main = "Point Pattern Quadrat Analysis: N = 240")
# plot(qcount, add = T, col = "red")
# dev.off()

#turn qcount into a dataframe
qcount.df <- as.data.frame(qcount)

##Second, count the number of quadrats with a distinct number of points.
qcount.df <- plyr::count(qcount.df,'Freq')

##Change the column names so that x=number of points per cell and f=frequency of quadrats with x point.
colnames(qcount.df) <- c("x","f")

sum.f.x2 <- sum(qcount.df$f * qcount.df$x ^ 2)
M <- sum(qcount.df$f) #number of cells
N <- sum(qcount.df$f * qcount.df$x) #Number of points in a dataset
sum.fx.2 <- sum(qcount.df$f * qcount.df$x) ^ 2

# VAR <- (sum.f.x2 - ((sum.fx.2)/M)) / (M - 1)
# MEAN <- N / M
# VMR <- VAR / MEAN
# 
# ##Finally, perform the test statistic to test for the existence of a random spatial pattern.
# chi.square = VMR * (M - 1)
# p.val.chi = 1 - pchisq(chi.square, (M - 1))
# 
# p.val.chi

##Nearest Neighbour Distance
###NEAREST NEIGHBOUR
nearestNeighbour <- nndist(spSample.ppp)

##Convert the nearestNeighbor object into a dataframe.
nearestNeighbour=as.data.frame(as.numeric(nearestNeighbour))
##Change the column name to "Distance"
colnames(nearestNeighbour) = "Distance"

##Calculate the nearest neighbor statistic to test for a random spatial distribution.
#mean nearest neighbour
nnd = sum(nearestNeighbour$Distance) / N

#mean nearest neighbour for random spatial distribution
#area of city boundary in m^2
studyArea <- gArea(spgeom = income.tracts, byid = FALSE)
pointDensity <- N / studyArea

r.nnd = 1 / (2 * sqrt(pointDensity))
d.nnd = 1.07453 / sqrt(pointDensity)
R = nnd / r.nnd
SE.NND <- 0.26136 / sqrt(N * pointDensity)
z.gwr = (nnd - r.nnd) / SE.NND

NNDResult <- data.frame(NND = nnd, NNDr = r.nnd, NNDd = d.nnd, R = R, Z = z.gwr)
#View(NNDResult)
