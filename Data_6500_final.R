#######################
#LOAD/PREPARE THE DATA#

#load required packages#
library(maptools)
library(spatstat)
library("rgdal")
library("sp")
library(splines)

#create directories for the data files#
datadir = "/Users/seanchapman/Downloads/6500_final/Traffic_Collisions_(ASR-T-TBL-001)"
KSIdir = "/Users/seanchapman/Downloads/6500_final/KSI"
bound = "/Users/seanchapman/Downloads/6500_final/torontoBoundary_wgs84" #bound - Toronto boundary polygon 

#read in the data using readOGR, collisions, KSI - spatial points dataframe, bound - spatial polygon
collisions = readOGR(datadir,"Traffic_Collisions_(ASR-T-TBL-001)") #shapefiles, that type of thing sp
KSI = readOGR(KSIdir,"KSI") #error cannot open layer
bound = readOGR(bound, "citygcs_regional_mun_wgs84")

#create function to extract coordinates from bound, required to get bound to work as owin for ppp object#
extractCoords <- function(bound)
{
  results <- list()
  for(i in 1:length(bound@polygons[[1]]@Polygons))
  {
    results[[i]] <- bound@polygons[[1]]@Polygons[[i]]@coords
  }
  results <- Reduce(rbind, results)
  results
}

#convert bound to a obervation window object called window#
coords <- extractCoords(bound)
poly <- sp::Polygon(coords)
poly1 <- sp::Polygons(list(poly),ID="A")
spoly <- sp::SpatialPolygons(list(poly1))
window = as.owin(spoly)
plot(window, axes = TRUE, main = "Traffic Collisions - KSI - Toronto", xlab = "Longitude", ylab = "Latitude")

#create ppp objects for analysis#

#KSI objects#
#df_ppp, ppp object contains all KSI information 
df = as.data.frame(KSI) #KSI
coordinates(df) = ~LONGITUDE+LATITUDE
plot(df,axes = TRUE)
df_ppp <- as(df,"ppp")
df_ppp$window <- as(bound,"owin")

#df2_ppp, ppp object contains KSI coordinates, no covariates, similar to japanesepines dataset
#used for F,K,L functions
df2 = df[,14:15]
View(df2)
coordinates(df2) = ~LONGITUDE+LATITUDE
df2_ppp <- as(df2,"ppp")
df2_ppp$window = window
plot(df2_ppp, add = TRUE)

#ppp_reg, ppp object contains KSI coordinates, binary fatal non-fatal covariate, for relative risk regression
dfreg = df2
dfreg$marks = ifelse(df$ACCLASS=="Fatal", "Fatal","Non-Fatal")
unique(dfreg$marks)
dfreg$marks = as.factor(dfreg$marks)
class(dfreg$marks)
coordinates(dfreg) = ~LONGITUDE+LATITUDE
View(dfreg$marks)
ppp_reg = as(dfreg,"ppp")
ppp_reg$window = window

#collisions objects#
#filter collisions with coordinates 0,0 for latitude/longitude
collisions_2 = subset(collisions, collisions$Longitude != 0)
plot(collisions_2) 

#df2_col_ppp, contains the coordinates of total collisions data set, similar to japanese pines, for F function
df_col = as.data.frame(collisions_2)
View(df2_col)
df2_col = df_col[,15:16]
coordinates(df2_col) = ~Longitude+Latitude
df2_col_ppp = as(df2_col,"ppp")
df2_col_ppp$window = window

#ppp_col, object containing coordinates and injury information for relative risk regression using full collisions dataset
dfreg_col = df2_col
dfreg_col$Injured = df_col$Injury_Col
dfreg_col$Injured = as.factor(dfreg_col$Injured)
class(dfreg_col$Injured)
coordinates(dfreg_col) = ~Longitude+Latitude
ppp_col <- as(dfreg_col,"ppp")
ppp_col$window = window

#ppp_sample, reduced version of ppp_col for relative risk regression that is able to work
?sample
sample = sample(1:452497, 20000, replace = FALSE)
df_sample <- dfreg_col[sample,]
View(df_sample)
coordinates(df_sample) = ~Longitude+Latitude
ppp_sample <- as(df_sample,"ppp")
ppp_sample$window = window

##########################
#RELATIVE RISK REGRESSION#

#KSI relative risk regression using ppp_reg object#
#cases as fatal, control as non-fatal
View(ppp_reg$marks)
cases<-unmark(subset(ppp_reg, marks(ppp_reg) == "Fatal"))
controls<-unmark(subset(ppp_reg, marks(ppp_reg) !="Fatal"))
length(cases$x)
length(controls$x)

par(mfrow=c(1,1))
plot(density(cases),col=gray.colors(12,rev=T), main="Fatal")
plot(density(controls),col=gray.colors(12,rev=T),main="Non-Fatal")

rrbwl=bw.relrisk(ppp_reg, warn=F)
relative_risk=relrisk(ppp_reg, rrbwl, relative=TRUE) 
plot(relative_risk,col=gray.colors(12,rev=T))

#Collisions#
#Attempt to model full dataset, won't run#
cases2<-unmark(subset(ppp_col, marks(ppp_col) =="YES"))
controls2<-unmark(subset(ppp_col, marks(ppp_col) =="NO"))

plot(density(cases2),col=gray.colors(12,rev=T), main="Injury YES")
plot(density(controls2),col=gray.colors(12,rev=T),main="injury NO")

rrbwl2=bw.relrisk(ppp_col, warn=F) #this is the line that not compute
relative_risk2=relrisk(ppp_col, rrbwl2, relative=TRUE) 
plot(relative_risk,col=gray.colors(12,rev=T))

#model of the reduced dataset using sample(), will run#
cases3<-unmark(subset(ppp_sample, marks(ppp_sample) =="YES"))
controls3<-unmark(subset(ppp_sample, marks(ppp_sample) =="NO"))

plot(density(cases3),col=gray.colors(12,rev=T), main="Injury YES")
plot(density(controls3),col=gray.colors(12,rev=T),main="Injury NO")

rrbwl3=bw.relrisk(ppp_sample, warn=F) #this line will take some time
relative_risk3=relrisk(ppp_sample, rrbwl3, relative=TRUE) 
plot(relative_risk3,col=gray.colors(12,rev=T))

############################
#SPATIAL DISTRIBUTION PLOTS#

set.seed(10) #set seed to repeat this procedure.

#K estimate envelope plot
envKSI_K=envelope(df2_ppp,fun=Kest,global=T,nrank=5,nsim=99)
envKSI_K 
plot(envKSI_K, main="KSI K-function")

#F estimate envelope plot
envKSI_F=envelope(df2_ppp,fun=Fest,global=T,nrank=5,nsim=99)
envKSI_F 
plot(envKSI_F, main="KSI F-function")

#L estimate envelope plot
envKSI_L=envelope(df2_ppp,fun=Linhom,global=T,nrank=5,nsim=99)
envKSI_L 
plot(envKSI_L, main="KSI L-Function")

#F inhom estimate envelope plot, takes much longer
envKSI_F=envelope(df2_ppp,fun=Fest,global=T,nrank=5,nsim=99)
envKSI_F 
plot(envKSI_F, main="KSI F-function")

#F inhom estimate plot, no envelope
hold2=Finhom(df2_ppp) 
names(hold2)
plot(hold2)
View(hold2)

#COLLISIONS#

#F estimate envelope plot
envCOL_F=envelope(df2_col_ppp,fun=Fest,global=T,nrank=5,nsim=99)
envCOL_F 
plot(envCOL_F, main="Collisions F-function")

#These estimates take too long, no completed sims after ~20mins
#L estimate envelope plot, takes too long to compute
envCOL_L=envelope(df2_ppp,fun=Lest,global=T,nrank=5,nsim=99)
envCOL_L 
plot(envCOL_L, main="Collisions L-Function")

#K estimate envelope plot, takes too long to compute
envCOL_K=envelope(df2_ppp,fun=Kest,global=T,nrank=5,nsim=99)
envCOL_K 
plot(envCOL_K, main="Collisions K-Function")







