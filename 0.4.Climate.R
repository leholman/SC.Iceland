## Climate code

library("mgcv")
library("dplR")
library("zoo")
# some functions

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}


### lets input data


climate <- read.csv("rawdata/climate/combined.csv")
sampledates <- read.csv("metadata/AgeOut.csv")
climatedat <- data.frame("Sample"=sampledates$ID,"Core"=sampledates$Core,"Date_bp_mean"=sampledates$mean,"Date_CE_mean"=1950-sampledates$mean)
unique(climate$record)


## Lets go! 


## delta 13C
subset <- "Arcticad13C"
climate$record==subset

value <- climate$value[climate$record==subset]
year <- climate$year[climate$record==subset]

valueSpline <- detrend.series(value,nyrs = 100,return.info=TRUE)
plot(year,value,col="grey")
points(year,valueSpline$curves$Spline,type="l",col="yellow",lwd=3)

gam1 <- gam(value ~ s(year,k=50), method = "REML")
plot(gam1)

prediction <- data.frame("year"=953:2000)
prediction <- cbind(prediction,predict(gam1,newdata = prediction,se.fit = TRUE))
prediction$uppCI <- prediction$fit+prediction$se.fit*1.96
prediction$lwrCI <- prediction$fit-prediction$se.fit*1.96

prediction2 <- data.frame("year"=climatedat$Date_CE_mean)
prediction2$year[prediction2$year>max(year)|prediction2$year<min(year)] <- NA
prediction2 <- cbind(prediction2,predict(gam1,newdata = prediction2,se.fit = TRUE))
prediction2$uppCI <- prediction2$fit+prediction2$se.fit*1.96
prediction2$lwrCI <- prediction2$fit-prediction2$se.fit*1.96

prediction2$lowrezspline <- valueSpline$curves$Spline[match(prediction2$year,year)]

prediction$fit
points(year,prediction$fit,type="l")

pdf("figures/climate/Artica_d13C.pdf",height=6,width =11)
plot(climate$year[climate$record==subset],climate$value[climate$record==subset],pch=16,
     xlab="YearCE",
     ylab="Arctica d13C",
     col="grey30")
polygon(c(prediction$year, rev(prediction$year)), c(prediction$uppCI, rev(prediction$lwrCI)), col=add.alpha('darkgreen',0.3), border=NA)
points(prediction$year,prediction$fit,lwd=3,col="darkgreen",type='l')
points(year,valueSpline$curves$Spline,type="l",col="yellow",lwd=3)
dev.off()

colnames(prediction2) <- c("year","d13C.GAM.fit","d13C.GAM.fit.se","d13C.GAM.fit.uppCI","d13C.GAM.fit.lwrCI","d13C.100yrspline")
climatedat <- cbind(climatedat,prediction2)


## delta 18O

subset <- "Arcticad18O"
climate$record==subset

value <- climate$value[climate$record==subset]
year <- climate$year[climate$record==subset]

valueSpline <- detrend.series(value,nyrs = 100,return.info=TRUE)
plot(year,value,col="grey")
points(year,valueSpline$curves$Spline,type="l",col="yellow",lwd=3)


gam1 <- gam(value ~ s(year,k=50), method = "REML")
plot(gam1)

prediction <- data.frame("year"=953:2000)
prediction <- cbind(prediction,predict(gam1,newdata = prediction,se.fit = TRUE))
prediction$uppCI <- prediction$fit+prediction$se.fit*1.96
prediction$lwrCI <- prediction$fit-prediction$se.fit*1.96


prediction2 <- data.frame("year"=climatedat$Date_CE_mean)
prediction2$year[prediction2$year>max(year)|prediction2$year<min(year)] <- NA
prediction2 <- cbind(prediction2,predict(gam1,newdata = prediction2,se.fit = TRUE))
prediction2$uppCI <- prediction2$fit+prediction2$se.fit*1.96
prediction2$lwrCI <- prediction2$fit-prediction2$se.fit*1.96

prediction2$lowrezspline <- valueSpline$curves$Spline[match(prediction2$year,year)]


pdf("figures/climate/Artica_d18O.pdf",height=6,width =11)
plot(climate$year[climate$record==subset],climate$value[climate$record==subset],pch=16,
     xlab="YearCE",
     ylab="Arctica d18O",
     col="grey30")
polygon(c(prediction$year, rev(prediction$year)), c(prediction$uppCI, rev(prediction$lwrCI)), col=add.alpha('darkgreen',0.3), border=NA)
points(prediction$year,prediction$fit,lwd=3,col="darkgreen",type='l')
points(year,valueSpline$curves$Spline,type="l",col="yellow",lwd=3)
dev.off()

colnames(prediction2) <- c("year","d18O.GAM.fit","d18O.GAM.fit.se","d18O.GAM.fit.uppCI","d18O.GAM.fit.lwrCI","d18O.100yrspline")
climatedat <- cbind(climatedat,prediction2)


## Sea ice index
subset <- "SeaIceIndex"
climate$record==subset

value <- climate$value[climate$record==subset]
year <- climate$year[climate$record==subset]

valueSpline <- detrend.series(value,nyrs = 100,return.info=TRUE)
plot(year,value,col="grey")
points(year,valueSpline$curves$Spline,type="l",col="yellow",lwd=3)


gam1 <- gam(value ~ s(year,k=50), method = "REML")
plot(gam1)

max(year)
min(year)

prediction <- data.frame("year"=1600:2007)
prediction <- cbind(prediction,predict(gam1,newdata = prediction,se.fit = TRUE))
prediction$uppCI <- prediction$fit+prediction$se.fit*1.96
prediction$lwrCI <- prediction$fit-prediction$se.fit*1.96

prediction2 <- data.frame("year"=climatedat$Date_CE_mean)
prediction2$year[prediction2$year>max(year)|prediction2$year<min(year)] <- NA
prediction2 <- cbind(prediction2,predict(gam1,newdata = prediction2,se.fit = TRUE))
prediction2$uppCI <- prediction2$fit+prediction2$se.fit*1.96
prediction2$lwrCI <- prediction2$fit-prediction2$se.fit*1.96
prediction2$lowrezspline <- valueSpline$curves$Spline[match(prediction2$year,year)]

pdf("figures/climate/SeaIceIndex.pdf",height=6,width =11)
plot(climate$year[climate$record==subset],climate$value[climate$record==subset],pch=16,
     xlab="YearCE",
     ylab="Seaice Index",
     col="grey30")
polygon(c(prediction$year, rev(prediction$year)), c(prediction$uppCI, rev(prediction$lwrCI)), col=add.alpha('darkgreen',0.3), border=NA)
points(prediction$year,prediction$fit,lwd=3,col="darkgreen",type='l')
points(year,valueSpline$curves$Spline,type="l",col="yellow",lwd=3)
dev.off()

colnames(prediction2) <- c("year","SeaIceIndex.GAM.fit","SeaIceIndex.GAM.fit.se","SeaIceIndex.GAM.fit.uppCI","SeaIceIndex.GAM.fit.lwrCI","SeaIceIndex.100yrspline")
climatedat <- cbind(climatedat,prediction2)



##  "MD99-2275sstAlkenone" 
subset <- "MD99-2275sstAlkenone"
climate$record==subset

value <- climate$value[climate$record==subset]
year <- climate$year[climate$record==subset]

## interpolate and expand data 
seq(min(year), max(year))
interpol <- merge(data.frame("year"= seq(min(year), max(year))),data.frame("year" = year,"value"=value), by = "year", all = TRUE)
interpol$value <- na.approx(interpol$value,x=interpol$year,na.rm = FALSE)                  

valueSpline <- detrend.series(interpol$value,nyrs = 100,return.info=TRUE)
plot(interpol$year,interpol$value,col="grey")
points(interpol$year,valueSpline$curves$Friedman,type="l",col="yellow",lwd=3)

gam1 <- gam(value ~ s(year,k=100), method = "REML")
plot(gam1)

max(year)
min(year)

prediction <- data.frame("year"=-51:9598)
prediction <- cbind(prediction,predict(gam1,newdata = prediction,se.fit = TRUE))
prediction$uppCI <- prediction$fit+prediction$se.fit*1.96
prediction$lwrCI <- prediction$fit-prediction$se.fit*1.96

prediction2 <- data.frame("yearCE"=climatedat$Date_CE_mean,"year"=climatedat$Date_bp_mean)
prediction2$year[prediction2$year>max(year)|prediction2$year<min(year)] <- NA
prediction2 <- cbind(prediction2,predict(gam1,newdata = prediction2,se.fit = TRUE))
prediction2$uppCI <- prediction2$fit+prediction2$se.fit*1.96
prediction2$lwrCI <- prediction2$fit-prediction2$se.fit*1.96
prediction2$lowrezspline <- valueSpline$curves$Friedman[match(prediction2$year,interpol$year)]


pdf("figures/climate/MD99-2275sstAlkenone.pdf",height=6,width =11)
plot(climate$year[climate$record==subset],climate$value[climate$record==subset],pch=16,cex=0.4,
     xlim=c(max(year),min(year)),
     xlab="YearBP",
     ylab="MD99-2275sstAlkenone",
     col="grey30")
polygon(c(prediction$year, rev(prediction$year)), c(prediction$uppCI, rev(prediction$lwrCI)), col=add.alpha('darkgreen',0.3), border=NA)
points(prediction$year,prediction$fit,lwd=1,col="darkgreen",type='l')
points(interpol$year,valueSpline$curves$Friedman,type="l",col="yellow",lwd=3)
dev.off()

colnames(prediction2) <- c("yearCE","yearBP","MD99-2275sstAlkenone.GAM.fit","MD99-2275sstAlkenone.GAM.fit.se","MD99-2275sstAlkenone.GAM.fit.uppCI","MD99-2275sstAlkenone.GAM.fit.lwrCI","MD99-2275sstAlkenone.100yrspline")
climatedat <- cbind(climatedat,prediction2)


##  "MD99-2275sstDiatom"   
subset <- "MD99-2275sstDiatom"
climate$record==subset

value <- climate$value[climate$record==subset]
year <- climate$year[climate$record==subset]

## interpolate and expand data 
seq(min(year), max(year))
interpol <- merge(data.frame("year"= seq(min(year), max(year))),data.frame("year" = year,"value"=value), by = "year", all = TRUE)
interpol$value <- na.approx(interpol$value,x=interpol$year,na.rm = FALSE)                  

valueSpline <- detrend.series(interpol$value,nyrs = 100,return.info=TRUE)
plot(interpol$year,interpol$value,col="grey")
points(interpol$year,valueSpline$curves$Friedman,type="l",col="yellow",lwd=3)

gam1 <- gam(value ~ s(year,k=100), method = "REML")
plot(gam1)

max(year)
min(year)

prediction <- data.frame("year"=75:9271)
prediction <- cbind(prediction,predict(gam1,newdata = prediction,se.fit = TRUE))
prediction$uppCI <- prediction$fit+prediction$se.fit*1.96
prediction$lwrCI <- prediction$fit-prediction$se.fit*1.96

prediction2 <- data.frame("yearCE"=climatedat$Date_CE_mean,"year"=climatedat$Date_bp_mean)
prediction2$year[prediction2$year>max(year)|prediction2$year<min(year)] <- NA
prediction2 <- cbind(prediction2,predict(gam1,newdata = prediction2,se.fit = TRUE))
prediction2$uppCI <- prediction2$fit+prediction2$se.fit*1.96
prediction2$lwrCI <- prediction2$fit-prediction2$se.fit*1.96
prediction2$lowrezspline <- valueSpline$curves$Friedman[match(prediction2$year,interpol$year)]


pdf("figures/climate/MD99-2275sstDiatom.pdf",height=6,width =11)
plot(climate$year[climate$record==subset],climate$value[climate$record==subset],pch=16,
     xlim=c(max(year),min(year)),
     xlab="YearBP",
     ylab="MD99-2275sstDiatom",
     col="grey30")
polygon(c(prediction$year, rev(prediction$year)), c(prediction$uppCI, rev(prediction$lwrCI)), col=add.alpha('darkgreen',0.3), border=NA)
points(prediction$year,prediction$fit,lwd=3,col="darkgreen",type='l')
points(interpol$year,valueSpline$curves$Friedman,type="l",col="yellow",lwd=3)
dev.off()

colnames(prediction2) <- c("yearCE","yearBP","MD99-2275sstDiatom.GAM.fit","MD99-2275sstDiatom.GAM.fit.se","MD99-2275sstDiatom.GAM.fit.uppCI","MD99-2275sstDiatom.GAM.fit.lwrCI","MD99-2275sstDiatom.100yrspline")
climatedat <- cbind(climatedat,prediction2)


## "MD99-2275IP25" 
subset <- "MD99-2275sstIP25" 
climate$record==subset

value <- climate$value[climate$record==subset]
year <- climate$year[climate$record==subset]

## interpolate and expand data ))
interpol <- merge(data.frame("year"= seq(min(year), max(year))),data.frame("year" = year,"value"=value), by = "year", all = TRUE)
interpol$value <- na.approx(interpol$value,x=interpol$year,na.rm = FALSE)                  

valueSpline <- detrend.series(interpol$value,nyrs = 100,return.info=TRUE)
plot(interpol$year,interpol$value,col="grey")
points(interpol$year,valueSpline$curves$Spline,type="l",col="yellow",lwd=3)

gam1 <- gam(value ~ s(year,k=100), method = "REML")
plot(gam1)

max(year)
min(year)

prediction <- data.frame("year"=952:1894)
prediction <- cbind(prediction,predict(gam1,newdata = prediction,se.fit = TRUE))
prediction$uppCI <- prediction$fit+prediction$se.fit*1.96
prediction$lwrCI <- prediction$fit-prediction$se.fit*1.96

prediction2 <- data.frame("year"=climatedat$Date_CE_mean,"yearBP"=climatedat$Date_bp_mean)
prediction2$year[prediction2$year>max(year)|prediction2$year<min(year)] <- NA
prediction2 <- cbind(prediction2,predict(gam1,newdata = prediction2,se.fit = TRUE))
prediction2$uppCI <- prediction2$fit+prediction2$se.fit*1.96
prediction2$lwrCI <- prediction2$fit-prediction2$se.fit*1.96
prediction2$lowrezspline <- valueSpline$curves$Friedman[match(prediction2$year,interpol$year)]


pdf("figures/climate/MD99-2275IP25.pdf",height=6,width =11)
plot(climate$year[climate$record==subset],climate$value[climate$record==subset],pch=16,
     xlim=c(min(year),max(year)),
     xlab="YearCE",
     ylab="MD99-2275IP25",
     col="grey30")
polygon(c(prediction$year, rev(prediction$year)), c(prediction$uppCI, rev(prediction$lwrCI)), col=add.alpha('darkgreen',0.3), border=NA)
points(prediction$year,prediction$fit,lwd=3,col="darkgreen",type='l')
points(interpol$year,valueSpline$curves$Spline,type="l",col="yellow",lwd=3)
dev.off()

colnames(prediction2) <- c("yearCE","yearBP","MD99-2275IP25.GAM.fit","MD99-2275IP25.GAM.fit.se","MD99-2275IP25.GAM.fit.uppCI","MD99-2275IP25.GAM.fit.lwrCI","MD99-2275IP25.100yrspline")
climatedat <- cbind(climatedat,prediction2)

## output climate data 

write.csv(climatedat,"metadata/climate.csv")






