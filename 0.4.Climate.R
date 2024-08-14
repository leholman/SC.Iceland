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

pdf("figures/fig2/d13C.pdf",height=4,width =10)
par(mar=c(4.1,4.1,2.1,6.1))
plot(climate$year[climate$record==subset],climate$value[climate$record==subset],pch=16,
     xlab="",
     col=add.alpha('lightgreen',0.3),
     xlim=c(-1550,2000),
     ylim = rev(range(climate$value[climate$record==subset])),
     xaxt = 'n', bty = 'n',
     yaxt ="n",ylab="")
axis(4,at=c(1.5,2,2.5,3),las=2)
mtext("Arctica d13C", side = 4, line = 3,col="green4")
points(year,valueSpline$curves$Spline,type="l",col="green4",lwd=3)
dev.off()


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


pdf("figures/fig2/d18O.pdf",height=4,width =10)
par(mar=c(4.1,4.1,2.1,6.1))
plot(climate$year[climate$record==subset],climate$value[climate$record==subset],pch=16,
     xlab="",
     col=add.alpha('green3',0.3),
     xlim=c(-1550,2000),
     ylim = rev(range(climate$value[climate$record==subset])),
     xaxt = 'n', bty = 'n',
     yaxt ="n",ylab="")
axis(4,at=c(3,3.5,4),las=2)
mtext("Arctica d18O", side = 4, line = 3,col="green3")
points(year,valueSpline$curves$Spline,type="l",col="seagreen",lwd=3)
dev.off()



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

pdf("figures/fig2/MD99Alkenone.pdf",height=4,width =10)
par(mar=c(4.1,4.1,2.1,6.1))
plot(1950-climate$year[climate$record==subset],climate$value[climate$record==subset],pch=16,
     xlab="",
     col=add.alpha('darkolivegreen',0.3),
     xlim=c(-1550,2000),
     ylim=c(6.5,10),
     xaxt = 'n',yaxt = 'n', bty = 'n',
     ylab ="MD99-2275Alkenone SST",
     col.lab = "darkgreen")
axis(2,seq(7,10,1),las=2)
points(1950-interpol$year,valueSpline$curves$Friedman,type="l",col="darkgreen",lwd=3)
dev.off()





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

pdf("figures/fig2/MD99DiatomTF.pdf",height=4,width =10)
par(mar=c(4.1,4.1,2.1,6.1))
plot(1950-climate$year[climate$record==subset],climate$value[climate$record==subset],pch=16,
     xlab="",
     col=add.alpha('green',0.3),
     xlim=c(-1550,2000),
     ylim=c(6.5,10),
     xaxt = 'n',yaxt = 'n', bty = 'n',
     ylab ="")
abline(h=8,lty=2,col="grey")
axis(2,seq(7,10,1),las=2)
mtext("MD99-2275DiatomTF SST", side = 2, line = 2,col="springgreen2")
points(1950-interpol$year,valueSpline$curves$Friedman,type="l",col="springgreen2",lwd=3)
dev.off()






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



pdf("figures/fig2/MD99-2275IP25.pdf",height=4,width=10)
par(mar=c(4.1,4.1,2.1,6.1))
plot(climate$year[climate$record==subset],climate$value[climate$record==subset],pch=16,
     xlab="",
     col=add.alpha('yellowgreen',0.3),
     xlim=c(-1550,2000),
     ylim=c(0,3.2),
     xaxt = 'n', bty = 'n',
     yaxt="n",
     ylab ="")
axis(4,at=c(0,1,2,3),las=2)
mtext("MD99-2275 IP25 SeaIce", side = 4, line = 3,col="olivedrab")
points(interpol$year,valueSpline$curves$Friedman,type="l",col="olivedrab",lwd=3)
dev.off()


## output climate data 

write.csv(climatedat,"metadata/climate.csv")



#### now lets look at people 


pop <- read.csv("rawdata/human/settlement.csv")

       

pdf("figures/fig2/HumanPOP.pdf",height=5,width=10)
par(mar=c(4.1,4.1,2.1,6.1))
plot(pop$Date,pop$Population,
     ylim=c(0,310000),
     xlim=c(-1550,2000),
      bty = 'n',xaxt='n',yaxt="n",ylab="",xlab="",
     type="l",lwd=3,col="orange")
mtext("Iceland Human Population", side = 4, line = 4,col="orange")
mtext("    vs                  Arch Finds        ", side = 4, line = 3,col="black")
mtext("Marine                                                         ", side = 4, line = 3,col="#AED7E5")
mtext("                     Terrestrial                                        ", side = 4, line = 3,col="#8B5B13")
axis(4,at=c(0,50000,100000,150000,200000,250000,300000),label=paste0(c(0,50000,100000,150000,200000,250000,300000)/1000,c("",rep("k",6))),las=2,cex.axis=0.8)
axis(1,at=seq(-1500,2000,500),labels=paste0(sqrt(seq(-1500,2000,500)^2),c("BCE","BCE","BCE","","CE","CE","CE","CE")),lwd.ticks = 2,cex=2)
dev.off()


pdf("figures/fig1/axis.pdf",height=5,width=10)
par(mar=c(4.1,4.1,2.1,6.1))
plot(NA,ylim=c(0,310000),
     xlim=c(-1500,1500),
     bty = 'n',xaxt='n',yaxt="n",ylab="",xlab="",
     type="l",lwd=3,col="orange")

axis(1,at=seq(-1500,1500,500),labels=paste0(sqrt(seq(-1500,1200,500)^2),c("BCE","BCE","BCE","","CE","CE","CE")),lwd.ticks = 2,cex=2,lwd=0)
dev.off()




## How about he NISP data

nisp <- read.csv("rawdata/human/NISP.csv")

nisp.out <- data.frame("cent"=unique(nisp$Century),
                       "NISPter"=rep(0,length(unique(nisp$Century))),
                       "NISPmar"=rep(0,length(unique(nisp$Century))),
                       "meanMarProp"=rep(0,length(unique(nisp$Century))))

## first lets try try a way where each site is an independent draw of the proportion. 

for (cent in nisp.out$cent){
  loopdat <- nisp[nisp$Century==cent,]
  nisp.out$NISPter[nisp.out$cent==cent] <- sum(loopdat$Terrestrial.Mammals.NISP*loopdat$Spread)
  nisp.out$NISPmar[nisp.out$cent==cent] <- sum(loopdat$Marine.Fish.NISP*loopdat$Spread)
  nisp.out$meanMarProp[nisp.out$cent==cent] <- sum(loopdat$Marine.Fish*loopdat$Spread) / sum(loopdat$Spread)
}

nisp.out$meanTarProp <- 100-nisp.out$meanMarProp

pdf("figures/Fig1/HumDiet.pdf",height=2,width=4.5)
#par(mfrow=c(2,1))
par(mar=c(1.1,1.1,2.1,3.1))
barplot(t(as.matrix(nisp.out[,4:5])),col=c("lightblue","orange4"),names.arg = nisp.out$cent,border=NA,yaxt="n",main="")
dev.off()

barplot(t(as.matrix(cbind(nisp.out[,3]/(nisp.out[,3]+nisp.out[,2]),1-(nisp.out[,3]/(nisp.out[,3]+nisp.out[,2]))))),col=c("darkorange3","orange4"),names.arg = nisp.out$cent,border=NA,yaxt="n",main="Raw.NISP.prop")



#axis(4,at=seq(0,100,20),label=seq(0,100,20),las=2)
dev.off()



### what about hcnage sint he trophic level?

troph <- read.csv("rawdata/human/FishIsotopes.csv")
troph <- troph[troph$Species=="cod",]

value <- troph$N
year <- troph$MiddleDate

gam1 <- gam(value ~ s(year,k=7), method = "REML")
plot(gam1)

max(year)
min(year)

prediction <- data.frame("year"=970:2021)
prediction <- cbind(prediction,predict(gam1,newdata = prediction,se.fit = TRUE))
prediction$uppCI <- prediction$fit+prediction$se.fit*1.96
prediction$lwrCI <- prediction$fit-prediction$se.fit*1.96


pdf("figures/fig1/codN.pdf",height=3,width=10)
par(mar=c(4.1,4.1,2.1,6.1))
plot(troph$MiddleDate,troph$N,pch=16,xlim=c(-1550,2000),
     col="lightblue",
     bty = 'n',xaxt='n',yaxt="n",ylab="",xlab="",
     ylim=c(12,15))
axis(4,at=c(12,13,14,15))
mtext(bquote("Atlantic Cod " ~ delta^15 * N), side = 4, line = 3, col = "cadetblue4")
polygon(c(prediction$year, rev(prediction$year)), c(prediction$uppCI, rev(prediction$lwrCI)), col=add.alpha('lightblue4',0.3), border=NA)
points(prediction$year,prediction$fit,type="l",col="cadetblue4",lwd=3)
dev.off()

