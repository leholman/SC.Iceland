## Climate code

library("mgcv")

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


unique(climate$record)


## Lets go! 


## delta 13C
subset <- "Arcticad13C"
climate$record==subset

value <- climate$value[climate$record==subset]
year <- climate$year[climate$record==subset]

gam1 <- gam(value ~ s(year,k=50), method = "REML")
plot(gam1)

prediction <- data.frame("year"=953:2000)
prediction <- cbind(prediction,predict(gam1,newdata = prediction,se.fit = TRUE))
prediction$uppCI <- prediction$fit+prediction$se.fit*1.96
prediction$lwrCI <- prediction$fit-prediction$se.fit*1.96

pdf("figures/climate/Artica_d13C.pdf",height=6,width =11)
plot(climate$year[climate$record==subset],climate$value[climate$record==subset],pch=16,
     xlab="YearCE",
     ylab="Arctica d13C",
     col="grey30")
polygon(c(prediction$year, rev(prediction$year)), c(prediction$uppCI, rev(prediction$lwrCI)), col=add.alpha('darkgreen',0.3), border=NA)
points(prediction$year,prediction$fit,lwd=3,col="darkgreen",type='l')
dev.off()

## delta 18O

subset <- "Arcticad18O"
climate$record==subset

value <- climate$value[climate$record==subset]
year <- climate$year[climate$record==subset]


gam1 <- gam(value ~ s(year,k=50), method = "REML")
plot(gam1)

prediction <- data.frame("year"=953:2000)
prediction <- cbind(prediction,predict(gam1,newdata = prediction,se.fit = TRUE))
prediction$uppCI <- prediction$fit+prediction$se.fit*1.96
prediction$lwrCI <- prediction$fit-prediction$se.fit*1.96

pdf("figures/climate/Artica_d18O.pdf",height=6,width =11)
plot(climate$year[climate$record==subset],climate$value[climate$record==subset],pch=16,
     xlab="YearCE",
     ylab="Arctica d18O",
     col="grey30")
polygon(c(prediction$year, rev(prediction$year)), c(prediction$uppCI, rev(prediction$lwrCI)), col=add.alpha('darkgreen',0.3), border=NA)
points(prediction$year,prediction$fit,lwd=3,col="darkgreen",type='l')
dev.off()

## Sea ice index
subset <- "SeaIceIndex"
climate$record==subset

value <- climate$value[climate$record==subset]
year <- climate$year[climate$record==subset]


gam1 <- gam(value ~ s(year,k=50), method = "REML")
plot(gam1)

max(year)
min(year)

prediction <- data.frame("year"=1600:2007)
prediction <- cbind(prediction,predict(gam1,newdata = prediction,se.fit = TRUE))
prediction$uppCI <- prediction$fit+prediction$se.fit*1.96
prediction$lwrCI <- prediction$fit-prediction$se.fit*1.96

pdf("figures/climate/SeaIceIndex.pdf",height=6,width =11)
plot(climate$year[climate$record==subset],climate$value[climate$record==subset],pch=16,
     xlab="YearCE",
     ylab="Seaice Index",
     col="grey30")
polygon(c(prediction$year, rev(prediction$year)), c(prediction$uppCI, rev(prediction$lwrCI)), col=add.alpha('darkgreen',0.3), border=NA)
points(prediction$year,prediction$fit,lwd=3,col="darkgreen",type='l')
dev.off()

##  "MD99-2275sstAlkenone" "MD99-2275sstDiatom"   "MD99-2275sstIP25" 
subset <- "MD99-2275sstAlkenone"
climate$record==subset

value <- climate$value[climate$record==subset]
year <- climate$year[climate$record==subset]


gam1 <- gam(value ~ s(year,k=100), method = "REML")
plot(gam1)

max(year)
min(year)

prediction <- data.frame("year"=-51:9598)
prediction <- cbind(prediction,predict(gam1,newdata = prediction,se.fit = TRUE))
prediction$uppCI <- prediction$fit+prediction$se.fit*1.96
prediction$lwrCI <- prediction$fit-prediction$se.fit*1.96

pdf("figures/climate/MD99-2275sstAlkenone.pdf",height=6,width =11)
plot(climate$year[climate$record==subset],climate$value[climate$record==subset],pch=16,cex=0.4,
     xlim=c(max(year),min(year)),
     xlab="YearBP",
     ylab="MD99-2275sstAlkenone",
     col="grey30")
polygon(c(prediction$year, rev(prediction$year)), c(prediction$uppCI, rev(prediction$lwrCI)), col=add.alpha('darkgreen',0.3), border=NA)
points(prediction$year,prediction$fit,lwd=1,col="darkgreen",type='l')
dev.off()



##  "MD99-2275sstDiatom"   "MD99-2275sstIP25" 
subset <- "MD99-2275sstDiatom"
climate$record==subset

value <- climate$value[climate$record==subset]
year <- climate$year[climate$record==subset]


gam1 <- gam(value ~ s(year,k=100), method = "REML")
plot(gam1)

max(year)
min(year)

prediction <- data.frame("year"=75:9271)
prediction <- cbind(prediction,predict(gam1,newdata = prediction,se.fit = TRUE))
prediction$uppCI <- prediction$fit+prediction$se.fit*1.96
prediction$lwrCI <- prediction$fit-prediction$se.fit*1.96

pdf("figures/climate/MD99-2275sstDiatom.pdf",height=6,width =11)
plot(climate$year[climate$record==subset],climate$value[climate$record==subset],pch=16,
     xlim=c(max(year),min(year)),
     xlab="YearBP",
     ylab="MD99-2275sstDiatom",
     col="grey30")
polygon(c(prediction$year, rev(prediction$year)), c(prediction$uppCI, rev(prediction$lwrCI)), col=add.alpha('darkgreen',0.3), border=NA)
points(prediction$year,prediction$fit,lwd=3,col="darkgreen",type='l')
dev.off()

## "MD99-2275IP25" 
subset <- "MD99-2275sstIP25" 
climate$record==subset

value <- climate$value[climate$record==subset]
year <- climate$year[climate$record==subset]


gam1 <- gam(value ~ s(year,k=100), method = "REML")
plot(gam1)

max(year)
min(year)

prediction <- data.frame("year"=952:1894)
prediction <- cbind(prediction,predict(gam1,newdata = prediction,se.fit = TRUE))
prediction$uppCI <- prediction$fit+prediction$se.fit*1.96
prediction$lwrCI <- prediction$fit-prediction$se.fit*1.96

pdf("figures/climate/MD99-2275IP25.pdf",height=6,width =11)
plot(climate$year[climate$record==subset],climate$value[climate$record==subset],pch=16,
     xlim=c(max(year),min(year)),
     xlab="YearBP",
     ylab="MD99-2275IP25",
     col="grey30")
polygon(c(prediction$year, rev(prediction$year)), c(prediction$uppCI, rev(prediction$lwrCI)), col=add.alpha('darkgreen',0.3), border=NA)
points(prediction$year,prediction$fit,lwd=3,col="darkgreen",type='l')
dev.off()



