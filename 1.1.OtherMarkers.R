










### now some whale plots

##RIZ
RIZtax.h <- read.csv("taxonomy/byHand/RIZtax_assigned2105.csv")

unique(RIZtax.h$ID[RIZtax.h$Level=="Family"])
unique(RIZtax.h$ID[RIZtax.h$Level=="Genus"])
unique(RIZtax.h$ID[RIZtax.h$Level=="Species" & RIZtax.h$X.1.base.in.difference=="Y"])


for (taxa in unique(RIZtax.h$ID[RIZtax.h$Level=="Species" & RIZtax.h$X.1.base.in.difference=="Y"])){
  loopASVs <- RIZtax.h$OTU[RIZtax.h$ID==taxa]
  loopData1 <- colSums(RIZ.P19[loopASVs,])
  loopData1 <- tapply(loopData1 > 0, gsub("(.*)_[0-9]$","\\1",names(loopData1)), sum)
  loopData2 <- colSums(RIZ.GC1[loopASVs,])
  loopData2 <- tapply(loopData2 > 0, gsub("(.*)_[0-9]$","\\1",names(loopData2)), sum)
  pdf(paste0("figures/SpecificTaxa/RIZ.",taxa,".pdf"),height=7,width = 9)
  par(mfrow=c(2,1),mar=c(4.1, 4.1, 1.1, 1.1))
  plot(1950-ages$mean[match(names(loopData1),ages$ID2)],jitter(as.numeric(loopData1)),main=paste0("PC19 Species:",taxa),pch=16,xlab="Year (CE)",ylab="replicates +ive",xlim=c(-1550,1700))
  plot(1950-ages$mean[match(names(loopData2),ages$ID2)],jitter(as.numeric(loopData2)),main=paste0("GC01 Species:",taxa),pch=16,xlab="Year (CE)",ylab="replicates +ive",xlim=c(-1550,1700))
  dev.off()
}


for (taxa in unique(RIZtax.h$ID[RIZtax.h$Level=="Genus"])){
  loopASVs <- RIZtax.h$OTU[RIZtax.h$ID==taxa]
  loopData1 <- colSums(RIZ.P19[loopASVs,])
  loopData1 <- tapply(loopData1 > 0, gsub("(.*)_[0-9]$","\\1",names(loopData1)), sum)
  loopData2 <- colSums(RIZ.GC1[loopASVs,])
  loopData2 <- tapply(loopData2 > 0, gsub("(.*)_[0-9]$","\\1",names(loopData2)), sum)
  pdf(paste0("figures/SpecificTaxa/RIZ.",taxa,".pdf"),height=7,width = 9)
  par(mfrow=c(2,1),mar=c(4.1, 4.1, 1.1, 1.1))
  plot(1950-ages$mean[match(names(loopData1),ages$ID2)],jitter(as.numeric(loopData1)),main=paste0("PC19 Genus:",taxa),pch=16,xlab="Year (CE)",ylab="replicates +ive",xlim=c(-1550,1700))
  plot(1950-ages$mean[match(names(loopData2),ages$ID2)],jitter(as.numeric(loopData2)),main=paste0("GC01 Genus:",taxa),pch=16,xlab="Year (CE)",ylab="replicates +ive",xlim=c(-1550,1700))
  dev.off()
}


for (taxa in unique(RIZtax.h$Assignment[RIZtax.h$Level=="Family"])){
  loopASVs <- RIZtax.h$OTU[RIZtax.h$ID==taxa]
  loopData1 <- colSums(RIZ.P19[loopASVs,])
  loopData1 <- tapply(loopData1 > 0, gsub("(.*)_[0-9]$","\\1",names(loopData1)), sum)
  loopData2 <- colSums(RIZ.GC1[loopASVs,])
  loopData2 <- tapply(loopData2 > 0, gsub("(.*)_[0-9]$","\\1",names(loopData2)), sum)
  pdf(paste0("figures/SpecificTaxa/RIZ.",taxa,".pdf"),height=7,width = 9)
  par(mfrow=c(2,1),mar=c(4.1, 4.1, 1.1, 1.1))
  plot(1950-ages$mean[match(names(loopData1),ages$ID2)],jitter(as.numeric(loopData1)),main=paste0("PC19 Family:",taxa),pch=16,xlab="Year (CE)",ylab="replicates +ive",xlim=c(-1550,1700))
  plot(1950-ages$mean[match(names(loopData2),ages$ID2)],jitter(as.numeric(loopData2)),main=paste0("GC01 Family:",taxa),pch=16,xlab="Year (CE)",ylab="replicates +ive",xlim=c(-1550,1700))
  dev.off()
}

##MAM

MAMtax.h <- read.csv("taxonomy/byHand/MAMtax_assigned2105.csv")

unique(MAMtax.h$ID[MAMtax.h$Level=="Family"])
unique(MAMtax.h$ID[MAMtax.h$Level=="Genus"])
unique(MAMtax.h$ID[MAMtax.h$Level=="Species" & MAMtax.h$X.1.base.in.difference=="Y"])
taxa <- "Cystophora cristata"

for (taxa in unique(MAMtax.h$ID[MAMtax.h$Level=="Species" & MAMtax.h$X.1.base.in.difference=="Y"])){
  loopASVs <- MAMtax.h$OTU[MAMtax.h$ID==taxa]
  loopData1 <- colSums(MAM.P19[loopASVs,])
  loopData1 <- tapply(loopData1 > 0, gsub("(.*)_[0-9]$","\\1",names(loopData1)), sum)
  loopData2 <- colSums(MAM.GC1[loopASVs,])
  loopData2 <- tapply(loopData2 > 0, gsub("(.*)_[0-9]$","\\1",names(loopData2)), sum)
  pdf(paste0("figures/SpecificTaxa/MAM.",taxa,".pdf"),height=7,width = 9)
  par(mfrow=c(2,1),mar=c(4.1, 4.1, 1.1, 1.1))
  plot(1950-ages$mean[match(names(loopData1),ages$ID2)],jitter(as.numeric(loopData1)),main=paste0("PC19 Species:",taxa),pch=16,xlab="Year (CE)",ylab="replicates +ive",xlim=c(-1550,1700))
  plot(1950-ages$mean[match(names(loopData2),ages$ID2)],jitter(as.numeric(loopData2)),main=paste0("GC01 Species:",taxa),pch=16,xlab="Year (CE)",ylab="replicates +ive",xlim=c(-1550,1700))
  dev.off()
}


for (taxa in unique(MAMtax.h$ID[MAMtax.h$Level=="Genus"])){
  loopASVs <- MAMtax.h$OTU[MAMtax.h$ID==taxa]
  loopData1 <- colSums(MAM.P19[loopASVs,])
  loopData1 <- tapply(loopData1 > 0, gsub("(.*)_[0-9]$","\\1",names(loopData1)), sum)
  loopData2 <- colSums(MAM.GC1[loopASVs,])
  loopData2 <- tapply(loopData2 > 0, gsub("(.*)_[0-9]$","\\1",names(loopData2)), sum)
  pdf(paste0("figures/SpecificTaxa/MAM.",taxa,".pdf"),height=7,width = 9)
  par(mfrow=c(2,1),mar=c(4.1, 4.1, 1.1, 1.1))
  plot(1950-ages$mean[match(names(loopData1),ages$ID2)],jitter(as.numeric(loopData1)),main=paste0("PC19 Genus:",taxa),pch=16,xlab="Year (CE)",ylab="replicates +ive",xlim=c(-1550,1700))
  plot(1950-ages$mean[match(names(loopData2),ages$ID2)],jitter(as.numeric(loopData2)),main=paste0("GC01 Genus:",taxa),pch=16,xlab="Year (CE)",ylab="replicates +ive",xlim=c(-1550,1700))
  dev.off()
}


for (taxa in unique(MAMtax.h$ID[MAMtax.h$Level=="Family"])){
  loopASVs <- MAMtax.h$OTU[MAMtax.h$ID==taxa]
  loopData1 <- colSums(MAM.P19[loopASVs,])
  loopData1 <- tapply(loopData1 > 0, gsub("(.*)_[0-9]$","\\1",names(loopData1)), sum)
  loopData2 <- colSums(MAM.GC1[loopASVs,])
  loopData2 <- tapply(loopData2 > 0, gsub("(.*)_[0-9]$","\\1",names(loopData2)), sum)
  pdf(paste0("figures/SpecificTaxa/MAM.",taxa,".pdf"),height=7,width = 9)
  par(mfrow=c(2,1),mar=c(4.1, 4.1, 1.1, 1.1))
  plot(1950-ages$mean[match(names(loopData1),ages$ID2)],jitter(as.numeric(loopData1)),main=paste0("PC19 Family:",taxa),pch=16,xlab="Year (CE)",ylab="replicates +ive",xlim=c(-1550,1700))
  plot(1950-ages$mean[match(names(loopData2),ages$ID2)],jitter(as.numeric(loopData2)),main=paste0("GC01 Family:",taxa),pch=16,xlab="Year (CE)",ylab="replicates +ive",xlim=c(-1550,1700))
  dev.off()
}