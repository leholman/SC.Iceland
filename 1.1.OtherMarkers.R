







### MAM datasets
MAM.P19 <- read.csv("cleaneddata/combinedcoredata/MAM.PC019.csv",row.names = 1)
MAM.P19.nREPS <- NrepsMaker(MAM.P19,gsub("(.*)_[0-9]$","\\1",colnames(MAM.P19)))
MAM.GC1 <- read.csv("cleaneddata/combinedcoredata/MAM.GC01.csv",row.names = 1)
MAM.GC1.nREPS <- NrepsMaker(MAM.GC1,gsub("(.*)_[0-9]$","\\1",colnames(MAM.GC1)))

### RIZ datasets
RIZ.P19 <- read.csv("cleaneddata/combinedcoredata/RIZ.PC019.csv",row.names = 1)
RIZ.P19.nREPS <- NrepsMaker(RIZ.P19,gsub("(.*)_[0-9]$","\\1",colnames(RIZ.P19)))
RIZ.GC1 <- read.csv("cleaneddata/combinedcoredata/RIZ.GC01.csv",row.names = 1)
RIZ.GC1.nREPS <- NrepsMaker(RIZ.GC1,gsub("(.*)_[0-9]$","\\1",colnames(RIZ.GC1)))



###Make taxonomic files

MAMtax <- read.csv("taxonomy/MAM.combined.parsed.csv")
MAMasv <- readDNAStringSet("cleaneddata/ASVs/MAM.cleaned.fasta")
MAMtax$ASVseq <- as.character(MAMasv)[match(MAMtax$OTU,names(MAMasv))]

write.csv(MAMtax,file = "taxonomy/byHand/MAMtax.csv")

RIZtax <- read.csv("taxonomy/RIZ.combined.parsed.csv")
RIZasv <- readDNAStringSet("cleaneddata/ASVs/RIZ.cleaned.fasta")
RIZtax$ASVseq <- as.character(RIZasv)[match(RIZtax$OTU,names(RIZasv))]

write.csv(RIZtax,file = "taxonomy/byHand/RIZtax.csv")




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






#### workshop 

taxa <- "Cystophora cristata"
#for (taxa in unique(MAMtax.h$ID[MAMtax.h$Level=="Family"])){
#for (taxa in unique(MAMtax.h$ID[MAMtax.h$Level=="Genus"])){
for (taxa in unique(MAMtax.h$ID[MAMtax.h$Level=="Species" & MAMtax.h$X.1.base.in.difference=="Y"])){
  loopASVs <- MAMtax.h$OTU[MAMtax.h$ID==taxa]
  loopData1 <- colSums(MAM.P19[loopASVs,])
  loopData1 <- tapply(loopData1 > 0, gsub("(.*)_[0-9]$","\\1",names(loopData1)), sum)
  loopData2 <- colSums(MAM.GC1[loopASVs,])
  loopData2 <- tapply(loopData2 > 0, gsub("(.*)_[0-9]$","\\1",names(loopData2)), sum)
  pdf(paste0("figures/SpecificTaxa/test/MAM.",taxa,".pdf"),height=7,width = 9)
  par(mfrow=c(2,1),mar=c(4.1, 4.1, 1.1, 1.1))
  plot(1950-ages$mean[match(names(loopData1),ages$ID2)],jitter(as.numeric(loopData1)),main=paste0("PC19 Species:",taxa),pch=16,xlab="Year (CE)",ylab="replicates +ive",xlim=c(-1550,1700))
  plot(1950-ages$mean[match(names(loopData2),ages$ID2)],jitter(as.numeric(loopData2)),main=paste0("GC01 Species:",taxa),pch=16,xlab="Year (CE)",ylab="replicates +ive",xlim=c(-1550,1700))
  dev.off()
}

tapply(loopData1 < 0, names(loopData1), sum)



### whale biz

plot(ages$mean[match(colnames(MAM.P19.nREPS[245,]),ages$ID2)],jitter(as.numeric(MAM.P19.nREPS[245,])),pch=16,xlab="Cal Yr BP",ylab="replicates +ive")

plot(ages$mean[match(colnames(MAM.P19.nREPS[3,]),ages$ID2)],jitter(as.numeric(MAM.P19.nREPS[43,])),pch=16,xlab="Cal Yr BP",ylab="replicates +ive")


#right whale ASV number
c(10,169,230,293,336,392,490,491,502,504,514,517,534,565,576,594,612)

plot(ages$mean[match(colnames(MAM.P19.nREPS[c(10,169,230,293,336,392,490,491,502,504,514,517,534,565,576,594,612),]),ages$ID2)],
     jitter(as.numeric(colSums(MAM.P19.nREPS[c(10,169,230,293,336,392,490,491,502,504,514,517,534,565,576,594,612),],na.rm = TRUE))),pch=16,xlab="Cal Yr BP",ylab="replicates +ive")


summary(test)
plot(test)

ASV <-838

plot(1950-ages$mean[match(colnames(EUK.P19.nREPS[ASV,]),ages$ID2)],jitter(as.numeric(EUK.P19.nREPS[ASV,])),pch=16,xlab="Cal Yr BP",ylab="replicates +ive",xlim=c(-1550,2000))
plot(1950-ages$mean[match(colnames(EUK.GC1.nREPS[ASV,]),ages$ID2)],jitter(as.numeric(EUK.GC1.nREPS[ASV,])),pch=16,xlab="Cal Yr BP",ylab="replicates +ive",xlim=c(-1550,2000))



plot(1950-ages$mean[match(colnames(EUK.P19.nREPS[168,]),ages$ID2)],jitter(as.numeric(EUK.P19.nREPS[168,])),pch=16,xlab="Cal Yr BP",ylab="replicates +ive",xlim=c(-1550,2000))
plot(1950-ages$mean[match(colnames(EUK.P19.nREPS[776,]),ages$ID2)],jitter(as.numeric(EUK.P19.nREPS[776,])),pch=16,xlab="Cal Yr BP",ylab="replicates +ive",xlim=c(-1550,2000))
plot(1950-ages$mean[match(colnames(EUK.P19.nREPS[838,]),ages$ID2)],jitter(as.numeric(EUK.P19.nREPS[838,])),pch=16,xlab="Cal Yr BP",ylab="replicates +ive",xlim=c(-1550,2000))
plot(1950-ages$mean[match(colnames(EUK.P19.nREPS[1022,]),ages$ID2)],jitter(as.numeric(EUK.P19.nREPS[1022,])),pch=16,xlab="Cal Yr BP",ylab="replicates +ive",xlim=c(-1550,2000))
plot(1950-ages$mean[match(colnames(EUK.P19.nREPS[3101,]),ages$ID2)],jitter(as.numeric(EUK.P19.nREPS[3101,])),pch=16,xlab="Cal Yr BP",ylab="replicates +ive",xlim=c(-1550,2000))
plot(1950-ages$mean[match(colnames(EUK.P19.nREPS[4881,]),ages$ID2)],jitter(as.numeric(EUK.P19.nREPS[4881,])),pch=16,xlab="Cal Yr BP",ylab="replicates +ive",xlim=c(-1550,2000))



pdf("figures/MAM.ASV2.pdf",width = 6,height = 6)
plot(ages$mean[match(colnames(MAM.P19.nREPS[2,]),ages$ID2)],jitter(as.numeric(MAM.P19.nREPS[2,])),pch=16,xlab="Cal Yr BP",ylab="replicates +ive")
dev.off()

pdf("figures/MAM.ASV3.pdf",width = 6,height = 6)
plot(ages$mean[match(colnames(MAM.P19.nREPS[3,]),ages$ID2)],jitter(as.numeric(MAM.P19.nREPS[3,])),pch=16,xlab="Cal Yr BP",ylab="replicates +ive")
dev.off()

pdf("figures/MAM.ASV25.pdf",width = 6,height = 6)
plot(ages$mean[match(colnames(MAM.P19.nREPS[25,]),ages$ID2)],jitter(as.numeric(MAM.P19.nREPS[4,])),pch=16,xlab="Cal Yr BP",ylab="replicates +ive")
dev.off()

ASV_838


## whale GLM

year <- ages$mean[match(colnames(MAM.P19.nREPS[3,]),ages$ID2)]
detections <- as.numeric(MAM.P19.nREPS[3,])

test <- glm(detections ~year,poisson(link = "log"))
summary(test)

prediction <- data.frame("year"=200:3000)
prediction <- cbind(prediction,predict(test,prediction,se.fit=TRUE))
plot(prediction$year,prediction$fit)





