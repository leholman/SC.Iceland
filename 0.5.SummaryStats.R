################################################################
####======     Summary Statistics.                =======#######
####==== Luke E. Holman====26.07.2024====================#######
################################################################

# read in data 
mtDat <- readxl::read_excel("metadata/AllSampleData.xlsx")

GC01 <- read.csv("cleaneddata/combinedcoredata/EUK.GC01.csv",row.names = 1)
PC19 <- read.csv("cleaneddata/combinedcoredata/EUK.PC019.csv",row.names = 1)

neg.B5 <- read.csv("cleaneddata/negativedata/neg.B05.EUK.csv",row.names = 1)
neg.B5.dat <- neg.B5[,2:127]
neg.L1 <- read.csv("cleaneddata/negativedata/neg.LN1.EUK.csv",row.names = 1)
neg.L1.dat <- neg.L1[,2:72]
neg.L2 <- read.csv("cleaneddata/negativedata/neg.LN2.EUK.csv",row.names = 1)
neg.L2.dat <- neg.L2[,2:87]

##### READ COUNTS #####

## make a dataframe to collect results 
output <- data.frame("ExperimentalMean"=rep(0,3),"ExperimentalSD"=rep(0,3),"NegativeControlMean"=rep(0,3),"NegativeControlSD"=rep(0,3),"PositiveControlMean"=rep(0,3),"PositiveControlSD"=rep(0,3),"SummedTotalExperimental"=rep(0,3),"SummedTotalNegative"=rep(0,3),"SummedTotalPositive"=rep(0,3),"SummedTotals"=rep(0,3))
rownames(output) <- c("InitialDemulitplex","RevPrimerStripped","PostDADA2")

# combine sense and antisense counts
mtDat$ReadCountTotal <- mtDat$SenseDemultiplex+mtDat$AntiSenseDemultiplex
mtDat$ReadCountTotalS <- mtDat$SenseStripped+mtDat$AntSenseStripped
# How many reads were demultiplexed? 
output["InitialDemulitplex","SummedTotals"] <- sum(mtDat$ReadCountTotal)
output["RevPrimerStripped","SummedTotals"] <- sum(mtDat$ReadCountTotalS)


# how many reads per each type of sample?
# experimental
output["InitialDemulitplex","SummedTotalExperimental"] <- sum(mtDat$ReadCountTotal[mtDat$SampleType=="Experimental"])
output["InitialDemulitplex","ExperimentalMean"] <- mean(mtDat$ReadCountTotal[mtDat$SampleType=="Experimental"])
output["InitialDemulitplex","ExperimentalSD"] <- sd(mtDat$ReadCountTotal[mtDat$SampleType=="Experimental"])
output["RevPrimerStripped","SummedTotalExperimental"] <- sum(mtDat$ReadCountTotalS[mtDat$SampleType=="Experimental"])
output["RevPrimerStripped","ExperimentalMean"] <- mean(mtDat$ReadCountTotalS[mtDat$SampleType=="Experimental"])
output["RevPrimerStripped","ExperimentalSD"] <- sd(mtDat$ReadCountTotalS[mtDat$SampleType=="Experimental"])

# negative controls 
output["InitialDemulitplex","SummedTotalNegative"] <- sum(mtDat$ReadCountTotal[mtDat$SampleType=="ControlN"])
output["InitialDemulitplex","NegativeControlMean"] <- mean(mtDat$ReadCountTotal[mtDat$SampleType=="ControlN"])
output["InitialDemulitplex","NegativeControlSD"] <- sd(mtDat$ReadCountTotal[mtDat$SampleType=="ControlN"])
output["RevPrimerStripped","SummedTotalNegative"] <- sum(mtDat$ReadCountTotalS[mtDat$SampleType=="ControlN"])
output["RevPrimerStripped","NegativeControlMean"] <- mean(mtDat$ReadCountTotalS[mtDat$SampleType=="ControlN"])
output["RevPrimerStripped","NegativeControlSD"] <- sd(mtDat$ReadCountTotalS[mtDat$SampleType=="ControlN"])

# Positive controls
output["InitialDemulitplex","SummedTotalPositive"] <- sum(mtDat$ReadCountTotal[mtDat$SampleType=="ControlP"])
output["InitialDemulitplex","PositiveControlMean"] <- mean(mtDat$ReadCountTotal[mtDat$SampleType=="ControlP"])
output["InitialDemulitplex","PositiveControlSD"] <- sd(mtDat$ReadCountTotal[mtDat$SampleType=="ControlP"])
output["RevPrimerStripped","SummedTotalPositive"] <- sum(mtDat$ReadCountTotalS[mtDat$SampleType=="ControlP"])
output["RevPrimerStripped","PositiveControlMean"] <- mean(mtDat$ReadCountTotalS[mtDat$SampleType=="ControlP"])
output["RevPrimerStripped","PositiveControlSD"] <- sd(mtDat$ReadCountTotalS[mtDat$SampleType=="ControlP"])


## How many reads remaining after dada2
#Experimental
output["PostDADA2","SummedTotalExperimental"]  <- sum(c(colSums(GC01),colSums(PC19)))
output["PostDADA2","ExperimentalMean"] <- mean(c(colSums(GC01),colSums(PC19)))
output["PostDADA2","ExperimentalSD"] <- sd(c(colSums(GC01),colSums(PC19)))

#Negative 
#NOTE becuase some negatives have been removed from the data as they are empty we need to add them back in to the caluclations of mean and sd - which is why they look complicated 
output["PostDADA2","SummedTotalNegative"] <- sum(c(colSums(neg.B5.dat),colSums(neg.L1.dat),colSums(neg.L2.dat)))
output["PostDADA2","NegativeControlMean"] <- sum(c(colSums(neg.B5.dat),colSums(neg.L1.dat),colSums(neg.L2.dat)))/320
output["PostDADA2","NegativeControlSD"] <- sqrt((sum((mean(c(colSums(neg.B5.dat),colSums(neg.L1.dat),colSums(neg.L2.dat)))-c(colSums(neg.B5.dat),colSums(neg.L1.dat),colSums(neg.L2.dat)))^2))/320-1)


hist(log10(c(colSums(neg.B5.dat),colSums(neg.L1.dat),colSums(neg.L2.dat))),breaks=100,xlab="log10 reads",main="")
hist(c(colSums(neg.B5.dat),colSums(neg.L1.dat),colSums(neg.L2.dat)),breaks=100,xlab="reads",main="")

## what if we get rid of the samples with large numbers of reads?
negs <- c(colSums(neg.B5.dat),colSums(neg.L1.dat),colSums(neg.L2.dat))
negs2 <- negs[negs<3000]
sum(negs2)/318
sqrt((sum((mean(negs2)-negs2)^2))/318-1)


dim(GC01)[2]
dim(PC19)[2]




##### ASV counts #####

#GC01
dim(GC01[rowSums(GC01)!=0,])[1]
#PC19
dim(PC19[rowSums(PC19)!=0,])[1]
#BOTH
length(unique(c(row.names(GC01[rowSums(GC01)!=0,]),row.names(GC01[rowSums(PC19)!=0,]))))
# overlap 
table(rownames(GC01[rowSums(GC01)!=0,]) %in% rownames(PC19[rowSums(PC19)!=0,]))



##### Taxonomic Assignments #####

tax <- read.csv("taxonomy/EUK.cleaned.PR2.csv",row.names = 1)
tax.subset <- tax[match(unique(c(row.names(GC01[rowSums(GC01)!=0,]),row.names(GC01[rowSums(PC19)!=0,]))),tax$X.1),]
tax.subset.euk <- tax.subset[tax.subset$Domain=="Eukaryota",]

#taxonomic domain 
sum(table(tax.subset$Domain[tax.subset$Domain.1>80]))

#Class 
length(na.omit(unique(tax.subset.euk$Class[tax.subset.euk$Class.1>80])))
#Order
length(na.omit(unique(tax.subset.euk$Order[tax.subset.euk$Order.1>80])))
#Family
length(na.omit(unique(tax.subset.euk$Family[tax.subset.euk$Family.1>80])))
#Genus
length(na.omit(unique(tax.subset.euk$Genus[tax.subset.euk$Genus.1>80])))

### NCBI assignments