### Here we are prepping data for a quick first look ###

library("seqinr")
library(dada2)


test <- collapseNoMismatch(t(read.csv("~/Desktop/WP4.LN2/6.mappings/OTUtabs/EUK.LN2.raw.csv",row.names = 1)))

# bring in metadata
ages<-read.csv("metadata/AgeOut.csv") 
metadata.B05 <- read.csv("metadata/B05-2006.metadata.csv")
metadata.LN1 <- read.csv("metadata/LN1.metadata.csv")
metadata.LN2 <- read.csv("metadata/LN2.metadata.csv")

#metadata$rep <- gsub(".*-([0-9])$","\\1",metadata$SampleID)

#B05
B05.EUK.data<- read.csv("rawData/EUK.B05-2006.raw.names.csv",row.names = 1)
B05.EUK.tax <- read.csv("taxonomy/EUK.B05-2006.parsed.csv",row.names = 1)
B05.EUK.asvs <- read.fasta("rawdata/ASVs/EUK.B05-2006.DADA2.ASVs.fasta",as.string = TRUE)
write.csv(file="cleaneddata/testdata/B05.EUK.test.csv",cbind(B05.EUK.data,unlist(B05.EUK.asvs),B05.EUK.tax[match(rownames(B05.EUK.data),B05.EUK.tax$OTU),],nchar(B05.EUK.asvs)))
B05.RIZ.data<- read.csv("rawData/RIZ.B05-2006.raw.names.csv",row.names = 1)
B05.RIZ.tax <- read.csv("taxonomy/RIZ.B05-2006.parsed.csv",row.names = 1)
B05.RIZ.asvs <- read.fasta("rawdata/ASVs/RIZ.B05-2006.DADA2.ASVs.fasta",as.string = TRUE)
write.csv(file="cleaneddata/testdata/B05.RIZ.test.csv",cbind(B05.RIZ.data,unlist(B05.RIZ.asvs),B05.RIZ.tax[match(rownames(B05.RIZ.data),B05.RIZ.tax$OTU),],nchar(B05.RIZ.asvs)))
B05.MAM.data<- read.csv("rawData/MAM.B05-2006.raw.names.csv",row.names = 1)
B05.MAM.tax <- read.csv("taxonomy/MAM.B05-2006.parsed.csv",row.names = 1)
B05.MAM.asvs <- read.fasta("rawdata/ASVs/MAM.B05-2006.DADA2.ASVs.fasta",as.string = TRUE)
write.csv(file="cleaneddata/testdata/B05.MAM.test.csv",cbind(B05.MAM.data,unlist(B05.MAM.asvs),B05.MAM.tax[match(rownames(B05.MAM.data),B05.MAM.tax$OTU),],nchar(B05.MAM.asvs)))


#LN1
LN1.EUK.data<- read.csv("rawData/EUK.LN1.raw.names.csv",row.names = 1)
LN1.EUK.tax <- read.csv("taxonomy/EUK.LN1.parsed.csv",row.names = 1)
LN1.EUK.asvs <- read.fasta("rawdata/ASVs/EUK.LN1.DADA2.ASVs.fasta",as.string = TRUE)
write.csv(file="cleaneddata/testdata/LN1.EUK.test.csv",cbind(LN1.EUK.data,unlist(LN1.EUK.asvs),LN1.EUK.tax[match(rownames(LN1.EUK.data),LN1.EUK.tax$OTU),],nchar(LN1.EUK.asvs)))

LN1.EUK.lulu.data<- read.csv("rawData/EUK.LN1.lulu.csv",row.names = 1)
write.csv(file="cleaneddata/testdata/LN1.EUK.lulu.test.csv",cbind(LN1.EUK.lulu.data,
                                                             unlist(LN1.EUK.asvs)[match(rownames(LN1.EUK.lulu.data),names(LN1.EUK.asvs))],
                                                             LN1.EUK.tax[match(rownames(LN1.EUK.lulu.data),LN1.EUK.tax$OTU),],
                                                             nchar(LN1.EUK.asvs)[match(rownames(LN1.EUK.lulu.data),names(LN1.EUK.asvs))]))


LN1.RIZ.data<- read.csv("rawData/RIZ.LN1.raw.names.csv",row.names = 1)
LN1.RIZ.tax <- read.csv("taxonomy/RIZ.LN1.parsed.csv",row.names = 1)
LN1.RIZ.asvs <- read.fasta("rawdata/ASVs/RIZ.LN1.DADA2.ASVs.fasta",as.string = TRUE)
write.csv(file="cleaneddata/testdata/LN1.RIZ.test.csv",cbind(LN1.RIZ.data,unlist(LN1.RIZ.asvs),LN1.RIZ.tax[match(rownames(LN1.RIZ.data),LN1.RIZ.tax$OTU),],nchar(LN1.RIZ.asvs)))
LN1.MAM.data<- read.csv("rawData/MAM.LN1.raw.names.csv",row.names = 1)
LN1.MAM.tax <- read.csv("taxonomy/MAM.LN1.parsed.csv",row.names = 1)
LN1.MAM.asvs <- read.fasta("rawdata/ASVs/MAM.LN1.DADA2.ASVs.fasta",as.string = TRUE)
write.csv(file="cleaneddata/testdata/LN1.MAM.test.csv",cbind(LN1.MAM.data,unlist(LN1.MAM.asvs),LN1.MAM.tax[match(rownames(LN1.MAM.data),LN1.MAM.tax$OTU),],nchar(LN1.MAM.asvs)))

#LN2
LN2.EUK.data<- read.csv("rawData/EUK.LN2.raw.names.csv",row.names = 1)
LN2.EUK.tax <- read.csv("taxonomy/EUK.LN2.parsed.csv",row.names = 1)
LN2.EUK.asvs <- read.fasta("rawdata/ASVs/EUK.LN2.DADA2.ASVs.fasta",as.string = TRUE)
write.csv(file="cleaneddata/testdata/LN2.EUK.test.csv",cbind(LN2.EUK.data,unlist(LN2.EUK.asvs),LN2.EUK.tax[match(rownames(LN2.EUK.data),LN2.EUK.tax$OTU),],nchar(LN2.EUK.asvs)))
LN2.RIZ.data<- read.csv("rawData/RIZ.LN2.raw.names.csv",row.names = 1)
LN2.RIZ.tax <- read.csv("taxonomy/RIZ.LN2.parsed.csv",row.names = 1)
LN2.RIZ.asvs <- read.fasta("rawdata/ASVs/RIZ.LN2.DADA2.ASVs.fasta",as.string = TRUE)
write.csv(file="cleaneddata/testdata/LN2.RIZ.test.csv",cbind(LN2.RIZ.data,unlist(LN2.RIZ.asvs),LN2.RIZ.tax[match(rownames(LN2.RIZ.data),LN2.RIZ.tax$OTU),],nchar(LN2.RIZ.asvs)))
LN2.MAM.data<- read.csv("rawData/MAM.LN2.raw.names.csv",row.names = 1)
LN2.MAM.tax <- read.csv("taxonomy/MAM.LN2.parsed.csv",row.names = 1)
LN2.MAM.asvs <- read.fasta("rawdata/ASVs/MAM.LN2.DADA2.ASVs.fasta",as.string = TRUE)
write.csv(file="cleaneddata/testdata/LN2.MAM.test.csv",cbind(LN2.MAM.data,unlist(LN2.MAM.asvs),LN2.MAM.tax[match(rownames(LN2.MAM.data),LN2.MAM.tax$OTU),],nchar(LN2.MAM.asvs)))



