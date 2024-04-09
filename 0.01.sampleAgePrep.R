##subsetting data 


lane1 <- readxl::read_excel("DataSheet.xlsx",sheet = "Lane1")
samplesGC1 <- sort(unique(lane1$SampleID))[grep(pattern = "GC1",sort(unique(lane1$SampleID)))]
samplesGC1 <- gsub("IS-","",samplesGC1)


laneA <- readxl::read_excel("DataSheet.xlsx",sheet = "LaneA")
samplesPC019<- sort(unique(laneA$SampleID))[grep(pattern = "PC019",sort(unique(laneA$SampleID)))]
samplesPC019 <- samplesPC019[-grep("Blank",samplesPC019)]
samplesPC019 <- sort(gsub("_","-",samplesPC019 ))


laneB <- readxl::read_excel("DataSheet.xlsx",sheet = "LaneB")
samplesPC019.2 <- sort(unique(laneB$SampleID))[grep(pattern = "PC019",sort(unique(laneB$SampleID)))]
samplesPC019.2 <- samplesPC019.2[-grep("Blank",samplesPC019.2)]
samplesPC019.2 <- sort(gsub("_","-",samplesPC019.2))
samplesPC019.2 <- samplesPC019.2[-grep("B6",samplesPC019.2)]

samplesPC019all <- sort(c(samplesPC019,samplesPC019.2))
unique(samplesPC019all)

samplesPC022 <- sort(unique(laneB$SampleID))[grep(pattern = "PC022",sort(unique(laneB$SampleID)))]
samplesPC022 <- samplesPC022[-grep("B[0-9]",samplesPC022)]
samplesPC022 <- sort(gsub("_","-",samplesPC022))

allSamples <- data.frame("ID"=c(samplesGC1,samplesPC019all,samplesPC022),
                         "Core"=c(sapply(strsplit(samplesGC1,"-"), `[`, 1),sapply(strsplit(samplesPC019all,"-"), `[`, 1),sapply(strsplit(samplesPC022,"-"), `[`, 1)),
                         "Depth"=as.numeric(c(sapply(strsplit(samplesGC1,"-"), `[`, 2),sapply(strsplit(samplesPC019all,"-"), `[`, 2),sapply(strsplit(samplesPC022,"-"), `[`, 2))))
write.csv(allSamples,"Age.csv")
