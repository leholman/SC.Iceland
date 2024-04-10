### Here we are prepping data for a quick first look ###

library("seqinr")
library(dada2)


# bring in metadata
ages<-read.csv("metadata/AgeOut.csv") 
metadata.B05 <- read.csv("metadata/B05-2006.metadata.csv")
metadata.LN1 <- read.csv("metadata/LN1.metadata.csv")
metadata.LN2 <- read.csv("metadata/LN2.metadata.csv")

#metadata$rep <- gsub(".*-([0-9])$","\\1",metadata$SampleID)


## Cleaning data 

# separate positive / negatives


## some cleaning functions

minimumReads <- function(data,minreads) {
  # Ensure the input is either a dataframe or matrix and is numeric
  if (!is.data.frame(data) && !is.matrix(data)) {stop("Input must be a dataframe or matrix")}
  if (!all(sapply(data, is.numeric))) {stop("Error: All values must be numeric or integer.")}
  
  # Set data below mininmum as zero
  data[] <- lapply(data, function(x) ifelse(x <= minreads, 0, x))

  # Get rid of empty rows
  outdata <- data[rowSums(data) > 0,]
  
  # Return data 
  return(outdata)
}



minimumReps <- function(data,min_obs) {
  # Check if data is a dataframe or matrix
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Input must be a dataframe or matrix.")
  }
  
  # Check for NA values
  if (any(is.na(data))) {
    stop("Data contains NA values, which are not allowed.")
  }
  
  # Determine the input type to return the same type
  input_type <- if (is.data.frame(data)) "data.frame" else "matrix"
  
  # Calculate the number of non-zero observations per row
  total_obs <- apply(data, 1, function(row) sum(row > 0, na.rm = TRUE))
  
  # Filter rows based on the user-specified minimum number of non-zero observations
  filtered_data <- data[total_obs >= min_obs, ]

  # Ensure the output is of the same type as the input
  if (input_type == "matrix") {
    return(as.matrix(filtered_data))
  } else {
    return(as.data.frame(filtered_data))
  }
  
}

dataCleanBy <- function(inputdata, cleaningdata, method) {
  # Check if inputs are matrices or dataframes and convert dataframes to matrices
  if (is.data.frame(inputdata)) inputdata <- as.matrix(inputdata)
  if (is.data.frame(cleaningdata)) cleaningdata <- as.matrix(cleaningdata)
  
  # Check both contain only numeric values
  if (!all(sapply(inputdata, is.numeric)) || !all(sapply(cleaningdata, is.numeric))) {
    stop("Both inputs must contain only numeric values.")
  }
  
  # Check they have the same number of rows
  if (nrow(inputdata) != nrow(cleaningdata)) {
    stop("Both inputs must have the same number of rows.")
  }
  
  # Calculate metric for each row in the second matrix/dataframe
  metric_values <- apply(cleaningdata, 1, function(row) {
    if (method == "max") {
      return(max(row))
    } else if (method == "min") {
      return(min(row))
    } else if (method == "avr") {
      return(mean(row))
    } else {
      stop("Method must be 'max', 'min', or 'avr'.")
    }
  })
  
  # Update the first dataset based on the metric
  for (i in 1:nrow(inputdata)) {
    inputdata[i, inputdata[i, ] >= metric_values[i]] <- 0
  }
  
  # Remove rows with only zero values
  inputdata <- inputdata[rowSums(inputdata != 0) > 0, ]
  
  return(inputdata)
}





for (dataset in list.files("rawdata",pattern="....raw.csv")){
  
  datasetname <- substr(dataset,1,3)  
  
  indata <- read.csv(paste0("rawdata/",datasetname,".raw.csv"),row.names = 1)
  
  expSamples <- indata[,na.omit(match(gsub(",|-",".",sort(metadata$SampleID[metadata$SampleType=="Experimental"])),substring(colnames(indata),5)))]
  
  ctlSamples <-indata[,na.omit(match(gsub(",|-",".",sort(metadata$SampleID[metadata$SampleType=="ControlN"])),substring(colnames(indata),5)))]
  
  #Filter 1 - minimum number of reads for any ID
  expSamples[expSamples< minreads] <- 0
  expSamples <- expSamples[rowSums(expSamples) > 0,]
  
  
  
  #Filter 2 - within samples OTU must appear in more than one sample (this works because there are lots of reps per site and sample)
  filtersam <- expSamples
  filtersam[filtersam>0 ] <- 1
  filtersam <-filtersam[rowSums(filtersam) > 1,]
  expSamples <- expSamples[rownames(expSamples) %in% rownames(filtersam),]
  
  #Filter 3 -Maximum value in neg = 0 value in samples
  
  controlsCONTAM <- ctlSamples[rowSums(ctlSamples) > 0,]
  for (contamOTU in 1:length(controlsCONTAM[,1])){
    loopOTU <- row.names(controlsCONTAM[contamOTU,])
    loopMax <- max(as.numeric(controlsCONTAM[contamOTU,]))
    #loopSum <- sum(as.numeric(controlsCONTAM[contamOTU,]))
    if (any(is.na(expSamples[loopOTU,]))){next}
    expSamples[loopOTU,expSamples[loopOTU,]<loopMax] <- 0
    print(paste("Cleaning contaminants",contamOTU))
  }
  
  















#### CODE basement 

## making test datasets to look at 

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



