### Here we are prepping data for a quick first look ###

library("seqinr")
library(dada2)
library(Biostrings)

set.seed(1234)


# bring in metadata
ages<-read.csv("metadata/AgeOut.csv") 
metadata.B05 <- read.csv("metadata/B05-2006.metadata.csv")
metadata.LN1 <- read.csv("metadata/LN1.metadata.csv")
metadata.LN2 <- read.csv("metadata/LN2.metadata.csv")

## some cleaning functions
minimumReads <- function(data, minreads) {
  # Ensure the input is either a dataframe or matrix and is numeric
  if (!is.data.frame(data) && !is.matrix(data)) {stop("Input must be a dataframe or matrix")}
  if (!all(sapply(data, is.numeric))) {stop("Error: All values must be numeric or integer.")}
  # Check for NA values
  if (any(is.na(data))) {stop("Data contains NA values, which are not allowed.")}
  
  # Store original row names
  original_row_names <- row.names(data)
  
  # Apply thresholding operation based on the structure of the data
  if (is.matrix(data)) {
    # For matrices, use vectorized operation
    data[data <= minreads] <- 0
  } else if (is.data.frame(data)) {
    # For data frames, apply operation column-wise
    data <- as.data.frame(lapply(data, function(x) ifelse(x <= minreads, 0, x)))
  }
  
  
  # Get rid of empty rows
  # For both data frames and matrices, rowSums works and subsetting by rowSums > 0 removes empty rows
  outdata <- data[rowSums(data > 0) > 0, ]
  row.names(outdata) <- original_row_names[rowSums(data > 0) > 0]
  
  
  # Return the cleaned data
  return(outdata)
}
minimumReps <- function(data,min_obs) {
  # Ensure the input is either a dataframe or matrix and is numeric
  if (!is.data.frame(data) && !is.matrix(data)) {stop("Input must be a dataframe or matrix")}
  if (!all(sapply(data, is.numeric))) {stop("Error: All values must be numeric or integer.")}
  # Check for NA values
  if (any(is.na(data))) {stop("Data contains NA values, which are not allowed.")}
  
  # Store original row names
  original_row_names <- row.names(data)
  
  # Determine the input type to return the same type
  input_type <- ifelse(is.data.frame(data), "data.frame", "matrix")
  
  # Calculate the number of non-zero observations per row
  total_obs <- apply(data, 1, function(row) sum(row > 0, na.rm = TRUE))
 
   # Filter rows based on the user-specified minimum number of non-zero observations
  filtered_data <- data[total_obs >= min_obs, ]

  # Reapply the original row names to the filtered data
  row.names(filtered_data) <- original_row_names[total_obs >= min_obs]
  
  # Ensure the output is of the same type as the input
  if (input_type == "matrix") {
    return(as.matrix(filtered_data))
  } else {
    return(as.data.frame(filtered_data))
  }
}
dataCleanBy <- function(inputdata, cleaningdata, method) {
  # Ensure the input is either a dataframe or matrix and is numeric
  if (!is.data.frame(inputdata) && !is.matrix(inputdata)) {stop("Input must be a dataframe or matrix")}
  if (!all(sapply(inputdata, is.numeric))) {stop("Error: All values must be numeric or integer.")}
  # Check for NA values
  if (any(is.na(inputdata))) {stop("Data contains NA values, which are not allowed.")}
  # Validate method parameter
  if (!method %in% c("max", "min", "avr")) {
    stop("Method must be 'max', 'min', or 'avr'.")
  }
  
  # Preserve original format of inputdata
  inputIsDataFrame <- is.data.frame(inputdata)
  
  # Convert dataframes to matrices for uniform processing
  inputdata <- as.matrix(inputdata)
  cleaningdata <- as.matrix(cleaningdata)
  
  # Ensure both inputs contain only numeric values
  if (!all(is.numeric(inputdata)) || !all(is.numeric(cleaningdata))) {
    stop("Both inputs must contain only numeric values.")
  }
  
  # Ensure both inputs have the same number of rows
  if (nrow(inputdata) != nrow(cleaningdata)) {
    stop("Both inputs must have the same number of rows.")
  }
  
  # Calculate a metric for each row in cleaningdata based on the specified method
  metric_values <- switch(method,
                          "max" = apply(cleaningdata, 1, max),
                          "min" = apply(cleaningdata, 1, min),
                          "avr" = apply(cleaningdata, 1, mean))
  
  # Vectorized update of inputdata based on the calculated metrics
  metric_matrix <- matrix(rep(metric_values, each = ncol(inputdata)), 
                          nrow = nrow(inputdata), 
                          byrow = TRUE)
  inputdata[inputdata <= metric_matrix] <- 0
  
  # Remove rows in inputdata with only zero values
  inputdata2 <- inputdata[rowSums(inputdata != 0) > 0, ]
  message(paste0("Removed ",dim(inputdata)[1]-dim(inputdata2)[1]," ASV/OTUs from samples"))
  
  # Return inputdata in its original format
  if (inputIsDataFrame) {
    return(as.data.frame(inputdata2))
  } else {
    return(inputdata2)
  }
}

### Global settings
min_reads <- 2
min_reps <- 2
comp_method <- "avr"

## Lets start with B05 
# EUK
B05.EUK <- read.csv("rawdata/EUK.B05.lulu.csv")
B05.EUK.asv <- read.fasta("rawdata/ASVs/EUK.B05-2006.DADA2.ASVs.fasta",as.string = TRUE)
B05.EUK.tax <- read.csv("taxonomy/EUK.B05-2006.parsed.csv",row.names = 1)
# set aside negative and experimental samples 
neg <- na.omit(colnames(B05.EUK)[metadata.B05$SampleType[match(gsub("EUK.","",colnames(B05.EUK)), gsub("-",".",metadata.B05$UniqueID))] == "ControlN"])
exp <- na.omit(colnames(B05.EUK)[metadata.B05$SampleType[match(gsub("EUK.","",colnames(B05.EUK)), gsub("-",".",metadata.B05$UniqueID))] == "Experimental"])
# create data subsets
B05.EUK.neg <- B05.EUK[,match(neg,colnames(B05.EUK))]
B05.EUK.exp <- B05.EUK[,match(exp,colnames(B05.EUK))]
# clean the data using cleaning functions 
B05.EUK.exp.1 <- dataCleanBy(B05.EUK.exp,B05.EUK.neg,method = "avr")
B05.EUK.exp.2 <- minimumReads(B05.EUK.exp.1,min_reads)
B05.EUK.exp.3 <- minimumReps(B05.EUK.exp.2,min_reps)

B05.EUK.neg.1 <- minimumReads(B05.EUK.neg,1)
#write out the data
write.csv(file="cleaneddata/B05.EUK.csv",cbind(unlist(B05.EUK.asv)[match(rownames(B05.EUK.exp.3),names(unlist(B05.EUK.asv)))],
                                               B05.EUK.exp.3))
write.csv(file="cleaneddata/negativedata/neg.B05.EUK.csv",cbind(unlist(B05.EUK.asv)[match(rownames(B05.EUK.neg.1),names(unlist(B05.EUK.asv)))],
                                               B05.EUK.neg.1,
                                               B05.EUK.tax[match(rownames(B05.EUK.neg.1),B05.EUK.tax$OTU),]))

# MAM
B05.MAM <- read.csv("rawdata/MAM.B05-2006.raw.names.csv")
B05.MAM.asv <- read.fasta("rawdata/ASVs/MAM.B05-2006.DADA2.ASVs.fasta",as.string = TRUE)
B05.MAM.tax <- read.csv("taxonomy/MAM.B05-2006.parsed.csv",row.names = 1)
# set aside negative and experimental samples 
neg <- na.omit(colnames(B05.MAM)[metadata.B05$SampleType[match(gsub("MAM.","",colnames(B05.MAM)), gsub("-",".",metadata.B05$UniqueID))] == "ControlN"])
exp <- na.omit(colnames(B05.MAM)[metadata.B05$SampleType[match(gsub("MAM.","",colnames(B05.MAM)), gsub("-",".",metadata.B05$UniqueID))] == "Experimental"])
# create data subsets
B05.MAM.neg <- B05.MAM[,match(neg,colnames(B05.MAM))]
B05.MAM.exp <- B05.MAM[,match(exp,colnames(B05.MAM))]
# clean the data using cleaning functions 
B05.MAM.exp.1 <- dataCleanBy(B05.MAM.exp,B05.MAM.neg,method = "avr")
B05.MAM.exp.2 <- minimumReads(B05.MAM.exp.1,min_reads)
B05.MAM.exp.3 <- minimumReps(B05.MAM.exp.2,min_reps)

B05.MAM.neg.1 <- minimumReads(B05.MAM.neg,1)
#write out the data
write.csv(file="cleaneddata/B05.MAM.csv",cbind(unlist(B05.MAM.asv)[match(rownames(B05.MAM.exp.3),names(unlist(B05.MAM.asv)))],
                                               B05.MAM.exp.3))
write.csv(file="cleaneddata/negativedata/neg.B05.MAM.csv",cbind(unlist(B05.MAM.asv)[match(rownames(B05.MAM.neg.1),names(unlist(B05.MAM.asv)))],
                                                                B05.MAM.neg.1,
                                                                B05.MAM.tax[match(rownames(B05.MAM.neg.1),B05.MAM.tax$OTU),]))

# RIZ
B05.RIZ <- read.csv("rawdata/RIZ.B05-2006.raw.names.csv")
B05.RIZ.asv <- read.fasta("rawdata/ASVs/RIZ.B05-2006.DADA2.ASVs.fasta",as.string = TRUE)
B05.RIZ.tax <- read.csv("taxonomy/RIZ.B05-2006.parsed.csv",row.names = 1)
# set aside negative and experimental samples 
neg <- na.omit(colnames(B05.RIZ)[metadata.B05$SampleType[match(gsub("RIZ.","",colnames(B05.RIZ)), gsub("-",".",metadata.B05$UniqueID))] == "ControlN"])
exp <- na.omit(colnames(B05.RIZ)[metadata.B05$SampleType[match(gsub("RIZ.","",colnames(B05.RIZ)), gsub("-",".",metadata.B05$UniqueID))] == "Experimental"])
# create data subsets
B05.RIZ.neg <- B05.RIZ[,match(neg,colnames(B05.RIZ))]
B05.RIZ.exp <- B05.RIZ[,match(exp,colnames(B05.RIZ))]
# clean the data using cleaning functions 
B05.RIZ.exp.1 <- dataCleanBy(B05.RIZ.exp,B05.RIZ.neg,method = "avr")
B05.RIZ.exp.2 <- minimumReads(B05.RIZ.exp.1,min_reads)
B05.RIZ.exp.3 <- minimumReps(B05.RIZ.exp.2,min_reps)

B05.RIZ.neg.1 <- minimumReads(B05.RIZ.neg,1)
#write out the data
write.csv(file="cleaneddata/B05.RIZ.csv",cbind(unlist(B05.RIZ.asv)[match(rownames(B05.RIZ.exp.3),names(unlist(B05.RIZ.asv)))],
                                               B05.RIZ.exp.3))
write.csv(file="cleaneddata/negativedata/neg.B05.RIZ.csv",cbind(unlist(B05.RIZ.asv)[match(rownames(B05.RIZ.neg.1),names(unlist(B05.RIZ.asv)))],
                                                                B05.RIZ.neg.1,
                                                                B05.RIZ.tax[match(rownames(B05.RIZ.neg.1),B05.RIZ.tax$OTU),]))


## Now LN1
# EUK
LN1.EUK <- read.csv("rawdata/EUK.LN1.lulu.csv")
LN1.EUK.asv <- read.fasta("rawdata/ASVs/EUK.LN1.DADA2.ASVs.fasta",as.string = TRUE)
LN1.EUK.tax <- read.csv("taxonomy/EUK.LN1.parsed.csv",row.names = 1)
# set aside negative and experimental samples 
neg <- na.omit(colnames(LN1.EUK)[metadata.LN1$SampleType[match(gsub("\\.","_",colnames(LN1.EUK)), gsub("-","_",gsub("_[0-9]-[0-9]","",metadata.LN1$UniqueID)))] == "ControlN"])
exp <- na.omit(colnames(LN1.EUK)[metadata.LN1$SampleType[match(gsub("\\.","_",colnames(LN1.EUK)), gsub("-","_",gsub("_[0-9]-[0-9]","",metadata.LN1$UniqueID)))] == "Experimental"])
# create data subsets
LN1.EUK.neg <- LN1.EUK[,match(neg,colnames(LN1.EUK))]
LN1.EUK.exp <- LN1.EUK[,match(exp,colnames(LN1.EUK))]
# clean the data using cleaning functions 
LN1.EUK.exp.1 <- dataCleanBy(LN1.EUK.exp,LN1.EUK.neg,method = "avr")
LN1.EUK.exp.2 <- minimumReads(LN1.EUK.exp.1,min_reads)
LN1.EUK.exp.3 <- minimumReps(LN1.EUK.exp.2,min_reps)

LN1.EUK.neg.1 <- minimumReads(LN1.EUK.neg,1)
#write out the data
write.csv(file="cleaneddata/LN1.EUK.csv",cbind(unlist(LN1.EUK.asv)[match(rownames(LN1.EUK.exp.3),names(unlist(LN1.EUK.asv)))],
                                               LN1.EUK.exp.3))
write.csv(file="cleaneddata/negativedata/neg.LN1.EUK.csv",cbind(unlist(LN1.EUK.asv)[match(rownames(LN1.EUK.neg.1),names(unlist(LN1.EUK.asv)))],
                                                                LN1.EUK.neg.1,
                                                                LN1.EUK.tax[match(rownames(LN1.EUK.neg.1),LN1.EUK.tax$OTU),]))

# MAM
LN1.MAM <- read.csv("rawdata/MAM.LN1.raw.names.csv")
LN1.MAM.asv <- read.fasta("rawdata/ASVs/MAM.LN1.DADA2.ASVs.fasta",as.string = TRUE)
LN1.MAM.tax <- read.csv("taxonomy/MAM.LN1.parsed.csv",row.names = 1)
# set aside negative and experimental samples 
neg <- na.omit(colnames(LN1.MAM)[metadata.LN1$SampleType[match(gsub("\\.","_",colnames(LN1.MAM)), gsub("-","_",gsub("_[0-9]-[0-9]","",metadata.LN1$UniqueID)))] == "ControlN"])
exp <- na.omit(colnames(LN1.MAM)[metadata.LN1$SampleType[match(gsub("\\.","_",colnames(LN1.MAM)), gsub("-","_",gsub("_[0-9]-[0-9]","",metadata.LN1$UniqueID)))] == "Experimental"])
# create data subsets
LN1.MAM.neg <- LN1.MAM[,match(neg,colnames(LN1.MAM))]
LN1.MAM.exp <- LN1.MAM[,match(exp,colnames(LN1.MAM))]
# clean the data using cleaning functions 
LN1.MAM.exp.1 <- dataCleanBy(LN1.MAM.exp,LN1.MAM.neg,method = "avr")
LN1.MAM.exp.2 <- minimumReads(LN1.MAM.exp.1,min_reads)
LN1.MAM.exp.3 <- minimumReps(LN1.MAM.exp.2,min_reps)

LN1.MAM.neg.1 <- minimumReads(LN1.MAM.neg,1)
#write out the data
write.csv(file="cleaneddata/LN1.MAM.csv",cbind(unlist(LN1.MAM.asv)[match(rownames(LN1.MAM.exp.3),names(unlist(LN1.MAM.asv)))],
                                               LN1.MAM.exp.3))
write.csv(file="cleaneddata/negativedata/neg.LN1.MAM.csv",cbind(unlist(LN1.MAM.asv)[match(rownames(LN1.MAM.neg.1),names(unlist(LN1.MAM.asv)))],
                                                                LN1.MAM.neg.1,
                                                                LN1.MAM.tax[match(rownames(LN1.MAM.neg.1),LN1.MAM.tax$OTU),]))

# RIZ
LN1.RIZ <- read.csv("rawdata/RIZ.LN1.raw.names.csv")
LN1.RIZ.asv <- read.fasta("rawdata/ASVs/RIZ.LN1.DADA2.ASVs.fasta",as.string = TRUE)
LN1.RIZ.tax <- read.csv("taxonomy/RIZ.LN1.parsed.csv",row.names = 1)
# set aside negative and experimental samples 
neg <- na.omit(colnames(LN1.RIZ)[metadata.LN1$SampleType[match(gsub("\\.","_",colnames(LN1.RIZ)), gsub("-","_",gsub("_[0-9]-[0-9]","",metadata.LN1$UniqueID)))] == "ControlN"])
exp <- na.omit(colnames(LN1.RIZ)[metadata.LN1$SampleType[match(gsub("\\.","_",colnames(LN1.RIZ)), gsub("-","_",gsub("_[0-9]-[0-9]","",metadata.LN1$UniqueID)))] == "Experimental"])
# create data subsets
LN1.RIZ.neg <- LN1.RIZ[,match(neg,colnames(LN1.RIZ))]
LN1.RIZ.exp <- LN1.RIZ[,match(exp,colnames(LN1.RIZ))]
# clean the data using cleaning functions 
LN1.RIZ.exp.1 <- dataCleanBy(LN1.RIZ.exp,LN1.RIZ.neg,method = "avr")
LN1.RIZ.exp.2 <- minimumReads(LN1.RIZ.exp.1,min_reads)
LN1.RIZ.exp.3 <- minimumReps(LN1.RIZ.exp.2,min_reps)

LN1.RIZ.neg.1 <- minimumReads(LN1.RIZ.neg,1)
#write out the data
write.csv(file="cleaneddata/LN1.RIZ.csv",cbind(unlist(LN1.RIZ.asv)[match(rownames(LN1.RIZ.exp.3),names(unlist(LN1.RIZ.asv)))],
                                               LN1.RIZ.exp.3))
write.csv(file="cleaneddata/negativedata/neg.LN1.RIZ.csv",cbind(unlist(LN1.RIZ.asv)[match(rownames(LN1.RIZ.neg.1),names(unlist(LN1.RIZ.asv)))],
                                                                LN1.RIZ.neg.1,
                                                                LN1.RIZ.tax[match(rownames(LN1.RIZ.neg.1),LN1.RIZ.tax$OTU),]))

## Now LN2
# EUK
LN2.EUK <- read.csv("rawdata/EUK.LN2.lulu.csv")
LN2.EUK.asv <- read.fasta("rawdata/ASVs/EUK.LN2.DADA2.ASVs.fasta",as.string = TRUE)
LN2.EUK.tax <- read.csv("taxonomy/EUK.LN2.parsed.csv",row.names = 1)
# set aside negative and experimental samples 
neg <- na.omit(colnames(LN2.EUK)[metadata.LN2$SampleType[match(gsub("\\.","_",colnames(LN2.EUK)), gsub("-","_",gsub("_[0-9]-[0-9]|_B-[0-9]","",metadata.LN2$UniqueID)))] == "ControlN"])
exp <- na.omit(colnames(LN2.EUK)[metadata.LN2$SampleType[match(gsub("\\.","_",colnames(LN2.EUK)), gsub("-","_",gsub("_[0-9]-[0-9]|_B-[0-9]","",metadata.LN2$UniqueID)))] == "Experimental"])
# create data subsets
LN2.EUK.neg <- LN2.EUK[,match(neg,colnames(LN2.EUK))]
LN2.EUK.exp <- LN2.EUK[,match(exp,colnames(LN2.EUK))]
# clean the data using cleaning functions 
LN2.EUK.exp.1 <- dataCleanBy(LN2.EUK.exp,LN2.EUK.neg,method = "avr")
LN2.EUK.exp.2 <- minimumReads(LN2.EUK.exp.1,min_reads)
LN2.EUK.exp.3 <- minimumReps(LN2.EUK.exp.2,min_reps)

LN2.EUK.neg.1 <- minimumReads(LN2.EUK.neg,1)
#write out the data
write.csv(file="cleaneddata/LN2.EUK.csv",cbind(unlist(LN2.EUK.asv)[match(rownames(LN2.EUK.exp.3),names(unlist(LN2.EUK.asv)))],
                                               LN2.EUK.exp.3))
write.csv(file="cleaneddata/negativedata/neg.LN2.EUK.csv",cbind(unlist(LN2.EUK.asv)[match(rownames(LN2.EUK.neg.1),names(unlist(LN2.EUK.asv)))],
                                                                LN2.EUK.neg.1,
                                                                LN2.EUK.tax[match(rownames(LN2.EUK.neg.1),LN2.EUK.tax$OTU),]))

# MAM
LN2.MAM <- read.csv("rawdata/MAM.LN2.raw.names.csv")
LN2.MAM.asv <- read.fasta("rawdata/ASVs/MAM.LN2.DADA2.ASVs.fasta",as.string = TRUE)
LN2.MAM.tax <- read.csv("taxonomy/MAM.LN2.parsed.csv",row.names = 1)
# set aside negative and experimental samples 
neg <- na.omit(colnames(LN2.MAM)[metadata.LN2$SampleType[match(gsub("\\.","_",colnames(LN2.MAM)), gsub("-","_",gsub("_[0-9]-[0-9]|_B-[0-9]","",metadata.LN2$UniqueID)))] == "ControlN"])
exp <- na.omit(colnames(LN2.MAM)[metadata.LN2$SampleType[match(gsub("\\.","_",colnames(LN2.MAM)), gsub("-","_",gsub("_[0-9]-[0-9]|_B-[0-9]","",metadata.LN2$UniqueID)))] == "Experimental"])
# create data subsets
LN2.MAM.neg <- LN2.MAM[,match(neg,colnames(LN2.MAM))]
LN2.MAM.exp <- LN2.MAM[,match(exp,colnames(LN2.MAM))]
# clean the data using cleaning functions 
LN2.MAM.exp.1 <- dataCleanBy(LN2.MAM.exp,LN2.MAM.neg,method = "avr")
LN2.MAM.exp.2 <- minimumReads(LN2.MAM.exp.1,min_reads)
LN2.MAM.exp.3 <- minimumReps(LN2.MAM.exp.2,min_reps)

LN2.MAM.neg.1 <- minimumReads(LN2.MAM.neg,1)
#write out the data
write.csv(file="cleaneddata/LN2.MAM.csv",cbind(unlist(LN2.MAM.asv)[match(rownames(LN2.MAM.exp.3),names(unlist(LN2.MAM.asv)))],
                                               LN2.MAM.exp.3))
write.csv(file="cleaneddata/negativedata/neg.LN2.MAM.csv",cbind(unlist(LN2.MAM.asv)[match(rownames(LN2.MAM.neg.1),names(unlist(LN2.MAM.asv)))],
                                                                LN2.MAM.neg.1,
                                                                LN2.MAM.tax[match(rownames(LN2.MAM.neg.1),LN2.MAM.tax$OTU),]))

# RIZ
LN2.RIZ <- read.csv("rawdata/RIZ.LN2.raw.names.csv")
LN2.RIZ.asv <- read.fasta("rawdata/ASVs/RIZ.LN2.DADA2.ASVs.fasta",as.string = TRUE)
LN2.RIZ.tax <- read.csv("taxonomy/RIZ.LN2.parsed.csv",row.names = 1)
# set aside negative and experimental samples 
neg <- na.omit(colnames(LN2.RIZ)[metadata.LN2$SampleType[match(gsub("\\.","_",colnames(LN2.RIZ)), gsub("-","_",gsub("_[0-9]-[0-9]|_B-[0-9]","",metadata.LN2$UniqueID)))] == "ControlN"])
exp <- na.omit(colnames(LN2.RIZ)[metadata.LN2$SampleType[match(gsub("\\.","_",colnames(LN2.RIZ)), gsub("-","_",gsub("_[0-9]-[0-9]|_B-[0-9]","",metadata.LN2$UniqueID)))] == "Experimental"])
# create data subsets
LN2.RIZ.neg <- LN2.RIZ[,match(neg,colnames(LN2.RIZ))]
LN2.RIZ.exp <- LN2.RIZ[,match(exp,colnames(LN2.RIZ))]
# clean the data using cleaning functions 
LN2.RIZ.exp.1 <- dataCleanBy(LN2.RIZ.exp,LN2.RIZ.neg,method = "avr")
LN2.RIZ.exp.2 <- minimumReads(LN2.RIZ.exp.1,min_reads)
LN2.RIZ.exp.3 <- minimumReps(LN2.RIZ.exp.2,min_reps)

LN2.RIZ.neg.1 <- minimumReads(LN2.RIZ.neg,1)
#write out the data
write.csv(file="cleaneddata/LN2.RIZ.csv",cbind(unlist(LN2.RIZ.asv)[match(rownames(LN2.RIZ.exp.3),names(unlist(LN2.RIZ.asv)))],
                                               LN2.RIZ.exp.3))
write.csv(file="cleaneddata/negativedata/neg.LN2.RIZ.csv",cbind(unlist(LN2.RIZ.asv)[match(rownames(LN2.RIZ.neg.1),names(unlist(LN2.RIZ.asv)))],
                                                                LN2.RIZ.neg.1,
                                                                LN2.RIZ.tax[match(rownames(LN2.RIZ.neg.1),LN2.RIZ.tax$OTU),]))





## Now lets combine into cores 
EUK.B05.clean <- read.csv("cleaneddata/B05.EUK.csv",row.names = 2)[,-1]
colnames(EUK.B05.clean) <- gsub("\\.", "_", gsub("^EUK\\.IS\\.(GC[0-9]+\\.[0-9]+\\.[0-9]+)\\..*$", "\\1", colnames(EUK.B05.clean))) 
EUK.LN1.clean <- read.csv("cleaneddata/LN1.EUK.csv",row.names = 2)[,-1]
colnames(EUK.LN1.clean)  <- gsub("\\.","_",colnames(EUK.LN1.clean))
EUK.LN2.clean <- read.csv("cleaneddata/LN2.EUK.csv",row.names = 2)[,-1]
colnames(EUK.LN2.clean)  <- gsub("\\.","_",colnames(EUK.LN2.clean))

#combine all tables
EUK.combined <- as.data.frame(t(mergeSequenceTables(t(as.matrix(EUK.B05.clean)),
                    t(as.matrix(EUK.LN1.clean)),
                    t(as.matrix(EUK.LN2.clean)),
                    tryRC = TRUE)))

#get rid of short and long seqs
EUK.combined.2 <- EUK.combined[nchar(row.names(EUK.combined))>70 & nchar(row.names(EUK.combined))<170,]
#output ASVs
ASVs <- DNAStringSet(rownames(EUK.combined.2))
names(ASVs) <- paste0("ASV_",1:length(ASVs))
writeXStringSet(ASVs,"cleaneddata/ASVs/EUK.cleaned.fasta", append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")
EUK.combined.3 <- EUK.combined.2[,order(colnames(EUK.combined.2))]
rownames(EUK.combined.3) <- paste0("ASV_",1:length(ASVs))

#reassign using PR2 and save output 
PR2assign <- assignTaxonomy( unlist(read.fasta("cleaneddata/ASVs/EUK.cleaned.fasta",as.string = T)),refFasta = "taxonomy/pr2_version_5.0.0_SSU_dada2.fasta.gz",tryRC = TRUE,multithread = TRUE,taxLevels = c("Domain","Supergroup","Division","Subdivision","Class","Order","Family","Genus","Species"),outputBootstraps = TRUE)
PR2master <- cbind(names(ASVs),PR2assign$tax,PR2assign$boot)
write.csv(PR2master,"taxonomy/EUK.cleaned.PR2.csv")

#subset by core
EUK.GC01 <- EUK.combined.3[,grep("GC1",colnames(EUK.combined.3))]
EUK.GC06 <- EUK.combined.3[,grep("GC6",colnames(EUK.combined.3))]
EUK.PC019 <- EUK.combined.3[,grep("PC019",colnames(EUK.combined.3))]
EUK.PC022 <- EUK.combined.3[,grep("PC022",colnames(EUK.combined.3))]

## write core data 
write.csv(EUK.GC01,"cleaneddata/combinedcoredata/EUK.GC01.csv")
write.csv(EUK.GC06,"cleaneddata/combinedcoredata/EUK.GC06.csv")
write.csv(EUK.PC019,"cleaneddata/combinedcoredata/EUK.PC019.csv")
write.csv(EUK.PC022,"cleaneddata/combinedcoredata/EUK.PC022.csv")



# now MAM


## Now lets combine into cores 
MAM.B05.clean <- read.csv("cleaneddata/B05.MAM.csv",row.names = 2)[,-1]
colnames(MAM.B05.clean) <- gsub("\\.", "_", gsub("^MAM\\.IS\\.(GC[0-9]+\\.[0-9]+\\.[0-9]+)\\..*$", "\\1", colnames(MAM.B05.clean))) 
MAM.LN1.clean <- read.csv("cleaneddata/LN1.MAM.csv",row.names = 2)[,-1]
colnames(MAM.LN1.clean)  <- gsub("\\.","_",colnames(MAM.LN1.clean))
MAM.LN2.clean <- read.csv("cleaneddata/LN2.MAM.csv",row.names = 2)[,-1]
colnames(MAM.LN2.clean)  <- gsub("\\.","_",colnames(MAM.LN2.clean))

#combine all tables
MAM.combined <- as.data.frame(t(mergeSequenceTables(t(as.matrix(MAM.B05.clean)),
                                                    t(as.matrix(MAM.LN1.clean)),
                                                    t(as.matrix(MAM.LN2.clean)),
                                                    tryRC = TRUE)))

#get rid of short and long seqs
MAM.combined.2 <- MAM.combined[nchar(row.names(MAM.combined))>80 & nchar(row.names(MAM.combined))<120,]
#output ASVs
ASVs <- DNAStringSet(rownames(MAM.combined.2))
names(ASVs) <- paste0("ASV_",1:length(ASVs))
writeXStringSet(ASVs,"cleaneddata/ASVs/MAM.cleaned.fasta", append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")
MAM.combined.3 <- MAM.combined.2[,order(colnames(MAM.combined.2))]
rownames(MAM.combined.3) <- paste0("ASV_",1:length(ASVs))

#subset by core
MAM.GC01 <- MAM.combined.3[,grep("GC1",colnames(MAM.combined.3))]
MAM.GC06 <- MAM.combined.3[,grep("GC6",colnames(MAM.combined.3))]
MAM.PC019 <- MAM.combined.3[,grep("PC019",colnames(MAM.combined.3))]
MAM.PC022 <- MAM.combined.3[,grep("PC022",colnames(MAM.combined.3))]

## write core data 
write.csv(MAM.GC01,"cleaneddata/combinedcoredata/MAM.GC01.csv")
write.csv(MAM.GC06,"cleaneddata/combinedcoredata/MAM.GC06.csv")
write.csv(MAM.PC019,"cleaneddata/combinedcoredata/MAM.PC019.csv")
write.csv(MAM.PC022,"cleaneddata/combinedcoredata/MAM.PC022.csv")



# now RIZ


## Now lets combine into cores 
RIZ.B05.clean <- read.csv("cleaneddata/B05.RIZ.csv",row.names = 2)[,-1]
colnames(RIZ.B05.clean) <- gsub("\\.", "_", gsub("^RIZ\\.IS\\.(GC[0-9]+\\.[0-9]+\\.[0-9]+)\\..*$", "\\1", colnames(RIZ.B05.clean))) 
RIZ.LN1.clean <- read.csv("cleaneddata/LN1.RIZ.csv",row.names = 2)[,-1]
colnames(RIZ.LN1.clean)  <- gsub("\\.","_",colnames(RIZ.LN1.clean))
RIZ.LN2.clean <- read.csv("cleaneddata/LN2.RIZ.csv",row.names = 2)[,-1]
colnames(RIZ.LN2.clean)  <- gsub("\\.","_",colnames(RIZ.LN2.clean))

#combine all tables
RIZ.combined <- as.data.frame(t(mergeSequenceTables(t(as.matrix(RIZ.B05.clean)),
                                                    t(as.matrix(RIZ.LN1.clean)),
                                                    t(as.matrix(RIZ.LN2.clean)),
                                                    tryRC = TRUE)))

#get rid of short and long seqs
RIZ.combined.2 <- RIZ.combined[nchar(row.names(RIZ.combined))>80 & nchar(row.names(RIZ.combined))<150,]
#output ASVs
ASVs <- DNAStringSet(rownames(RIZ.combined.2))
names(ASVs) <- paste0("ASV_",1:length(ASVs))
writeXStringSet(ASVs,"cleaneddata/ASVs/RIZ.cleaned.fasta", append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")
RIZ.combined.3 <- RIZ.combined.2[,order(colnames(RIZ.combined.2))]
rownames(RIZ.combined.3) <- paste0("ASV_",1:length(ASVs))

#subset by core
RIZ.GC01 <- RIZ.combined.3[,grep("GC1",colnames(RIZ.combined.3))]
RIZ.GC06 <- RIZ.combined.3[,grep("GC6",colnames(RIZ.combined.3))]
RIZ.PC019 <- RIZ.combined.3[,grep("PC019",colnames(RIZ.combined.3))]
RIZ.PC022 <- RIZ.combined.3[,grep("PC022",colnames(RIZ.combined.3))]

## write core data 
write.csv(RIZ.GC01,"cleaneddata/combinedcoredata/RIZ.GC01.csv")
write.csv(RIZ.GC06,"cleaneddata/combinedcoredata/RIZ.GC06.csv")
write.csv(RIZ.PC019,"cleaneddata/combinedcoredata/RIZ.PC019.csv")
write.csv(RIZ.PC022,"cleaneddata/combinedcoredata/RIZ.PC022.csv")




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



### function tests
library(MASS)
library(dplyr)



# Generate example dataset
rows <- 1000
columns <- 80
data_matrix <- matrix(rlnorm(rows * columns, meanlog = log(1), sdlog = 0.6), nrow = rows, ncol = columns)


# Normalize and scale the data to fit between 1 and 250000, peaking around 1-5
scaled_matrix <- (data_matrix / max(data_matrix)) * 250000

# Apply floor to avoid decimals and ensure the range starts at 1
scaled_matrix <- floor(scaled_matrix)
hist(scaled_matrix,breaks=1000)


generate_data_poisson <- function(rows, columns, lambda = 1, zero_inflation = 0.3) {
  # Generate data using a Poisson distribution
  data_matrix <- matrix(rpois(rows * columns, lambda = lambda), nrow = rows, ncol = columns)
  
  # Introduce zero inflation
  set.seed(42)  # Ensure reproducibility
  zero_indices <- sample(length(data_matrix), size = zero_inflation * length(data_matrix))
  data_matrix[zero_indices] <- 0
  
  # Scaling - a naive approach might not work as Poisson distribution might not give us the high values needed
  # Instead, we might simulate the higher values separately or adjust lambda dynamically (not typical for Poisson)
  
  # Directly modifying some values to simulate the extended range up to 250,000
  # This approach is somewhat artificial but demonstrates how you could inject larger values
  high_value_indices <- sample((1:length(data_matrix))[-zero_indices], size = 0.01 * length(data_matrix))
  data_matrix[high_value_indices] <- sample(200:250, length(high_value_indices), replace = TRUE) * 1000
  
  return(data_matrix)
}

# Generate example dataset
rows <- 1000
columns <- 80
example_dataset <- generate_data_poisson(rows, columns, lambda = 1.5, zero_inflation = 0.1)
example_dataset2 <- generate_data_poisson(rows, columns, lambda = 1.5, zero_inflation = 0.5)


example_dataset.c1 <- minimumReads(as.data.frame(example_dataset),minreads = 5)
example_dataset.c2 <- minimumReps(example_dataset.c1,min_obs = 60)

example_dataset.c3 <- dataCleanBy(example_dataset.c1,example_dataset2[,1:10],method = "max")





