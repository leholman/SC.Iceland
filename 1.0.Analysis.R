################################################################
####======SeaChange Icelandic Dataset Analysis ==========#######
####==== Luke E. Holman====15.04.2024====================#######
################################################################

####====0.1 Packages====####
library(vegan)
library(RColorBrewer)
library(breakaway)
library("mgcv")
library("Biostrings")
library("seqinr")

#### WAVELET STUFF
#xwt {biwavelet}





####====0.2 Functions====####
NrepsMaker <- function(INdataframe,vector){
  ##write these checks
  #check the dataframe is a dataframe
  if(!is.data.frame(INdataframe)){stop("Input dataframe doesn't look like a dataframe")}
  #check the vector is a vector
  if(!is.vector(vector)){stop("Input vector doesn't look like a vector")}
  #check the dataframe contains the vector
  ## TO DO
  #make a new dataframe to captuire the output
  newDataFrame <- data.frame(matrix(0,nrow = length(INdataframe[,1]),ncol=length(unique(vector))))
  #name stuff
  colnames(newDataFrame) <- unique(vector)
  rownames(newDataFrame) <- rownames(INdataframe)
  #make it binary
  INdataframe[INdataframe<1] <- 0
  INdataframe[INdataframe>0] <- 1
  #loop over all the samples with replicates, summing according to the vector
  for (column in 1:length(newDataFrame[1,])){
    # this if statement checks in case there is only one replicate remaining (this sometimes happens in controls)
    if(is.vector(INdataframe[,vector %in% colnames(newDataFrame)[column]])){
      newDataFrame[,column]  <- INdataframe[,vector %in% colnames(newDataFrame)[column]]}else{
        newDataFrame[,column] <- rowSums(INdataframe[,vector %in% colnames(newDataFrame)[column]])}
  }
  return(newDataFrame)
}

CountTable <- function(in.taxonomy,in.data,output="Count",some.unassigned=T){
  if(length(in.taxonomy)!=length(in.data[,1])){stop("Dataframe and corresponding taxonomy are not the same length")}
  in.taxonomy[is.na(in.taxonomy)] <- ""
  out.dat <- as.data.frame(matrix(ncol=length(in.data[1,]),nrow=length(unique(in.taxonomy))))
  rownames(out.dat) <- sort(unique(in.taxonomy))
  colnames(out.dat) <- colnames(in.data)    
  out.dat.abundance <- out.dat
  for (sample in 1:length(in.data[1,])){
    out.dat[,sample] <- table(in.taxonomy[in.data[,sample]>0])[match(sort(unique(in.taxonomy)),names(table(in.taxonomy[in.data[,sample]>0])))]
    out.dat.abundance[,sample] <- aggregate(in.data[,sample], by=list(Category=in.taxonomy), FUN=sum)[,2]
  }
  out.dat[is.na(out.dat)] <- 0
  if(some.unassigned==T){rownames(out.dat)[1] <- "Unassigned"}
  if(output=="Count"){return(out.dat)}else if(
    output=="Abundance"){return(out.dat.abundance)}
}

minAbundance <- function(inputtable = NA, minAbun = 0.01) {
  others <- rep(0, ncol(inputtable))
  
  for (col in 1:ncol(inputtable)) {
    threshold = sum(inputtable[, col]) * minAbun
    below_threshold_indices = inputtable[, col] < threshold
    others[col] = sum(inputtable[below_threshold_indices, col])
    inputtable[below_threshold_indices, col] = 0
  }
  
  inputtable <- rbind(inputtable, others)
  rownames(inputtable)[nrow(inputtable)] = "Others"
  
  # Remove rows that sum to zero
  inputtable <- inputtable[rowSums(inputtable) != 0, ]
  
  return(inputtable)
}


make_binary <- function (df,threshold){
  df <- sapply(df, function(x) ifelse(is.numeric(x) & x < threshold, 0, 1))
  return(as.data.frame(df))
}

average_by_reps <- function (df,sampleIndex){
  # Check if the dataframe and cols have compatible sizes
  if (ncol(df) != length(sampleIndex)) {
    stop("The number of columns in the dataframe must match the length of the vector")
  }

  # Get unique sample names
  unique_samples <- unique(sampleIndex)
  
  # Initialize the output dataframe
  avg_df <- data.frame(matrix(ncol = length(unique_samples), nrow = nrow(df)))
  colnames(avg_df) <- unique_samples
  
  # Compute averages for each sample
  for (sample in unique_samples) {
    indices <- which(sampleIndex == sample)  # Get indices for columns corresponding to this sample
    if (length(indices) == 1) {
      # If only one column, use it directly
      avg_df[, sample] <- df[, indices, drop = FALSE]
    } else {
      # Compute row means for multiple columns
      avg_df[, sample] <- rowMeans(df[, indices, drop = FALSE], na.rm = TRUE)
    }
  }
  
  return(avg_df)
}

map_values_to_colors <- function(values, color1, color2) {
  # Create a palette function with the two colors
  palette <- colorRampPalette(c(color1, color2))
  
  # Compute the number of unique values
  num_colors <- length(unique(values))
  
  # Generate the colors from the palette
  colors <- palette(num_colors)
  
  # Normalize the values to range from 1 to num_colors
  scaled_values <- cut(values, breaks = pretty(range(values), num_colors), labels = FALSE)
  
  # Assign colors based on scaled values
  assigned_colors <- colors[scaled_values]
  
  return(assigned_colors)
}


plot_with_gradient_legend <- function(values, color1, color2) {
  # Use the mapping function to get assigned colors
  assigned_colors <- map_values_to_colors(values, color1, color2)
  
  
  # Define the gradient legend parameters
  num_colors <- length(unique(values))
  colors <- colorRampPalette(c(color1, color2))(num_colors)  # Recreate palette for the legend
  legend_height <- (par("usr")[4] - par("usr")[3]) / 8  # 1/8 of the plot window height
  legend_width <- 0.1 * (par("usr")[2] - par("usr")[1])  # Relative width of the legend
  
  # Position of the legend at the bottom left
  legend_x <- par("usr")[1] + 0.02 * (par("usr")[2] - par("usr")[1])  # Slightly offset from the left
  legend_y <- par("usr")[3] + 0.02 * (par("usr")[4] - par("usr")[3])  # Slightly offset from the bottom
  
  # Draw the gradient legend
  for (i in 1:num_colors) {
    rect_x <- legend_x
    rect_y <- legend_y + (i-1) * (legend_height / num_colors)
    rect(rect_x, rect_y, rect_x + legend_width, rect_y + (legend_height / num_colors), col = colors[i], border = NA)
  }
  
  # Labeling the min and max values at the bottom and top of the legend
  text(legend_x + legend_width / 2, legend_y - 0.1 * legend_height, labels = min(values), cex = 0.8, adj = 0.5)
  text(legend_x + legend_width / 2, legend_y + legend_height + 0.05 * legend_height, labels = max(values), cex = 0.8, adj = 0.5)
}


add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}



####====0.3 Data Prep====####
ages <- read.csv("metadata/AgeOut.csv")
ages$ID2 <- gsub("-","_",ages$ID)
EUK.tax.PR2 <- read.csv("taxonomy/EUK.cleaned.PR2.csv",row.names = 1)
climate <- read.csv("metadata/climate.csv",row.names = 1)



## EUK datasets
EUK.GC1 <- read.csv("cleaneddata/combinedcoredata/EUK.GC01.csv",row.names = 1)
EUK.GC1.nREPS <- NrepsMaker(EUK.GC1,gsub("(.*)_[0-9]$","\\1",colnames(EUK.GC1)))
EUK.GC1.avr <- average_by_reps(prop.table(as.matrix(EUK.GC1),2),gsub("(.*)_[0-9]$","\\1",colnames(EUK.GC1)))
EUK.GC1.bin <- make_binary(EUK.GC1,1)
EUK.GC1.ages <- ages$median[match(colnames(EUK.GC1.avr),ages$ID2)]

EUK.P19 <- read.csv("cleaneddata/combinedcoredata/EUK.PC019.csv",row.names = 1)
## make a subset of the data gettign rid of high resolution volcano samples 
ages$ID2[184:204]


EUK.P19.nREPS <- NrepsMaker(EUK.P19,gsub("(.*)_[0-9]$","\\1",colnames(EUK.P19)))
EUK.P19.avr <- average_by_reps(prop.table(as.matrix(EUK.P19),2),gsub("(.*)_[0-9]$","\\1",colnames(EUK.P19)))
EUK.P19.bin <- make_binary(EUK.P19,1)
EUK.P19.ages <- ages$median[match(colnames(EUK.P19.avr),ages$ID2)]

EUK.p19.avr.bac <- EUK.P19.avr[match(EUK.tax.PR2$X.1[EUK.tax.PR2$Domain=="Bacteria" & EUK.tax.PR2$Domain.1>80],EUK.tax.PR2$X.1),]
row.names(EUK.p19.avr.bac) <- EUK.tax.PR2$X.1[EUK.tax.PR2$Domain=="Bacteria" & EUK.tax.PR2$Domain.1>80]


EUK.p19.avr.met <- EUK.P19.avr[match(na.omit(EUK.tax.PR2$X.1[EUK.tax.PR2$Subdivision=="Metazoa" & EUK.tax.PR2$Subdivision>80]),EUK.tax.PR2$X.1),]
row.names(EUK.p19.avr.met) <- na.omit(EUK.tax.PR2$X.1[EUK.tax.PR2$Subdivision=="Metazoa" & EUK.tax.PR2$Subdivision>80])


EUK.p19.avr.pro <- EUK.P19.avr[match(na.omit(EUK.tax.PR2$X.1[EUK.tax.PR2$Domain!="Bacteria" & EUK.tax.PR2$Subdivision!="Metazoa" & EUK.tax.PR2$Subdivision!="Fungi" & EUK.tax.PR2$Subdivision>80]),EUK.tax.PR2$X.1),]
row.names(EUK.p19.avr.pro) <- na.omit(EUK.tax.PR2$X.1[EUK.tax.PR2$Domain!="Bacteria" & EUK.tax.PR2$Subdivision!="Metazoa" & EUK.tax.PR2$Subdivision!="Fungi" & EUK.tax.PR2$Subdivision>80])




### MAM datasets
MAM.P19 <- read.csv("cleaneddata/combinedcoredata/MAM.PC019.csv",row.names = 1)
MAM.P19.nREPS <- NrepsMaker(MAM.P19,gsub("(.*)_[0-9]$","\\1",colnames(MAM.P19)))


### RIZ datasets
RIZ.P19 <- read.csv("cleaneddata/combinedcoredata/RIZ.PC019.csv",row.names = 1)
RIZ.P19.nREPS <- NrepsMaker(RIZ.P19,gsub("(.*)_[0-9]$","\\1",colnames(RIZ.P19)))


###Make taxonomic files
MAMtax <- read.csv("taxonomy/MAM.combined.parsed.csv")
MAMasv <- readDNAStringSet("cleaneddata/ASVs/MAM.cleaned.fasta")
MAMtax$ASVseq <- as.character(MAMasv)[match(MAMtax$OTU,names(MAMasv))]

write.csv(MAMtax,file = "taxonomy/byHand/MAMtax.csv")

RIZtax <- read.csv("taxonomy/RIZ.combined.parsed.csv")
RIZasv <- readDNAStringSet("cleaneddata/ASVs/RIZ.cleaned.fasta")
RIZtax$ASVseq <- as.character(RIZasv)[match(RIZtax$OTU,names(RIZasv))]

write.csv(RIZtax,file = "taxonomy/byHand/RIZtax.csv")



####====1.0 Taxonomic overview====####

getPalette = colorRampPalette(brewer.pal(9, "Set3"))


#GC1
EUK.GC1.tax.a <- minAbundance(CountTable(as.character(EUK.tax.PR2$Family),EUK.GC1.avr,output = "Abundance"),minAbun=0.01)
EUK.GC1.tax.c <- minAbundance(CountTable(as.character(EUK.tax.PR2$Family),EUK.GC1.avr,output = "Count"),minAbun=0.01)
row.names(EUK.GC1.tax.a)[1] <- "Unknown"
row.names(EUK.GC1.tax.c)[1] <- "Unknown"
EUK.GC1.tax.a <- as.matrix(prop.table(as.matrix(EUK.GC1.tax.a[,order(ages$median[match(colnames(EUK.GC1.tax.a),ages$ID2)])]),margin = 2))
EUK.GC1.tax.c <- as.matrix(prop.table(as.matrix(EUK.GC1.tax.c[,order(ages$median[match(colnames(EUK.GC1.tax.c),ages$ID2)])]),margin = 2))
barplot(EUK.GC1.tax.a,col=rev(getPalette(dim(EUK.GC1.tax.a)[1])))
barplot(EUK.GC1.tax.c,col=rev(getPalette(dim(EUK.GC1.tax.c)[1])))

pdf("figures/EUK.GC1.tax.Family.pdf",width = 12,height = 9)
par(mfrow=c(2,1),mar=c(5.1, 4.1, 1.1, 6.1),xpd=TRUE)
barplot(EUK.GC1.tax.a,las=2,cex.names=0.6,col=rev(getPalette(dim(EUK.GC1.tax.a)[1])),ylab="Read Abundance",names.arg=ages$median[match(colnames(EUK.GC1.tax.a),ages$ID2)])
legend(56,1,rev(rownames(EUK.GC1.tax.a)),fill=getPalette(dim(EUK.GC1.tax.a)[1]),cex=0.5,bty = "n",y.intersp=0.75)
barplot(EUK.GC1.tax.c,las=2,cex.names=0.6,col=rev(getPalette(dim(EUK.GC1.tax.c)[1])),ylab="ASV Counts",names.arg=ages$median[match(colnames(EUK.GC1.tax.c),ages$ID2)])
legend(56,1,rev(rownames(EUK.GC1.tax.c)),fill=getPalette(dim(EUK.GC1.tax.c)[1]),cex=0.5,bty = "n",y.intersp=0.75)
dev.off()



##PC19
EUK.P19.tax.a <- minAbundance(CountTable(as.character(EUK.tax.PR2$Family),EUK.P19.avr,output = "Abundance"),minAbun=0.01)
EUK.P19.tax.c <- minAbundance(CountTable(as.character(EUK.tax.PR2$Family),EUK.P19.avr,output = "Count"),minAbun=0.01)
row.names(EUK.P19.tax.a)[1] <- "Unknown"
row.names(EUK.P19.tax.c)[1] <- "Unknown"
EUK.P19.tax.a <- as.matrix(prop.table(as.matrix(EUK.P19.tax.a[,order(ages$median[match(colnames(EUK.P19.tax.a),ages$ID2)])]),margin = 2))
EUK.P19.tax.c <- as.matrix(prop.table(as.matrix(EUK.P19.tax.c[,order(ages$median[match(colnames(EUK.P19.tax.c),ages$ID2)])]),margin = 2))

pdf("figures/EUK.P19.tax.Family.pdf",width = 12,height = 9)
par(mfrow=c(2,1),mar=c(5.1, 4.1, 1.1, 6.1),xpd=TRUE)
barplot(EUK.P19.tax.a,las=2,cex.names=0.6,col=rev(getPalette(dim(EUK.P19.tax.a)[1])),ylab="Read Abundance",names.arg=ages$median[match(colnames(EUK.P19.tax.a),ages$ID2)])
legend(190,1,rev(rownames(EUK.P19.tax.a)),fill=getPalette(dim(EUK.P19.tax.a)[1]),cex=0.5,bty = "n",y.intersp=0.75)
barplot(EUK.P19.tax.c,las=2,cex.names=0.6,col=rev(getPalette(dim(EUK.P19.tax.c)[1])),ylab="ASV Counts",names.arg=ages$median[match(colnames(EUK.P19.tax.c),ages$ID2)])
legend(190,1,rev(rownames(EUK.P19.tax.c)),fill=getPalette(dim(EUK.P19.tax.c)[1]),cex=0.5,bty = "n",y.intersp=0.75)
dev.off()


####====2.0 Alpha Diversity====####
## GC1
## calc raw richness 
EUK.GC1.raw.richness <- colSums(make_binary(EUK.GC1,2))

## calc breakaway richness
EUK.richEst <- breakaway(EUK.GC1)
EUK.GC1.breakaway.richEstimate <- unlist(lapply(EUK.richEst,FUN = function(x){x[["estimate"]]}))

##cacl rarefy rich
EUK.GC1.rare.richness.1 <- colSums(make_binary(as.data.frame(t(rrarefy(t(EUK.GC1),7000))),2))
EUK.GC1.rare.richness.2 <- colSums(make_binary(as.data.frame(t(rrarefy(t(EUK.GC1),14000))),2))


## compare values
plot(EUK.GC1.raw.richness,EUK.GC1.breakaway.richEstimate)
plot(EUK.GC1.raw.richness,EUK.GC1.rare.richness.1)
plot(EUK.GC1.raw.richness,EUK.GC1.rare.richness.2)

EUK.GC1.mean.rich <- tapply(EUK.GC1.raw.richness,FUN=mean,INDEX = gsub("(.*)_[0-9]$","\\1",names(EUK.GC1.raw.richness)))

## plot vs time
pdf("figures/EUK.GC1.richness.pdf",width = 9,height = 6.5)
plot(ages$median[match(gsub("(.*)_[0-9]$","\\1",colnames(EUK.GC1)),ages$ID2)],
     colSums(make_binary(EUK.GC1,2)),
     pch=16,
     xlab="Cal yr BP",
     ylab="ASV Richness",
     col="grey")
points(ages$median[match(names(EUK.GC1.mean.rich),ages$ID2)],
       EUK.GC1.mean.rich,
       pch=16,
       col="black",
       cex=2)
dev.off()

plot(ages$median[match(gsub("(.*)_[0-9]$","\\1",names(EUK.GC1)),ages$ID2)],
     colSums(make_binary(as.data.frame(t(rrarefy(t(EUK.GC1),14000))),2)),
     pch=16,
     xlab="Cal yr BP",
     ylab="ASV Richness")

## with GAM

GC1.rich.gam.dat  <- data.frame("year"=ages$median[match(names(EUK.GC1.mean.rich),ages$ID2)],
                           "richness"=EUK.GC1.mean.rich)

m <- gam(richness ~ s(year, k = 20), data = GC1.rich.gam.dat, method = "REML")
prediction <- data.frame("year"=280:3600)
prediction <- cbind(prediction,predict(m,prediction,se.fit=TRUE))
deriv <- gratia::derivatives(m)
prediction$uppCI <- prediction$fit+prediction$se.fit*1.96
prediction$lwrCI <- prediction$fit-prediction$se.fit*1.96


pdf("figures/EUK.GC1.richness.gam.pdf",width = 9,height = 4.5)
par(mar=c(5.1, 4.1, 1.1, 1.1))
plot(ages$median[match(gsub("(.*)_[0-9]$","\\1",colnames(EUK.GC1)),ages$ID2)],
     colSums(make_binary(EUK.GC1,2)),
     pch=16,
     xlab="Cal yr BP",
     ylab="ASV Richness",
      col="grey")
polygon(c(prediction$year, rev(prediction$year)), c(prediction$uppCI, rev(prediction$lwrCI)), col=add.alpha('darkgrey',0.3), border=NA)
points(prediction$year,prediction$fit,lwd=5,col="black",type='l')
dev.off()


## P19
## calc raw richness 
EUK.P19.raw.richness <- colSums(make_binary(EUK.P19,2))

## calc breakaway richness
EUK.richEst <- breakaway(EUK.P19)
EUK.P19.breakaway.richEstimate <- unlist(lapply(EUK.richEst,FUN = function(x){x[["estimate"]]}))

##cacl rarefy rich
EUK.P19.rare.richness.1 <- colSums(make_binary(as.data.frame(t(rrarefy(t(EUK.P19),7000))),2))
EUK.P19.rare.richness.2 <- colSums(make_binary(as.data.frame(t(rrarefy(t(EUK.P19),14000))),2))


## compare values
plot(EUK.P19.raw.richness,EUK.P19.breakaway.richEstimate)
plot(EUK.P19.raw.richness,EUK.P19.rare.richness.1)
plot(EUK.P19.raw.richness,EUK.P19.rare.richness.2)

EUK.P19.mean.rich <- tapply(EUK.P19.raw.richness,FUN=mean,INDEX = gsub("(.*)_[0-9]$","\\1",names(EUK.P19.raw.richness)))

## plot vs time
pdf("figures/EUK.P19.richness.pdf",width = 9,height = 6.5)
plot(ages$median[match(gsub("(.*)_[0-9]$","\\1",colnames(EUK.P19)),ages$ID2)],
     colSums(make_binary(EUK.P19,2)),
     pch=16,
     xlab="Cal yr BP",
     ylab="ASV Richness",
     col="grey")
points(ages$median[match(names(EUK.P19.mean.rich),ages$ID2)],
       EUK.P19.mean.rich,
       pch=16,
       col="black",
       cex=2)
dev.off()

plot(ages$median[match(gsub("(.*)_[0-9]$","\\1",names(EUK.P19)),ages$ID2)],
     colSums(make_binary(as.data.frame(t(rrarefy(t(EUK.P19),14000))),2)),
     pch=16,
     xlab="Cal yr BP",
     ylab="ASV Richness")



## with GAM

P19.rich.gam.dat  <- data.frame("year"=ages$median[match(names(EUK.P19.mean.rich),ages$ID2)],
                                "richness"=EUK.P19.mean.rich)

m <- gam(richness ~ s(year, k = 20), data = P19.rich.gam.dat, method = "REML")
prediction <- data.frame("year"=200:3100)
prediction <- cbind(prediction,predict(m,prediction,se.fit=TRUE))
deriv <- gratia::derivatives(m)
prediction$uppCI <- prediction$fit+prediction$se.fit*1.96
prediction$lwrCI <- prediction$fit-prediction$se.fit*1.96


pdf("figures/EUK.P19.richness.gam.pdf",width = 9,height = 4.5)
par(mar=c(5.1, 4.1, 1.1, 1.1))
plot(ages$median[match(gsub("(.*)_[0-9]$","\\1",colnames(EUK.P19)),ages$ID2)],
     colSums(make_binary(EUK.P19,2)),
     pch=16,
     xlab="Cal yr BP",
     ylab="ASV Richness",
     col="grey")
     #ylim=c(500,1500))
polygon(c(prediction$year, rev(prediction$year)), c(prediction$uppCI, rev(prediction$lwrCI)), col=add.alpha('darkgrey',0.3), border=NA)
points(prediction$year,prediction$fit,lwd=5,col="black",type='l')
dev.off()



### Just metazoa

EUK.p19.avr.met.b <- make_binary((EUK.p19.avr.met)*1000000000000000,1)
EUK.p19.avr.met.tax <- EUK.tax.PR2[na.omit(match(rownames(EUK.p19.avr.met),EUK.tax.PR2$X.1)),]
EUK.p19.avr.met.tax$Genus[EUK.p19.avr.met.tax$Genus.1<80] <- NA
genus_richness_per_site <- apply(EUK.p19.avr.met.b, 2, function (x) length(unique(EUK.p19.avr.met.tax$Genus[which(x == 1)])))


P19.rich.gam.met  <- data.frame("year"=ages$median[match(colnames(EUK.p19.avr.met),ages$ID2)],
                                "richness"=colSums(make_binary((EUK.p19.avr.met)*1000000000000000,1)),
                                "genera.richness"=genus_richness_per_site)

m <- gam(richness ~ s(year, k = 20), data = P19.rich.gam.met, method = "REML")
prediction <- data.frame("year"=200:3100)
prediction <- cbind(prediction,predict(m,prediction,se.fit=TRUE))
prediction$uppCI <- prediction$fit+prediction$se.fit*1.96
prediction$lwrCI <- prediction$fit-prediction$se.fit*1.96

m <- gam(genera.richness ~ s(year, k = 20), data = P19.rich.gam.met, method = "REML")
prediction2 <- data.frame("year"=200:3100)
prediction2 <- cbind(prediction2,predict(m,prediction2,se.fit=TRUE))
prediction2$uppCI <- prediction2$fit+prediction2$se.fit*1.96
prediction2$lwrCI <- prediction2$fit-prediction2$se.fit*1.96

  ## all
pdf("figures/EUK.P19.richness.gam.met.pdf",width = 9,height = 4.5)
plot(ages$median[match(colnames(EUK.p19.avr.met),ages$ID2)],
     colSums(make_binary((EUK.p19.avr.met)*1000000000000000,1)),
     pch=16,
     xlab="Cal yr BP",
     ylab="ASV Richness",
     col="grey")
polygon(c(prediction$year, rev(prediction$year)), c(prediction$uppCI, rev(prediction$lwrCI)), col=add.alpha('darkgrey',0.3), border=NA)
points(prediction$year,prediction$fit,lwd=5,col="black",type='l')
dev.off()

  ## genus
pdf("figures/EUK.P19.richness.gam.met.genus.pdf",width = 9,height = 4.5)
plot(ages$median[match(names(genus_richness_per_site),ages$ID2)],
     genus_richness_per_site,
     pch=16,
     xlab="Cal yr BP",
     ylab="ASV Richness",
     col="grey")
polygon(c(prediction2$year, rev(prediction2$year)), c(prediction2$uppCI, rev(prediction2$lwrCI)), col=add.alpha('darkgrey',0.3), border=NA)
points(prediction2$year,prediction2$fit,lwd=5,col="black",type='l')
dev.off()

### protists

EUK.p19.avr.pro.b <- make_binary((EUK.p19.avr.pro)*1000000000000000,1)
EUK.p19.avr.pro.tax <- EUK.tax.PR2[na.omit(match(rownames(EUK.p19.avr.pro),EUK.tax.PR2$X.1)),]
EUK.p19.avr.pro.tax$Genus[EUK.p19.avr.pro.tax$Genus.1<80] <- NA
genus_richness_per_site <- apply(EUK.p19.avr.pro.b, 2, function (x) length(unique(EUK.p19.avr.pro.tax$Genus[which(x == 1)])))


P19.rich.gam.pro  <- data.frame("year"=ages$median[match(colnames(EUK.p19.avr.pro),ages$ID2)],
                                "richness"=colSums(make_binary((EUK.p19.avr.pro)*1000000000000000,1)),
                                "genera.richness"=genus_richness_per_site)

m <- gam(richness ~ s(year, k = 20), data = P19.rich.gam.pro, method = "REML")
prediction <- data.frame("year"=200:3100)
prediction <- cbind(prediction,predict(m,prediction,se.fit=TRUE))
prediction$uppCI <- prediction$fit+prediction$se.fit*1.96
prediction$lwrCI <- prediction$fit-prediction$se.fit*1.96

m <- gam(genera.richness ~ s(year, k = 20), data = P19.rich.gam.pro, method = "REML")
prediction2 <- data.frame("year"=200:3100)
prediction2 <- cbind(prediction2,predict(m,prediction2,se.fit=TRUE))
prediction2$uppCI <- prediction2$fit+prediction2$se.fit*1.96
prediction2$lwrCI <- prediction2$fit-prediction2$se.fit*1.96

## all
pdf("figures/EUK.P19.richness.gam.pro.pdf",width = 9,height = 4.5)
plot(ages$median[match(colnames(EUK.p19.avr.pro),ages$ID2)],
     colSums(make_binary((EUK.p19.avr.pro)*1000000000000000,1)),
     pch=16,
     xlab="Cal yr BP",
     ylab="ASV Richness",
     col="grey")
polygon(c(prediction$year, rev(prediction$year)), c(prediction$uppCI, rev(prediction$lwrCI)), col=add.alpha('darkgrey',0.3), border=NA)
points(prediction$year,prediction$fit,lwd=5,col="black",type='l')
dev.off()

## genus
pdf("figures/EUK.P19.richness.gam.pro.genus.pdf",width = 9,height = 4.5)
plot(ages$median[match(names(genus_richness_per_site),ages$ID2)],
     genus_richness_per_site,
     pch=16,
     xlab="Cal yr BP",
     ylab="ASV Richness",
     col="grey")
polygon(c(prediction2$year, rev(prediction2$year)), c(prediction2$uppCI, rev(prediction2$lwrCI)), col=add.alpha('darkgrey',0.3), border=NA)
points(prediction2$year,prediction2$fit,lwd=5,col="black",type='l')
dev.off()

### bacteria

EUK.p19.avr.bac.b <- make_binary((EUK.p19.avr.bac)*1000000000000000,1)
EUK.p19.avr.bac.tax <- EUK.tax.PR2[na.omit(match(rownames(EUK.p19.avr.bac),EUK.tax.PR2$X.1)),]
EUK.p19.avr.bac.tax$Genus[EUK.p19.avr.bac.tax$Genus.1<80] <- NA
genus_richness_per_site <- apply(EUK.p19.avr.bac.b, 2, function (x) length(unique(EUK.p19.avr.bac.tax$Genus[which(x == 1)])))


P19.rich.gam.bac  <- data.frame("year"=ages$median[match(colnames(EUK.p19.avr.bac),ages$ID2)],
                                "richness"=colSums(make_binary((EUK.p19.avr.bac)*1000000000000000,1)),
                                "genera.richness"=genus_richness_per_site)

m <- gam(richness ~ s(year, k = 20), data = P19.rich.gam.bac, method = "REML")
prediction <- data.frame("year"=200:3100)
prediction <- cbind(prediction,predict(m,prediction,se.fit=TRUE))
prediction$uppCI <- prediction$fit+prediction$se.fit*1.96
prediction$lwrCI <- prediction$fit-prediction$se.fit*1.96

m <- gam(genera.richness ~ s(year, k = 20), data = P19.rich.gam.bac, method = "REML")
prediction2 <- data.frame("year"=200:3100)
prediction2 <- cbind(prediction2,predict(m,prediction2,se.fit=TRUE))
prediction2$uppCI <- prediction2$fit+prediction2$se.fit*1.96
prediction2$lwrCI <- prediction2$fit-prediction2$se.fit*1.96

## all
pdf("figures/EUK.P19.richness.gam.bac.pdf",width = 9,height = 4.5)
plot(ages$mean[match(colnames(EUK.p19.avr.bac),ages$ID2)],
     colSums(make_binary((EUK.p19.avr.bac)*1000000000000000,1)),
     pch=16,
     xlab="Cal yr BP",
     ylab="ASV Richness",
     col="grey")
polygon(c(prediction$year, rev(prediction$year)), c(prediction$uppCI, rev(prediction$lwrCI)), col=add.alpha('darkgrey',0.3), border=NA)
points(prediction$year,prediction$fit,lwd=5,col="black",type='l')
dev.off()

## genus
pdf("figures/EUK.P19.richness.gam.bac.genus.pdf",width = 9,height = 4.5)
plot(ages$mean[match(names(genus_richness_per_site),ages$ID2)],
     genus_richness_per_site,
     pch=16,
     xlab="Cal yr BP",
     ylab="ASV Richness",
     col="grey")
polygon(c(prediction2$year, rev(prediction2$year)), c(prediction2$uppCI, rev(prediction2$lwrCI)), col=add.alpha('darkgrey',0.3), border=NA)
points(prediction2$year,prediction2$fit,lwd=5,col="black",type='l')
dev.off()



#### now lets try and regress the data!

#subset all observations between 225bp and 987bp
datain <- colSums(make_binary((EUK.p19.avr.met)*1000000000000000,1))
datain <- datain[ages$mean[match(names(datain),ages$ID2)] < 987 & ages$mean[match(names(datain),ages$ID2)] > 225]
climate_subset <- climate[match(names(datain),gsub("-","_",climate$Sample)),]

### make a massive GLM 

datasetin <- data.frame(
  richness = datain,
  dateBP = climate_subset$Date_bp_mean,
  d18O_GAM_fit = climate_subset$d18O.GAM.fit,
  d18O_100yrspline = climate_subset$d18O.100yrspline,
  d13C_GAM_fit = climate_subset$d13C.GAM.fit,
  d13C_100yrspline = climate_subset$d13C.100yrspline,
  MD99_2275sstAlkenone_GAM_fit = climate_subset$MD99.2275sstAlkenone.GAM.fit,
  MD99_2275sstAlkenone_100yrspline = climate_subset$MD99.2275sstAlkenone.100yrspline,
  MD99_2275sstDiatom_GAM_fit = climate_subset$MD99.2275sstDiatom.GAM.fit,
  MD99_2275sstDiatom_100yrspline = climate_subset$MD99.2275sstDiatom.100yrspline,
  MD99_2275IP25_GAM_fit = climate_subset$MD99.2275IP25.GAM.fit,
  MD99_2275IP25_100yrspline = climate_subset$MD99.2275IP25.100yrspline
)

datasetin.mini <- datasetin[,c(1,2,4,6,7,9,11)]

pdf("figures/climate/comparison.met.pdf",width = 9.5,height = 7)
plot(datasetin.mini,pch=16)
dev.off()

m1 <- glm(richness~dateBP+
            d18O_GAM_fit+
            d18O_100yrspline+
            d13C_GAM_fit+
            d13C_100yrspline+
            MD99_2275sstAlkenone_GAM_fit+
            MD99_2275sstAlkenone_100yrspline+
            MD99_2275sstDiatom_GAM_fit+
            MD99_2275sstDiatom_100yrspline+
            MD99_2275IP25_GAM_fit+
            MD99_2275IP25_100yrspline,
          data=datasetin)

plot(m1)
summary(m1)

summary(lm(richness~dateBP,data=datasetin))
summary(lm(richness~d18O_100yrspline,data=datasetin))
summary(lm(richness~d13C_100yrspline,data=datasetin))
summary(lm(richness~MD99_2275sstAlkenone_GAM_fit,data=datasetin))
summary(lm(richness~MD99_2275sstDiatom_GAM_fit,data=datasetin))
summary(lm(richness~MD99_2275IP25_GAM_fit,data=datasetin))


#Make in intelligent/interesting GLM


m2 <- glm(richness ~ dateBP + 
            d18O_100yrspline + 
            d13C_100yrspline+
            MD99_2275sstAlkenone_GAM_fit + 
            MD99_2275sstDiatom_GAM_fit+
            MD99_2275IP25_GAM_fit, 
          data = datasetin)


m3 <- glm(richness ~ dateBP + 
            d18O_100yrspline + 
            MD99_2275sstAlkenone_GAM_fit,
          data = datasetin)


summary(m2)
summary(m3)

pdf("figures/climate/GLM.met.pdf",width = 9.5,height = 6.5)
par(mfrow=c(3,1))
plot(allEffects(m3))
dev.off()

 #climate_subset <- climate[climate$Core=="PC019" & climate$Date_bp_mean > 225 & climate$Date_bp_mean < 987,]

### let's plot some climate lines
pdf("figures/climate/Pc019overview.pdf",width = 5,height = 9)
par(mfrow=c(5,1),mar=c(0.5, 4.1, 0.5, 0.5))
plot(climate_subset$Date_CE_mean[order(climate_subset$Date_CE_mean)],climate_subset$d13C.100yrspline[order(climate_subset$Date_CE_mean)],lwd=3,type="l",bty="n",xlab="Year CE",ylab="d13C",xaxt='n')
plot(climate_subset$Date_CE_mean[order(climate_subset$Date_CE_mean)],climate_subset$d18O.100yrspline[order(climate_subset$Date_CE_mean)],lwd=3,type="l",bty="n",xlab="Year CE",ylab="d18O",xaxt='n')
plot(climate_subset$Date_CE_mean[order(climate_subset$Date_CE_mean)],climate_subset$MD99.2275sstAlkenone.GAM.fit[order(climate_subset$Date_CE_mean)],lwd=3,type="l",bty="n",xlab="Year CE",ylab="Alkenone",xaxt='n')
plot(climate_subset$Date_CE_mean[order(climate_subset$Date_CE_mean)],climate_subset$MD99.2275sstDiatom.GAM.fit[order(climate_subset$Date_CE_mean)],lwd=3,type="l",bty="n",xlab="Year CE",ylab="DiatomTF",xaxt='n')
par(mar=c(3, 4.1, 0.5, 0.5))
plot(climate_subset$Date_CE_mean[order(climate_subset$Date_CE_mean)],climate_subset$MD99.2275IP25.GAM.fit[order(climate_subset$Date_CE_mean)],lwd=3,type="l",bty="n",xlab="Year CE",ylab="IP25",xaxt='n')
axis(1,seq(1000,1700,100),cex=3)
dev.off()

climate <- read.csv("rawdata/climate/combined.csv")

pdf("figures/climate/Pc019overlap.pdf",width = 8,height = 4.5)
par(mar=c(5.1, 10.1, 0.5, 0.5))
plot(1950-ages$mean[match(colnames(EUK.p19.avr.bac),ages$ID2)],
     rep(1,length(1950-ages$mean[match(colnames(EUK.p19.avr.bac),ages$ID2)])),pch=16,
     ylim=c(0,7),xlim=c(-8000,1900),yaxt='n',ylab="",xlab="Year CE")
polygon(c(1725,1725,963,963),c(-1,8,8,-1),col="grey",border = NA)
points(1950-ages$mean[match(colnames(EUK.p19.avr.bac),ages$ID2)],rep(1,length(1950-ages$mean[match(colnames(EUK.p19.avr.bac),ages$ID2)])),pch=16)
points(climate$year[climate$record=="Arcticad13C" ],rep(2,length(climate$year[climate$record=="Arcticad13C" ])),pch=16)
points(climate$year[climate$record=="Arcticad18O"],rep(3,length(climate$year[climate$record=="Arcticad18O"])),pch=16)
points(climate$year[climate$record=="SeaIceIndex"],rep(4,length(climate$year[climate$record=="SeaIceIndex"])),pch=16)
points(1950-climate$year[climate$record=="MD99-2275sstAlkenone"],rep(5,length(climate$year[climate$record=="MD99-2275sstAlkenone"])),pch=16)
points(1950-climate$year[climate$record=="MD99-2275sstDiatom"],rep(6,length(climate$year[climate$record=="MD99-2275sstDiatom"])),pch=16)
points(climate$year[climate$record=="MD99-2275sstIP25" ],rep(7,length(climate$year[climate$record=="MD99-2275sstIP25" ])),pch=16)
axis(2,1:7,labels=c("eDNA samples","Arcticad13C","Arcticad18O","SeaIceIndex","MD99-2275sstAlkenone","MD99-2275sstDiatom","MD99-2275sstIP25"),las=1)
dev.off()

####====3.0 Beta diversity ====####

## first lets identify outliers


## GC1
blacklist <- c("GC1_000_4","GC1_044_1","GC1_064_5","GC1_108_2","GC1_100_1","GC1_100_6","GC1_104_1","GC1_164_3","GC1_168_2","GC1_160_6","GC1_060_3","GC1_036_2","GC1_092_2","GC1_120_6")
blacklist <- c("GC1_148_1","GC1_180_7","GC1_168_2","GC1_104_1","GC1_100_6","GC1_160_6","GC1_164_3","GC1_108_2","GC1_044_1","GC1_116_1","GC1_116_2","GC1_116_3","GC1_116_4","GC1_116_5","GC1_116_6","GC1_116_7","GC1_116_8","GC1_064_5","GC1_000_4","GC1_140_1","GC1_140_2","GC1_140_3","GC1_140_4","GC1_140_5","GC1_140_6","GC1_140_7","GC1_140_8")

test <- EUK.GC1[,!colnames(EUK.GC1) %in% blacklist]
test <- EUK.GC1
EUK.GC1.nMDS.b <- metaMDS(t(prop.table(as.matrix(test),2)),k=3,trymax = 200,parallel=8)
EUK.GC1.nMDS.b <- metaMDS(t(prop.table(as.matrix(EUK.GC1[,!colnames(EUK.GC1) %in% blacklist]),2)),k=3,trymax = 200,parallel=8)

## here we use the average dataset 
EUK.GC1.MDS.b <- wcmdscale(vegdist(t(EUK.GC1.avr[,!colnames(EUK.GC1.avr) %in% c("GC1_116","GC1_140")])),eig = TRUE)
EUK.GC1.MDS.j <- wcmdscale(vegdist(t(EUK.GC1.avr[,!colnames(EUK.GC1.avr) %in% c("GC1_116","GC1_140")]),method = "jac",binary=TRUE),eig = TRUE)




pdf("figures/EUK.GC1.bray.nMDS.pdf",width = 9,height = 3.5)
par(mfrow=c(1,3),mar=c(5.1, 4.1, 1.1, 1.1),xpd=TRUE)
plot(EUK.GC1.nMDS.b$points[,1],EUK.GC1.nMDS.b$points[,2],xlab="nMDS1",ylab="nMDS2",pch=16,col=map_values_to_colors(ages$median[match(gsub("(.*)_[0-9]$","\\1",rownames(EUK.GC1.nMDS.b$points)),ages$ID2)], "dodgerblue", "darkred"))
plot_with_gradient_legend(ages$median[match(gsub("(.*)_[0-9]$","\\1",rownames(EUK.GC1.nMDS.b$points)),ages$ID2)], "dodgerblue", "darkred")
plot(EUK.GC1.nMDS.b$points[,1],EUK.GC1.nMDS.b$points[,3],xlab="nMDS1",ylab="nMDS3",pch=16,col=map_values_to_colors(ages$median[match(gsub("(.*)_[0-9]$","\\1",rownames(EUK.GC1.nMDS.b$points)),ages$ID2)], "dodgerblue", "darkred"))
plot(EUK.GC1.nMDS.b$points[,2],EUK.GC1.nMDS.b$points[,3],xlab="nMDS2",ylab="nMDS3",pch=16,col=map_values_to_colors(ages$median[match(gsub("(.*)_[0-9]$","\\1",rownames(EUK.GC1.nMDS.b$points)),ages$ID2)], "dodgerblue", "darkred"))
dev.off()

pdf("figures/EUK.GC1.bray.nMDS.linear.pdf",width = 8,height = 8)
par(mfrow=c(3,1),mar=c(5.1, 4.1, 1.1, 1.1),xpd=TRUE)
plot(ages$median[match(gsub("(.*)_[0-9]$","\\1",rownames(EUK.GC1.nMDS.b$points)),ages$ID2)],EUK.GC1.nMDS.b$points[,1],pch=16,ylab="nMDS 1",xlab="Yr Cal BP")
plot(ages$median[match(gsub("(.*)_[0-9]$","\\1",rownames(EUK.GC1.nMDS.b$points)),ages$ID2)],EUK.GC1.nMDS.b$points[,2],pch=16,ylab="nMDS 2",xlab="Yr Cal BP")
plot(ages$median[match(gsub("(.*)_[0-9]$","\\1",rownames(EUK.GC1.nMDS.b$points)),ages$ID2)],EUK.GC1.nMDS.b$points[,3],pch=16,ylab="nMDS 3",xlab="Yr Cal BP")
dev.off()



## Lets make some GAMS

GC1.MDS.b.gam.dat  <- data.frame("year"=ages$median[match(rownames(EUK.GC1.MDS.b$points),ages$ID2)],
                                "MDS1"=EUK.GC1.MDS.b$points[,1],
                                "MDS2"=EUK.GC1.MDS.b$points[,2],
                                "MDS3"=EUK.GC1.MDS.b$points[,3],
                                "MDS4"=EUK.GC1.MDS.b$points[,4])

m1 <- gam(MDS1 ~ s(year), data = GC1.MDS.b.gam.dat, method = "REML")
m2 <- gam(MDS2 ~ s(year, k = 20), data = GC1.MDS.b.gam.dat, method = "REML")
m3 <- gam(MDS3 ~ s(year, k = 20), data = GC1.MDS.b.gam.dat, method = "REML")
m4 <- gam(MDS4 ~ s(year, k = 20), data = GC1.MDS.b.gam.dat, method = "REML")

prediction <- data.frame("year"=280:3600)
prediction1 <- cbind(prediction,predict(m1,prediction,se.fit=TRUE))
prediction2 <- cbind(prediction,predict(m2,prediction,se.fit=TRUE))
prediction3 <- cbind(prediction,predict(m3,prediction,se.fit=TRUE))
prediction4 <- cbind(prediction,predict(m4,prediction,se.fit=TRUE))
prediction1$uppCI <- prediction1$fit+prediction1$se.fit*1.96
prediction1$lwrCI <- prediction1$fit-prediction1$se.fit*1.96
prediction2$uppCI <- prediction2$fit+prediction2$se.fit*1.96
prediction2$lwrCI <- prediction2$fit-prediction2$se.fit*1.96
prediction3$uppCI <- prediction3$fit+prediction3$se.fit*1.96
prediction3$lwrCI <- prediction3$fit-prediction3$se.fit*1.96
prediction4$uppCI <- prediction4$fit+prediction4$se.fit*1.96
prediction4$lwrCI <- prediction4$fit-prediction4$se.fit*1.96

deriv <- gratia::derivatives(m)


## plot



pdf("figures/EUK.GC1.bray.MDS.linear.pdf",width = 8,height = 8)
par(mfrow=c(4,1),mar=c(2, 4.1, 1.1, 1.1))
plot(ages$median[match(rownames(EUK.GC1.MDS.b$points),ages$ID2)],EUK.GC1.MDS.b$points[,1],pch=16,ylab="MDS 1",xlab="Yr Cal BP")
polygon(c(prediction1$year, rev(prediction1$year)), c(prediction1$uppCI, rev(prediction1$lwrCI)), col=add.alpha('dodgerblue',0.3), border=NA)
points(prediction1$year,prediction1$fit,lwd=2,col="dodgerblue",type='l')
plot(ages$median[match(rownames(EUK.GC1.MDS.b$points),ages$ID2)],EUK.GC1.MDS.b$points[,2],pch=16,ylab="MDS 2",xlab="Yr Cal BP")
polygon(c(prediction2$year, rev(prediction2$year)), c(prediction2$uppCI, rev(prediction2$lwrCI)), col=add.alpha('pink3',0.3), border=NA)
points(prediction2$year,prediction2$fit,lwd=2,col="pink3",type='l')
plot(ages$median[match(rownames(EUK.GC1.MDS.b$points),ages$ID2)],EUK.GC1.MDS.b$points[,3],pch=16,ylab="MDS 3",xlab="Yr Cal BP")
polygon(c(prediction3$year, rev(prediction3$year)), c(prediction3$uppCI, rev(prediction3$lwrCI)), col=add.alpha('darkgreen',0.3), border=NA)
points(prediction3$year,prediction3$fit,lwd=2,col="darkgreen",type='l')
plot(ages$median[match(rownames(EUK.GC1.MDS.b$points),ages$ID2)],EUK.GC1.MDS.b$points[,4],pch=16,ylab="MDS 4",xlab="Yr Cal BP")
polygon(c(prediction4$year, rev(prediction4$year)), c(prediction4$uppCI, rev(prediction4$lwrCI)), col=add.alpha('gold',0.3), border=NA)
points(prediction4$year,prediction4$fit,lwd=2,col="gold",type='l')
dev.off()



#GAMS
GC1.MDS.j.gam.dat  <- data.frame("year"=ages$median[match(rownames(EUK.GC1.MDS.j$points),ages$ID2)],
                                 "MDS1"=EUK.GC1.MDS.j$points[,1],
                                 "MDS2"=EUK.GC1.MDS.j$points[,2],
                                 "MDS3"=EUK.GC1.MDS.j$points[,3],
                                 "MDS4"=EUK.GC1.MDS.j$points[,4])

m1 <- gam(MDS1 ~ s(year, k = 20), data = GC1.MDS.j.gam.dat, method = "REML")
m2 <- gam(MDS2 ~ s(year, k = 20), data = GC1.MDS.j.gam.dat, method = "REML")
m3 <- gam(MDS3 ~ s(year, k = 20), data = GC1.MDS.j.gam.dat, method = "REML")
m4 <- gam(MDS4 ~ s(year, k = 20), data = GC1.MDS.j.gam.dat, method = "REML")

prediction <- data.frame("year"=280:3600)
prediction1 <- cbind(prediction,predict(m1,prediction,se.fit=TRUE))
prediction2 <- cbind(prediction,predict(m2,prediction,se.fit=TRUE))
prediction3 <- cbind(prediction,predict(m3,prediction,se.fit=TRUE))
prediction4 <- cbind(prediction,predict(m4,prediction,se.fit=TRUE))
prediction1$uppCI <- prediction1$fit+prediction1$se.fit*1.96
prediction1$lwrCI <- prediction1$fit-prediction1$se.fit*1.96
prediction2$uppCI <- prediction2$fit+prediction2$se.fit*1.96
prediction2$lwrCI <- prediction2$fit-prediction2$se.fit*1.96
prediction3$uppCI <- prediction3$fit+prediction3$se.fit*1.96
prediction3$lwrCI <- prediction3$fit-prediction3$se.fit*1.96
prediction4$uppCI <- prediction4$fit+prediction4$se.fit*1.96
prediction4$lwrCI <- prediction4$fit-prediction4$se.fit*1.96

pdf("figures/EUK.GC1.jacc.MDS.linear.pdf",width = 8,height = 8)
par(mfrow=c(4,1),mar=c(2, 4.1, 1.1, 1.1))
plot(ages$median[match(rownames(EUK.GC1.MDS.j$points),ages$ID2)],EUK.GC1.MDS.j$points[,1],pch=16,ylab="MDS 1",xlab="Yr Cal BP")
polygon(c(prediction1$year, rev(prediction1$year)), c(prediction1$uppCI, rev(prediction1$lwrCI)), col=add.alpha('dodgerblue',0.3), border=NA)
points(prediction1$year,prediction1$fit,lwd=2,col="dodgerblue",type='l')
plot(ages$median[match(rownames(EUK.GC1.MDS.j$points),ages$ID2)],EUK.GC1.MDS.j$points[,2],pch=16,ylab="MDS 2",xlab="Yr Cal BP")
polygon(c(prediction2$year, rev(prediction2$year)), c(prediction2$uppCI, rev(prediction2$lwrCI)), col=add.alpha('pink3',0.3), border=NA)
points(prediction2$year,prediction2$fit,lwd=2,col="pink3",type='l')
plot(ages$median[match(rownames(EUK.GC1.MDS.j$points),ages$ID2)],EUK.GC1.MDS.j$points[,3],pch=16,ylab="MDS 3",xlab="Yr Cal BP")
polygon(c(prediction3$year, rev(prediction3$year)), c(prediction3$uppCI, rev(prediction3$lwrCI)), col=add.alpha('darkgreen',0.3), border=NA)
points(prediction3$year,prediction3$fit,lwd=2,col="darkgreen",type='l')
plot(ages$median[match(rownames(EUK.GC1.MDS.j$points),ages$ID2)],EUK.GC1.MDS.j$points[,4],pch=16,ylab="MDS 4",xlab="Yr Cal BP")
polygon(c(prediction4$year, rev(prediction4$year)), c(prediction4$uppCI, rev(prediction4$lwrCI)), col=add.alpha('gold',0.3), border=NA)
points(prediction4$year,prediction4$fit,lwd=2,col="gold",type='l')
dev.off()



pdf("figures/EUK.GC1.bray.MDS.linear.prop.pdf",width = 5,height = 4)
par(mar=c(1, 3.1, 1.1, 1.1))
barplot(round((EUK.GC1.MDS.b$eig[EUK.GC1.MDS.b$eig>0]/sum(EUK.GC1.MDS.b$eig[EUK.GC1.MDS.b$eig>0]))*100,1),col=c("dodgerblue","pink3","darkgreen","gold",rep("grey",length(EUK.GC1.MDS.b$eig)-4)))
dev.off()
pdf("figures/EUK.GC1.jacc.MDS.linear.prop.pdf",width = 5,height = 4)
par(mar=c(1, 3.1, 1.1, 1.1))
barplot(round((EUK.GC1.MDS.j$eig[EUK.GC1.MDS.j$eig>0]/sum(EUK.GC1.MDS.j$eig[EUK.GC1.MDS.j$eig>0]))*100,1),col=c("dodgerblue","pink3","darkgreen","gold",rep("grey",length(EUK.GC1.MDS.j$eig)-4)))
dev.off()






## P19
#blacklist <- c("P19_000_4","P19_044_1","P19_064_5","P19_108_2","P19_100_1","P19_100_6","P19_104_1","P19_164_3","P19_168_2","P19_160_6","P19_060_3","P19_036_2","P19_092_2","P19_120_6")
#test <- EUK.P19[,!colnames(EUK.P19) %in% blacklist]

EUK.P19.2 <- EUK.P19[colSums(EUK.P19)!=0]

EUK.P19.nMDS.b <- metaMDS(t(prop.table(as.matrix(EUK.P19.2),2)),k=3,trymax = 200,parallel=8)
#EUK.P19.nMDS.b <- metaMDS(t(prop.table(as.matrix(EUK.P19[,!colnames(EUK.P19) %in% blacklist]),2)),k=3,trymax = 20)
EUK.P19.nMDS.b$points

## here we use the average dataset 
EUK.P19.MDS.b <- wcmdscale(vegdist(t(EUK.P19.avr)),eig = TRUE)
EUK.P19.MDS.j <- wcmdscale(vegdist(t(EUK.P19.avr),method = "jac",binary=TRUE),eig = TRUE)







pdf("figures/EUK.P19.bray.nMDS.pdf",width = 9,height = 3.5)
par(mfrow=c(1,3),mar=c(5.1, 4.1, 1.1, 1.1),xpd=TRUE)
plot(EUK.P19.nMDS.b$points[,1],EUK.P19.nMDS.b$points[,2],xlab="nMDS1",ylab="nMDS2",pch=16,col=map_values_to_colors(ages$median[match(gsub("(.*)_[0-9]$","\\1",rownames(EUK.P19.nMDS.b$points)),ages$ID2)], "dodgerblue", "darkred"))
plot_with_gradient_legend(ages$median[match(gsub("(.*)_[0-9]$","\\1",rownames(EUK.P19.nMDS.b$points)),ages$ID2)], "dodgerblue", "darkred")
plot(EUK.P19.nMDS.b$points[,1],EUK.P19.nMDS.b$points[,3],xlab="nMDS1",ylab="nMDS3",pch=16,col=map_values_to_colors(ages$median[match(gsub("(.*)_[0-9]$","\\1",rownames(EUK.P19.nMDS.b$points)),ages$ID2)], "dodgerblue", "darkred"))
plot(EUK.P19.nMDS.b$points[,2],EUK.P19.nMDS.b$points[,3],xlab="nMDS2",ylab="nMDS3",pch=16,col=map_values_to_colors(ages$median[match(gsub("(.*)_[0-9]$","\\1",rownames(EUK.P19.nMDS.b$points)),ages$ID2)], "dodgerblue", "darkred"))
dev.off()

pdf("figures/EUK.P19.bray.nMDS.linear.pdf",width = 8,height = 8)
par(mfrow=c(3,1),mar=c(5.1, 4.1, 1.1, 1.1),xpd=TRUE)
plot(ages$median[match(gsub("(.*)_[0-9]$","\\1",rownames(EUK.P19.nMDS.b$points)),ages$ID2)],EUK.P19.nMDS.b$points[,1],pch=16,ylab="nMDS 1",xlab="Yr Cal BP")
plot(ages$median[match(gsub("(.*)_[0-9]$","\\1",rownames(EUK.P19.nMDS.b$points)),ages$ID2)],EUK.P19.nMDS.b$points[,2],pch=16,ylab="nMDS 2",xlab="Yr Cal BP")
plot(ages$median[match(gsub("(.*)_[0-9]$","\\1",rownames(EUK.P19.nMDS.b$points)),ages$ID2)],EUK.P19.nMDS.b$points[,3],pch=16,ylab="nMDS 3",xlab="Yr Cal BP")
dev.off()

pdf("figures/EUK.P19.bray.nMDS.linear.short.pdf",width = 8,height = 8)
par(mfrow=c(3,1),mar=c(5.1, 4.1, 1.1, 1.1),xpd=TRUE)
plot(ages$median[match(gsub("(.*)_[0-9]$","\\1",rownames(EUK.P19.nMDS.b$points)),ages$ID2)],EUK.P19.nMDS.b$points[,1],pch=16,xlim=c(200,1850),ylab="nMDS 1",xlab="Yr Cal BP")
plot(ages$median[match(gsub("(.*)_[0-9]$","\\1",rownames(EUK.P19.nMDS.b$points)),ages$ID2)],EUK.P19.nMDS.b$points[,2],pch=16,xlim=c(200,1850),ylab="nMDS 2",xlab="Yr Cal BP")
plot(ages$median[match(gsub("(.*)_[0-9]$","\\1",rownames(EUK.P19.nMDS.b$points)),ages$ID2)],EUK.P19.nMDS.b$points[,3],pch=16,xlim=c(200,1850),ylab="nMDS 3",xlab="Yr Cal BP")
dev.off()


##### Example GAM 

gsub("(.*)_[0-9]$","\\1",rownames(EUK.P19.nMDS.b$points))

test <- tapply(ages$median[match(gsub("(.*)_[0-9]$","\\1",rownames(EUK.P19.nMDS.b$points)),ages$ID2)],
       INDEX=gsub("(.*)_[0-9]$","\\1",rownames(EUK.P19.nMDS.b$points)),FUN=mean)

example_dat  <- data.frame("example_x"=ages$median[match(gsub("(.*)_[0-9]$","\\1",rownames(EUK.P19.nMDS.b$points)),ages$ID2)],
                           "example_y"=EUK.P19.nMDS.b$points[,3])

example_dat_mean  <- data.frame("example_x"=tapply(ages$median[match(gsub("(.*)_[0-9]$","\\1",rownames(EUK.P19.nMDS.b$points)),ages$ID2)],INDEX=gsub("(.*)_[0-9]$","\\1",rownames(EUK.P19.nMDS.b$points)),FUN=mean),
                           "example_y"=tapply(EUK.P19.nMDS.b$points[,3],INDEX=gsub("(.*)_[0-9]$","\\1",rownames(EUK.P19.nMDS.b$points)),FUN=mean),
                           "example_y2"=tapply(EUK.P19.nMDS.b$points[,2],INDEX=gsub("(.*)_[0-9]$","\\1",rownames(EUK.P19.nMDS.b$points)),FUN=mean),
                           "example_y3"=tapply(EUK.P19.nMDS.b$points[,1],INDEX=gsub("(.*)_[0-9]$","\\1",rownames(EUK.P19.nMDS.b$points)),FUN=mean))

m <- gam(example_y ~ s(example_x, k = 20), data = example_dat, method = "REML")
prediction <- data.frame("example_x"=200:3000)
prediction <- cbind(prediction,predict(m,prediction,se.fit=TRUE))
deriv <- gratia::derivatives(m)


pdf("figures/EUK.GAM.all.pdf",width = 9,height = 7)
par(mfrow=c(2,1),mar=c(5.1, 4.1, 1.1, 1.1))

plot(example_dat$example_x,example_dat$example_y,pch=16,ylab="nMDS 3",xlab="Cal Yr BP",col="grey",cex=0.8)
polygon(c(prediction$example_x, rev(prediction$example_x)), c(prediction$fit+prediction$se.fit*1.96, rev(prediction$fit)+rev(prediction$se.fit*-1.96)), col=add.alpha('dodgerblue',0.5), border=NA)
points(prediction$example_x,prediction$fit,lwd=2,col="black",type='l')


plot(deriv$data,deriv$derivative,cex=0,ylab="Slope",xlab="Cal Yr BP")
polygon(c(deriv$data, rev(deriv$data)), c(deriv$upper, rev(deriv$lower)), col='gray', border=NA)
points(deriv$data,deriv$derivative,lwd=2,col="black",type='l')
abline(h=0,lty=2)
dev.off()



m <- gam(example_y ~ s(example_x, k = 30), data = example_dat_mean, method = "REML")
prediction <- data.frame("example_x"=200:3000)
prediction <- cbind(prediction,predict(m,prediction,se.fit=TRUE))
deriv <- gratia::derivatives(m)

m2 <- gam(example_y2 ~ s(example_x, k = 30), data = example_dat_mean, method = "REML")
prediction2 <- data.frame("example_x"=200:3000)
prediction2 <- cbind(prediction2,predict(m2,prediction2,se.fit=TRUE))
deriv2 <- gratia::derivatives(m2)


m3 <- gam(example_y3 ~ s(example_x, k = 30), data = example_dat_mean, method = "REML")
prediction3 <- data.frame("example_x"=200:3000)
prediction3 <- cbind(prediction3,predict(m3,prediction3,se.fit=TRUE))
deriv3 <- gratia::derivatives(m3)



pdf("figures/EUK.GAM.avr.pdf",width = 9,height = 7)
par(mfrow=c(2,1),mar=c(5.1, 4.1, 1.1, 1.1))

plot(example_dat_mean$example_x,example_dat_mean$example_y,pch=16,ylab="nMDS 3",xlab="Cal Yr BP",col="grey",cex=0.8)
polygon(c(prediction$example_x, rev(prediction$example_x)), c(prediction$fit+prediction$se.fit*1.96, rev(prediction$fit)+rev(prediction$se.fit*-1.96)), col=add.alpha('dodgerblue',0.5), border=NA)
points(prediction$example_x,prediction$fit,lwd=2,col="black",type='l')


plot(deriv$data,deriv$derivative,cex=0,ylab="Slope",xlab="Cal Yr BP")
polygon(c(deriv$data, rev(deriv$data)), c(deriv$upper, rev(deriv$lower)), col='gray', border=NA)
points(deriv$data,deriv$derivative,lwd=2,col="black",type='l')
abline(h=0,lty=2)
dev.off()


# all three nMDSs

m1 <- gam(example_y3 ~ s(example_x, k = 30), data = example_dat_mean, method = "REML")
prediction <- data.frame("example_x"=200:3000)
prediction <- cbind(prediction,predict(m1,prediction,se.fit=TRUE))
deriv <- gratia::derivatives(m1)

m2 <- gam(example_y2 ~ s(example_x, k = 30), data = example_dat_mean, method = "REML")
prediction2 <- data.frame("example_x"=200:3000)
prediction2 <- cbind(prediction2,predict(m2,prediction2,se.fit=TRUE))
deriv2 <- gratia::derivatives(m2)


m3 <- gam(example_y ~ s(example_x, k = 30), data = example_dat_mean, method = "REML")
prediction3 <- data.frame("example_x"=200:3000)
prediction3 <- cbind(prediction3,predict(m3,prediction3,se.fit=TRUE))
deriv3 <- gratia::derivatives(m3)


pdf("figures/EUK.GAM.slopeall.pdf",width = 9,heigh=6)
plot(deriv$data,deriv$derivative,lwd=2,col="dodgerblue",type='l',xlab="Cal Yr BP",ylab="Slope of GAM modelled component")
points(deriv2$data,deriv2$derivative,lwd=2,col="red",type='l')
points(deriv3$data,deriv3$derivative,lwd=2,col="darkgreen",type='l')


polygon(c(deriv$data, rev(deriv$data)), c(deriv$upper, rev(deriv$lower)), col=add.alpha('grey',0.3), border=NA)
polygon(c(deriv2$data, rev(deriv2$data)), c(deriv2$upper, rev(deriv2$lower)),col=add.alpha('grey',0.3), border=NA)
polygon(c(deriv3$data, rev(deriv3$data)), c(deriv3$upper, rev(deriv3$lower)), col=add.alpha('grey',0.3), border=NA)

dev.off()

## new stuff



## Lets make some GAMS

P19.MDS.b.gam.dat  <- data.frame("year"=ages$median[match(rownames(EUK.P19.MDS.b$points),ages$ID2)],
                                 "MDS1"=EUK.P19.MDS.b$points[,1],
                                 "MDS2"=EUK.P19.MDS.b$points[,2],
                                 "MDS3"=EUK.P19.MDS.b$points[,3],
                                 "MDS4"=EUK.P19.MDS.b$points[,4])

m1 <- gam(MDS1 ~ s(year,k=20), data = P19.MDS.b.gam.dat, method = "REML")
m2 <- gam(MDS2 ~ s(year, k = 20), data = P19.MDS.b.gam.dat, method = "REML")
m3 <- gam(MDS3 ~ s(year, k = 20), data = P19.MDS.b.gam.dat, method = "REML")
m4 <- gam(MDS4 ~ s(year, k = 20), data = P19.MDS.b.gam.dat, method = "REML")

prediction <- data.frame("year"=240:3100)
prediction1 <- cbind(prediction,predict(m1,prediction,se.fit=TRUE))
prediction2 <- cbind(prediction,predict(m2,prediction,se.fit=TRUE))
prediction3 <- cbind(prediction,predict(m3,prediction,se.fit=TRUE))
prediction4 <- cbind(prediction,predict(m4,prediction,se.fit=TRUE))
prediction1$uppCI <- prediction1$fit+prediction1$se.fit*1.96
prediction1$lwrCI <- prediction1$fit-prediction1$se.fit*1.96
prediction2$uppCI <- prediction2$fit+prediction2$se.fit*1.96
prediction2$lwrCI <- prediction2$fit-prediction2$se.fit*1.96
prediction3$uppCI <- prediction3$fit+prediction3$se.fit*1.96
prediction3$lwrCI <- prediction3$fit-prediction3$se.fit*1.96
prediction4$uppCI <- prediction4$fit+prediction4$se.fit*1.96
prediction4$lwrCI <- prediction4$fit-prediction4$se.fit*1.96



## plot



pdf("figures/EUK.P19.bray.MDS.linear.pdf",width = 8,height = 8)
par(mfrow=c(4,1),mar=c(2, 4.1, 1.1, 1.1))
plot(ages$median[match(rownames(EUK.P19.MDS.b$points),ages$ID2)],EUK.P19.MDS.b$points[,1],pch=16,ylab="MDS 1",xlab="Yr Cal BP")
polygon(c(prediction1$year, rev(prediction1$year)), c(prediction1$uppCI, rev(prediction1$lwrCI)), col=add.alpha('dodgerblue',0.3), border=NA)
points(prediction1$year,prediction1$fit,lwd=2,col="dodgerblue",type='l')
plot(ages$median[match(rownames(EUK.P19.MDS.b$points),ages$ID2)],EUK.P19.MDS.b$points[,2],pch=16,ylab="MDS 2",xlab="Yr Cal BP")
polygon(c(prediction2$year, rev(prediction2$year)), c(prediction2$uppCI, rev(prediction2$lwrCI)), col=add.alpha('pink3',0.3), border=NA)
points(prediction2$year,prediction2$fit,lwd=2,col="pink3",type='l')
plot(ages$median[match(rownames(EUK.P19.MDS.b$points),ages$ID2)],EUK.P19.MDS.b$points[,3],pch=16,ylab="MDS 3",xlab="Yr Cal BP")
polygon(c(prediction3$year, rev(prediction3$year)), c(prediction3$uppCI, rev(prediction3$lwrCI)), col=add.alpha('darkgreen',0.3), border=NA)
points(prediction3$year,prediction3$fit,lwd=2,col="darkgreen",type='l')
plot(ages$median[match(rownames(EUK.P19.MDS.b$points),ages$ID2)],EUK.P19.MDS.b$points[,4],pch=16,ylab="MDS 4",xlab="Yr Cal BP")
polygon(c(prediction4$year, rev(prediction4$year)), c(prediction4$uppCI, rev(prediction4$lwrCI)), col=add.alpha('gold',0.3), border=NA)
points(prediction4$year,prediction4$fit,lwd=2,col="gold",type='l')
dev.off()

##



## Lets make some GAMS

P19.MDS.j.gam.dat  <- data.frame("year"=ages$median[match(rownames(EUK.P19.MDS.j$points),ages$ID2)],
                                 "MDS1"=EUK.P19.MDS.j$points[,1],
                                 "MDS2"=EUK.P19.MDS.j$points[,2],
                                 "MDS3"=EUK.P19.MDS.j$points[,3],
                                 "MDS4"=EUK.P19.MDS.j$points[,4])

m1 <- gam(MDS1 ~ s(year,k=20), data = P19.MDS.j.gam.dat, method = "REML")
m2 <- gam(MDS2 ~ s(year, k = 20), data = P19.MDS.j.gam.dat, method = "REML")
m3 <- gam(MDS3 ~ s(year, k = 20), data = P19.MDS.j.gam.dat, method = "REML")
m4 <- gam(MDS4 ~ s(year, k = 20), data = P19.MDS.j.gam.dat, method = "REML")

prediction <- data.frame("year"=240:3100)
prediction1 <- cbind(prediction,predict(m1,prediction,se.fit=TRUE))
prediction2 <- cbind(prediction,predict(m2,prediction,se.fit=TRUE))
prediction3 <- cbind(prediction,predict(m3,prediction,se.fit=TRUE))
prediction4 <- cbind(prediction,predict(m4,prediction,se.fit=TRUE))
prediction1$uppCI <- prediction1$fit+prediction1$se.fit*1.96
prediction1$lwrCI <- prediction1$fit-prediction1$se.fit*1.96
prediction2$uppCI <- prediction2$fit+prediction2$se.fit*1.96
prediction2$lwrCI <- prediction2$fit-prediction2$se.fit*1.96
prediction3$uppCI <- prediction3$fit+prediction3$se.fit*1.96
prediction3$lwrCI <- prediction3$fit-prediction3$se.fit*1.96
prediction4$uppCI <- prediction4$fit+prediction4$se.fit*1.96
prediction4$lwrCI <- prediction4$fit-prediction4$se.fit*1.96



## plot



pdf("figures/EUK.P19.jacc.MDS.linear.pdf",width = 8,height = 8)
par(mfrow=c(4,1),mar=c(2, 4.1, 1.1, 1.1))
plot(ages$median[match(rownames(EUK.P19.MDS.j$points),ages$ID2)],EUK.P19.MDS.j$points[,1],pch=16,ylab="MDS 1",xlab="Yr Cal BP")
polygon(c(prediction1$year, rev(prediction1$year)), c(prediction1$uppCI, rev(prediction1$lwrCI)), col=add.alpha('dodgerblue',0.3), border=NA)
points(prediction1$year,prediction1$fit,lwd=2,col="dodgerblue",type='l')
plot(ages$median[match(rownames(EUK.P19.MDS.j$points),ages$ID2)],EUK.P19.MDS.j$points[,2],pch=16,ylab="MDS 2",xlab="Yr Cal BP")
polygon(c(prediction2$year, rev(prediction2$year)), c(prediction2$uppCI, rev(prediction2$lwrCI)), col=add.alpha('pink3',0.3), border=NA)
points(prediction2$year,prediction2$fit,lwd=2,col="pink3",type='l')
plot(ages$median[match(rownames(EUK.P19.MDS.j$points),ages$ID2)],EUK.P19.MDS.j$points[,3],pch=16,ylab="MDS 3",xlab="Yr Cal BP")
polygon(c(prediction3$year, rev(prediction3$year)), c(prediction3$uppCI, rev(prediction3$lwrCI)), col=add.alpha('darkgreen',0.3), border=NA)
points(prediction3$year,prediction3$fit,lwd=2,col="darkgreen",type='l')
plot(ages$median[match(rownames(EUK.P19.MDS.j$points),ages$ID2)],EUK.P19.MDS.j$points[,4],pch=16,ylab="MDS 4",xlab="Yr Cal BP")
polygon(c(prediction4$year, rev(prediction4$year)), c(prediction4$uppCI, rev(prediction4$lwrCI)), col=add.alpha('gold',0.3), border=NA)
points(prediction4$year,prediction4$fit,lwd=2,col="gold",type='l')
dev.off()



pdf("figures/EUK.P19.bray.MDS.linear.prop.pdf",width = 5,height = 4)
par(mar=c(1, 3.1, 1.1, 1.1))
barplot(round((EUK.P19.MDS.b$eig[EUK.P19.MDS.b$eig>0]/sum(EUK.P19.MDS.b$eig[EUK.P19.MDS.b$eig>0]))*100,1),col=c("dodgerblue","pink3","darkgreen","gold",rep("grey",length(EUK.P19.MDS.b$eig)-4)))
dev.off()
pdf("figures/EUK.P19.jacc.MDS.linear.prop.pdf",width = 5,height = 4)
par(mar=c(1, 3.1, 1.1, 1.1))
barplot(round((EUK.P19.MDS.j$eig[EUK.P19.MDS.j$eig>0]/sum(EUK.P19.MDS.j$eig[EUK.P19.MDS.j$eig>0]))*100,1),col=c("dodgerblue","pink3","darkgreen","gold",rep("grey",length(EUK.P19.MDS.j$eig)-4)))
dev.off()



#### here let's try some dbRDA


vegdist(t(EUK.p19.avr.clim))

ages$mean

colnames(EUK.P19.avr)
dbrda()


EUK.p19.avr.clim <- EUK.P19.avr[,ages$mean[match(colnames(EUK.P19.avr),ages$ID2)]<995 & ages$mean[match(colnames(EUK.P19.avr),ages$ID2)]>222]

colnames(EUK.p19.avr.clim)

PC019.climate <- climate[match(colnames(EUK.p19.avr.clim),gsub("-","_",climate$Sample)),]


Pc019.dbRDA.b <- dbrda(vegdist(t(EUK.p19.avr.clim))~d13C.GAM.fit+d13C.100yrspline+d18O.GAM.fit+d18O.100yrspline,PC019.climate)
Pc019.dbRDA.j <- dbrda(vegdist(t(EUK.p19.avr.clim),method="jaccard",binary = TRUE)~d13C.GAM.fit+d13C.100yrspline+d18O.GAM.fit+d18O.100yrspline,PC019.climate)



plot(Pc019.dbRDA.b)
plot(Pc019.dbRDA.j)
anova(Pc019.dbRDA.b,by="margin",permutations = 10000)
anova(Pc019.dbRDA.j,by="margin",permutations = 10000)


test <- metaMDS(vegdist(t(EUK.p19.avr.clim),method="jaccard",binary = TRUE),trymax = 200)
out.test <- envfit(test,PC019.climate)
plot(test)
ordisurf(test,PC019.climate$d18O.100yrspline)


test <- metaMDS(vegdist(t(EUK.p19.avr.clim)),trymax = 200)
out.test <- envfit(test,PC019.climate)
plot(test)
ordisurf(test,PC019.climate$d18O.100yrspline)


### here lets try some beta diversity plots of data subsets (metazoa/protist/bacteria)

#met
EUK.P19.nMDS.met.b <- metaMDS(vegdist(t(EUK.p19.avr.met[,match(gsub("-","_",climate_subset$Sample),colnames(EUK.p19.avr.met))])),k=2,trymax = 1000,parallel=8)
EUK.P19.nMDS.met.j <- metaMDS(vegdist(t(EUK.p19.avr.met[,match(gsub("-","_",climate_subset$Sample),colnames(EUK.p19.avr.met))]),binary = TRUE,method="jaccard"),k=2,trymax = 1000,parallel=8)
EUK.P19.MDS.met.b <- wcmdscale(vegdist(t(EUK.p19.avr.met[,match(gsub("-","_",climate_subset$Sample),colnames(EUK.p19.avr.met))])),eig = TRUE)
EUK.P19.MDS.met.j <- wcmdscale(vegdist(t(EUK.p19.avr.met[,match(gsub("-","_",climate_subset$Sample),colnames(EUK.p19.avr.met))]),binary = TRUE,method="jaccard"),eig = TRUE)

pdf("figures/climate/EUK.P19.climate.ord.met.pdf",width = 9 ,height = 9)
par(mfrow=c(2,2))
plot(EUK.P19.nMDS.met.b$points[,1],EUK.P19.nMDS.met.b$points[,2],col=map_values_to_colors(ages$mean[match(rownames(EUK.P19.nMDS.met.b$points),ages$ID2)], "dodgerblue", "darkred"),pch=16)
plot(EUK.P19.nMDS.met.j$points[,1],EUK.P19.nMDS.met.j$points[,2],col=map_values_to_colors(ages$mean[match(rownames(EUK.P19.nMDS.met.j$points),ages$ID2)], "dodgerblue", "darkred"),pch=16)
plot(EUK.P19.MDS.met.b$points[,1],EUK.P19.MDS.met.b$points[,2],col=map_values_to_colors(ages$mean[match(rownames(EUK.P19.MDS.met.b$points),ages$ID2)], "dodgerblue", "darkred"),pch=16)
plot(EUK.P19.MDS.met.j$points[,1],EUK.P19.MDS.met.j$points[,2],col=map_values_to_colors(ages$mean[match(rownames(EUK.P19.MDS.met.j$points),ages$ID2)], "dodgerblue", "darkred"),pch=16)
dev.off()

envfit(EUK.P19.nMDS.met.b,datasetin)
envfit(EUK.P19.nMDS.met.j,datasetin)
envfit(EUK.P19.MDS.met.b ,datasetin)
envfit(EUK.P19.MDS.met.j ,datasetin)

#pro
EUK.P19.nMDS.pro.b <- metaMDS(vegdist(t(EUK.p19.avr.pro[,match(gsub("-","_",climate_subset$Sample),colnames(EUK.p19.avr.pro))])),k=2,trymax = 1000,parallel=8)
EUK.P19.nMDS.pro.j <- metaMDS(vegdist(t(EUK.p19.avr.pro[,match(gsub("-","_",climate_subset$Sample),colnames(EUK.p19.avr.pro))]),binary = TRUE,method="jaccard"),k=2,trymax = 1000,parallel=8)
EUK.P19.MDS.pro.b <- wcmdscale(vegdist(t(EUK.p19.avr.pro[,match(gsub("-","_",climate_subset$Sample),colnames(EUK.p19.avr.pro))])),eig = TRUE)
EUK.P19.MDS.pro.j <- wcmdscale(vegdist(t(EUK.p19.avr.pro[,match(gsub("-","_",climate_subset$Sample),colnames(EUK.p19.avr.pro))]),binary = TRUE,method="jaccard"),eig = TRUE)

pdf("figures/climate/EUK.P19.climate.ord.pro.pdf",width = 9 ,height = 9)
par(mfrow=c(2,2))
plot(EUK.P19.nMDS.pro.b$points[,1],EUK.P19.nMDS.pro.b$points[,2],col=map_values_to_colors(ages$mean[match(rownames(EUK.P19.nMDS.pro.b$points),ages$ID2)], "dodgerblue", "darkred"),pch=16)
plot(EUK.P19.nMDS.pro.j$points[,1],EUK.P19.nMDS.pro.j$points[,2],col=map_values_to_colors(ages$mean[match(rownames(EUK.P19.nMDS.pro.j$points),ages$ID2)], "dodgerblue", "darkred"),pch=16)
plot(EUK.P19.MDS.pro.b$points[,1],EUK.P19.MDS.pro.b$points[,2],col=map_values_to_colors(ages$mean[match(rownames(EUK.P19.MDS.pro.b$points),ages$ID2)], "dodgerblue", "darkred"),pch=16)
plot(EUK.P19.MDS.pro.j$points[,1],EUK.P19.MDS.pro.j$points[,2],col=map_values_to_colors(ages$mean[match(rownames(EUK.P19.MDS.pro.j$points),ages$ID2)], "dodgerblue", "darkred"),pch=16)
dev.off()

envfit(EUK.P19.nMDS.pro.b,datasetin)
envfit(EUK.P19.nMDS.pro.j,datasetin)
envfit(EUK.P19.MDS.pro.b ,datasetin)
envfit(EUK.P19.MDS.pro.j ,datasetin)

#bac
EUK.P19.nMDS.bac.b <- metaMDS(vegdist(t(EUK.p19.avr.bac[,match(gsub("-","_",climate_subset$Sample),colnames(EUK.p19.avr.bac))])),k=2,trymax = 1000,parallel=8)
EUK.P19.nMDS.bac.j <- metaMDS(vegdist(t(EUK.p19.avr.bac[,match(gsub("-","_",climate_subset$Sample),colnames(EUK.p19.avr.bac))]),binary = TRUE,method="jaccard"),k=2,trymax = 1000,parallel=8)
EUK.P19.MDS.bac.b <- wcmdscale(vegdist(t(EUK.p19.avr.bac[,match(gsub("-","_",climate_subset$Sample),colnames(EUK.p19.avr.bac))])),eig = TRUE)
EUK.P19.MDS.bac.j <- wcmdscale(vegdist(t(EUK.p19.avr.bac[,match(gsub("-","_",climate_subset$Sample),colnames(EUK.p19.avr.bac))]),binary = TRUE,method="jaccard"),eig = TRUE)

pdf("figures/climate/EUK.P19.climate.ord.bac.pdf",width = 9 ,height = 9)
par(mfrow=c(2,2))
plot(EUK.P19.nMDS.bac.b$points[,1],EUK.P19.nMDS.bac.b$points[,2],col=map_values_to_colors(ages$mean[match(rownames(EUK.P19.nMDS.bac.b$points),ages$ID2)], "dodgerblue", "darkred"),pch=16)
plot(EUK.P19.nMDS.bac.j$points[,1],EUK.P19.nMDS.bac.j$points[,2],col=map_values_to_colors(ages$mean[match(rownames(EUK.P19.nMDS.bac.j$points),ages$ID2)], "dodgerblue", "darkred"),pch=16)
plot(EUK.P19.MDS.bac.b$points[,1],EUK.P19.MDS.bac.b$points[,2],col=map_values_to_colors(ages$mean[match(rownames(EUK.P19.MDS.bac.b$points),ages$ID2)], "dodgerblue", "darkred"),pch=16)
plot(EUK.P19.MDS.bac.j$points[,1],EUK.P19.MDS.bac.j$points[,2],col=map_values_to_colors(ages$mean[match(rownames(EUK.P19.MDS.bac.j$points),ages$ID2)], "dodgerblue", "darkred"),pch=16)
dev.off()

envfit(EUK.P19.nMDS.bac.b,datasetin)
envfit(EUK.P19.nMDS.bac.j,datasetin)
envfit(EUK.P19.MDS.bac.b ,datasetin)
envfit(EUK.P19.MDS.bac.j ,datasetin)





### whale biz

plot(ages$median[match(colnames(MAM.P19.nREPS[245,]),ages$ID2)],jitter(as.numeric(MAM.P19.nREPS[245,])),pch=16,xlab="Cal Yr BP",ylab="replicates +ive")

plot(ages$median[match(colnames(RIZ.P19.nREPS[3,]),ages$ID2)],jitter(as.numeric(RIZ.P19.nREPS[43,])),pch=16,xlab="Cal Yr BP",ylab="replicates +ive")


#right whale ASV number
c(10,169,230,293,336,392,490,491,502,504,514,517,534,565,576,594,612)

plot(ages$median[match(colnames(MAM.P19.nREPS[c(10,169,230,293,336,392,490,491,502,504,514,517,534,565,576,594,612),]),ages$ID2)],
     jitter(as.numeric(colSums(MAM.P19.nREPS[c(10,169,230,293,336,392,490,491,502,504,514,517,534,565,576,594,612),],na.rm = TRUE))),pch=16,xlab="Cal Yr BP",ylab="replicates +ive")


summary(test)
plot(test)


pdf("figures/RIZ.ASV2.pdf",width = 6,height = 6)
plot(ages$median[match(colnames(MAM.P19.nREPS[2,]),ages$ID2)],jitter(as.numeric(MAM.P19.nREPS[2,])),pch=16,xlab="Cal Yr BP",ylab="replicates +ive")
dev.off()

pdf("figures/RIZ.ASV3.pdf",width = 6,height = 6)
plot(ages$median[match(colnames(MAM.P19.nREPS[3,]),ages$ID2)],jitter(as.numeric(MAM.P19.nREPS[3,])),pch=16,xlab="Cal Yr BP",ylab="replicates +ive")
dev.off()

pdf("figures/RIZ.ASV25.pdf",width = 6,height = 6)
plot(ages$median[match(colnames(MAM.P19.nREPS[25,]),ages$ID2)],jitter(as.numeric(MAM.P19.nREPS[4,])),pch=16,xlab="Cal Yr BP",ylab="replicates +ive")
dev.off()


## whale GLM

year <- ages$median[match(colnames(MAM.P19.nREPS[3,]),ages$ID2)]
detections <- as.numeric(MAM.P19.nREPS[3,])

test <- glm(detections ~year,poisson(link = "log"))
summary(test)

prediction <- data.frame("year"=200:3000)
prediction <- cbind(prediction,predict(test,prediction,se.fit=TRUE))
plot(prediction$year,prediction$fit)


