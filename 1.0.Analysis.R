################################################################
####======SeaChange Icelandic Dataset Analysis ==========#######
####==== Luke E. Holman====15.04.2024====================#######
################################################################

####====0.1 Packages====####
library(vegan)
library(RColorBrewer)
library(breakaway)


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

minAbundance <- function(inputtable=NA,minAbun= 0.01){
  inputtable <- rbind(inputtable,rep(0, ncol(inputtable)))
  rownames(inputtable)[dim(inputtable)[1]] <- "Others"
  for (row in 1:dim(inputtable)[2]){
    min <- sum(inputtable[,row])*minAbun
    others <- sum(inputtable[inputtable[,row]<min,row])
    inputtable[inputtable[,row]<min,row] <- 0
    inputtable["Others",row] <- others
    inputtable <- inputtable[rowSums(inputtable)>1,]
  }
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

####====0.3 Data Prep====####
ages <- read.csv("metadata/AgeOut.csv")
ages$ID2 <- gsub("-","_",ages$ID)
EUK.tax.PR2 <- read.csv("taxonomy/EUK.cleaned.PR2.csv",row.names = 1)

## EUK datasets
EUK.GC1 <- read.csv("cleaneddata/combinedcodarkredata/EUK.GC01.csv",row.names = 1)
EUK.GC1.nREPS <- NrepsMaker(EUK.GC1,gsub("(.*)_[0-9]$","\\1",colnames(EUK.GC1)))
EUK.GC1.avr <- average_by_reps(EUK.GC1,gsub("(.*)_[0-9]$","\\1",colnames(EUK.GC1)))
EUK.GC1.bin <- make_binary(EUK.GC1,1)
EUK.GC1.ages <- ages$median[match(colnames(EUK.GC1.avr),ages$ID2)]

EUK.P19 <- read.csv("cleaneddata/combinedcodarkredata/EUK.PC019.csv",row.names = 1)
## make a subset of the data gettign rid of high resolution volcano samples 
ages$ID2[184:204]


EUK.P19.nREPS <- NrepsMaker(EUK.P19,gsub("(.*)_[0-9]$","\\1",colnames(EUK.P19)))
EUK.P19.avr <- average_by_reps(EUK.P19,gsub("(.*)_[0-9]$","\\1",colnames(EUK.P19)))
EUK.P19.bin <- make_binary(EUK.P19,1)
EUK.P19.ages <- ages$median[match(colnames(EUK.P19.avr),ages$ID2)]

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


## attempts here to change the spacing according to time
ages$median[match(colnames(EUK.GC1.tax.a),ages$ID2)]
axis(1, at=test, labels=ages$median[match(colnames(EUK.GC1.tax.a),ages$ID2)])

ages$median[match(colnames(EUK.GC1.tax.a),ages$ID2)]
positions <- c(0, diff(ages$median[match(colnames(EUK.GC1.tax.a),ages$ID2)]))
bar_widths <- 5
barplot(EUK.GC1.tax.a,
        space = spaces/10,
        width = bar_widths,
        xlim= c(1,860))
        #xlim = c(0, max(ages$median[match(colnames(EUK.GC1.tax.a),ages$ID2)]))*bar_widths)


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


####====3.0 Beta diversity ====####

colSums(prop.table(as.matrix(EUK.GC1),2))


## GC1
blacklist <- c("GC1_000_4","GC1_044_1","GC1_064_5","GC1_108_2","GC1_100_1","GC1_100_6","GC1_104_1","GC1_164_3","GC1_168_2","GC1_160_6","GC1_060_3","GC1_036_2","GC1_092_2","GC1_120_6",)
blacklist <- c("GC1_148_1","GC1_180_7","GC1_168_2","GC1_104_1","GC1_100_6","GC1_160_6","GC1_164_3","GC1_108_2","GC1_044_1","GC1_116_1","GC1_116_2","GC1_116_3","GC1_116_4","GC1_116_5","GC1_116_6","GC1_116_7","GC1_116_8","GC1_064_5","GC1_000_4","GC1_140_1","GC1_140_2","GC1_140_3","GC1_140_4","GC1_140_5","GC1_140_6","GC1_140_7","GC1_140_8")

test <- EUK.GC1[,!colnames(EUK.GC1) %in% blacklist]
EUK.GC1.nMDS.b <- metaMDS(t(prop.table(as.matrix(test),2)),k=3,trymax = 20)
EUK.GC1.nMDS.b <- metaMDS(t(prop.table(as.matrix(EUK.GC1[,!colnames(EUK.GC1) %in% blacklist]),2)),k=3,trymax = 20)
EUK.GC1.nMDS.b$points


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



## P19
blacklist <- c("P19_000_4","P19_044_1","P19_064_5","P19_108_2","P19_100_1","P19_100_6","P19_104_1","P19_164_3","P19_168_2","P19_160_6","P19_060_3","P19_036_2","P19_092_2","P19_120_6")

test <- EUK.P19[,!colnames(EUK.P19) %in% blacklist]
EUK.P19.nMDS.b <- metaMDS(t(prop.table(as.matrix(EUK.P19.avr),2)),k=3,trymax = 20)
EUK.P19.nMDS.b <- metaMDS(t(prop.table(as.matrix(EUK.P19[,!colnames(EUK.P19) %in% blacklist]),2)),k=3,trymax = 20)
EUK.P19.nMDS.b$points


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

