#!/usr/bin/R
########################################################################
# Alexey Pindyurin, Anton Ivankin, September 12, 2014, DAM_count_statistics.R
# Updated for brevity by Artem Ilyin, Jan-Dec 2016
# DESCRIPTION:
#   
#
# DATA:
#   To work with this script you need to:
#     a) name your folder for output via 'prefixDir' variable;
#     b) set the location of your Rdata files via 'sourceDir' variable. You can specify the highest folder
#        as it is possible. Searching runs recursively;
#     c) set location of your damid_description.csv file which lists
#        all the previous datasets for DamID-seq via 'damIdLocation' variable;
#     d) choose if you want to use all data or only edge reads via boolean variable 'onlyEdge';
#   
#
# OUTPUT:
#   
#
# VERSIONS:
#   140912: First revision!
#   210116: Added combining inner + edge; added MakeSamplesListFile function
#   260116: Added generation of single protein wigs an gffs for inner + edge experiments
#   280116: Reworked to generate bedgraphs instead of wigs.
#   xx1116: Many cosmetic changes with the use of dplyr package, added option to use heterochromatin,
#           note that your gff file for htseq-count must include heterochromatin if you want to use this
#           option.
#   131217: Scatter Plots via ggally package. Very slow but quite pretty.
#
#     
########################################################################
rm(list=ls())
library(gplots)
library(dplyr)
library(ggplot2)
library(GGally)

# Declare variables
###################
prefixDir <- "YourDir" # directory for other experiments
workDir <- getwd()	# working directory (WD)
foldersNames <- c("Bedgraph", "Statistics", "CSV") # names of folders for different objects created
outputDirs <- sapply (foldersNames, function(x) file.path(workDir, prefixDir, x))
binsFile <- paste0(workDir, "/bins300het.txt")	# location of your GATCs/bins file
sourceDir <- "YourRData" # location of your RData files. You can specify the highest folder as it is possible. Searching runs recursively.
damIdLocation <- "Your_damid_description.csv" # location of your DamID-Description file
#
# Some clarification about format of the damid_description file
# It has to be formatted like this:
#     TISSUE.PROTEIN.conditions.#_of_replica\tname of the fastq.gz file, which was used for mapping and htseq-counting
# e.g. BRAIN.PIWI.vasa(-).1   P155_CGGATG_Piwi1.fastq.gz
#      BRAIN.PIWI.vasa(-).2   P155_ATTGCC_Piwi2.fastq.gz
#      TEST.DAM.wt.1  D1_AGGTTA_LR1_R001.fastq.gz
# and so on
onlyEdge <- F
usePseudoCounts <- T	# do you want to add pseudo counts into source data (T or F)? Default is "T"
pseudoCounts <- c(1)		# the vector of pseudo counts, default "c(1)"
heatmapColors <- greenred(200)	# color scheme for heatmap

bins <- read.delim(binsFile)
lastcol <- ncol(bins)

lapply(outputDirs, dir.create, showWarnings = FALSE, recursive = T) # 


################ LOAD FUNCTIONS #################
source("DamID_count.functions.R")

# Make samples list file
MakeSamplesListFileIE(sourceDir, damIdLocation)


# Load GATC counts in data frame
################################

if (onlyEdge == T) {
    SamplesList <- samplesList[grep("edge" ,SamplesList$id), ]
  } else {
    modS <- SamplesList[grep("edge", SamplesList$id), ]
    modS$id <- sub("(.+)(edge)", paste0("\\1", "all"), modS$id, perl=T)
    
}



bins <- cbind(bins, matrix(data=NA, nrow=nrow(bins), ncol=nrow(SamplesList)))
for (i in 1:nrow(SamplesList)){
  colnames(bins)[lastcol+i] <- SamplesList$id[i]
  load(file=SamplesList$path[i])
  if (all(bins$ID == reads2GATC$ID)) bins[, lastcol + i] <- reads2GATC$count
}
rm(i)

if (!onlyEdge) {
    modG <- bins[, c(1:lastcol, grep("edge", names(bins)))]
    names(modG)[(lastcol + 1):ncol(modG)] <- gsub("(.+)(edge)", paste0("\\1", "all"),
                                                  names(modG)[(lastcol + 1):ncol(modG)], perl=T)
    for (enzyme in names(modG[(lastcol + 1):ncol(modG)])) {
      S <- gsub("(.+)(all)", paste0("\\1", "edge"), enzyme, perl=T)
      E <- gsub("(.+)(all)", paste0("\\1", "inner"), enzyme, perl=T)
      modG[[enzyme]] <- bins[[S]] + bins[[E]]
    }
    bins <- modG
    samplesList <- modS
    # rm(modS, modG)
} else {
  names(bins)[(lastcol + 1):ncol(bins)] <- sub("(.+\\.[0-9]?)$", "\\1_all",
                                               names(bins)[(lastcol + 1):ncol(bins)], perl=T)
  samplesList <- modS
}

save.gatc.df <- "01.Raw.Counts.csv"
write.table(bins, file=file.path(prefixDir, "CSV", save.gatc.df), sep=";",
            row.names=F, col.names=T, quote=F, dec=".", append=F)



DATA <- bins

# Count statistics
###################
  chrs <- unique(DATA$chr)
  DATA.only <- DATA[, (lastcol + 1):ncol(DATA)]
  stat <- as.data.frame(matrix(data=NA, nrow=length(chrs), ncol=ncol(DATA.only) + lastcol,
                               byrow=F, dimnames=NULL))
  names(stat) <- c("chr", "bins.number", "chr.length.bp", "chr.length.proportion", colnames(DATA.only)[1:(ncol(DATA.only))])
  stat$chr <- chrs
  if (any(grepl("Het", DATA$chr))){
    stat$chr.length.bp <- c(23011544, 21146708, 24543557, 27905053, 1351857, 22422827, 368872, 3288761, 2555491, 2517507, 204112, 347038)
  } else {
    stat$chr.length.bp <- c(23011544, 21146708, 24543557, 27905053, 22422827)
  }
  genome.length <- sum(stat$chr.length.bp)
  stat$chr.length.proportion <- round(100 * stat$chr.length.bp / genome.length, digits=2)
  for (j in 1:(ncol(DATA.only))){
     for (i in 1:length(chrs)){
        Data.only.chr <- DATA.only[(DATA$chr == chrs[i]), j]
        if (j == 1) stat$bins.number[i] <- length(Data.only.chr)
        stat[i, 4+j] <- sum(Data.only.chr)
        rm(Data.only.chr)
     }
     rm(i)
  }
  rm(j)
  
  statistics.a <- "02.Raw.Counts.Statistics.csv"
  
  write.table(stat, file=file.path(prefixDir, "CSV", statistics.a), sep=";",
              row.names=F, col.names=T, quote=F, dec=".", append=F)
	
  for (j in 1:(ncol(DATA.only))){
     totalCounts <- sum(stat[, 4+j])
     for (i in 1:length(chrs)){
        stat[i, 4+j] <- round(100 * stat[i, 4+j] / totalCounts, digits=2)
     }
     rm(i)
     rm(totalCounts)
  }
  rm(j)
  statistics.b <- "03.Statistics. Proportions.csv"
  write.table(stat, file=file.path(prefixDir, "CSV", statistics.b), sep=";",
              row.names=F, col.names=T, quote=F, dec=".", append=F)

# Add Pseudo counts
###################
 DATAs <- list(DATA=DATA)
if (usePseudoCounts == T) {
 for ( i in pseudoCounts) {
  num <- sub("^([0-1]*)(.?)([0-1]*$)", "\\1\\3", i)
  DATA.pseudo <- assign(paste("pseudo", num, sep=""), DATA)
  DATA.pseudo[, (lastcol + 1):ncol(DATA.pseudo)] <- DATA[, (lastcol + 1):ncol(DATA)] + i
  pseudo.filename <- assign(paste("pseudo.fn", num, sep=""), paste0("04.Pseudo.", num, "_Added.csv"))
  DATA.pseudo.strname <- assign(paste("pseudo", num, sep=""), paste("pseudo", num, sep=""))
  write.table(DATA.pseudo, file=file.path(prefixDir, "CSV", pseudo.filename), sep=";",
              row.names=F, col.names=T, quote=F, dec=".", append=F)
  
  DATAs[[DATA.pseudo.strname]] <- DATA.pseudo
 }
rm(i)
}

# Correlation on Counts
#######################
MainCorrelations(data=DATA, createPDF=T, corr.desc = "01.Raw_data")
print("Run calculate reads per million")

# Declare variables
DATAs.rpm <- DATAs
DATAs.norm <- DATAs
DATAs.norm.ave <- DATAs
####################################

for (name in names(DATAs.rpm)) {


# Calculation reads per million
###############################
for (i in (lastcol + 1):(ncol(DATAs.rpm[[name]]))){
    column.sum <- sum(DATAs.rpm[[name]][, i])
    DATAs.rpm[[name]][, i] <- DATAs.rpm[[name]][, i] / column.sum * 10^6
    rm(column.sum)
  }
  rm(i)
  calc.rpm.file <- paste0("05.RPMs_", name, ".csv")
  write.table(DATAs.rpm[[name]], file=file.path(prefixDir, "CSV", calc.rpm.file), sep=";",
              row.names=F, col.names=T, quote=F, dec=".", append=F)

# Correlation on Channels
#########################
MainCorrelations(data=DATAs.rpm[[name]], corr.desc = paste0("02.RPMs.", name), createPDF=T)

# Plot boxplots on RPMs
###########################
  bmp(filename=file.path(prefixDir, "Statistics", paste0("03.RPMs.Boxplot.", name, ".bmp")),
      width=2000, height=1000, units="px")
  par(mar=c(12, 8, 0.5, 0.5))
  boxplot(DATAs.rpm[[name]][, (lastcol + 1):(ncol(DATAs.rpm[[name]]))],
          names=colnames(DATAs.rpm[[name]])[(lastcol + 1):(ncol(DATAs.rpm[[name]]))], las=2, ylab="RPM")
  dev.off()
  
# Generate BedGraphs from RPM data
####################################
  MakeBedGraphFromDATA(data=DATAs.rpm[[name]], descr = paste0("single.", name))

# DAM Normalization
###################

DATAs.norm[[name]] <- DATAs.norm[[name]][, -c((lastcol + 1):ncol(DATAs.norm[[name]]))]
listNorm <- samplesList
listNorm$normalization <- paste0(listNorm$tissue, listNorm$cond, listNorm$rep, sep=".")
uniqueSamples <- unique(listNorm$normalization)
  for (sample in uniqueSamples) {
    tissue.id <- subset(subset(listNorm, normalization == sample), protein != "DAM")$id
    dam.id <- subset(subset(listNorm, normalization == sample), protein == "DAM")$id
    for (protein in tissue.id) {
      tissue.norm <- paste(protein, ".norm", sep="")
      DATAs.norm[[name]][[tissue.norm]] <- log2(DATAs.rpm[[name]][[protein]] / DATAs.rpm[[name]][[dam.id]])
    }
  }
  for (i in (lastcol + 1):(ncol(DATAs.norm[[name]]))){
    nan.index <- is.nan(DATAs.norm[[name]][, i])
    inf.index <- is.infinite(DATAs.norm[[name]][, i])
    DATAs.norm[[name]][nan.index, i] <- NA
    DATAs.norm[[name]][inf.index, i] <- NA
    rm(nan.index)
    rm(inf.index)
  }
  rm(i)
dam.norm <- paste0("06.Dam.Normalized.", name, ".csv")
write.table(DATAs.norm[[name]], file=file.path(prefixDir, "CSV", dam.norm), sep=";",
            row.names=F, col.names=T, quote=F, dec=".", append=F)

} 
rm(name)


for (name in names(DATAs.norm)) {
# Correlation on Normalized data
################################
MainCorrelations(data=DATAs.norm[[name]], use.opt="pairwise.complete.obs",
                 corr.desc = paste0("04.Dam.Normalized.", name), createPDF=T)

# Scatter Plots on Normalized data
##################################
ScatCor(data=DATAs.norm[[name]], pref=paste0("05.Dam.Normalized.", name))

# Averaging Replicates Only for two replicates
######################
DATAs.norm.ave[[name]] <- DATAs.norm.ave[[name]][, 1:lastcol]
listNormAve <- listNorm[!(listNorm$protein == "DAM"), ]  # remove rows with DAM
listNormAve$normalizationAve <- paste(listNormAve$tissue, listNormAve$protein, listNormAve$cond, sep=".")
uniqueAveSamples <- unique(listNormAve$normalizationAve)
	for (item in uniqueAveSamples) {
		item.norm.ave <- paste(item, ".norm.ave", sep="")
		first.item.id <- paste(subset(subset(listNormAve, normalizationAve == item), rep == 1)$id, ".norm", sep="")
    second.item.id <- paste(subset(subset(listNormAve, normalizationAve == item), rep != 1)$id, ".norm", sep="")
		DATAs.norm.ave[[name]][[item.norm.ave]] <- log2Mean(DATAs.norm[[name]][[first.item.id]], DATAs.norm[[name]][[second.item.id]])
	}
dam.norm.ave <- paste0("07.Normalized.Mean", name, ".csv")
write.table(DATAs.norm.ave[[name]], file=file.path(prefixDir, "CSV", dam.norm.ave), sep=";",
            row.names=F, col.names=T, quote=F, dec=".", append=F)

 


# ACF plots on Averaged
#######################
# AcfOnData(dataSet=DATAs.norm.ave[[name]], labelAcf=labelAcf, method="acf", suffixPDF="14_ACF_Averaged", ylab.val="ACF on rpms", na.data=T)

# Density on Averaged
#######################
# AcfOnData(dataSet=DATAs.norm.ave[[name]], labelAcf=labelAcf, method="density", suffixPDF="15_Density_Averaged", ylab.val="density", na.data=T)

# Generate BedGraphs on Averaged data
###################################
MakeBedGraphFromDATA(data=DATAs.norm.ave[[name]], descr = name)

############################################################################################
if ((ncol(DATAs.norm.ave[[name]])-lastcol) != 1) {
# Correlations on Averaged NormData
###################################
MainCorrelations(data=DATAs.norm.ave[[name]], use.opt="pairwise.complete.obs",
                 corr.desc = paste0("06.Dam.Normalized.Mean.", name), createPDF=T)
# Scatter Plots on Averaged data
################################
ScatCor(data=DATAs.norm.ave[[name]], pref=paste0("07.Dam.Normalized.Mean.", name))
}


}
print("Congratulations!!!")



