# Preparation of Inputs for BioHMM
library("dplyr")
library("snapCGH")
library("GenomicRanges")
# Empty workspace
rm(list=ls(all=T))
# VARIABLES
heterochrom = T
lastcol = 4 # last column of file with genomic intervals (e.g. gatcs.txt, bins.txt) that 
            # were used in Dam_count script

# Set directory with file from Dam_count_statistics as work directory
# Use file obtained from Dam_count_statistics script
# /CSV/06.Dam.Normalized.DATA.csv
setwd("~/IMG/DamID/Neurons_hp1_lam_pc/19.12.16_v2.0/")

# Load data
DATA <- as_data_frame(read.delim("CSV/06.Dam.Normalized.DATA.csv", header=T, as.is=T, dec=".", sep = ';'))

# Take only the chromosomes "2L", "2R", "3L", "3R" and "X" for the subsequent analysis
# If heterochrom is True take also heterochromatic data
# Prepare list with lists for each experiment with dataframes broke down into chromosomes

chroms = if (heterochrom){ c("2L", "2LHet", "2R", "2RHet", "3L", "3LHet", "3R", "3RHet", "X", "XHet", "4", "YHet")
}else{c("2L", "2R", "3L", "3R", "X")}

DATA.chrs <- lapply(DATA[(lastcol + 1):ncol(DATA)], function(x){
  names(x) = lapply(chroms, function(y) filter(cbind(DATA[1:lastcol], "DamID.value" = x), chr == y, !is.na(DamID.value)))
})


# run BioHMM

dir.create("BioHMM")
setwd("BioHMM")

# Functions
runBioHMM <- function (mval, datainfo, useCloneDists = TRUE, criteria = "AIC", 
                       delta = NA, var.fixed = FALSE, epsilon = 1e-06, numiter = 30000) 
{
  crit = TRUE
  if (criteria == "AIC") {
    aic = TRUE
    bic = FALSE
  } else if (criteria == "BIC") {
    bic = TRUE
    aic = FALSE
  } else crit = FALSE
  if ((crit == 1) || (crit == 2)) {
    if (criteria == "BIC") {
      if (is.na(delta)) {
        delta <- c(1)
      }
    }
    res <- try(fit.model(
      obs = mval, datainfo = datainfo, useCloneDists = useCloneDists, 
      aic = aic, bic = bic, 
      delta = delta, var.fixed = var.fixed, epsilon = epsilon, 
      numiter = numiter
    ))$out.list$state
  }
  else {
    cat("You must enter AIC or BIC for the criteria argument\n")
  }
}


fit.model <- function (obs, datainfo = NULL, useCloneDists = TRUE, 
                       aic = TRUE, bic = FALSE, delta = 1, var.fixed = FALSE, 
                       epsilon = 1e-06, numiter = 30000) 
{
  library(cluster)
  kb <- datainfo$start
  if (useCloneDists) {
    dists.pre = kb[2:length(kb)] - kb[1:(length(kb) - 1)]
    dists = dists.pre/(max(dists.pre))
  } else {
    dists <- rep(1, length(kb))
  }
  covars <- as.matrix(dists)
  obs.ord <- obs[order(kb)]
  kb.ord <- kb[order(kb)]
  ind.nonna <- which(!is.na(obs.ord))
  data <- obs.ord[ind.nonna]
  kb <- kb.ord[ind.nonna]
  numobs <- length(data)
  if (numobs > 5) {
    temp2 <- clara(data, 2)
    init.mean.two <- temp2$medoids
    init.var.two <- vector()
    if (var.fixed == FALSE) {
      for (i in 1:2) {
        if (length(temp2$data[temp2$clustering == i]) > 
            1) 
          init.var.two[i] <- log(sqrt(var(temp2$data[temp2$clustering == 
                                                       i])))
        else init.var.two[i] <- log(0.5)
      }
    } else {
      init.var.two[1:2] <- log(sqrt(var(data)))
    }
    z2.init <- c(init.mean.two[, 1], init.var.two, -1, -3.6, 
                 -3.6, 0)
    z.pre <- run.nelder(numobs, z2.init, data, 
                        covars, var.fixed, epsilon, numiter, i)
    if (!is.nan(z.pre$x[1])) {
      z2 <- find.param.two(z.pre, var.fixed)
    } else {
      z2 <- NULL
    }
    if (aic) {
      factor <- 2
    } else if (bic) {
      factor <- log(numobs) * delta
    } else {
      stop("No criteria selected")
    }
    z <- z2
    nstates <- 2
    trans.mat <- list()
    for (j in 1:(length(data) - 1)) {
      trans.mat[[j]] <- z$LH.trans + exp(-(covars[j, 
                                                  1]^(z$rate1)) * prod(covars[j, -1])) * z$RH.trans
    }
    Vit.seg <- Viterbi.two(data, 
                           z, trans.mat)
    maxstate.unique <- unique(Vit.seg)
    mean <- rep(0, length(data))
    var <- rep(0, length(data))
    for (m in 1:length(maxstate.unique)) {
      mean[Vit.seg == maxstate.unique[m]] <- mean(data[Vit.seg == 
                                                         maxstate.unique[m]])
      var[Vit.seg == maxstate.unique[m]] <- var(data[Vit.seg == 
                                                       maxstate.unique[m]])
    }
    out <- cbind(matrix(Vit.seg, ncol = 1), matrix(mean, 
                                                   ncol = 1), matrix(var, ncol = 1))
    out.all <- matrix(NA, nrow = length(kb.ord), ncol = 4)
    out.all[ind.nonna, 1:3] <- out
    out.all[, 4] <- obs.ord
    out.all <- as.data.frame(out.all)
    dimnames(out.all)[[2]] <- c("state", "mean", "var", "obs")
    numstates <- length(unique(Vit.seg))
  } else {
    out.all <- matrix(NA, nrow = length(kb.ord), ncol = 4)
    out.all[ind.nonna, 1] <- c(rep(1, numobs))
    out.all[ind.nonna, 2] <- c(rep(mean(obs.ord), numobs))
    out.all[ind.nonna, 3] <- c(rep(var(obs.ord), numobs))
    out.all[ind.nonna, 4] <- obs.ord
    out.all <- as.data.frame(out.all)
    dimnames(out.all)[[2]] <- c("state", "mean", "var", "obs")
    numstates = 1
  }
  list(out.list = out.all, nstates.list = numstates)
}

# Makes intersection of replicas for more stringent domain identification
IntersectRep <- function(dataBed){
  # USE ONLY IN CASE OF TWO REPLICAS!!
  proteins <- unique(sub("([^.]*)\\.([^.]*)\\..*", "\\1", names(dataBed)))
  intsc_list <- lapply(proteins, function(prot){
    inds <- grep(prot, names(dataBed))
    bedranges1 <- dataBed[[inds[1]]]
    # this use of GenomicRanges is to compress resulting bed files and unite extended domains
    gr1 <- GRanges(seqnames = Rle(bedranges1$chr),  ranges = IRanges(start = bedranges1$start, end = bedranges1$end))
    
    bedranges2 <- dataBed[[inds[2]]]
    gr2 <- GRanges(seqnames = Rle(bedranges2$chr),  ranges = IRanges(start = bedranges2$start, end = bedranges2$end))
    intersect.gr <- intersect(gr1, gr2)
    df_intersected <- data.frame(chr = seqnames(intersect.gr), start = start(intersect.gr), end = end(intersect.gr))
    # Remove gaps of size of one bin
    df_no_gap <- df_intersected[1,]
    for (rowN in 2:(nrow(df_intersected))){
      if (df_intersected[rowN, 2] - df_intersected[(rowN - 1), 3] == 300){
        df_no_gap[nrow(df_no_gap), 3] <- df_intersected[rowN, 3]
      } else{ df_no_gap <- rbind(df_no_gap, df_intersected[rowN,])}
    }
    df_no_gap$chr <- paste0("chr", df_no_gap$chr)
    return(df_no_gap)
  })
  names_intsc <- sapply(proteins, function(prot) {
    inds <- grep(prot, names(dataBed))
    sub("(.*)\\.\\d.*", "\\1.intersected", names(dataBed)[inds[1]])
  })
  names(intsc_list) <- names_intsc
  return(intsc_list)
  
}

DATA.HMM <- lapply(DATA.chrs, function(exp_list){
  lapply(exp_list, function(chromdf){
    BioHMM.output <- runBioHMM(mval = chromdf$DamID.value, datainfo = chromdf, useCloneDists = T)
    classifier <- aggregate(x=chromdf$DamID.value, by=list(BioHMM.output), FUN=mean)
    print(classifier)
    zero.prox <- classifier$x > 0
    if (length(zero.prox) == 2 & xor(zero.prox[1], zero.prox[2])){
      signals <- BioHMM.output == classifier[,1][zero.prox]
      nonsignals <- BioHMM.output == classifier[,1][!zero.prox]
      BioHMM.output[signals] <- 1
      BioHMM.output[nonsignals] <- 0
      chromdf$domain <- BioHMM.output
    } else {
      print("There's a problem with classification")
      return()
    }
    return(chromdf)
  })
})

# Remove Null elements from lists, though I think that GenomicRanges::intersect deals with it
DATA.HMM2 <- lapply(DATA.HMM, function(lll) Filter(Negate(function(x) is.null(unlist(x))), lll))

DATA.bed <- lapply(DATA.HMM2, function(x){
  bedranges <- do.call("rbind", x) %>% 
    # leave only domains
  filter(domain == 1) %>% 
    # add 'chr' to chromosomes names
  mutate(chrom = paste0("chr", chr))
  # this use of GenomicRanges is to compress resulting bed files and unite extended domains
  gr <- GRanges(seqnames = Rle(bedranges$chr),  ranges = IRanges(start = bedranges$start, end = bedranges$end)) %>% 
    reduce() 
  # return dataframe which is totally ready to be written in.bed
  data.frame(chr = seqnames(gr), start = start(gr), end = end(gr))

})

# Save as bed files
lapply(seq_along(DATA.bed), function(x){
  tissue = sub("([^.]*)\\.([^.]*)\\..*", "\\1", names(DATA.bed)[x])
  protein = sub("([^.]*)\\.([^.]*)\\..*", "\\2", names(DATA.bed)[x])
  number = sub(".*\\.(\\d).norm", "\\1", names(DATA.bed)[x])
  bed.name <- paste0(protein, '.', tissue, '.', number, '.', 'domains.bed') 
  track_name = paste0('track name="', protein, '.', number, ' ', tissue, ' HMM"')
  desc = paste0('description="', protein, ' HMM domains for ', tissue, '"')
  write.table(paste(track_name, desc), file=bed.name, sep="\t", row.names=F, col.names=F, quote=F, dec=".", append=F)
  write.table(DATA.bed[[x]], file=bed.name, sep="\t", row.names=F, col.names=F, quote=F, dec=".", append=T)
})

# Intersect domains between replicas
DATA.bed.intrsct <- IntersectRep(DATA.bed)

# Save intersected domains 
lapply(seq_along(DATA.bed.intrsct), function(x){
  tissue = sub("([^.]*)\\.([^.]*)\\..*", "\\1", names(DATA.bed.intrsct)[x])
  protein = sub("([^.]*)\\.([^.]*)\\..*", "\\2", names(DATA.bed.intrsct)[x])
  bed.name <- paste0(protein, '.', tissue, '.', 'domains.intersected.bed') 
  track_name = paste0('track name="', protein, ' ', tissue, ' HMM"')
  desc = paste0('description="', protein, ' HMM domains intersected between 2 replicas for ', tissue, '"')
  write.table(paste(track_name, desc), file=bed.name, sep="\t", row.names=F, col.names=F, quote=F, dec=".", append=F)
  write.table(DATA.bed.intrsct[[x]], file=bed.name, sep="\t", row.names=F, col.names=F, quote=F, dec=".", append=T)
})



