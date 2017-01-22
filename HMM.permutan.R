library(ggplot2)

pcentFun <- function(x, y) {
  100 * (x / y)
}

FeatureCalls.to.GFF.like <- function(start.coordinate,
                                     end.coordinate,
                                     feature.type) {
  e.ind <- cumsum(rle(feature.type)$lengths);
  s.ind <- c(1, e.ind[1:(length(e.ind)-1)]+1);
  gff <- data.frame(start=start.coordinate[s.ind],
                    end=end.coordinate[e.ind],
                    value=feature.type[s.ind]);
  invisible(gff);
}


statSet <- "NRN.HP1.mf_min"
origSet <- DOMAIN.data.filt$NRN.HP1.mf_min.domains #протяженные домены
kindOfbins <- "300 бины"
attempts <- 10000

statDF <- rbind.fill(HMM.data[[statSet]]) # выход с нмм 0 и 1
statDF <- DATA.rb[[1]]
HMM.stat <- list()
AllGenome <- 119029689 # размер генома
AllGenomewHet <- 129663327
if (exists("DOMAIN.statistics.df") == T) rm(DOMAIN.statistics.df) # ?

for (i in c(1:attempts)) {
	relyDF <- data.frame("chr" = as.factor(statDF$chr), "start" = statDF$start, "end" = statDF$end, "target" = sample(statDF$domain, nrow(statDF)))
	HMM.stat[[statSet]] <- split(relyDF, relyDF$chr) # разделить по хромосам
	names(HMM.stat[[statSet]]) <- paste0("chr_", names(HMM.stat[[statSet]]))

	if (exists("d.stat.full") == T) rm(d.stat.full)
	for (selectChr in names(HMM.stat[[statSet]])) {
		d.stat <- FeatureCalls.to.GFF.like(start.coordinate=HMM.stat[[statSet]][[selectChr]]$start, end.coordinate=HMM.stat[[statSet]][[selectChr]]$end, feature.type=HMM.stat[[statSet]][[selectChr]]$target)
		d.stat <- cbind(chr=HMM.stat[[statSet]][[selectChr]]$chr[1], d.stat, stringsAsFactors=F)
		d.stat <- d.stat[d.stat$value==1,]
		if (exists("d.stat.full") == F) {d.stat.full <- d.stat;} else {d.stat.full <- rbind(d.stat.full, d.stat);}
	}

	DomainCount <- nrow(d.stat.full)
	DomainLength <- d.stat.full$end - d.stat.full$start
	DomainSize <- sum(DomainLength)
	if (exists("DOMAIN.statistics.df") == F) {
		DOMAIN.statistics.df <- data.frame("Count" = DomainCount, MeanLength = mean(DomainLength), MedianLength = median(DomainLength), GenomeCoverage = round(pcentFun(DomainSize, AllGenome), digits=3))
	}	else {
		DOMAIN.statistics.df <- rbind(DOMAIN.statistics.df, data.frame("Count" = DomainCount, MeanLength = mean(DomainLength), MedianLength = median(DomainLength), GenomeCoverage = round(pcentFun(DomainSize, AllGenome), digits=3)))
	}
}

# xoxoxo
dat <- DATA.rb[[1]]
dat <- FeatureCalls.to.GFF.like(dat$start, dat$end, dat$domain)
dat <- dat[dat$value == 1,]
failsafe <- DOMAIN.statistics.df
# xoxoxo
DOMAIN.statistics.df <- rbind(data.frame("Count" = nrow(origSet), MeanLength = mean(origSet$end - origSet$start), MedianLength = median(origSet$end - origSet$start), GenomeCoverage = round(pcentFun(sum(origSet$end - origSet$start), AllGenomewHet), digits=3)), DOMAIN.statistics.df)
rownames(DOMAIN.statistics.df)[1] <- "original"
# xoxoxo
DOMAIN.statistics.df <- rbind(data.frame("Count" = nrow(dat), MeanLength = mean(dat$end - dat$start),
                                         MedianLength = median(dat$end - dat$start), GenomeCoverage = round(pcentFun(sum(dat$end - dat$start), AllGenome), digits=3)), DOMAIN.statistics.df)
rownames(DOMAIN.statistics.df)[1] <- "original"
DOMAIN.statistics.df$col <- "permutated"
DOMAIN.statistics.df$col[1] <- "original"

DOMAIN.statistics.df$col <- as.factor(DOMAIN.statistics.df$col)


tempstat <- ggplot(
  DOMAIN.statistics.df, aes(x=GenomeCoverage, fill = col)) + 
  geom_bar(width = 0.1)
pdf(file=file.path(prefixDir, paste0(statSet, "_reliability_on_", kindOfbins, "_", attempts, "_repeats", ".pdf")), width=14, height=7)
print(tempstat)
dev.off()