#MOLOC - Vehicle
library(moloc)
library(foreach)
library(doParallel)

#Import molQTLs and moloc bed file
	setwd('/scratch/msoliai/moloc/qtl')
	eqtl <- read.table('rv.eqtl.moloc.txt', header=TRUE)
	meqtl <- read.table('rv.meqtl.moloc.txt', header=TRUE)
	bed.data <- read.table('rv.bed.moloc.txt', header=TRUE)

#Import TAGC asthma GWAS
	setwd('/scratch/msoliai/moloc/gwas')
	gwas <- read.table('tagc.gwas.moloc.txt', header=TRUE)

#Identify genomic ranges for bed file (Independent signals from GARFIELD, +/- 1MB)
	setwd('/scratch/msoliai/moloc/garfield.moloc')
	ranges.gar <- read.table('asthma.tagc.range.moloc.txt', h=T)
	ranges.gar$down <- as.numeric(ranges.gar$pos) - 1e6
	ranges.gar$up <- as.numeric(ranges.gar$pos) + 1e6

#Determine which eQTL region from bed file falls within the molQTL enriched GWAS regions (GARFIELD)
	y  <- NULL

	for (i in 1:nrow(ranges.gar)){
	
		tmp <- bed.data[bed.data$CHR == as.numeric(ranges.gar[i,1]) & bed.data$START >= as.numeric(ranges.gar[i,6]) & bed.data$START <= as.numeric(ranges.gar[i,7]),]
  		y <- rbind(y, tmp)
  		
 	}

	bed.data <- y
	bed.data <- split( bed.data , f = bed.data$CHR )
	rm(y)

##############################################
##moloc - Genome-wide/multiple loci analysis##
##############################################
options(scipen = 1, digits = 2)

  gwas.moloc <- gwas[gwas$CHR == as.numeric(bed.data[[11]][1,3]),]
  eqtl.moloc <- eqtl[eqtl$CHR == as.numeric(bed.data[[11]][1,3]),]
  meqtl.moloc <- meqtl[meqtl$CHR == as.numeric(bed.data[[11]][1,3]),]
  data.all <- list(gwas.moloc, eqtl.moloc, meqtl.moloc)
  tmp <- coloc.genome(data.all, bed.data[[11]], prefix = "pref", save.SNP.info=FALSE, cores=2, have_alleles=FALSE, bychrpos=TRUE, prior_var="default", priors=c(1e-04, 1e-06, 1e-07), min_nsnps = 10, takelog = TRUE, write=FALSE, outfolder = "test", forcePVAL = FALSE)

	setwd('/scratch/msoliai/moloc/moloc.output/tagc')
	write.table(tmp, 'rv.moloc.asthma.tagc.11.txt', sep='\t', quote=FALSE, row.names=FALSE)
