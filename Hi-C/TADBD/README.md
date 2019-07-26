TADBD: a fast and sensitive tool for detection of TAD boundaries

Abstract<br>
 　  Hi-C technology allows for genome-wide profiling of chromatin interactions in space. Topological Associated Domain (TAD) is a self-interacting genomic block, which is conserved across species in mammalian genomes and in association with regulation of biological functions. The detection of TAD boundaries on Hi-C contact matrix is one of the most important issues in the analysis of 3D genome architecture at TAD level. Here, we present TADBD, a fast and sensitive computational tool for detection of TAD boundaries on Hi-C contact matrix. It implements a novel Haar-based algorithm considering Haar diagonal template, the acceleration via a compact integrogram, multi-scale aggregation at template size, and statistical filtering. The comparison results show that TADBD is superior in speed, and competitive in terms of accuracy, reproducibility and user-friendliness.<br>

1> Datasets<br>
  Simulated data:<br>
 　   Yu, W., He, B. and Tan, K. (2017) Identifying topologically associating domains and ubdomains by Gaussian Mixture model And 
      Proportion test, Nat Commun, 8, 535.https://bitbucket.org/mforcato/hictoolscompare/src/master/<br> 
  Real data:<br>
 　   Rao, S.S.P., et al. (2014) A 3D Map of the Human Genome at Kilobase Resolution Reveals Principles of Chromatin Looping, Cell, 159,         1665-1680.<br> 
		
2> TADBD<br>		
  Description:<br>
 　   An R package for detection of TAD boundaries<br> 
  Depends:<br>
 　   R (>= 3.5.1)<br>
  Example:
 
    rm(list=ls())
	library(TADBD)
	species <- "hg19"
	chr <- "chr18"
	resolution <- 50000
	options(scipen = 999)
	data(hicdata)
	hicmat <- DataLoad(hicdata, bsparse = F, species, chr, resolution)
	df_result <- TADBD(hicmat)
	Output(df_result, species, chr, resolution)
         
