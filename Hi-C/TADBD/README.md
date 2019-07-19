TADBD: a fast and accurate tool for detection of TAD boundaries

Abstract<br>
Motivation:<br> 
        Hi-C technology allows for genome-wide profiling of chromatin interactions in space. Topological Associated Domain (TAD) is a self-interacting genomic block, which is c-
		onserved across species in mammalian genomes and in association with regulation of 
		biological functions. The detection of TAD boundaries on Hi-C contact matrix is one
		of the most important issues in the analysis of 3D genome architecture at TAD level.
Results:<br> 
        Here, we present TADBD, a fast and accurate computational tool for detection of TAD
		boundaries on Hi-C contact matrix. It implements a novel Haar-based algorithm cons-
		idering Haar diagonal template, the acceleration via a compact integrogram, multis-
		cale aggregation at tem-plate size, and statistical filtering. The comparison resu-
		lts show that TADBD is superior in speed, and competitive in terms of accuracy,reproducibility and user-friendliness.

1> Dataset<br>
 Simulated data:<br>
        Yu, W., He, B. and Tan, K. (2017) Identifying topologically associating domains and subdomains by Gaussian Mixture model And Proportion test, Nat Commun, 8, 535.<br> 
 　　　 https://bitbucket.org/mforcato/hictoolscompare/src/master/<br> 
 Real data:<br>
		Rao, S.S.P., et al. (2014) A 3D Map of the Human Genome at Kilobase Resolution Reveals Principles of Chromatin Looping, Cell, 159, 1665-1680.<br> 
		
2> TADBD.R<br>		
 Function:<br>
 　　　 Detection of TAD boundaries<br> 
 Depends:<br>
 　　　 R (>= 3.5.1)<br>
 　　　 TADBD(R package, )<br>
 Example:
 
        rm(list=ls())
		library(TADBD)
		species <- "hg19"
		chr <- "chr18"
		resolution <- 50000
		options(scipen = 999)
		hicdata_path = "./data/001_chr18_50Kb-ICE.mat"
		hicmat <- DataLoad(hicdata_path, bsparse = F, species, chr, resolution)
		result <- TADBD(hicmat)
		Output(chr, resolution, hicmat, result)
