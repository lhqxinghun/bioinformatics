１. Source code for detecting TADs boundaries methods. The code has been tested under Ubuntu14.04 LTS.

1> HMDTB.R<br>
 Function:<br>
 　　　　Detecting TADs boundaries method<br> 
 Depends:<br>
 　　　　./HMDTB.R (provided)<br>
 　　　　R (>= 3.5.1)<br>
 　　　　quantmod (R package, 0.4-14)<br>
 Example:
		 
		 rm(list=ls())
		 setwd(path_folder)
		 source("./HMDTB.R")
		 species <- "hg19"
		 chr <- "chr14"
		 resolution <- 50000
		 hicdata_path = "./input/002_chr14_50Kb-ICE.txt"
		 hicdata<-read.table(hicdata_path, head=FALSE)
		 hicmat<-as.matrix(hicdata)
		 dir.create("output", showWarnings = FALSE, recursive = FALSE, mode = "0777")
		 result <- HMDTB(hicmat, window.size = NULL, statFilter = T)