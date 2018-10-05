normchromoRfun <- function(dirin, dirout, species, chr, resolution) 
{
  require(chromoR)
  
  hicobjlist <- gethicobjlist(dirin, species, chr, resolution, normethod="chromoR")
  chromoRformat <- function(i)
  {
    hicnorm <- correctCIM(hicobjlist$matrixlist[[i]], hicobjlist$seg)
    write.table(round(as.matrix(hicnorm$mCorrected), digits=6), paste0(dirout, "/", names(hicobjlist$matrixlist)[i],"-chromoR.txt"), row.names = FALSE, col.names = FALSE, sep = "\t")
  }
  lapply(seq_along(hicobjlist$matrixlist), chromoRformat)
}