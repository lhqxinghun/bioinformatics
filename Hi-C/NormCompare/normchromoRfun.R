#function name:normchromoRfun
#function:read data from disks and use the method chromoR to normalize data
#parameters:
#             dirin:the path of input data
#             dirout:the path of output result
#             species:hg19 or mm9 
#             chr:the numberID of the chromosome
#             resolution:the resolution of input data
normchromoRfun <- function(dirin, dirout, species, chr, resolution) 
{
  require(chromoR)
  #read data from disks and transform data to the data format of chromoR
  hicobjlist <- gethicobjlist(dirin, species, chr, resolution, normethod="chromoR")
  #function:use the method of chromoR to normalize data
  chromoRformat <- function(i)
  {
    hicnorm <- correctCIM(hicobjlist$matrixlist[[i]], hicobjlist$seg)
    write.table(round(as.matrix(hicnorm$mCorrected), digits=6), paste0(dirout, "/", names(hicobjlist$matrixlist)[i],"-chromoR.txt"), row.names = FALSE, col.names = FALSE, sep = "\t")
  }
  lapply(seq_along(hicobjlist$matrixlist), chromoRformat)
}
