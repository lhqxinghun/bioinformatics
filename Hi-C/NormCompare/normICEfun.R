#function name:normICEfun
#function:read data from disks and use the method ICE to normalize data
#parameters:
#             dirin:the path of input data
#             dirout:the path of output result
#             dirtemp:the path of temporary data
#             species:hg19 or mm9 
#             chr:the numberID of the chromosome
#             resolution:the resolution of input data
normICEfun <- function(dirin, dirout, dirtemp, species, chr, resolution) 
{
  require(HiTC)
  #read data from disks and transform data to the data format of ICE
  hicobjlist <- gethicobjlist(dirin, species, chr, resolution, normethod="ICE")
  #function:use the method of ICE to normalize data
  ICEformat <- function(i)
  {
    write.table(hicobjlist[[i]], paste0(dirtemp,"/",names(hicobjlist)[i],"-ICE",".txt"), quote=FALSE, sep = "\t")
    datalist <- sapply(list.files(dirtemp, pattern=paste0(names(hicobjlist)[i],"-ICE",".txt"), full.names=TRUE), import.my5C)
    hiC <- HTClist(datalist)
    hiC <- hiC[isIntraChrom(hiC)]
    hicnorm <- HTClist(lapply(hiC, normICE, max_iter=50))
    write.table(round(as.matrix(hicnorm[[1]]@intdata), digits=6), paste0(dirout, "/", names(hicobjlist)[i],"-ICE.txt"), row.names = FALSE, col.names = FALSE, sep = "\t")
  }
  lapply(seq_along(hicobjlist), ICEformat)
}
