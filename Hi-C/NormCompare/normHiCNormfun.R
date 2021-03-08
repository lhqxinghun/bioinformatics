#function name:normHiCNormfun
#function:read data from disks and use the method HiCNorm to normalize data
#parameters:
#             dirin:the path of input data
#             dirout:the path of output result
#             dirtemp:the path of temporary data
#             species:hg19 or mm9 
#             chr:the numberID of the chromosome
#             resolution:the resolution of input data
#             enzymes:according to the restriction enzyme, the sequence of restriction site can be set
normHiCNormfun <- function(dirin, dirout, dirtemp, species, chr, resolution, enzymes) 
{
  
  require(rtracklayer)
  require(HiTC)
  
  #read data from disks and transform data to the data format of HiCNorm
  hicobjlist <- gethicobjlist(dirin, species, chr, resolution, normethod="HiCNorm")
  #function:use the method of HiCNorm to normalize data
  HiCNormformat <- function(i)
  {
    
    write.table(hicobjlist[[i]], paste0(dirtemp,"/",names(hicobjlist)[i],"-HiCNorm",".txt"), quote=FALSE, sep = "\t")
    datalist <- sapply(list.files(dirtemp, pattern=paste0(names(hicobjlist)[i],"-HiCNorm",".txt"), full.names=TRUE), import.my5C)
    hiC <- HTClist(datalist)
    hiC <- hiC[isIntraChrom(hiC)]
    #according to the restriction enzyme, the sequence of restriction site can be set
    switch(enzymes[i],
      HindIII=
      {
        cutSites <- "AAGCTT"
      },
      MboI=
      {
        cutSites <- "GATC"
      },
      DpnII=
      {
        cutSites <- "GATC"
      }
    )
    #species:hg19 or mm9
    switch(species,
      hg19=
      {
        map <- import(system.file("extdata", "wgEncodeCrgMapabilityAlign100mer-hg19.bigWig", package = "HiTC"), format="BigWig")
        cutSites <- getAnnotatedRestrictionSites(resSite=cutSites, overhangs5=1, chromosomes=seqlevels(hiC), genomePack="BSgenome.Hsapiens.UCSC.hg19", wingc=200, mappability=map, winmap=500)
      },
      mm9=
      {
        map <- import(system.file("extdata", "wgEncodeCrgMapabilityAlign100mer-mm9.bigWig", package = "HiTC"), format="BigWig")
        cutSites <- getAnnotatedRestrictionSites(resSite=cutSites, overhangs5=1, chromosomes=seqlevels(hiC), genomePack="BSgenome.Mmusculus.UCSC.mm9", wingc=200, mappability=map, winmap=500)
      }
    )
    hiC_annot <- HTClist(lapply(hiC, setGenomicFeatures, cutSites))
    hicnorm <- HTClist(lapply(hiC_annot, normLGF))
    write.table(round(as.matrix(hicnorm[[1]]@intdata), digits=6), paste0(dirout, "/", names(hicobjlist)[i],"-HiCNorm.txt"), row.names = FALSE, col.names = FALSE, sep = "\t")
  }
  lapply(seq_along(hicobjlist), HiCNormformat)
}
