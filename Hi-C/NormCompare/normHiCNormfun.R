normHiCNormfun <- function(dirin, dirout, dirtemp, species, chr, resolution, enzymes) 
{
  #require(BSgenome.Hsapiens.UCSC.hg19)
  require(rtracklayer)
  require(HiTC)
  
  hicobjlist <- gethicobjlist(dirin, species, chr, resolution, normethod="HiCNorm")
  HiCNormformat <- function(i)
  {
    # filtedobj <- zerofilter(as.matrix(hicobjlist[[i]]))
    # #cat(paste0("#GMBothall.0.maq.",basename(dir),".hm.newtracks12forBryan.heatmap.matrix.tab\n\t"),file=paste0(dirname(dir),"/temp/",names(hicobjlist)[i]))
    # #write.table(hicobjlist[[i]], paste0(dirname(dir),"/temp/",names(hicobjlist)[i]), append=TRUE, quote=FALSE, sep = "\t")
    # write.table(filtedobj$matrix, paste0(dirname(dir),"/temp/",names(hicobjlist)[i],"-HiCNorm",".txt"), quote=FALSE, sep = "\t")
    # datalist <- sapply(list.files(paste0(dirname(dir),"/temp"), pattern=paste0(names(hicobjlist)[i],"-HiCNorm",".txt"), full.names=TRUE), import.my5C)
    # hiC <- HTClist(datalist)
    # hiC <- hiC[isIntraChrom(hiC)]
    # map_hg19<- import(system.file("extdata", "wgEncodeCrgMapabilityAlign100mer.bigWig", package = "HiTC"), format="BigWig")
    # cutSites <- getAnnotatedRestrictionSites(resSite="AAGCTT", overhangs5=1, chromosomes=seqlevels(hiC), genomePack="BSgenome.Hsapiens.UCSC.hg19", wingc=200, mappability=map_hg19, winmap=500)
    # hiC_annot <- HTClist(lapply(hiC, setGenomicFeatures, cutSites))
    # hicnorm <- HTClist(lapply(hiC_annot, normLGF))
    # recoveredhicnorm <- zerorecover(as.matrix(hicnorm[[1]]@intdata), filtedobj)
    # write.table(round(recoveredhicnorm, digits=6), paste0(dirname(dir),"/result/normalized/",names(hicobjlist)[i],"-HiCNorm.txt"), row.names = FALSE, col.names = FALSE, sep = "\t")

    #cat(paste0("#GMBothall.0.maq.",basename(dir),".hm.newtracks12forBryan.heatmap.matrix.tab\n\t"),file=paste0(dirname(dir),"/temp/",names(hicobjlist)[i]))
    #write.table(hicobjlist[[i]], paste0(dirname(dir),"/temp/",names(hicobjlist)[i]), append=TRUE, quote=FALSE, sep = "\t")
    write.table(hicobjlist[[i]], paste0(dirtemp,"/",names(hicobjlist)[i],"-HiCNorm",".txt"), quote=FALSE, sep = "\t")
    datalist <- sapply(list.files(dirtemp, pattern=paste0(names(hicobjlist)[i],"-HiCNorm",".txt"), full.names=TRUE), import.my5C)
    hiC <- HTClist(datalist)
    hiC <- hiC[isIntraChrom(hiC)]
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
  #unlink(paste0(dirname(dir),"/temp"), recursive=TRUE)
  #dir.create(paste0(dirname(dir),"/temp"), recursive=TRUE)
  lapply(seq_along(hicobjlist), HiCNormformat)
}