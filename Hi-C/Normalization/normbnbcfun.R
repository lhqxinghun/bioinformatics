normbnbcfun <- function(dirin1, dirin2, dirin3, dirout1, dirout2, dirout3, species, chr, resolution, batch) 
{
  require(bnbc)

  hicobjlist1 <- gethicobjlist(dirin1, species, chr, resolution, normethod="bnbc")
  hicobjlist2 <- gethicobjlist(dirin2, species, chr, resolution, normethod="bnbc")
  hicobjlist3 <- gethicobjlist(dirin3, species, chr, resolution, normethod="bnbc")
  groups=c(rep(1,length(hicobjlist1$matrixlist)),rep(2,length(hicobjlist2$matrixlist)),rep(3,length(hicobjlist3$matrixlist)))
  names <-c(names(hicobjlist1$matrixlist), names(hicobjlist2$matrixlist), names(hicobjlist3$matrixlist))
  rowdata <- makeGRangesFromDataFrame(hicobjlist1$indexdf)
  hicobjlist1 <- c(hicobjlist1$matrixlist,hicobjlist2$matrixlist,hicobjlist3$matrixlist)
  rm(hicobjlist2, hicobjlist3)
  coldata <- DataFrame(Batch = batch, row.names = names(hicobjlist1))
  hicobjlist1 <- ContactGroup(rowdata, contacts=hicobjlist1, coldata)
  # if(resolution == 5000)
  # {
  #   if((nrow(hicobjlist1@contacts[[1]])-1) >= 2)
  #     nbands <- 2
  #   else
  #     nbands <- nrow(hicobjlist1@contacts[[1]])-1
  # }
  # else if(resolution == 10000)
  # {
  #   if((nrow(hicobjlist1@contacts[[1]])-1) >= 2)
  #     nbands <- 2
  #   else
  #     nbands <- nrow(hicobjlist1@contacts[[1]])-1    
  # }
  # else if(resolution == 25000)
  # {
  #   if((nrow(hicobjlist1@contacts[[1]])-1) >= 4)
  #     nbands <- 4
  #   else
  #     nbands <- nrow(hicobjlist1@contacts[[1]])-1    
  # }    
  # else if(resolution == 50000)
  # {
  #   if((nrow(hicobjlist1@contacts[[1]])-1) >= 4)
  #     nbands <- 4
  #   else
  #     nbands <- nrow(hicobjlist1@contacts[[1]])-1    
  # }
  # else if(resolution == 100000)
  # {
  #   if((nrow(hicobjlist1@contacts[[1]])-1) >= 6)
  #     nbands <- 6
  #   else
  #     nbands <- nrow(hicobjlist1@contacts[[1]])-1    
  # }
  # else if(resolution == 250000)
  # {
  #   if((nrow(hicobjlist1@contacts[[1]])-1) >= 6)
  #     nbands <- 6
  #   else
  #     nbands <- nrow(hicobjlist1@contacts[[1]])-1    
  # }
  # else if(resolution == 500000)
  # {
  #   if((nrow(hicobjlist1@contacts[[1]])-1) >= 13)
  #     nbands <- 13
  #   else
  #     nbands <- nrow(hicobjlist1@contacts[[1]])-1    
  # }
  # else if(resolution == 1000000)
  # {
  #   if((nrow(hicobjlist1@contacts[[1]])-1) >= 25)
  #     nbands <- 25
  #   else
  #     nbands <- nrow(hicobjlist1@contacts[[1]])-1    
  # }
  # else if(resolution == 2500000)
  # {
  #   if((nrow(hicobjlist1@contacts[[1]])-1) >= 63)
  #     nbands <- 63
  #   else
  #     nbands <- nrow(hicobjlist1@contacts[[1]])-1    
  # }
  hicobjlist1.cpm <- logCPM(hicobjlist1)
  hicobjlist1.smooth <- boxSmoother(hicobjlist1.cpm, h=5, mc.cores=1)
  hicobjlist1 <- bnbc(hicobjlist1.smooth, colData(hicobjlist1)$Batch, threshold=2500000, step=resolution, bstart=2, nbands=nrow(hicobjlist1@contacts[[1]])-1) #nbands=nrow(hicobjlist1@contacts[[1]])-1 / nbands=nbands
  # hicobjlist1 <- bnbc(hicobjlist1, colData(hicobjlist1)$Batch, threshold=2500000, step=resolution, bstart=2, nbands=(nrow(hicobjlist1@contacts[[1]])-1)) #nbands=(nrow(hicobjlist1@contacts[[1]])-1)
  bnbcformat <- function(i)
  {
    if(groups[i] == 1)
      dirout <- dirout1
    else if(groups[i] == 2)
      dirout <- dirout2
    else if(groups[i] == 3)
      dirout <- dirout3
    hicmatrix <- as.matrix(hicobjlist1@contacts[[i]])
    hicmatrix <- exp(hicmatrix)
    write.table(round(hicmatrix, digits=6), paste0(dirout, "/", names[i],"-bnbc.txt"), row.names = FALSE, col.names = FALSE, sep = "\t")
  }
  lapply(seq_along(hicobjlist1@contacts), bnbcformat)
}