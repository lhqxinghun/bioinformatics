R2vsiteration <- function(filein.y, filein.p, repnum, selycon, selyrep, selpcon, conum, piternum, fileout)
{
  library(caret)
  library(lattice)
  library(ggplot2)
  library(data.table)
  
  yvalue <- as.matrix(fread(filein.y, header=F))
  pvalue <- as.matrix(fread(filein.p, header=F))
  if(selyrep != 0)
  {
    selyrow <- (selycon-1)*repnum+selyrep
    RMSE <- c()
    R2 <- c()
    for(i in 1:piternum)
    {
      selprow <- (i-1)*conum+selpcon
      RMSER2 <- postResample(pvalue[selprow,],yvalue[selyrow,])
      #RMSER2 <- postResample(rep(pvalue[selprow,], repnum), as.vector(t(yvalue)))
      RMSE[i] <- RMSER2[1]
      R2[i] <- RMSER2[2]
    }
    #R2 VS iteration
    fileout.R2vsiter <- paste0(dirname(fileout), "/", substr(basename(fileout),1,nchar(basename(fileout))-5),"-rep",selyrep,".TIFF")
    tiff(fileout.R2vsiter, width = 4.15, height = 3.28, units = "cm", pointsize=1, res=350)
    par(mar=c(11, 11, 0, 0), mgp=c(8,2,0), lwd=0.2,lty =2, cex.axis=5, cex.lab=5)
    plot(1:piternum, R2, type="o", col="black", bg="black", lwd=0.4, pch=21, axes = FALSE,ann = FALSE, cex=5, lty = 1)
    arrows(2,R2[2],0,R2[2],col="green",length = 0,lty = 2,lwd=0.6)
    arrows(2,R2[2],2,min(R2),col="green",length = 0,lty = 2,lwd=0.6)
    points(2, R2[2], col="green", bg="green", pch = 21, cex=5)
    axis(side = 1, lwd = 0.2, mgp = c(8, 3, 0))
    axis(side = 2, lwd = 0.2)
    mtext("Iterative number", side = 1, line=7, cex=5)
    mtext("R2", side = 2, line=7, cex=5)
    dev.off()
  }
  else
  {
    for (selyrep in 1:repnum) {
      selyrow <- (selycon-1)*repnum+selyrep
      RMSE <- c()
      R2 <- c()
      for(i in 1:piternum)
      {
        selprow <- (i-1)*conum+selpcon
        RMSER2 <- postResample(pvalue[selprow,],yvalue[selyrow,])
        RMSE[i] <- RMSER2[1]
        R2[i] <- RMSER2[2]
      }
      #R2 VS iteration
      fileout.R2vsiter <- paste0(dirname(fileout), "/", substr(basename(fileout),1,nchar(basename(fileout))-5),"-rep",selyrep,".TIFF")
      fileout.R2vsiter_out <- paste0(dirname(fileout), "/", substr(basename(fileout),1,nchar(basename(fileout))-5),"-rep",selyrep,".txt")
      write.table(R2, file = fileout.R2vsiter_out, quote = FALSE, row.names = F, col.names = F)
      tiff(fileout.R2vsiter, width = 4.15, height = 3.28, units = "cm", pointsize=1, res=350)
      par(mar=c(11, 11, 0, 0), mgp=c(8,2,0), lwd=0.2,lty =2, cex.axis=5, cex.lab=5)
      plot(1:piternum, R2, type="o", col="black", bg="black", lwd=0.4, pch=21, axes = FALSE,ann = FALSE, cex=5, lty = 1)
      arrows(2,R2[2],0,R2[2],col="green",length = 0,lty = 2,lwd=0.6)
      arrows(2,R2[2],2,min(R2),col="green",length = 0,lty = 2,lwd=0.6)
      points(2, R2[2], col="green", bg="green", pch = 21, cex=5)
      axis(side = 1, lwd = 0.2, mgp = c(8, 3, 0))
      axis(side = 2, lwd = 0.2)
      mtext("Iterative number", side = 1, line=7, cex=5)
      mtext("R2", side = 2, line=7, cex=5)
      dev.off()
    }
  }
}
