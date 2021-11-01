regevaluate <- function(filein.y, filein.p, repnum, selycon, selyrep, selpcon, conum, selpiter, fileout)
{
  library(caret)
  library(lattice)
  library(ggplot2)
  require(data.table)
  
  yvalue <- as.matrix(fread(filein.y, header=F))
  pvalue <- as.matrix(fread(filein.p, header=F))
  selyrow <- (selycon-1)*repnum+selyrep
  selprow <- (selpiter-1)*conum+selpcon
  if(selyrep != 0)
  {
    #Data for four evaluation plots
    originalval <- c(yvalue[selyrow,])
    fittedval <-  c(pvalue[selprow,])
    residuals <- originalval-fittedval
    stdresiduals <- residuals / sd(residuals)
  }
  else
  {
    #Data for four evaluation plots
    originalval <- as.vector(t(yvalue))
    fittedval <-  rep(pvalue[selprow,], repnum)
    residuals <- originalval-fittedval
    stdresiduals <- residuals / sd(residuals)
  }
    #Residuals VS Fitted
    fileout.resvsfit <- paste0(dirname(fileout), "/", substr(basename(fileout),1,nchar(basename(fileout))-5),"-rep",selyrep,"-iter",selpiter,"-resvsfit.TIFF")
    tiff(fileout.resvsfit, width = 4.15, height = 3.28, units = "cm", pointsize=1, res=350)
    par(mar=c(9, 11, 0.2, 0), mgp=c(6,3,0), lwd=0.2,lty =2, cex.axis=5, cex.lab=5)
    plot(fittedval, residuals, pch=20, lty = 2, lwd=0.6, axes = FALSE,ann = FALSE, ylim = c(-1.2,1.2))
    lfmodel <- loess(residuals~fittedval)
    idx <- order(fittedval)
    lines(fittedval[idx],lfmodel$fitted[idx],col="red", lwd = 0.8, lty=2, cex=5)
    abline(h=0, col="blue", lwd = 0.8, lty=2, cex=5)
    axis(side = 1, lwd = 0.2)
    axis(side = 2, lwd = 0.2)
    mtext("Fitted values", side = 1, line=7, cex=5)
    mtext("Residuals", side = 2, line=7, cex=5)
    dev.off()
    #Scale-Location
    fileout.scavsfit <- paste0(dirname(fileout), "/", substr(basename(fileout),1,nchar(basename(fileout))-5),"-rep",selyrep,"-iter",selpiter, "-scavsfit.TIFF")
    tiff(fileout.scavsfit, width = 4.15, height = 3.28, units = "cm", pointsize=1, res=350)
    par(mar=c(9, 11, 0.2, 0), mgp=c(6,3,0), lwd=0.2,lty =2, cex.axis=5, cex.lab=5)
    plot(fittedval, sqrt(abs(stdresiduals)), pch=20, lty = 2, lwd=0.6, axes = FALSE,ann = FALSE, ylim = c(0,3.6))
    lfmodel <- loess(sqrt(abs(stdresiduals))~fittedval)
    idx <- order(fittedval)
    lines(fittedval[idx],lfmodel$fitted[idx],col="red", lwd = 0.8, lty=2, cex=5)
    axis(side = 1, lwd = 0.2)
    axis(side = 2, lwd = 0.2)
    mtext("Fitted values", side = 1, line=7, cex=5)
    mtext("sqrt(|Std residuals|)", side = 2, line=7, cex=5)
    dev.off()
    #Normal Q-Q
    fileout.norQQ <- paste0(dirname(fileout), "/", substr(basename(fileout),1,nchar(basename(fileout))-5),"-rep",selyrep,"-iter",selpiter, "-norQQ.TIFF")
    tiff(fileout.norQQ, width = 4.15, height = 3.28, units = "cm", pointsize=1, res=350)
    par(mar=c(9, 11, 0.2, 0), mgp=c(6,3,0), lwd=0.2,lty =2, cex.axis=5, cex.lab=5)
    qqnorm(stdresiduals, pch=20, lty = 2, lwd=0.6, axes = FALSE,ann = FALSE)
    abline(0, 1, col="blue", lwd = 0.8, lty=2, cex=5)
    axis(side = 1, lwd = 0.2)
    axis(side = 2, lwd = 0.2)
    mtext("Theoretical quantiles(norm)", side = 1, line=7, cex=5)
    mtext("Std residuals", side = 2, line=7, cex=5)
    legend("topleft", legend=c("Actual", "Idea"),col=c("red", "blue"), bty="n", lwd = 0.8, lty=2, cex=5)
    dev.off()
    #Fitted VS Original
    fileout.fitvsori <- paste0(dirname(fileout), "/", substr(basename(fileout),1,nchar(basename(fileout))-5),"-rep",selyrep,"-iter",selpiter, "-fitvsori.TIFF")
    tiff(fileout.fitvsori, width = 4.15, height = 3.28, units = "cm", pointsize=1, res=350)
    par(mar=c(9, 11, 0.2, 0), mgp=c(6,3,0), lwd=0.2,lty =2, cex.axis=5, cex.lab=5)
    plot(originalval, fittedval, pch=20, lty = 2, lwd=0.6, axes = FALSE,ann = FALSE)
    lfmodel <- loess(fittedval~originalval)
    idx <- order(originalval)
    lines(originalval[idx],lfmodel$fitted[idx],col="red", lwd = 1, lty=2, cex=5)
    abline(0, 1, col="blue", lwd = 0.8, lty=2, cex=5)
    axis(side = 1, lwd = 0.2)
    axis(side = 2, lwd = 0.2)
    mtext("Original values", side = 1, line=7, cex=5)
    mtext("Fitted values", side = 2, line=7, cex=5)
    dev.off()
  
}
