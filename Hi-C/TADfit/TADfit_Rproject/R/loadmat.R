loadmat <- function(dirin1, dirin2 = NULL, isFile=FALSE)
{
matlist  <- list()
if(!isFile)
{
  
  if(is.null(dirin1) == FALSE)
  {
    listcon1 <- sapply(list.files(dirin1, pattern=".txt", full.names=TRUE), read.table, simplify = FALSE, USE.NAMES = TRUE)
    listcon1 <- lapply(listcon1, as.matrix)
    namesplit <- strsplit(basename(names(listcon1)),split=".",fixed=TRUE)
    names(listcon1) <- unlist(lapply(namesplit,head,1)) 
    matlist[["listcon1"]] <- listcon1
  }
  if(is.null(dirin2) == FALSE)
  {
    listcon2 <- sapply(list.files(dirin2, pattern=".txt", full.names=TRUE), read.table, simplify = FALSE, USE.NAMES = TRUE)
    listcon2 <- lapply(listcon2, as.matrix)
    namesplit <- strsplit(basename(names(listcon2)),split=".",fixed=TRUE)
    names(listcon2) <- unlist(lapply(namesplit,head,1))
    matlist[["listcon2"]] <- listcon2
  }
  
}
else
{
  if(is.null(dirin1) == FALSE)
  {
    mat1 <- NULL
    mat1 <- read.table(dirin1, stringsAsFactors = FALSE, header = FALSE)
    matlist[["listcon1"]] <- list(as.matrix(mat1))
  }
  if(is.null(dirin2) == FALSE)
  {
    mat2 <- NULL
    mat2 <- read.table(dirin2, stringsAsFactors = FALSE, header = FALSE)
    matlist[["listcon2"]] <- list(as.matrix(mat2))
  }
  
}

return(matlist)
  
  
}
