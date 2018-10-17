#function name:zerorecover
#function:recover the recorded full zero rows and columns of a matrix after using the function zerofilter
#parameters:
#            mat:the input matrix which need to be recovered
#            filtedobj:the recorded full zero rows or columns of a matrix
zerorecover <- function(mat, filtedobj) 
{
  if(is.null(filtedobj$lines) == FALSE)
  {
    lines <- filtedobj$lines
    rnames <- filtedobj$rownames
    cnames <- filtedobj$colnames
    for (i in 1:length(lines))
    {
      if(lines[i] == 1)
      {
        mat <- rbind(0, mat[1:nrow(mat),])
        mat <- cbind(0, mat[,1:ncol(mat)])
      }
      else if
      (lines[i] == nrow(mat)+1)
      {
        mat <- rbind(mat[1:nrow(mat),], 0)
        mat <- cbind(mat[,1:ncol(mat)], 0)
      }
      else if(lines[i] >1 && lines[i] <= nrow(mat))
      {
        mat <- rbind(mat[1:(lines[i]-1),], 0, mat[lines[i]:nrow(mat),])
        mat <- cbind(mat[,1:(lines[i]-1)], 0, mat[,lines[i]:ncol(mat)])
      }
      else
      {
        stop("unkonwn columns found")
      }
      rownames(mat)[lines[i]] <- rnames[i]
      colnames(mat)[lines[i]] <- cnames[i]
    }    
  }
  return (mat)
}
