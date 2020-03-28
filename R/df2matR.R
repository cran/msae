#' @title Transform Dataframe to Matrix R
#' @description This function transforms dataframe contains sampling variance to block diagonal matrix R
#' @param var.df dataframe of sampling variances of direct estimators.
#' @param r number of variables
#' @return Block diagonal matrix R
#' @examples NULL
#' @export df2matR
df2matR <- function(var.df , r){
  if(dim(var.df)[2] != sum(1:r))
    stop("Length of vardir is not appropiate with data")
  var.df <- as.data.frame(var.df)
  n <- nrow(var.df)
  R <- lapply(var.df,diag)

  R_1n <- matrix()
  for (i in 1:r){
    R.row <- R[[i]]
    for (j in i:r){
      if(i!=j)
        R.row <- cbind(R.row, R[[sum((r-i):r) + j - r]])
    }
    if(i == 1 ){
      R_1n <- R.row
    } else {
      tmp <- matrix(rep(0,n*n*(i-1)), n, n*(i-1))
      R.row <- cbind(tmp, R.row)
      R_1n <- rbind(R_1n, R.row)
    }
  }
  for(i in 1 : (r*n)){
    for(j in i : (r*n)){
      if(R_1n[j,i] != R_1n[i,j]){
        R_1n[j,i] <- R_1n[i,j]
      }
    }
  }
  return(R_1n)
}
