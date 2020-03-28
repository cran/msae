#' @title EBLUPs based on a Univariate Fay Herriot (Model 0)
#' @description This function gives the EBLUP and MSE based on a univariate Fay Herriot model (model 0)
#' @param data dataframe containing the variables named in \code{formula} and \code{vardir}
#' @param formula an object of class list of formula, describe the model to be fitted
#' @param vardir if data is available, it is vector containing name of sampling variances of direct estimators. if not, it is data frame of sampling variances of direct estimators. The order is : \code{var1, var2, . , var(k) , cov12, . cov1k, cov23, . , cov(k-1)(k)}
#' @param samevar logical input, true if variance of the data is same, Default: \code{FALSE}
#' @param MAXITER maximum number of iterations allowed in the Fisher-scoring algorithm, Default: \code{100}
#' @param PRECISION convergence tolerance limit for the Fisher-scoring algorithm, Default: \code{1e-4}
#' @return The function returns a list with the following objects:
#' \describe{
#'   \item{eblup}{a dataframe with the values of the EBLUP estimators}
#'   \item{MSE}{a dataframe with the estimated mean squared errors of the EBLUPs for the small domains}
#'   \item{randomEffect}{a dataframe with the values of the random effect estimators}
#'   \item{Rmatrix}{a block diagonal matrix composed of sampling errors}
#'   \item{fit}{a list containing the following objects:}
#'   \itemize{
#'     \item method : type of fitting method, named "REML"
#'     \item convergence : a logical value of convergence of Fisher Scoring algorithm
#'     \item iterations : number of iterations performed by Fisher-Scoring algorithm
#'     \item estcoef : a dataframe with the estimated model coefficient in the first column, their standard error in the second column, the t statistics in the third column, and the p-values of the significance of each coefficient in the last column
#'     \item refvar : a dataframe with the estimated random effect variance
#'     \item informationFisher : a matrix of information Fisher of Fisher-Scoring algorithm
#'   }
#' }
#' @examples
#' ## Load dataset
#' data(datasae1)
#'
#' # Compute EBLUP and MSE of Y1 Y2 and Y3  based on Model 0
#' # using auxiliary variables X1 and X2 for each dependent variable
#'
#' ## Using parameter 'data'
#' Fo <- list(f1=Y1~X1+X2,
#'            f2=Y2~X1+X2,
#'            f3=Y3~X1+X2)
#' vardir <- c("v1", "v2", "v3", "v12", "v13", "v23")
#' un <- eblupUFH(Fo, vardir, data=datasae1)
#'
#' ## Without parameter 'data'
#' Fo <- list(f1=datasae1$Y1~datasae1$X1+datasae1$X2,
#'            f2=datasae1$Y2~datasae1$X1+datasae1$X2,
#'            f3=datasae1$Y3~datasae1$X1+datasae1$X2)
#' vardir <- datasae1[,c("v1", "v2", "v3", "v12", "v13", "v23")]
#' un <- eblupMFH1(Fo, vardir)
#'
#' un$eblup   # see the EBLUP estimators
#' un$MSE   # see MSE of EBLUP estimators
#'
#' @export eblupUFH
eblupUFH <- function (formula, vardir, samevar = FALSE, MAXITER = 100, PRECISION = 1e-04, data){
  r <- length(formula)

  if(!missing(data)){
    formula.matrix <- lapply(formula, function(x){model.frame(x,na.action = na.omit,data)})
    y.matrix <- unlist(lapply(formula, function(x){model.frame(x,na.action = na.omit,data)}[[1]]))
    x.matrix <- Reduce(adiag,lapply(formula, function(x){model.matrix(x,data)}))
    n <- length(y.matrix)/r
    if(!all(vardir %in% names(data)))
      stop("Object vardir is not appropiate with data")
    if(length(vardir) != sum(1:r))
      stop("Length of vardir is not appropiate with data")
    if (any(is.na(data[,vardir])))
      stop("Object vardir contains NA values.")
    R <- df2matR(data[,vardir], r)
  } else {
    formula.matrix <- lapply(formula, function(x){model.frame(x,na.action = na.omit)})
    y.matrix <- unlist(lapply(formula, function(x){model.frame(x, na.action = na.omit)}[[1]]))
    x.matrix <- Reduce(adiag,lapply(formula, function(x){model.matrix(x)}))
    n <- length(y.matrix)/r
    if(dim(vardir)[2]!= sum(1:r)){
      stop("Object vardir is not appropiate with data")
    }
    if (any(is.na(vardir)))
      stop("Object vardir contains NA values.")
    R <- df2matR(vardir, r)
  }

  for(i in 1:r){
    if (attr(attributes(formula.matrix[[i]])$terms,"response")==1)
      textformula <- paste(formula[[i]][2],formula[[i]][1],formula[[i]][3])
    else
      textformula <- paste(formula[[i]][1],formula[[i]][2])

    if (length(na.action(formula.matrix[[i]]))>0){
      stop("Argument formula= ",textformula," contains NA values.")
    }
  }

  y.var <- sapply(formula, "[[", 2)

  R <- diag(diag(R))
  I=diag(n)
  Id <- diag(r)
  d.omega <- list()
  d.Omega <- list()
  for(i in 1:r){
    d.omega[[i]] <- matrix(0, r,r)
    d.omega[[i]][i,i] <- 1
    d.Omega[[i]] <- kronecker(d.omega[[i]], I)
  }
  convergence <- TRUE

  if(samevar){
    sigmau <- median(diag(R))
    k <- 0
    diff <- rep(PRECISION + 1,r)
    while (any(diff > PRECISION) & (k < MAXITER)) {
      k <- k + 1
      sigmau1=sigmau
      GI=kronecker(kronecker(sigmau1, Id),I)
      Omega <- solve(GI + R)
      Q <- solve(t(Omega %*% x.matrix) %*% x.matrix)
      P <- Omega - (Omega %*% x.matrix) %*% Q %*% t(Omega %*% x.matrix)
      Py <- P %*% y.matrix
      s <- (-0.5) * sum(diag(P)) + 0.5 * (t(Py) %*% Py)
      iF <- 0.5 * sum(diag(P%*%P))
      sigmau <- sigmau1 + solve(iF)%*%s
      diff <- abs((sigmau - sigmau1)/sigmau1)
    }
    sigmau <- as.vector(mapply(max, sigmau, rep(0,r)))
    names(sigmau) <- y.var
    if (k >= MAXITER && diff >= PRECISION) {
      convergence <- FALSE }

    GI=kronecker(diag(sigmau),I)
    Omega <- solve(GI + R)
    Qh <- solve(t(Omega %*% x.matrix) %*% x.matrix)
    P <- Omega - Omega %*% x.matrix %*% Qh %*% t(Omega %*% x.matrix)
    Py <- P %*% y.matrix

    beta=Qh%*%t(Omega %*% x.matrix)%*%y.matrix
    res=y.matrix-x.matrix%*%beta
    eblup=x.matrix%*%beta+GI%*%Omega%*%(y.matrix-x.matrix%*%beta)
    eblup.df <- data.frame(matrix(eblup, n, r))
    names(eblup.df) <- y.var
    se.b=sqrt(diag(Qh))
    t.val=beta/se.b
    pv <- 2 * pnorm(abs(t.val), lower.tail = FALSE)
    coef <- cbind(beta, se.b, t.val, pv)
    colnames(coef) <- c("beta", "std.error", "t.statistics", "p.value")

    d=kronecker(Id,I)-GI%*%Omega
    gg1=diag(GI%*%Omega%*%R)
    gg2=diag(d%*%x.matrix%*%Qh%*%t(x.matrix)%*%t(d))
    dg <- Omega - GI %*% Omega %*% Omega
    g3 <- (dg%*%(GI+R) %*% t(dg)) / iF
    gg3 <- diag(g3)
    mse = matrix(gg1+gg2+2*gg3,n)
    mse.df <- data.frame(matrix(0,n,r))
    names(mse.df) <- y.var
    for(i in 1:r)
      mse.df[,i] <- mse[((i-1)*n+1):(i*n)]
  } else {
    sigmau <- apply(matrix(diag(R), n, r), 2,  median)
    k <- 0
    diff <- rep(PRECISION + 1,r)
    while (any(diff > rep(PRECISION,r)) & (k < MAXITER)) {
      k <- k + 1
      sigmau1=sigmau
      if(r == 1){
        G <- sigmau1
      } else {
        G <- diag(as.vector(sigmau1))
      }
      GI=kronecker(G,I)
      Omega <- solve(GI + R)
      Q <- solve(t(Omega %*% x.matrix) %*% x.matrix)
      P <- Omega - Omega %*% x.matrix %*% Q %*% t(Omega %*% x.matrix)
      Py <- P %*% y.matrix
      s <- sapply(d.Omega, function(x) (-0.5) * sum(diag(P%*%x)) + 0.5 * (t(Py)%*%x %*% Py))
      iF <- matrix(unlist(lapply(d.Omega, function(x) lapply(d.Omega, function(y) 0.5 * sum(diag(P%*%x%*%P%*%y))))),r)
      sigmau <- sigmau1 + solve(iF)%*%s
      diff <- abs((sigmau - sigmau1)/sigmau1)
    }
    sigmau <- as.vector(mapply(max, sigmau, rep(0,r)))
    names(sigmau) <- y.var
    if (k >= MAXITER && diff >= PRECISION) {
      convergence <- FALSE}

    if(r == 1){
      G <- sigmau
    } else {
      G <- diag(as.vector(sigmau))
    }
    GI=kronecker(G,I)
    Omega <- solve(GI + R)
    Qh <- solve(t(Omega %*% x.matrix) %*% x.matrix)
    P <- Omega - Omega %*% x.matrix %*% Qh %*% t(Omega %*% x.matrix)
    Py <- P %*% y.matrix

    beta=Qh%*%t(Omega %*% x.matrix)%*%y.matrix
    res=y.matrix-x.matrix%*%beta
    eblup=x.matrix%*%beta+GI%*%Omega%*%(y.matrix-x.matrix%*%beta)
    eblup.df <- data.frame(matrix(eblup, n, r))
    names(eblup.df) <- y.var
    se.b=sqrt(diag(Qh))
    t.val=beta/se.b
    pv <- 2 * pnorm(abs(t.val), lower.tail = FALSE)
    coef <- cbind(beta, se.b, t.val, pv)
    colnames(coef) <- c("beta", "std.error", "t.statistics", "p.value")

    FI=solve(iF)
    d=kronecker(Id,I)-GI%*%Omega
    gg1=diag(GI%*%Omega%*%R)
    gg2=diag(d%*%x.matrix%*%Qh%*%t(x.matrix)%*%t(d))
    dg <- lapply(d.Omega, function(x) x %*% Omega - GI %*% Omega %*% x %*% Omega)
    g3 <- list()
    for (i in 1:r){
      for (j in 1:r){
        g3[[(i-1)*r+j]] <- FI[i,j]*(dg[[i]]%*%(GI+R) %*% t(dg[[j]]))
      }
    }
    gg3 <- diag(Reduce('+', g3))
    mse = matrix(gg1 +gg2 +2*gg3,n)
    mse.df <- data.frame(matrix(0,n,r))
    names(mse.df) <- y.var
    for(i in 1:r)
      mse.df[,i] <- mse[((i-1)*n+1):(i*n)]
  }
  u.cap <- GI %*% Omega %*% res
  u.cap.df <- as.data.frame(matrix(u.cap, n, r))
  names(u.cap.df) <- y.var

  result <- list(eblup = NA,
                 MSE = NA,
                 randomEffect=NA,
                 Rmatrix=NA,
                 fit = list(method = NA,
                            convergence = NA ,
                            iterations = NA,
                            estcoef= NA,
                            refvar=NA,
                            informationFisher=NA
                 ))
  result$eblup <- signif(eblup.df, digits=5)
  result$MSE <- signif(mse.df, digits=5)
  result$randomEffect <- signif(u.cap.df, digits=5)
  result$Rmatrix <- signif(R, digits=5)
  result$fit$method <- "REML"
  result$fit$convergence <- convergence
  result$fit$iterations <- k
  result$fit$refvar <- signif(data.frame(t(sigmau)), digits=5)
  result$fit$informationFisher <- signif(iF, digits=5)
  result$fit$estcoef <- signif(coef, digits=5)
  return(result)
}
