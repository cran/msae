#' @title EBLUPs based on a Autoregressive Multivariate Fay Herriot (Model 2)
#' @description This function gives the EBLUP and MSE based on a  autoregressive multivariate Fay-Herriot model (model 2).
#' @param data dataframe containing the variables named in \code{formula} and \code{vardir}
#' @param formula an object of class list of formula, describe the model to be fitted
#' @param vardir if data is available, it is vector containing name of sampling variances of direct estimators. if not, it is data frame of sampling variances of direct estimators. The order is : \code{var1, var2, . , var(k) , cov12, . cov1k, cov23, . , cov(k-1)(k)}
#' @param MAXITER maximum number of iterations allowed in the Fisher-scoring algoritm, Default: \code{100}
#' @param PRECISION convergence tolerance limit for the Fisher-scoring algorithm, Default: \code{1e-4}
#' @return The function returns a list with the following objects:
#' \describe{
#'   \item{eblup}{a dataframe with the values of the EBLUP estimators}
#'   \item{MSE}{a dataframe with the estimated mean squared errors of the EBLUPs for the small domains}
#'   \item{randomEffect}{a dataframe with the values of the random effect estimators}
#'   \item{fit}{a list containing the following objects:}
#'   \itemize{
#'     \item method : type of fitting method, named "REML"
#'     \item convergence : a logical value of convergence of Fisher Scoring algorithm
#'     \item iterations : number of iterations performed by Fisher-Scoring algorithm
#'     \item estcoef : a dataframe with the estimated model coefficient in the first column, their standart error in the second column, the t statistics in the third column, and the p-values of the significance of each coefficient in the last column
#'     \item refvar : a dataframe with the estimated random effect variance
#'     \item rho : a dataframe with the estimated rho of random effect variance and their rho parameter test based on Model 2
#'     \item informationFisher : a matrix of information Fisher of Fisher-Scoring algorithm
#'   }
#' }
#' @examples
#' data(datasae2)
#' Fo <- list(f1=Y1~X1+X2,
#' f2=Y2~X1+X2,
#' f3=Y3~X1+X2)
#' vardir <- c("v1", "v2", "v3", "v12", "v13", "v23")
#' eblupMFH2(Fo, vardir, data=datasae2)
#'
#' @export eblupMFH2
eblupMFH2 <- function (formula, vardir, MAXITER = 100, PRECISION = 1e-04, data){
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
    vardir <- data[,vardir]
    R <- df2matR(vardir, r)
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

  I=diag(n)
  Id <- diag(r)
  d.Omega <- list()
  sigmau <- c(mean(sapply(vardir[,1:r], median)), 1 - 1e-04)
  ki <- 0
  diff <- rep(PRECISION + 1,2)
  convergence <- TRUE

  while (any(diff > rep(PRECISION,2)) & (ki < MAXITER)) {
    ki <- ki + 1
    sigmau1=sigmau
    Omega_AR <- matrix(0, r,r)
    if(r == 1){
      G <- sigmau1[1]/(1-sigmau1[2]^2)
    } else {
      for(i in 1:r){
        for(j in 1:r){
          Omega_AR[i,j] <- sigmau1[2]^(abs(i-j))/(1-sigmau1[2]^2)
        }
      }
      G <- sigmau1[1] * Omega_AR
    }
    GI=kronecker(G,I)
    Omega <- solve(GI + R)
    Xto <- t(Omega %*% x.matrix)
    Q <- solve(Xto %*% x.matrix)
    P <- Omega - t(Xto) %*% Q %*% Xto
    d.Omega[[1]] <- kronecker(Omega_AR, I)
    d.Omega[[2]] <- matrix(NA, r,r)
    for(i in 1:r){
      for(j in 1:r){
        k <- abs(i-j)
        d.Omega[[2]][i,j] <- sigmau1[1]*(k*sigmau1[2]^(k-1) + (2-k)*(sigmau1[2]^(k+1)))/(1-sigmau1[2]^2)^2
      }
    }
    d.Omega[[2]] <- kronecker(d.Omega[[2]], I)

    Py <- P %*% y.matrix
    s <- sapply(d.Omega, function(x) (-0.5) * sum(diag(P%*%x)) + 0.5 * (t(Py)%*%x %*% Py))
    iF <- matrix(unlist(lapply(d.Omega, function(x) lapply(d.Omega, function(y) 0.5 * sum(diag(P%*%x%*%P%*%y))))),2)
    sigmau <- sigmau1 + solve(iF)%*%s
    if(abs(sigmau[2]) > 1){
      sigmau[2] <- sigmau1[2]
    }
    diff <- abs((sigmau - sigmau1)/sigmau1)
  }
  sigmau[1] <- max(sigmau[1], 0)
  if (ki >= MAXITER && diff >= PRECISION) {
    convergence <- FALSE }

  if(r == 1){
    G <- sigmau[1]/(1-sigmau[2]^2)
  } else {
    for(i in 1:r){
      for(j in 1:r){
        Omega_AR[i,j] <- sigmau[2]^(abs(i-j))/(1-sigmau[2]^2)
      }
    }
    G <- sigmau[1] * Omega_AR
  }
  GI=kronecker(G,I)
  Omega <- solve(GI + R)
  Xto <- t(Omega %*% x.matrix)
  Qh <- solve(Xto %*% x.matrix)
  P <- Omega - t(Xto) %*% Qh %*% Xto
  Py <- P %*% y.matrix
  d.Omega[[1]] <- kronecker(Omega_AR, I)
  d.Omega[[2]] <- matrix(NA, r,r)
  for(i in 1:r){
    for(j in 1:r){
      k <- abs(i-j)
      d.Omega[[2]][i,j] <- sigmau[1]*(k*sigmau[2]^(k-1) + (2-k)*(sigmau[2]^(k+1)))/((1-sigmau[2]^2)^2)
    }
  }
  d.Omega[[2]] <- kronecker(d.Omega[[2]], I)

  beta=Qh%*%Xto%*%y.matrix
  res=y.matrix-x.matrix%*%beta
  eblup=x.matrix%*%beta+GI%*%Omega%*%res
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
  for (i in 1:2){
    for (j in 1:2){
      g3[[(i-1)*2+j]] <- FI[i,j]*(dg[[i]]%*%(GI+R) %*% t(dg[[j]]))
    }
  }
  gg3 <- diag(Reduce('+', g3))
  mse = gg1 + gg2 + 2*gg3
  mse.df <- data.frame(matrix(0,n,r))
  names(mse.df) <- y.var
  for(i in 1:r)
    mse.df[,i] <- mse[((i-1)*n+1):(i*n)]

  u.cap <- GI %*% Omega %*% res
  u.cap.df <- as.data.frame(matrix(u.cap, n, r))
  names(u.cap.df) <- y.var

  F.inv <- solve(iF)
  T.test <- (sigmau[2,]) / (sqrt(F.inv[2,2]))
  if(T.test > 0){
    p.val <- 1-pnorm(T.test)
  } else {
    p.val <- pnorm(T.test)
  }

  rho.df <- data.frame(signif(data.frame(sigmau[2], T.test, p.val), digits = 5))
  names(rho.df) <- c("rho","T.test", "p-value")

  result <- list(eblup = NA,
                 MSE = NA,
                 randomEffect=NA,
                 Rmatrix=NA,
                 fit = list(method = NA,
                            convergence = NA ,
                            iterations = NA,
                            estcoef= NA,
                            refvar=NA,
                            rho=NA,
                            informationFisher=NA
                 ))
  result$eblup <- signif(eblup.df, digits = 5)
  result$MSE <- signif(mse.df, digits = 5)
  result$randomEffect <- signif(u.cap.df, digits = 5)
  result$Rmatrix <- signif(R, digits = 5)
  result$fit$method <- "REML"
  result$fit$convergence <- convergence
  result$fit$iterations <- ki
  result$fit$estcoef <- coef
  result$fit$refvar <- signif(sigmau[1,], digits = 5)
  result$fit$rho <- rho.df
  result$fit$informationFisher <- signif(iF, digits = 5)
  return(result)
}
