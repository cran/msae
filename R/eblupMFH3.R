#' @title EBLUPs based on a Heteroscedastic Autoregressive Multivariate Fay Herriot (Model 3)
#' @description This function gives the EBLUP and MSE based on a heteroscedastic autoregressive multivariate Fay-Herriot model (model 3).
#' @param data dataframe containing the variables named in \code{formula} and \code{vardir}
#' @param formula an object of class list of formula, describe the model to be fitted
#' @param vardir if data is available, it is vector containing name of sampling variances of direct estimators. if not, it is data frame of sampling variances of direct estimators. The order is : \code{var1, var2, . , var(k) , cov12, . cov1k, cov23, . , cov(k-1)(k)}
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
#'     \item refvarTest : homogeneity of random effect variance test based on Model 3
#'     \item rho : a dataframe with the estimated rho of random effect variance and their rho parameter test based on Model 2
#'     \item informationFisher : a matrix of information Fisher of Fisher-Scoring algorithm
#'   }
#' }
#' @examples
#' ## Load dataset
#' data(datasae3)
#'
#' # Compute EBLUP and MSE of Y1 Y2 and Y3  based on Model 3
#' # using auxiliary variables X1 and X2 for each dependent variable
#'
#' ## Using parameter 'data'
#' Fo <- list(f1=Y1~X1+X2,
#'            f2=Y2~X1+X2,
#'            f3=Y3~X1+X2)
#' vardir <- c("v1", "v2", "v3", "v12", "v13", "v23")
#' m3 <- eblupMFH3(Fo, vardir, data=datasae3)
#'
#' ## Without parameter 'data'
#' Fo <- list(f1=datasae3$Y1~datasae3$X1+datasae3$X2,
#'            f2=datasae3$Y2~datasae3$X1+datasae3$X2,
#'            f3=datasae3$Y3~datasae3$X1+datasae3$X2)
#' vardir <- datasae3[,c("v1", "v2", "v3", "v12", "v13", "v23")]
#' m3 <- eblupMFH3(Fo, vardir)
#'
#' m3$eblup   # see the EBLUP estimators
#' m3$MSE   # see MSE of EBLUP estimators
#'
#' @export eblupMFH3
eblupMFH3 <- function (formula, vardir, MAXITER = 100, PRECISION = 1e-04, data){
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

  I=diag(n)
  Id <- diag(r)
  d.Omega <- list()
  sigmau <- apply(matrix(diag(R), n, r), 2,  median)
  sigmau[r+1] <- 1e-04
  kit <- 0
  diff <- rep(PRECISION + 1,(r+1))
  convergence <- TRUE

  while (any(diff > rep(PRECISION,(r+1))) & (kit < MAXITER)) {
    kit <- kit + 1
    sigmau1=sigmau
    G <- matrix(NA, r,r)
    if(r == 1){
      G <- sigmau1[r] +  sigmau1[r+1]^2
    } else {
      for(i in 1:r){
        for(j in 1:r){
          if(i == j){
            G[i,j] <- sigmau1[i]
            for(l in 1:i){
              if(i == l){
                G[i,j] <- G[i,j] + sigmau1[r+1]^(2*l)
              } else {
                G[i,j] <- G[i,j] + sigmau1[r+1]^(2*l)*sigmau1[i-l]
              }
            }
          } else {
            ki <- abs(i-j)
            G[i,j] <- sigmau1[ki]  * sigmau1[r+1]^(ki)
            for(l in 1:ki){
              if(ki == l){
                G[i,j] <- G[i,j] + sigmau1[r+1]^(2*l+ki)
              } else {
                G[i,j] <- G[i,j] + sigmau1[r+1]^(2*l+ki)*sigmau1[ki-l]
              }
            }
          }
        }
      }
    }
    GI=kronecker(G,I)
    Omega <- solve(GI + R)
    Xto <- t(Omega %*% x.matrix)
    Q <- solve(Xto %*% x.matrix)
    P <- Omega - t(Xto) %*% Q %*% Xto
    for(ki in 1:r){
      d.Omega[[ki]] <- matrix(0,r,r)
      for(i in 1:r){
        for(j in 1:r){
          if(i==j){
            if(i<ki){
              d.Omega[[ki]][i,j] <- 0
            } else {
              d.Omega[[ki]][i,j] <- sigmau1[r+1]^(2*(i-ki))
            }
          } else {
            s <- abs(i-j)
            if(s<ki){
              d.Omega[[ki]][i,j] <- 0
            } else {
              d.Omega[[ki]][i,j] <- sigmau1[r+1]^(3*(s-ki)+ki)
            }
          }
        }
      }
      d.Omega[[ki]] <- kronecker(d.Omega[[ki]], I)
    }
    d.Omega[[r+1]] <- matrix(0,r,r)
    for(i in 1:r){
      for(j in 1:r){
        if(i == j){
          for(ki in 1:i){
            if(ki == i){
              d.Omega[[r+1]][i,j] <- d.Omega[[r+1]][i,j] + 2*ki*sigmau1[r+1]^(2*ki-1)
            } else {
              d.Omega[[r+1]][i,j] <-  d.Omega[[r+1]][i,j] + 2*ki*sigmau1[r+1]^(2*ki-1)  * sigmau1[i-ki]
            }
          }
        } else {
          n <- abs(i-j)
          d.Omega[[r+1]][i,j] <- n*sigmau1[r+1]^(n-1)*sigmau1[n]
          for(ki in 1:n){
            if(ki == n){
              d.Omega[[r+1]][i,j] <-  d.Omega[[r+1]][i,j] + (2*ki+n)*sigmau1[r+1]^(n+2*ki-1)
            } else {
              d.Omega[[r+1]][i,j] <- d.Omega[[r+1]][i,j] + (2*ki+n)*sigmau1[r+1]^(n+2*ki-1) * sigmau1[n-ki]
            }
          }
        }
      }
    }
    d.Omega[[r+1]] <- kronecker(d.Omega[[r+1]], I)
    Py <- P %*% y.matrix
    s <- sapply(d.Omega, function(x) (-0.5) * sum(diag(P%*%x)) + 0.5 * (t(Py)%*%x %*% Py))
    iF <- matrix(unlist(lapply(d.Omega, function(x) lapply(d.Omega, function(y) 0.5 * sum(diag(P%*%x%*%P%*%y))))),(r+1))
    sigmau <- sigmau1 + solve(iF)%*%s
    if(abs(sigmau[r+1]) > 1){
      sigmau[r+1] <- sigmau1[r+1]
    }
    diff <- abs((sigmau - sigmau1)/sigmau1)
  }
  sigmau[1:r] <- mapply(max, sigmau[1:r] , rep(0,r))
  if (kit >= MAXITER && diff >= PRECISION) {
    convergence <- FALSE }

  if(r == 1){
    G <- sigmau[r] +  sigmau[r+1]^2
  } else {
    for(i in 1:r){
      for(j in 1:r){
        if(i == j){
          G[i,j] <- sigmau[i]
          for(l in 1:i){
            if(i == l){
              G[i,j] <- G[i,j] + sigmau[r+1]^(2*l)
            } else {
              G[i,j] <- G[i,j] + sigmau[r+1]^(2*l)*sigmau[i-l]
            }
          }
        } else {
          k <- abs(i-j)
          G[i,j] <- sigmau[k]*sigmau[r+1]^(k)
          for(l in 1:k){
            if(k == l){
              G[i,j] <- G[i,j] + sigmau[r+1]^(2*l+k)
            } else {
              G[i,j] <- G[i,j] + sigmau[r+1]^(2*l+k)*sigmau[k-l]
            }
          }
        }
      }
    }
  }
  GI=kronecker(G,I)
  Omega <- solve(GI + R)
  Xto <- t(Omega %*% x.matrix)
  Qh <- solve(Xto %*% x.matrix)
  P <- Omega - t(Xto) %*% Qh %*% Xto
  Py <- P %*% y.matrix

  for(k in 1:r){
    d.Omega[[k]] <- matrix(0,r,r)
    for(i in 1:r){
      for(j in 1:r){
        if(i==j){
          if(i<k){
            d.Omega[[k]][i,j] <- 0
          } else {
            d.Omega[[k]][i,j] <- sigmau[r+1]^(2*(i-k))
          }
        } else {
          si <- abs(i-j)
          if(si<k){
            d.Omega[[k]][i,j] <- 0
          } else {
            d.Omega[[k]][i,j] <- sigmau[r+1]^(3*(si-k)+k)
          }
        }
      }
    }
    d.Omega[[k]] <- kronecker(d.Omega[[k]], I)
  }
  d.Omega[[r+1]] <- matrix(0,r,r)
  for(i in 1:r){
    for(j in 1:r){
      if(i == j){
        for(k in 1:i){
          if(k == i){
            d.Omega[[r+1]][i,j] <- d.Omega[[r+1]][i,j] + 2*k*sigmau[r+1]^(2*k-1)
          } else {
            d.Omega[[r+1]][i,j] <-  d.Omega[[r+1]][i,j] + 2*k*sigmau[r+1]^(2*k-1)  * sigmau[i-k]
          }
        }
      } else {
        n <- abs(i-j)
        d.Omega[[r+1]][i,j] <- n*sigmau[r+1]^(n-1)*sigmau[n]
        for(k in 1:n){
          if(k == n){
            d.Omega[[r+1]][i,j] <-  d.Omega[[r+1]][i,j] + (2*k+n)*sigmau[r+1]^(n+2*k-1)
          } else {
            d.Omega[[r+1]][i,j] <- d.Omega[[r+1]][i,j] + (2*k+n)*sigmau[r+1]^(n+2*k-1) * sigmau[n-k]
          }
        }
      }
    }
  }
  d.Omega[[r+1]] <- kronecker(d.Omega[[r+1]], I)

  n <- length(y.matrix)/r
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
  n <- length(y.matrix)/r
  d=kronecker(Id,I)-GI%*%Omega
  gg1=diag(GI%*%Omega%*%R)
  gg2=diag(d%*%x.matrix%*%Qh%*%t(x.matrix)%*%t(d))
  dg <- lapply(d.Omega, function(x) x %*% Omega - GI %*% Omega %*% x %*% Omega)
  g3 <- list()
  for (i in 1:(r+1)){
    for (j in 1:(r+1)){
      g3[[(i-1)*(r+1)+j]] <- FI[i,j]*(dg[[i]]%*%(GI+R) %*% t(dg[[j]]))
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

  T.test <- matrix(NA, r,r)
  p.val <- matrix(NA, r,r)
  for(i in 1:r){
    for(j in i:r){
      T.test[i,j] <- signif((sigmau[i] - sigmau[j]) / (sqrt(FI[i,i] + FI[j,j] - 2*FI[i,j])), digits = 5)
      if(is.na(T.test[i,j])){
        p.val[i,j] <- NA
      } else {
        p.val[i,j] <- signif(2 * pnorm(abs(T.test[i,j]), lower.tail = FALSE), digits = 5)
      }
    }
  }
  var <- NULL
  for(i in 1:(dim(p.val)[1]-1)){
    for(j in (i+1):dim(p.val)[2]){
      var[(dim(p.val)[1])*(j-1)+i] <-   paste(y.var[i] ," vs ", y.var[j])
    }
  }
  T.test.temp <- cbind(var, as.vector(T.test)[1:length(var)], as.vector(p.val)[1:length(var)])
  T.test.dt <- as.data.frame(na.omit(T.test.temp))
  colnames(T.test.dt) <- c("refvar", "T-test", "p-value")

  T.test <- (sigmau[r+1]) / (sqrt(FI[r+1,r+1]))
  if(T.test > 0){
    p.val <- 1-pnorm(T.test)
  } else {
    p.val <- pnorm(T.test)
  }
  rho.df <- data.frame(signif(data.frame(sigmau[r+1], T.test, p.val), digits = 5))
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
                            refvarTest=NA,
                            rho=NA,
                            informationFisher=NA
                 ))

  result$eblup <- signif(eblup.df, digits = 5)
  result$MSE <- signif(mse.df, digits = 5)
  result$randomEffect <- signif(u.cap.df, digits = 5)
  result$Rmatrix <- signif(R, digits = 5)
  result$fit$method <- "REML"
  result$fit$convergence <- convergence
  result$fit$iterations <- kit
  result$fit$estcoef <- signif(coef, digits = 5)
  result$fit$refvar <- signif(sigmau[1:r], digits = 5)
  names(result$fit$refvar) <- c(y.var)
  result$fit$refvarTest <- T.test.dt
  result$fit$rho <- rho.df
  result$fit$informationFisher <- signif(iF, digits = 5)
  return(result)
}
