#' msae : Multivariate Fay Herriot Models for Small Area Estimation
#'
#' Implements multivariate Fay-Herriot models for small area estimation. It uses empirical best linear unbiased prediction (EBLUP) estimator. Multivariate models consider the correlation of several target variable and borrow strength from auxiliary variables to improve the effectiveness of a domain sample size. Models which accomodated by this package are univariate model with several target variables (model 0), multivariate model (model 1), autoregressive multivariate model (model 2), and heteroscedastic autoregressive multivariate model (model 3). Functions provide EBLUP estimators and mean squared error (MSE) estimator for each model. These models were developed by Roberto Benavent and Domingo Morales (2015) <doi:10.1016/j.csda.2015.07.013>.
#'
#'
#' @section Author(s):
#' Novia Permatasari, Azka Ubaidillah
#'
#' \strong{Maintainer}: Novia Permatasari \email{16.9335@@stis.ac.id}
#'
#' @section Functions:
#' \describe{
#'   \item{\code{\link{eblupUFH}}}{Gives the EBLUPs and MSE of Univariate SAE (Model 0)}
#'   \item{\code{\link{eblupMFH1}}}{Gives the EBLUPs and MSE of Multivariate SAE (Model 1)}
#'   \item{\code{\link{eblupMFH2}}}{Gives the EBLUPs and MSE of Autoregressive Multivariate SAE (Model 2)}
#'   \item{\code{\link{eblupMFH3}}}{Gives the EBLUPs and MSE of Heteroscedastics Autoregressive Multivariate SAE (Model 3)}
#' }
#'
#' @section Reference:
#' \itemize{
#'   \item{Benavent, Roberto & Morales, Domingo. (2015). Multivariate Fay-Herriot models for small area estimation. Computational Statistics & Data Analysis. 100. 372-390. DOI: 10.1016/j.csda.2015.07.013.}
#'   \item{Rao, J.N.K & Molina. (2015). Small Area Estimation 2nd Edition. New York: John Wiley and Sons, Inc.}
#'   \item{Ubaidillah, Azka et al. (2019). Multivariate Fay-Herriot models for small area estimation with application to household consumption per capita expenditure in Indonesia. Journal of Applied Statistics. 46:15. 2845-2861. DOI: 10.1080/02664763.2019.1615420.}
#'   }
#'
#' @docType package
#' @name msae
#' @import magic
NULL
