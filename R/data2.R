#' Data generated based on Autoregressive Multivariate Fay Herriot Model (Model 2)
#'
#' This data is generated based on autoregressive multivariate Fay-Herriot model (model 2) by following these steps:
#' \enumerate{
#'   \item Generate sampling error \code{e}, random effect \code{u}, and auxiliary variables \code{X1 X2}.
#'   \cr For sampling error \code{e}, we take \eqn{\sigma_{e11}}{\sigmae11} = 0.1 , \eqn{\sigma_{e22}}{\sigmae22} = 0.2 ,  \eqn{\sigma_{e33}}{\sigmae33} = 0.3 , and \eqn{\rho_{e}}{\rhoe} = 0.5. Sampling variance-covariance is square of these sampling errors.
#'   \cr For random effect \code{u}, we take \eqn{\sigma_{u}}{\sigmau} = 0.4 , and  \eqn{\rho_{u}}{\rhou} = 0.5.
#'   \cr For auxiliary variables \code{X1 and X2}, we take \eqn{\mu_{X1}}{\muX1} = 5 , \eqn{\mu_{X2}}{\muX2} = 5 , \eqn{\sigma_{X11}}{\sigmaX11} = 1 ,  \eqn{\sigma_{X22}}{\sigmaX22} = 2 , and \eqn{\rho_{x}}{\rhox} = 0.5. The formula of auxiliary variables is following Roberto Benavent and Domingo Morales (2015) <doi:10.1016/j.csda.2015.07.013>.
#'   \item Calculate direct estimation \code{Y1 Y2 and Y3} , where \eqn{Y_{i}}{Yi} = \eqn{X * \beta + u_{i} + e_{i}}{X\beta+ui+ei}
#' }
#' Auxiliary variables \code{X1 X2}, direct estimation \code{Y1 Y2 Y3}, and sampling variance-covariance \code{v1 v2 v3 v12 v13 v23}  are combined into a dataframe called datasae2.
#'
#' @format A data frame with 30 rows and 11 variables:
#' \describe{
#'   \item{X1}{Auxiliary variable of X1}
#'   \item{X2}{Auxiliary variable of X2}
#'   \item{Y1}{Direct Estimation of Y1}
#'   \item{Y2}{Direct Estimation of Y2}
#'   \item{Y3}{Direct Estimation of Y3}
#'   \item{v1}{Sampling Variance of Y1}
#'   \item{v12}{Sampling Covariance of Y1 and Y2}
#'   \item{v13}{Sampling Covariance of Y1 and Y3}
#'   \item{v2}{Sampling Variance of Y2}
#'   \item{v23}{Sampling Covariance of Y2 and Y3}
#'   \item{v3}{Sampling Variance of Y3}
#' }
#' @section Reference: Benavent, Roberto & Morales, Domingo. (2015). Multivariate Fay-Herriot models for small area estimation. Computational Statistics & Data Analysis. 100. 372-390. DOI: 10.1016/j.csda.2015.07.013.
"datasae2"
