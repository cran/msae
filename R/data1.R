#' Data generated based on Multivariate Fay Herriot Model (Model 1)
#'
#' This data is generated based on multivariate Fay-Herriot model (model 1) by these following steps:
#' \enumerate{
#'   \item Generate sampling error \code{e}, random effect \code{u}, and auxiliary variables \code{X1 X2}.
#'   \itemize{
#'     \item For sampling error \code{e}, we set \eqn{e_{d}}{ed} ~ \eqn{N_{3}(0, V_{ed})}{N3(0, Ved)} , where \eqn{V_{ed} = (\sigma_{dij})_{i,j=1,2,3}}{Ved = (\sigmadij)i,j=1,2,3}, with  \eqn{\sigma_{11}}{\sigma11}  ~ \eqn{InvGamma(11, 1)}{InvGamma(11, 1)}, \eqn{\sigma_{22}}{\sigma22}  ~ \eqn{InvGamma(11, 2)}{InvGamma(11, 2)},  \eqn{\sigma_{33}}{\sigma33}  ~ \eqn{InvGamma(11, 3)}{InvGamma(11, 3)}, and \eqn{\rho_{e}}{\rhoe} = 0.5.
#'     \item For random effect \code{u}, we set \eqn{u}{u} ~ \eqn{N_{3}(0, V_{u})}{N2(0, Vu)} , where \eqn{\sigma_{u11}}{\sigmau11} = 0.2 , \eqn{\sigma_{u22}}{\sigmau22} = 0.4 , and \eqn{\sigma_{u33}}{\sigmau33} = 1.2.
#'     \item For auxiliary variables \code{X1 and X2}, we set \eqn{X1}{X1} ~ \eqn{N(5, 0.1)}{N(5, 0.1)} and \eqn{X2}{X2} ~ \eqn{N(10, 0.2)}{N(10, 0.2)}.
#'   }
#'   \item Calculate direct estimation \code{Y1 Y2 and Y3} , where \eqn{Y_{i}}{Yi} = \eqn{X * \beta + u_{i} + e_{i}}{X\beta+ui+ei}. We take \eqn{\beta_{1} = 5}{\beta1 = 5} and \eqn{\beta_{2} = 10}{\beta2 = 10}.
#' }
#' Auxiliary variables \code{X1 X2}, direct estimation \code{Y1 Y2 Y3}, and sampling variance-covariance \code{v1 v2 v3 v12 v13 v23} are combined into a dataframe called datasae1.
#'
#' @format A data frame with 50 rows and 11 variables:
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
"datasae1"
