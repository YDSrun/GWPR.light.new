#' Moran's I Test for Panel Regression
#'
#' @description Moran's I test for spatial autocorrelation in residuals from
#'              an estimated panel linear model (plm).
#'
#' @usage GWPR.moran.test(plm_model, SDF, bw, adaptive = FALSE, p = 2,
#'                        kernel = "bisquare", longlat = FALSE, alternative = "greater")
#'
#' @param plm_model     An object of class inheriting from "plm", see plm
#' @param SDF           Spatial*DataFrame on which is based the data, with the "ID" in the index
#' @param bw            The optimal bandwidth, either adaptive or fixed distance
#' @param adaptive      If TRUE, adaptive distance bandwidth is used, otherwise, fixed distance bandwidth.
#' @param p             The power of the Minkowski distance, default is 2, i.e. the Euclidean distance
#' @param kernel        bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise (default);
#'                            gaussian: wgt = exp(-.5*(vdist/bw)^2);
#'                            exponential: wgt = exp(-vdist/bw);
#'                            tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise;
#'                            boxcar: wgt=1 if dist < bw, wgt=0 otherwise
#' @param longlat       If TRUE, great circle distances will be calculated
#' @param alternative   A character string specifying the alternative hypothesis, must be one of greater (default), less or two.sided.
#'
#' @import dplyr
#' @import GWmodel
#' @importFrom sp merge coordinates
#' @importFrom plm pdim index
#' @importFrom methods is
#' @importFrom stats pnorm
#'
#' @return A list of result:
#' \describe{
#' \item{statistic}{the value of the standard deviate of Moran's I.}
#' \item{p.value}{the p-value of the test.}
#' \item{Estimated.I}{the value of the observed Moran's I.}
#' \item{Excepted.I}{the value of the expectation of Moran's I.}
#' \item{V2}{the value of the variance of Moran's I.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' }
#' @export
#'
#' @note: Current version of panel Moran's I test can only chech the balanced panel data.
#'
#' @author Chao Li <chaoli0394@gmail.com> Shunsuke Managi
#'
#' @references Beenstock, M., Felsenstein, D., 2019. The econometric analysis of non-stationary spatial panel data. Springer.
#'
#' @examples
#' data(TransAirPolCalif)
#' data(California)
#' formula.GWPR <- pm25 ~ co2_mean + Developed_Open_Space_perc + Developed_Low_Intensity_perc +
#'    Developed_Medium_Intensity_perc + Developed_High_Intensity_perc +
#'    Open_Water_perc + Woody_Wetlands_perc + Emergent_Herbaceous_Wetlands_perc +
#'    Deciduous_Forest_perc + Evergreen_Forest_perc + Mixed_Forest_perc +
#'    Shrub_perc + Grassland_perc + Pasture_perc + Cultivated_Crops_perc +
#'    pop_density + summer_tmmx + winter_tmmx + summer_rmax + winter_rmax
#'
#' pdata <- plm::pdata.frame(TransAirPolCalif, index = c("GEOID", "year"),
#'                           stringsAsFactors = FALSE)
#' moran.plm.model <- plm::plm(formula = formula.GWPR, data = pdata, model = "within")
#' summary(moran.plm.model)
#'
#' #precomputed bandwidth
#' bw.AIC.Fix <- 2.010529
#'
#' # moran's I test
#' GWPR.moran.test(moran.plm.model, SDF = California, bw = bw.AIC.Fix, kernel = "bisquare",
#'                  adaptive = FALSE, p = 2, longlat = FALSE, alternative = "greater")
GWPR.moran.test <- function(plm_model, SDF, bw, adaptive = FALSE, p = 2, kernel = "bisquare",
                             longlat = FALSE, alternative = "greater")
{
  if(!inherits(plm_model, "plm"))
  {
    stop("This test only accepts the object \"plm\".")
  }
  if(!plm::pdim(plm_model$model)$balanced)
  {
    stop("Current version only accepts balanced panel regression.")
  }
  if(!methods::is(SDF, "Spatial"))
  {
    stop("SDF should be the spatial data frame based on \"sp\".")
  }
  if(!(colnames(plm::index(plm_model$model))[1] %in% colnames(SDF@data)))
  {
    stop("Indexes in plm and SDP are not consistent.")
  }

  plm.resid <- as.data.frame(as.matrix(plm_model$residuals))
  n <- nrow(plm.resid)
  Ti <- ncol(plm.resid)

  id_col <- colnames(plm::index(plm_model$model))[1]
  resid_id <- rownames(plm.resid)
  if (is.null(resid_id))
  {
    stop("The residuals matrix must have row names for individuals.")
  }
  resid_id <- as.character(resid_id)
  sdf_id <- as.character(SDF@data[[id_col]])
  match_idx <- match(resid_id, sdf_id)
  if (anyNA(match_idx))
  {
    stop("Indexes in plm and SDP are not consistent.")
  }
  dp.locat <- as.matrix(sp::coordinates(SDF))[match_idx, , drop = FALSE]
  dMat <- GWmodel::gw.dist(dp.locat = dp.locat, rp.locat = dp.locat,
                           focus = 0, p = p, longlat=longlat)
  if(adaptive)
  {
    bw <- bw * Ti
  } # if adaptive bandwidth, the input number is the numbers of individuals rather than records.
  weight <- GWmodel::gw.weight(dMat, bw=bw, kernel=kernel, adaptive=adaptive)
  diag(weight) <- 0
  rs <- rowSums(weight)
  if(any(rs == 0))
  {
    stop("Some individuals have no neighbours (row sum of weights is zero); please increase bw.")
  }
  weight <- weight/rs

  resid_mat <- as.matrix(plm.resid)
  weighted_resid <- weight %*% resid_mat
  sum.weight.residuals <- colSums(resid_mat * weighted_resid)
  sum.residuals <- colSums(resid_mat^2)
  I.vector <- sum.weight.residuals / sum.residuals
  I.mean <- mean(I.vector)

  # V2 of average I
  V2.upper <-  n*sum(weight^2) + 3 * (sum(weight))^2 - n* sum((colSums(weight))^2)
  V2.lower <- Ti * (n^2 - 1) * (sum(weight))^2
  V2 <- V2.upper/V2.lower

  # Expected value of I
  E <- -1/(n - 1)

  ZI <- (I.mean - E) / sqrt(V2)

  if (alternative == "two.sided")
  {
    PrI <- 2 * stats::pnorm(abs(ZI), lower.tail=FALSE)
  }
  else
  {
    if (alternative == "greater")
    {
      PrI <- stats::pnorm(ZI, lower.tail=FALSE)
    }
    else
    {
      PrI <- stats::pnorm(ZI)
    }
  }

  res <- list(statistic = ZI, p.value=PrI, Estimated.I = I.mean, Excepted.I = E, V2 = V2,
              alternative=alternative)
  return(res)
}
