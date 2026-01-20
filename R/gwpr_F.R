#' Geographically Weighted Panel Regression Based on the Optimal Fixed Distance Bandwidth
#'
#' @param bw                The optimal adaptive bandwidth
#' @param data              The data.frame has been washed
#' @param SDF               Spatial*DataFrame on which is based the data, with the "ID" in the index
#' @param ID_list           The data.frame with individuals' ID
#' @param formula           The regression formula: : Y ~ X1 + ... + Xk
#' @param p                 The power of the Minkowski distance, default is 2, i.e. the Euclidean distance (see GWmodel::bw.gwr)
#' @param longlat           If TRUE, great circle distances will be calculated
#' @param adaptive          If TRUE, adaptive distance bandwidth is used, otherwise, fixed distance bandwidth.
#' @param model             Panel model transformation: (c("within", "random", "pooling"))
#' @param index             A vector for the indexes : (c("ID", "Time"))
#' @param kernel            bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise (default);
#'                           gaussian: wgt = exp(-.5*(vdist/bw)^2);
#'                           exponential: wgt = exp(-vdist/bw);
#'                           tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise;
#'                           boxcar: wgt=1 if dist < bw, wgt=0 otherwise
#' @param effect            The effects introduced in the model, one of "individual" (default) , "time", "twoways", or "nested"
#' @param random.method     Method of estimation for the variance components in the random effects model, one of "swar" (default), "amemiya", "walhus", or "nerlove"
#' @param huge_data_size    If TRUE, the "progress_bar" function will be launched
#'
#' @import dplyr
#' @import GWmodel
#' @importFrom  plm plm pdata.frame r.squared
#' @importFrom lmtest coeftest
#' @import sp
#'
#' @return A list of result:
#' \describe{
#' \item{GW.arguments}{a list class object including the model fitting parameters for generating the report file}
#' \item{R2}{global r2}
#' \item{index}{the index used in the result, Note: in order to avoid mistakes, we forced a rename of the individuals'ID as id.}
#' \item{plm.result}{an object of class inheriting from plm, see plm}
#' \item{raw.data}{the data.frame used in the regression}
#' \item{GWPR.residuals}{the data.frame includes Y, Y hat, and residuals from GWPR}
#' \item{SDF}{a Spatial*DataFrame (either Points or Polygons, see sp) integrated with fit.points,GWPR coefficient estimates,coefficient standard errors and t-values in its "data" slot.}
#' }
#'
#' @references Fotheringham, A. Stewart, Chris Brunsdon, and Martin Charlton. Geographically weighted regression: the analysis of spatially varying relationships. John Wiley & Sons, 2003.
#' @noRd
gwpr_F <- function(bw = bw, data, SDF, ID_list,
                   formula = formula, p = p, longlat = longlat, adaptive = FALSE,
                   model = model, index = index, kernel = kernel, effect = effect,
                   random.method = random.method, huge_data_size = FALSE)
{
  GW.arguments <- list(formula = formula, individual.number = nrow(ID_list), bw = bw,
                       kernel = kernel, adaptive = adaptive, p = p, longlat = longlat)
  message("************************ GWPR Begin *************************\n",
          "Formula: ", paste(as.character(formula)[2], " = ", as.character(formula)[3]), " -- Individuals: ", nrow(ID_list), "\n",
          "Bandwidth: ", bw, " ---- ", "Adaptive: ", adaptive, "\n",
          "Model: ", model, " ---- ", "Effect: ", effect, "\n")
  global_plm_data <- plm::pdata.frame(data, index = index, drop.index = FALSE, row.names = FALSE,
                                      stringsAsFactors = FALSE)
  global_plm <- plm::plm(formula=formula, model=model, data = global_plm_data,
                         effect = effect, index=index, random.method = random.method)
  ID_list_single <- as.vector(ID_list[[1]])
  varibale_name_in_equation <- all.vars(formula)
  coef_names <- colnames(stats::model.matrix(formula, data = data))
  if (model == "within")
  {
    coef_names <- coef_names[coef_names != "(Intercept)"]
  }
  varibale_name_in_equation_out <- gsub("^\\(Intercept\\)$", "Intercept", coef_names)
  coef_count <- length(varibale_name_in_equation_out)
  failed_ids <- c()
  output_rows <- vector("list", length(ID_list_single))
  resid_rows <- vector("list", length(ID_list_single))
  loop_times <- 1
  for (i in seq_along(ID_list_single))
  {
    ID_individual <- ID_list_single[i]
    data$aim[data$id == ID_individual] <- 1
    data$aim[data$id != ID_individual] <- 0
    subsample <- data
    subsample <- subsample[order(-subsample$aim),]
    dp_locat_subsample <- dplyr::select(subsample, dplyr::all_of(c("X", "Y")))
    dp_locat_subsample <- as.matrix(dp_locat_subsample)
    dMat <- GWmodel::gw.dist(dp.locat = dp_locat_subsample, rp.locat = dp_locat_subsample,
                             focus = 1, p=p, longlat=longlat)
    weight <- GWmodel::gw.weight(as.numeric(dMat), bw=bw, kernel=kernel, adaptive=adaptive)
    subsample$wgt <- as.vector(weight)
    subsample <- subsample[(subsample$wgt > 0),]
    Psubsample <- plm::pdata.frame(subsample, index = index, drop.index = FALSE, row.names = FALSE,
                                   stringsAsFactors = FALSE)
    wgt <- Psubsample$wgt
    plm_subsample <- tryCatch(
      plm::plm(formula=formula, model=model, data=Psubsample,
               effect = effect, index=index, weights = wgt,
               random.method = random.method),
      error = function(e) e
    )
    if(!inherits(plm_subsample, "error"))
    {
      coefMat <- lmtest::coeftest(plm_subsample)
      rn_orig <- rownames(coefMat)
      rn_norm <- gsub("^\\(Intercept\\)$", "Intercept", rn_orig)
      coef_hat <- rep(NA_real_, coef_count); names(coef_hat) <- varibale_name_in_equation_out
      coef_se  <- rep(NA_real_, coef_count); names(coef_se)  <- varibale_name_in_equation_out
      coef_t   <- rep(NA_real_, coef_count); names(coef_t)   <- varibale_name_in_equation_out
      idx <- match(varibale_name_in_equation_out, rn_norm, nomatch = 0)
      fill_pos <- which(idx > 0)
      if (length(fill_pos) > 0)
      {
        src <- idx[fill_pos]
        coef_hat[fill_pos] <- coefMat[src, 1]
        coef_se[fill_pos]  <- coefMat[src, 2]
        coef_t[fill_pos]   <- coefMat[src, 3]
      }
      local_r2 <- plm::r.squared(plm_subsample)
      result_line <- c(ID_individual, coef_hat, coef_se, coef_t, local_r2)
      output_rows[[i]] <- result_line
      dataset_add_resid <- cbind(Psubsample, plm_subsample$residuals)
      dataset_add_resid <- as.data.frame(dataset_add_resid)
      dataset_add_resid <- dplyr::select(dataset_add_resid, dplyr::all_of(index), dplyr::all_of(varibale_name_in_equation)[1],
                                         "plm_subsample$residuals")
      colnames(dataset_add_resid) <- c(index, "y", "resid")
      dataset_add_resid$yhat <- dataset_add_resid$y - dataset_add_resid$resid
      dataset_add_resid <- dplyr::filter(dataset_add_resid, id == ID_individual)
      resid_rows[[i]] <- dataset_add_resid
    }
    else
    {
      result_line <- c(ID_individual, rep(NA_real_, coef_count * 3 + 1))
      output_rows[[i]] <- result_line
      dataset_add_resid <- dplyr::select(Psubsample, dplyr::all_of(index),
                                         dplyr::all_of(varibale_name_in_equation)[1])
      dataset_add_resid <- as.data.frame(dataset_add_resid)
      colnames(dataset_add_resid) <- c(index, "y")
      dataset_add_resid$resid <- NA_real_
      dataset_add_resid$yhat <- NA_real_
      dataset_add_resid <- dplyr::filter(dataset_add_resid, id == ID_individual)
      resid_rows[[i]] <- dataset_add_resid
      failed_ids <- c(failed_ids, ID_individual)
    }
    if (huge_data_size)
    {
      progress_bar(loop_times = loop_times, nrow(ID_list))
      loop_times <- loop_times + 1
    }
  }
  output_result <- as.data.frame(do.call(rbind, output_rows))
  y_yhat_resid <- do.call(rbind, resid_rows)
  colnames(output_result) <- c("id", varibale_name_in_equation_out, paste0(varibale_name_in_equation_out,"_SE"),
                               paste0(varibale_name_in_equation_out,"_TVa"), "Local_R2")
  SDF <- sp::merge(SDF, output_result, by = "id")
  y_yhat_resid[,1] <- as.numeric(as.character(y_yhat_resid[,1]))
  y_yhat_resid[,2] <- as.numeric(as.character(y_yhat_resid[,2]))
  if(length(failed_ids) > 0)
  {
    warning("Some local regressions failed; coefficients and residuals contain NA.", call. = FALSE)
  }
  if(anyNA(y_yhat_resid$resid) || anyNA(y_yhat_resid$y))
  {
    r2 <- NA_real_
  }
  else
  {
    r2 <- 1 - sum(y_yhat_resid$resid^2)/(sum((y_yhat_resid$y - mean(y_yhat_resid$y))^2))
  }
  result_list <- list(GW.arguments = GW.arguments, R2 = r2, index = index, plm.result = global_plm,
                      raw.data = data, GWPR.residuals = y_yhat_resid, SDF = SDF)
  return(result_list)
}
