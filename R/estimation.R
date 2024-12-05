#' Compare GLMM Fit for True and False Families
#'
#' @param params A matrix or data frame of parameter values for GLMM data
#' generation.
#' @param ns Number of samples or observations per batch.
#' @param gen_family The family distribution used for data generation.
#' @param another_family An alternative family distribution to compare.
#' @param gen_link A link function for generating the data
#' @param fit_formula The model formula for GLMM fitting.
#' @param X A design matrix for the GLMM.
#' @param vpc_input_values A vector for which vpc is calculated
#' @param num_cores Number of cores to fit GLMM with
#' @param fit_transform logical wether a VST model should be fitter and its vpc
#' computed.
#' @param iter Number of iterations for data generation (default is 1).
#' @param seed An optional seed for reproducibility.
#'
#' @return A list containing the coefficients of the true family, false family,
#' and VST-transformed model (if applicable).
#'
#' @export
generatecompareEstimation <- function(params, ns, gen_family,
                              another_family, gen_link, fit_formula,
                              X = NULL, vpc_input_values = NULL,num_cores=1,
                              fit_transform=TRUE,iter = 1, seed = NULL){

  RNGversion("4.3.0")
  RNGkind(kind = "default", normal.kind = "default", sample.kind = "default")

  if (!is.null(seed)) {
    set.seed(seed)
  }

  vpc_true <- as.data.frame(sapply(vpc_input_values, function(x) {
    result <- glmmVpc::vpc_from_paramMat(family=gen_family, paramMat=params, x = x)
    as.numeric(result)
  }))

  colnames(vpc_true) <- paste0("vpc", vpc_input_values)

  datafrmMat <- glmmVpc::batchGLMMDataUsingMatrix(paramMat=params,
                                                   ns=ns, X=X,
                                                   family = gen_family,
                                                   link=gen_link, iter = iter,
                                                   parallel = TRUE)

  result <- compareModels(dataMat=datafrmMat, fit_formula=fit_formula,
                                family1=gen_family, family2=another_family,
                                fit_transform=fit_transform,
                                vpc_input_values = vpc_input_values,
                                num_cores=num_cores)

  result$true <- vpc_true
  return(result)
}


#' Compare GLMM Fit for Groups when handled independently
#'
#' @param params A matrix or data frame of parameter values for GLMM data
#' generation.
#' @param ns Number of samples or observations per batch.
#' @param family The family distribution used for data generation.
#' @param link A link function for generating the data
#' @param formula The model formula for GLMM fitting.
#' @param X A design matrix for the GLMM.
#' @param vpc_input_values A vector for which vpc is calculated
#' @param num_cores Number of cores to fit GLMM with
#' @param iter Number of iterations for data generation (default is 1).
#' @param seed An optional seed for reproducibility.
#'
#' @return A list containing the coefficients of the true vpcs, vpcs with
#' single mixed effect model, and vpcs for independent mixed models.
#'
#' @export
generatecompareGroups <- function(params, ns, family, link, formula,
                                      X = NULL, vpc_input_values = NULL,
                                      num_cores=1,iter = 1, seed = NULL){

  RNGversion("4.3.0")
  RNGkind(kind = "default", normal.kind = "default", sample.kind = "default")

  if (!is.null(seed)) {
    set.seed(seed)
  }

  vpc_true <- as.data.frame(sapply(vpc_input_values, function(x) {
    result <- glmmVpc::vpc_from_paramMat(family=family, paramMat=params, x = x)
    as.numeric(result)
  }))

  colnames(vpc_true) <- paste0("vpc", vpc_input_values)

  datafrmMat <- glmmVpc::batchGLMMDataUsingMatrix(paramMat=params,
                                                  ns=ns, X=X,
                                                  family = family,
                                                  link=link, iter = iter,
                                                  parallel = TRUE)

  fits <- glmmVpc::batchGLMMFit(formula = Feature ~ X + (X|cluster),
                                dataMat = datafrmMat, family = family,
                                num_cores = num_cores)

  vpc_mixed <- as.data.frame(sapply(vpc_input_values, function(x) {
    result <- glmmVpc::vpc(model_fit = fits, x = x)
    result <- sapply(result, function(res) res$vpc)
    result
  }))
  colnames(vpc_mixed) <- paste0("vpc", vpc_input_values)

  result <- compareGroups(datafrmMat, family, num_cores=num_cores)
  result$vpc_mixed <- vpc_mixed
  result$vpc_true <- vpc_true
  class(result) <- "vpcestgr"
  return(result)
}


compareModels <- function(dataMat, fit_formula,
                              family1, family2,
                              fit_transform=TRUE,
                              vpc_input_values,
                              num_cores=1, save_family2 = TRUE) {

  fit_family1 <- glmmVpc::batchGLMMFit(formula=fit_formula,
                                    dataMat=dataMat,
                                    family = family1,
                                    num_cores=num_cores)

  vpc_family1 <- as.data.frame(sapply(vpc_input_values, function(x) {
    result <- glmmVpc::vpc(model_fit = fit_family1, x = x)
    result <- sapply(result, function(res) res$vpc)
    result
  }))


  fit_family2 <- glmmVpc::batchGLMMFit(formula=fit_formula,
                                    dataMat=dataMat,
                                    family = family2,
                                    num_cores=num_cores)

  if(save_family2) {
    coeff <- stats::coef(fit_family2)
    utils::write.csv(coeff, file = paste0(family2,".csv"), row.names=F)
  }

  vpc_family2 <- as.data.frame(sapply(vpc_input_values, function(x) {
    result <- glmmVpc::vpc(model_fit = fit_family2, x = x)
    result <- sapply(result, function(res) res$vpc)
    result
  }))
  col_names <- paste0("vpc", vpc_input_values)
  colnames(vpc_family1) <- colnames(vpc_family2)  <- col_names
  result <- list()
  result[[family1]] <- vpc_family1
  result[[family2]] <- vpc_family2

  if(fit_transform) {
    ys <- glmmVpc::get_y(dataMat)
    X <- glmmVpc::get_x(dataMat)
    cluster <- glmmVpc::get_cluster(dataMat)
    vst_ys <- transform_counts(ys, glmmVpc::vstransform, num_cores)
    dataMat <- glmmVpc::makeDataMatrix(X,vst_ys,cluster)
    fit_vst <- glmmVpc::batchGLMMFit(formula=fit_formula,
                                       dataMat=dataMat,
                                       family = "gaussian",
                                       num_cores=num_cores)
    vpc_vst <- as.data.frame(sapply(vpc_input_values, function(x) {
      result <- glmmVpc::vpc(model_fit = fit_vst, x = x)
      result <- sapply(result, function(res) res$vpc)
      result
    }))
    colnames(vpc_vst) <- col_names
    result$gaussian <- vpc_vst
  }
  class(result) <- "vpcestmo"
  return(result)
}


compareGroups <- function(dataMat, family, num_cores=1) {
  data <- unclass(dataMat)
  dataMat0 <- data[, data["X",]==0][-1,]
  dataMat1 <- data[, data["X",]==1][-1,]
  attr(dataMat0, "num_covriate") <- attr(dataMat1, "num_covriate") <- 0
  attr(dataMat0, "num_feat") <- attr(dataMat1, "num_feat") <- attr(dataMat, "num_feat")
  class(dataMat0) <- class(dataMat1) <- class(dataMat)
  fit0 <- glmmVpc::batchGLMMFit(formula= Feature ~ 1 + (1|cluster),
                        dataMat=dataMat0,
                        family = family,
                        num_cores=num_cores)
  fit1 <- glmmVpc::batchGLMMFit(formula= Feature ~ 1 + (1|cluster),
                                dataMat=dataMat1,
                                family = family,
                                num_cores=num_cores)

  vpc0 <- as.numeric(glmmVpc::vpc(model_fit = fit0))
  vpc1 <- as.numeric(glmmVpc::vpc(model_fit = fit1))
  return(list("vpc0"=vpc0, "vpc1"=vpc1))
}


#' Compare Estimated VPC with True VPC
#'
#' This function compares the true Variance Partition Coefficient (VPC) values
#' calculated from parameter matrices with the estimated VPC values from GLMM
#' model fits.
#'
#' @param params A parameter matrix used for generating the data.
#'   Columns may include parameters like \code{b0}, \code{b1}, \code{sig11},
#'   \code{sig12}, \code{sig22}.
#' @param ns A vector specifying the sample sizes for each group in the data
#' generation.
#' @param X The covariate matrix for fixed effects. Default is \code{X = X}.
#' @param family A character string specifying the family for the model
#' (e.g., "negative_binomial", "tweedie").
#' @param link The link function to be used in the GLMM model.
#' @param formula The formula for fitting the GLMM model (e.g.,
#' \code{Feature ~ X + (X|cluster)}).
#' @param vpc_input_values A vector of input values at which to compute the VPC.
#' Default is \code{NULL}.
#' @param iter Number of iterations for generating data. Default is 1.
#' @param seed An optional random seed for reproducibility. Default is
#' \code{NULL}.
#' @param num_cores The number of CPU cores to use for parallel processing.
#' Default is 1.
#'
#' @return A list with two elements: \code{true} and \code{vpcest}.
#'   - \code{true}: The true VPC values based on the input parameters.
#'   - \code{vpcest}: The estimated VPC values from the fitted GLMM models.
#'   The function also stores the number of iterations as an attribute called
#'   \code{num_iter}.
#'
#' @export
comparePoints <- function(params, ns, X=X, family, link, formula,
                          vpc_input_values=NULL, iter=1, seed = NULL,
                          num_cores=1) {

  RNGversion("4.3.0")
  RNGkind(kind = "default", normal.kind = "default", sample.kind = "default")

  if (!is.null(seed)) {
    set.seed(seed)
  }

  vpc_true <- as.data.frame(sapply(vpc_input_values, function(x) {
    result <- glmmVpc::vpc_from_paramMat(family=family, paramMat=params, x = x)
    as.numeric(result)
  }))

  datafrmMat <- glmmVpc::batchGLMMDataUsingMatrix(paramMat=params,
                                                  ns=ns, X=X,
                                                  family = family,
                                                  link=link, iter = iter,
                                                  parallel = TRUE)
  colnames(vpc_true) <- paste0("vpc", vpc_input_values)

  fits <- glmmVpc::batchGLMMFit(formula = formula,
                                dataMat = datafrmMat, family = family,
                                num_cores = num_cores)

  vpcest <- as.data.frame(sapply(vpc_input_values, function(x) {
    result <- glmmVpc::vpc(model_fit = fits, x = x)
    result <- sapply(result, function(res) res$vpc)
    result
  }))

  colnames(vpcest) <- paste0("vpc", vpc_input_values)

  grouping <- as.integer(sub("Feature", "",sapply(fits,
                      function(fit) colnames(fit$modObj$frame)[1])))
  grouping_index <- factor((grouping - 1) %/% iter + 1)
  split_vpcest <- split(vpcest, grouping_index)
  # split_vpcest <- split(vpcest, rep(1:nrow(params), each = iter))
  result <- list("true"=vpc_true, "vpcest" = split_vpcest)
  class(result) <- "vpcestpt"
  return(result)
}


