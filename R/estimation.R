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

  result <- compareEstimation(dataMat=datafrmMat, fit_formula=fit_formula,
                                family1=gen_family, family2=another_family,
                                fit_transform=fit_transform,
                                vpc_input_values = vpc_input_values,
                                num_cores=num_cores)

  result$true <- vpc_true
  return(result)
}


compareEstimation <- function(dataMat, fit_formula,
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
    as.numeric(result)
  }))

  fit_family2 <- glmmVpc::batchGLMMFit(formula=fit_formula,
                                    dataMat=dataMat,
                                    family = family2,
                                    num_cores=num_cores)
  if(save_family2) {
    coeff <- coef(fit_family2)
    write.csv(coeff, file = paste0(family2,".csv"), row.names=F)
  }

  vpc_family2 <- as.data.frame(sapply(vpc_input_values, function(x) {
    result <- glmmVpc::vpc(model_fit = fit_family2, x = x)
    as.numeric(result)
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
      as.numeric(result)
    }))
    colnames(vpc_vst) <- col_names
    result$gaussian <- vpc_vst
  }
  class(result) <- "vpcest"
  return(result)
}

