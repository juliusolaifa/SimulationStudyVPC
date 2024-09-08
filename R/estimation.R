#' Compare GLMM Fit for True and False Families
#'
#' @param params A matrix or data frame of parameter values for GLMM data
#' generation.
#' @param covariates A vector or matrix of covariates for the GLMM.
#' @param X A design matrix for the GLMM.
#' @param ns Number of samples or observations per batch.
#' @param true_family The family distribution assumed to be correct.
#' @param false_family An alternative family distribution to compare.
#' @param formula The model formula for GLMM fitting.
#' @param iter Number of iterations for data generation (default is 1).
#' @param seed An optional seed for reproducibility.
#'
#' @return A list containing the coefficients of the true family, false family,
#' and VST-transformed model (if applicable).
#'
#' @export
compareGLMMFit <- function(params, covariates, X, ns, true_family,
                           false_family, formula, iter = 1, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  print(paste("Generating Counts from ", true_family))
  datafrmMat <- vpc::parallelbatchGLMMData(params, X, ns, family = true_family,
                                      iter = iter)

  ys <- datafrmMat[-1, ]
  X <- datafrmMat[1, ]
  group <- colnames(datafrmMat)

  print("Tranforming Counts with VST")
  vst_ys <- tryCatch({
    vpc::vstransform(ys, num_cores=4)
  }, error = function(e) {
    message("Error in VST transformation: ", e)
    return(NULL)
  })

  print(paste("Fitting GLMM to", true_family, "Data Generated from", true_family))
  fit_true <- vpc::batchGLMMFit(formula, ys, X, group, family = true_family,
                           cov_values = covariates)

  print(paste("Fitting GLMM to", false_family, "Data Generated from", true_family))
  fit_false <- vpc::batchGLMMFit(formula, ys, X, group, family = false_family,
                            cov_values = covariates)

  if (!is.null(vst_ys)) {
    print("Fitting GLMM to VST transformed data")
    fit_vst <- vpc::batchGLMMFit(formula, vst_ys, X, group, family = "ga",
                                 cov_values = covariates)
    print("Obtaining coef for VST and its Vpcs")
    coefs_vst <- stats::coef(fit_vst)
  } else {
    coefs_vst <- NULL
  }

  print(paste("Obtaining coef for", true_family, "and it's vpcs"))
  coefs_true <- stats::coef(fit_true)
  print(paste("Obtaining coef for", false_family, "and it's vpcs"))
  coefs_false <- stats::coef(fit_false)

  return(list(coefs_true = coefs_true, coefs_false = coefs_false, coefs_vst = coefs_vst))
}
