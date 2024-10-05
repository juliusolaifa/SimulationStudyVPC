#' Title: Test Variance Partition Coefficient (VPC)
#'
#' @description This function generates data using the specified parameters,
#' fits a GLMM,
#' computes the Variance Partition Coefficient (VPC), and performs a hypothesis
#' test on the VPC estimate.
#'
#' @param null A formula specifying the null hypothesis for the VPC test.
#' @param alt A formula specifying the alternative hypothesis.
#' @param beta A numeric vector of fixed effect coefficients.
#' @param ns A numeric vector specifying the sample sizes for each group.
#' @param Sigma A covariance matrix for the random effects. If NULL, an identity
#' matrix is used.
#' @param alpha A significance level for the hypothesis test (default is 0.05).
#' @param num Number of datasets to simulate (default is 1).
#' @param X Design matrix for fixed effects (optional).
#' @param covariate A vector or for which to compute the VPC.
#' @param family The family of the response variable (default is "gaussian").
#' @param link The link function to use (default is "identity").
#' @param num_cores Number of cores for parallel processing (default is 1).
#' @param type Type of hypothesis test to perform.
#' @param seed Seed for random number generation (optional).
#' @param ... Additional arguments passed to data generation or fitting functions.
#'
#' @return A list containing:
#' \item{vpc_estimate}{Estimated Variance Partition Coefficient.}
#' \item{hypothesis_test}{Results of the hypothesis test for the VPC.}
#'
#' @export
testVpc <- function(null, alt, beta, ns, Sigma = NULL, alpha = 0.05, num = 1,
                    X = NULL, covariate = NULL, family = "gaussian",
                    link = "identity", num_cores = 1, type,
                    seed = NULL, ...) {

  if (!is.null(seed)) {
    set.seed(seed)
  }
  data <- glmmVpc::batchGLMMData(beta = beta, ns = ns, Sigma = Sigma, num = num,
                                 X = X, family = family, link = link, ...)
  fits <- glmmVpc::batchGLMMFit(formula = alt, dataMat = data, family = family)
  vpc.est <- glmmVpc::vpc(fits, x = covariate)
  hyp.test <- glmmVpc::vpc.test(vpc.est, null_formula = null, type = type)
  return(hyp.test)
}
