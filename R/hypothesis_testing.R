testVpc <- function( null, alt, beta, ns, Sigma = NULL, alpha = 0.05, num = 1,
                     X = NULL,family = "gaussian", link = "identity",
                     num_cores=1, seed = NULL,...) {

    if(!is.null(seed)) {
      set.seed(seed)
    }

    data <- glmmVpc::batchGLMMData(beta=beta,ns=ns,Sigma=Sigma,num=num,
                                X=X,family=family,link=link,...)

    null_fit <- glmmVpc::batchGLMMFit(null, data, family, num_cores=num_cores)
    alt_fit <- glmmVpc::batchGLMMFit(alt, data, family, num_cores=num_cores)

    logLik_null <- lapply(null_fit, stats::logLik)
    logLik_alt <- lapply(alt_fit, stats::logLik)
    df_null <- unlist(lapply(logLik_null, attr, "df"))
    df_alt <- unlist(lapply(logLik_alt, attr, "df"))

    test_stat <- round(2*(as.numeric(logLik_alt) - as.numeric(logLik_null)),6)

    df_diff <- df_alt - df_null

    p_values <- t(sapply(1:length(df_diff), function(i) {
      result_c <- stats::pchisq(test_stat[i], df = df_diff[i], lower.tail = FALSE)
      result_a <- 0.5*stats::pchisq(test_stat[i], df = df_diff[i],
                                    lower.tail = FALSE) + 0.5*(test_stat[i]==0)
      return(c(result_c, result_a))
    }))

    colnames(p_values) <- c("classical", "adjusted")
    result <- colMeans(p_values < alpha)

    c("num.iterations"=num, result)

}

# beta <- 3
# Sigma <- 2
# theta = 0.3
# ns <- rep(10,10)
# null <- Feature ~ 1
# alt <- Feature ~ 1 + (1|cluster)
# link <- "log"
# family <- "negative_binomial"
# nums <- c(20, 50, 100, 200, 500, 1000, 10000)





# Define the function to run the simulations for different num values
run_simulations <- function(null, alt, beta, ns, Sigma = NULL, alpha = 0.05,
                            nums=1, X = NULL, family = "gaussian", link = "identity",
                            num_cores = 1, seed = NULL, ...) {

  results <- data.frame(num = integer(), power_classical = numeric(), power_adjusted = numeric())

  for (num in nums) {
    result <- testVpc(null = null, alt = alt, beta = beta, ns = ns,
                      Sigma = Sigma, num = num, X = X, family = family,
                      link = link, num_cores = num_cores, seed = seed, ...)

    results <- rbind(results, data.frame(num_iterations = num,
                                         classical = result["classical"],
                                         adjusted = result["adjusted"]))
  }
  rownames(results) <- NULL

  return(results)
}


# simulation_results <- run_simulations(null = null, alt = alt, beta = beta,
#                                       ns = ns, Sigma = Sigma, theta=theta,
#                                       nums = nums,
#                                       family = family, link = link,
#                                       num_cores = 8, seed = 123)
#
# # Print the results
# print(simulation_results)


