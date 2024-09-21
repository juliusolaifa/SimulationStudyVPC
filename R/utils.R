#' Transform Count Data
#'
#' This function takes a matrix of count data and applies a transformation
#' function (e.g., variance stabilizing transformation)
#' to it. The transformation is executed with error handling in case the
#' transformation fails.
#'
#' @param ys A matrix or data frame containing the count data to be transformed.
#' @param fun A function to be applied to the count data. This function should
#' accept the count data as its first argument
#'   and optionally additional parameters like `num_cores` for parallel execution.
#' @param num_cores number of PC cores
#'
#' @return The transformed count data as returned by the function \code{fun}.
#' If the transformation fails, the function
#'   returns \code{NULL} and an error message is printed.
#' @export
#'
#' @examples
#' ys <- matrix(rpois(20, lambda = 10), ncol = 5)
#' transform_counts(ys, vst_transform)
transform_counts <- function(ys, fun, num_cores) {
  print("Transforming Counts")
  vst_ys <- tryCatch({
    fun(ys, num_cores = num_cores)  # Adjust num_cores as needed
  }, error = function(e) {
    message("Error in Counts transformation: ", e)
    return(NULL)
  })
  return(vst_ys)
}

#' @export
#' @method plot vpcest
plot.vpcest <- function(x, ...) {
  result_names <- c("negative_binomial", "tweedie", "gaussian")

  result_order <- intersect(names(x), result_names)

  true_result <- x[["true"]]
  true_vpc0 <- true_result[,"vpc0"]
  true_vpc1 <- true_result[,"vpc1"]

  graphics::par(mfrow = c(3, 2), oma = c(0, 0, 3, 0))

  for (result_name in result_order) {
    result <- x[[result_name]]

    result_vpc0 <- result[,"vpc0"]
    result_vpc1 <- result[,"vpc1"]

    plot(true_vpc0, result_vpc0 - true_vpc0, ylim = c(-1, 1), xlim=c(0,1),
         pch = 16, col = "blue", ylab = "Bias", xlab = "True Vpc0",
         main = result_name)
    plot(true_vpc1, result_vpc1 - true_vpc1, ylim = c(-1, 1),  xlim=c(0,1),
         pch = 16, col = "red", ylab = "Bias", xlab = "True Vpc1",
         main = result_name)
  }

  graphics::mtext(paste("Data Generation Family:", result_order[1]),
                  side = 3, line = -1, outer = TRUE, cex = 1)

  graphics::par(mfrow = c(1, 1))
}





