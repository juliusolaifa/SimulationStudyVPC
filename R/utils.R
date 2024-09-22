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
#' @method plot vpcestmo
plot.vpcestmo <- function(x, ...) {
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
    graphics::abline(h=0)
    plot(true_vpc1, result_vpc1 - true_vpc1, ylim = c(-1, 1),  xlim=c(0,1),
         pch = 16, col = "red", ylab = "Bias", xlab = "True Vpc1",
         main = result_name)
    graphics::abline(h=0)
  }

  graphics::mtext(paste("Data Generation Family:", result_order[1]),
                  side = 3, line = -1, outer = TRUE, cex = 1)

  graphics::par(mfrow = c(1, 1))
}


' Plot Method for vpcestgr Objects
#'
#' This function creates a 2x2 grid of plots for visualizing various aspects of vpcestgr objects.
#'
#' @param x An object of class vpcestgr
#' @param ... Additional arguments passed to plot
#'
#' @return No return value, called for side effects (plotting)
#'
#' @export
#' @method plot vpcestgr
plot.vpcestgr <- function(x, ...) {
  # Set up the plotting area
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  graphics::par(mfrow = c(2, 2), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 4))

  # Common plotting function
  plot_vpc <- function(x0, x1, col, main) {
    plot(x0, x1, col = col, xlab = "Vpc0", ylab = "Vpc1",
         main = main, ylim = c(0, 1), xlim = c(0, 1), ...)
    graphics::abline(0, 1)
  }

  # Individual plots
  plot_vpc(x[["vpc0"]], x[["vpc1"]], "red", "Independent Fitting")
  plot_vpc(x[["vpc_true"]][["vpc0"]], x[["vpc_true"]][["vpc1"]], "blue", "True")
  plot_vpc(x[["vpc_mixed"]][["vpc0"]], x[["vpc_mixed"]][["vpc1"]], "green", "Combined")

  # Comparison plot
  plot_vpc(x[["vpc0"]], x[["vpc1"]], "red", "Comparison")
  graphics::points(x[["vpc_true"]][["vpc0"]], x[["vpc_true"]][["vpc1"]], col = "blue")
  graphics::points(x[["vpc_mixed"]][["vpc0"]], x[["vpc_mixed"]][["vpc1"]], col = "green")

  graphics::mtext("VPC Estimation Comparison", outer = TRUE, cex = 1.5)

  graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  graphics::legend("right", legend = c("Independent", "True", "Combined"),
                   col = c("red", "blue", "green"), pch = 1, bty = "n",
                   cex = 1, xpd = TRUE, inset = c(-0.1, 0))
}


#' @export
#' @method boxplot vpcestpt
#' @importFrom graphics boxplot
boxplot.vpcestpt <- function(x, ...) {
  boxplot(do.call(cbind, x[["vpcest"]][as.character(1:9)]),
          names = rep(c("vpc0", "vpc1"), 9),
          main = "Boxplots of VPCs",
          ylab = "VPC Values",
          col = rep(c("lightblue", "pink"), 9),
          las=2)

  for (i in 1:9) {
    graphics::points(2*i - 1, x[["true"]]$vpc0[i], pch = 19, cex = 1.2)
    graphics::points(2*i, x[["true"]]$vpc1[i], pch = 19, cex = 1.2)
  }

}





