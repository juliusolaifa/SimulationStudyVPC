#' Append VPC Values to Params Data Frame
#'
#' This function calculates Variance Partition Coefficient (VPC) values for
#' each row in the `params` data frame and appends the results as new columns.
#' Each value of `x` generates a corresponding VPC column.
#'
#' @param params A data frame where each row contains parameter values, including
#' the necessary variables like `b0`, `b1`, `sig11`, `sig12`, `sig22`, and `phi`.
#' @param vpcfunc A function that calculates the VPC, taking `beta`, `Sigma`, `phi`,
#' and `x` as arguments.
#' @param x A numeric vector. Each element in `x` will be used to generate
#' a VPC column in the `params` data frame.
#'
#' @return A modified version of the `params` data frame with new columns for
#' each VPC calculated based on the `x` values.
#'
#' @export
append_vpc_to_params <- function(params, vpcfunc, x) {
  for (i in seq_along(x)) {
    vpc_col_name <- paste0("vpc", x[i])

    params[[vpc_col_name]] <- apply(params, 1, function(row) {
      beta <- as.numeric(c(row["b0"], row["b1"]))
      Sigma <- matrix(as.numeric(c(row["sig11"], row["sig12"],
                        row["sig12"], row["sig22"])), 2, 2)
      phi <- as.numeric(row["phi"])
      vpcfunc(beta, Sigma, phi, x[i])
    })

  }
  return(params)
}

#' Generate Parameter Grid with VPC Values
#'
#' This function generates a parameter grid based on inputs of `b0s`, `b1s`,
#' `sigmas`, and `phis`, and then appends Variance Partition Coefficient (VPC)
#' values for each row based on specified values of `x`.
#'
#' @param b0s A vector of values for the intercept term `b0`.
#' @param b1s A vector of values for the slope term `b1`.
#' @param sigmas A list of sigma matrices (with each element containing sig11,
#' sig12, sig22).
#' @param phis A vector of `phi` values for the model.
#'
#' @return A data frame containing the parameter grid with VPC values appended.
#' @export
#'
#' @examples
#' b0s <- c(3, 7)
#' b1s <- c(-5, 5)
#' sigmas <- list(c(2, 1, 2), c(2, -1, 2), c(2, 0, 4))
#' phis <- seq(0.01, 1, 0.05)
#' params_with_vpc <- paramgridWithVPC(b0s, b1s, sigmas, phis)
paramgridWithVPC <- function(b0s, b1s, sigmas, phis) {
  params <- expand.grid(b0=b0s, b1=b1s, sigmas = sigmas, phi=phis)
  params$sig11 <- sapply(params$sigmas, `[`, 1)
  params$sig12 <- sapply(params$sigmas, `[`, 2)
  params$sig22 <- sapply(params$sigmas, `[`, 3)
  params$sigmas <- NULL
  params <- append_vpc_to_params(params, vpc::vpc.nb, x = c(0,1))
  return(params)
}


#' Generate bias plots for multiple models using various vpc columns
#'
#' This function generates bias plots for each model's vpc columns compared to the true values.
#' It allows customization of the color palette using `RColorBrewer`.
#'
#' @param true_params A data frame containing the true values for vpc columns. The columns should be named using the prefix 'vpc', for example, `vpc0`, `vpc1`, etc.
#' @param ... One or more data frames containing the model estimates for vpc columns to be compared against `true_params`.
#' @param palette_name A string specifying the name of the `RColorBrewer` palette to be used. The default is `"Set1"`. Use `RColorBrewer::display.brewer.all()` to see available palettes.
#'
#' @return Generates plots for each model's vpc columns. The function does not return any value.
#' @export
#'
#' @examples
percomPlot <- function(true_params, ... , palette_name = "Set1") {
  args <- list(...)
  # Capture the names of the arguments passed via ...
  modelnames <- sapply(substitute(list(...))[-1], deparse)
  vpcnames <- names(true_params)[grep("vpc", names(true_params))]
  estvpcnames <- unique(unlist(lapply(args, function(df) {
                            grep("vpc", names(df), value = TRUE)})))
  missing_vpcs <- setdiff(vpcnames, estvpcnames)
  if(length(missing_vpcs) > 0) {
    stop("The following vpc columns are not common: ",
         paste(missing_vpcs, collapse = ", "))
  }

  nrows <- length(args)
  ncols <- length(vpcnames)
  graphics::par(mfrow = c(nrows, ncols))

  colors <- RColorBrewer::brewer.pal(max(ncols, 3), palette_name)

  p <- lapply(seq_along(args), function(i) {
    df <- args[[i]]
    modelname <- modelnames[i]

    for (j in seq_along(vpcnames)) {
      vpc <- vpcnames[j]
      true <- true_params[[vpc]]
      bias <- true_params[[vpc]] - df[[vpc]]

      plot(true, bias, ylim = c(-1, 1),
           xlab = paste("True", vpc), ylab = paste("Bias ", vpc),
           main = paste("Model:", modelname),
           col = colors[j], pch = 19)
    }
  })

}


