custom_boxplot_stats <- function(x, coef = 1.5, do.conf = TRUE, do.out = TRUE) {
  nna <- !is.na(x)
  n <- sum(nna)

  # Compute the statistics: min, Q1, mean, Q3, max
  stats <- c(
    min(x, na.rm = TRUE),                     # Minimum
    stats::quantile(x, probs = 0.25, na.rm = TRUE), # First Quartile (Q1)
    mean(x, na.rm = TRUE),                    # Mean
    stats::quantile(x, probs = 0.75, na.rm = TRUE), # Third Quartile (Q3)
    max(x, na.rm = TRUE)                      # Maximum
  )

  # Calculate IQR and outliers
  iqr <- stats[4] - stats[2]

  out <- if (!is.na(iqr)) {
    x < (stats[2] - coef * iqr) | x > (stats[4] + coef * iqr)
  } else !is.finite(x)

  if (any(out[nna], na.rm = TRUE)) {
    stats[c(1, 5)] <- range(x[!out], na.rm = TRUE)
  }

  conf <- if (do.conf) {
    stats[3] + c(-1.58, 1.58) * iqr / sqrt(n)
  }

  list(stats = stats, n = n, conf = conf, out = if (do.out) x[out & nna] else numeric())
}

# Custom boxplot function using custom_boxplot_stats
custom_boxplot <- function(x, ..., range = 1.5, width = NULL, varwidth = FALSE,
                           notch = FALSE, outline = TRUE, names, plot = TRUE,
                           border = graphics::par("fg"), col = "lightgray", log = "",
                           pars = list(boxwex = 0.8, staplewex = 0.5, outwex = 0.5),
                           ann = !add, horizontal = FALSE, add = FALSE, at = NULL)
{
  args <- list(x, ...)
  namedargs <- if (!is.null(attributes(args)$names))
    attributes(args)$names != ""
  else rep_len(FALSE, length(args))
  groups <- if (is.list(x))
    x
  else args[!namedargs]

  if (0L == (n <- length(groups)))
    stop("invalid first argument")

  if (length(class(groups)))
    groups <- unclass(groups)

  if (!missing(names))
    attr(groups, "names") <- names
  else {
    if (is.null(attr(groups, "names")))
      attr(groups, "names") <- 1L:n
    names <- attr(groups, "names")
  }

  # Replace boxplot.stats with custom_boxplot_stats
  for (i in 1L:n) groups[i] <- list(custom_boxplot_stats(unclass(groups[[i]]),
                                                         range))

  stats <- matrix(0, nrow = 5L, ncol = n)
  conf <- matrix(0, nrow = 2L, ncol = n)
  ng <- out <- group <- numeric(0L)
  ct <- 1

  for (i in groups) {
    stats[, ct] <- i$stats
    conf[, ct] <- i$conf
    ng <- c(ng, i$n)
    if ((lo <- length(i$out))) {
      out <- c(out, i$out)
      group <- c(group, rep.int(ct, lo))
    }
    ct <- ct + 1
  }

  z <- list(stats = stats, n = ng, conf = conf, out = out,
            group = group, names = names)

  if (plot) {
    if (is.null(pars$boxfill) && is.null(args$boxfill))
      pars$boxfill <- col
    do.call(graphics::bxp, c(list(z, notch = notch, width = width,
                        varwidth = varwidth, log = log, border = border,
                        pars = pars, outline = outline, horizontal = horizontal,
                        add = add, ann = ann, at = at), args[namedargs]),
            quote = TRUE)

    invisible(z)
  } else z
}

#' Custom Boxplot for Variance Partition Coefficients (VPCs)
#'
#' This function generates boxplots of variance partition coefficients (VPCs) from a given
#' vpcestpt object. It overlays the true VPC values on the boxplots for comparison.
#'
#' @param x A vpcestpt object containing VPC estimates and true VPC values.
#' @param ... Additional arguments to be passed to the custom_boxplot function.
#'
#' @export
customboxplot <- function(x, ...) {
  # Set up the layout to accommodate the plot and legend
  # Adjust bottom margin to make room for the legend
  #par(mar = c(6, 4, 4, 2) + 0.1)  # Default margins plus small adjustment

  max_rows <- max(sapply(x[["vpcest"]], nrow))
  x[["vpcest"]] <- lapply(x[["vpcest"]], function(df) {
    if (nrow(df) < max_rows) {
      rbind(df, matrix(NA, nrow = max_rows - nrow(df), ncol = ncol(df)))
    } else {
      df
    }
  })

  num <- nrow(x[["true"]])
  args <- list(...)
  title <- ifelse(is.null(args$title), "Boxplots of VPCs", args$title)
  custom_boxplot(do.call(cbind, x[["vpcest"]][as.character(1:num)]),
                 cex.axis = 0.8, names = args$names,
                 main = title,
                 ylab = "VPC Values",
                 col = rep(c("lightblue", "pink"), num))

  # Add legend horizontally at the bottom
  graphics::legend("bottomright",
         legend = c("VPC0","VPC1"),
         col = c("lightblue", "pink"),
         pch = 15,
         horiz = F,#TRUE,  # Make legend horizontal
         xpd = F,#TRUE,    # Allow plotting outside the plot region
         #inset = c(0, -0.15) # Adjust position below plot
         )

  for (i in 1:num) {
    graphics::points(2*i - 1, x[["true"]]$vpc0[i], pch = 19, cex = 1.2)
    graphics::points(2*i, x[["true"]]$vpc1[i], pch = 19, cex = 1.2)
  }
}
