

#' Make a scatter plot of measured and predicted occupancy
#'
#' @param x x values (measured)
#' @param y y values (predicted)
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param title title of the figure
#' @param xlim range of x-axis values
#' @param ylim range of y-axis values
#' @param color color of the dots
#' @import ggplot2
#'
#' @export
#'
scatterplot_predictions <- function(x, y,
                                    xlab = 'measured', ylab = 'predicted',
                                    title = '',
                                    xlim = c(0,10),
                                    ylim = c(0,10),
                                    color = "black"){

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package \"ggplot2\" must be installed to use this function.",
      call. = FALSE
    )
  }

  df <- data.frame(x=x, y=y)

  # Make a scatter plot
  p <- ggplot(df, aes(x=x, y=y)) +
    geom_abline(intercept = 0, slope = 1, color="darkgray",
                size = 0.5) +
    geom_point(shape=19,      # Use solid circles
               alpha=0.3,     # opacity
               color = color,
               size = 0.5) +
    scale_x_continuous(breaks=seq(xlim[1],xlim[2],length.out = 5), limits = c(xlim[1],xlim[2])) +
    scale_y_continuous(breaks=seq(ylim[1],ylim[2],length.out = 5), limits = c(ylim[1],ylim[2])) +
    labs(x = xlab, y = ylab, title = title,
         subtitle = paste('R =', round(stats::cor(df[,1], df[,2]),3))) +
    theme_classic()

  return(p)
}

#' Simple scatter plot
#'
#' @param x x values (measured)
#' @param y y values (predicted)
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param title title of the figure
#' @param xlim range of x-axis values
#' @param ylim range of y-axis values
#' @param color color of the dots
#' @import ggplot2
#'
#' @export
#'
scatterplot <- function(x, y,
                        xlab = '',
                        ylab = '',
                        title = '',
                        xlim = NULL,
                        ylim = NULL,
                        color = "black"){

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package \"ggplot2\" must be installed to use this function.",
      call. = FALSE
    )
  }

  df <- data.frame(x=x, y=y)

  # Make a simple scatter plot
  p <- ggplot(df, aes(x=x, y=y)) +
    geom_point(shape=19,      # Use solid circles
               alpha=0.3,     # opacity
               color = color,
               size = 0.5) +
    labs(x = xlab, y = ylab, title = title) +
    theme_classic()

  if(!is.null(xlim)){
    p <- p + scale_x_continuous(breaks=seq(xlim[1],xlim[2],length.out = 5), limits = c(xlim[1],xlim[2]))
  }

  if(!is.null(ylim)){
    p <- p + scale_y_continuous(breaks=seq(ylim[1],ylim[2],length.out = 5), limits = c(ylim[1],ylim[2]))
  }

  return(p)
}
