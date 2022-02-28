
#' @title Scatter plot of measured and predicted occupancy
#' @description Make a scatter plot of measured and predicted occupancy,
#' with Pearson's correlation (R) between measured and predicted occupancy
#' @param x x-axis values of points in the plot (measured).
#' @param y y-axis values of points in the plot (predicted).
#' @param xlab a label for the x axis.
#' @param ylab a label for the y axis.
#' @param title a title for the plot.
#' @param xlim range of x-axis values.
#' @param ylim range of y-axis values.
#' @param color The plotting color.
#' @return A \code{ggplot} object for
#' the scatter plot of measured and predicted occupancy.
#' @import ggplot2
#' @export
#'
scatterplot_predictions <- function(x, y,
                                    xlab = 'measured',
                                    ylab = 'predicted',
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

#' @title Simple x-y scatter plot
#' @description Make a simple x-y scatter plot
#' @param x x-axis values of points in the plot.
#' @param y y-axis values of points in the plot.
#' @param xlab a label for the x axis.
#' @param ylab a label for the y axis.
#' @param title a title for the plot.
#' @param xlim range of x-axis values.
#' @param ylim range of y-axis values.
#' @param color The plotting color.
#' @return A \code{ggplot} object.
#' @import ggplot2
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
