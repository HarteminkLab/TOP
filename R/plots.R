

#' Make a scatter plot of measured and predicted occupancy
#'
#' @param x x values (measured)
#' @param y y values (predicted)
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param title title of the figure
#' @param transform transformation of the occupancy (asinh, log2 or sqrt)
#' @param xlim range of x-axis values
#' @param ylim range of y-axis values
#' @import ggplot2
#' @import cowplot
#' @export
#'
scatterplot_predictions <- function(x, y,
                           xlab = 'measured', ylab = 'predicted',
                           title = '',
                           transform = c('asinh', 'log2', 'sqrt'),
                           xlim = c(0,10),
                           ylim = c(0,10)){
  transform <- match.arg(transform)

  df <- data.frame(x=x, y=y)

  # Transform the occupancy
  if (transform == 'asinh') {
    df <- asinh(df)
  } else if (transform == 'log2') {
    df <- log2(df + 1)
  } else if (transform == 'sqrt') {
    df <- sqrt(df)
  }
  df[df < 0] <- 0

  # Make a scatter plot

  p <- ggplot(df, aes(x=x, y=y)) +
    geom_abline(intercept = 0, slope = 1, color="darkgray",
                size = 0.5) +
    geom_point(shape=19,      # Use solid circles
               alpha=0.3,     # opacity
               color = "blue",
               size = 0.5) +
    scale_x_continuous(breaks=seq(xlim[1],xlim[2],length.out = 5), limits = c(xlim[1],xlim[2])) +
    scale_y_continuous(breaks=seq(ylim[1],ylim[2],length.out = 5), limits = c(ylim[1],ylim[2])) +
    labs(x = xlab, y = ylab, title = title,
         subtitle = paste('R =', round(stats::cor(df[,1], df[,2]),3))) +
    theme_cowplot()

  return(p)
}
