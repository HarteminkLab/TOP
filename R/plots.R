
#' @title Scatter plot of measured and predicted occupancy
#' @description Makes a scatter plot of measured and predicted occupancy,
#' with Pearson's correlation (R) between measured and predicted occupancy.
#' @param x x-axis values of points in the plot (measured).
#' @param y y-axis values of points in the plot (predicted).
#' @param xlab Label for the x axis.
#' @param ylab Label for the y axis.
#' @param title Title for the plot.
#' @param xlim Range of x-axis values.
#' @param ylim Range of y-axis values.
#' @param color The plotting color.
#' @return A \code{ggplot} object for
#' the scatter plot of measured and predicted occupancy.
#' @import ggplot2
#' @export
#' @examples
#'
#' \dontrun{
#' scatterplot_predictions(x = asinh(chip),
#'                         y = asinh(predicted),
#'                         xlab = 'asinh(measured occupancy)',
#'                         ylab = 'asinh(predicted occupancy)')
#' }
scatterplot_predictions <- function(x, y,
                                    xlab = 'measured',
                                    ylab = 'predicted',
                                    title = '',
                                    xlim = c(0,10),
                                    ylim = c(0,10),
                                    color = 'black'){

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package \"ggplot2\" must be installed to use this function.",
      call. = FALSE
    )
  }

  df <- data.frame(x=x, y=y)

  p <- ggplot(df, aes(x=x, y=y)) +
    geom_abline(intercept = 0, slope = 1, color='darkgray',
                size = 0.5) +
    geom_point(shape=19,      # Use solid circles
               alpha=0.3,     # opacity
               color=color,
               size=0.5) +
    scale_x_continuous(breaks=seq(xlim[1],xlim[2],length.out = 5), limits = c(xlim[1],xlim[2])) +
    scale_y_continuous(breaks=seq(ylim[1],ylim[2],length.out = 5), limits = c(ylim[1],ylim[2])) +
    labs(x = xlab, y = ylab, title = title,
         subtitle = paste('R =', round(stats::cor(df[,1], df[,2]),3))) +
    theme_classic()

  return(p)
}

#' @title Plots DNase or ATAC profiles
#' @param cuts Matrix of DNase or ATAC cuts
#' @param mlen Motif length
#' @param title Title of the plot
#' @export
#' @examples
#' \dontrun{
#' plot_profile(cuts, title='Footprint profiles around binding sites')
#' }
plot_profile <- function(cuts, mlen=ncol(cuts)/2-200, title=''){
  profile <- colMeans(cuts)
  fwd_profile <- profile[1:(length(profile)/2)]
  fwd_profile <- fwd_profile/max(fwd_profile)
  rev_profile <- profile[(length(profile)/2+1):length(profile)]
  rev_profile <- rev_profile/max(rev_profile)
  plot(fwd_profile, type = 'l', col = 'darkblue', xaxt='n',
       xlab = 'Relative position', ylab = 'Normalized cuts', lwd = 2, ylim = c(0,1), main = title)
  lines(rev_profile, type = 'l', col = 'darkred', lwd = 2)
  axis(1, at = c(1, round(length(fwd_profile)/2-mlen/2), round(length(fwd_profile)/2+mlen/2), length(fwd_profile)), labels = c('-100bp', '', '','100bp'))
  legend('topright', legend = c('Forward strand', 'Reverse strand'), lty = 1, col = c('darkblue', 'darkred'), bty = 'n', lwd = 2)
}

#' @title Plots DNase or ATAC profiles by strands of motif matches
#' @param cuts Matrix of DNase or ATAC cuts
#' @param sites Data frame of candidate sites
#' @param mlen Motif length
#' @param title Title of the plot
#' @export
#' @examples
#' \dontrun{
#' plot_profile_strands(cuts, sites, title='Footprint profiles around binding sites')
#' }
plot_profile_strands <- function(cuts, sites, mlen=ncol(cuts)/2-200, title = '', strand = c('both', '+', '-')){
  if(nrow(cuts) != nrow(sites)){
    stop("Number of sites do not match between cuts and sites!")
  }

  strand <- match.arg(strand)
  cuts <- as.matrix(cuts)
  pos_profile <- colMeans(cuts[which(sites$strand == '+'), ])
  pos_profile <- pos_profile/max(pos_profile)
  neg_profile <- colMeans(cuts[which(sites$strand == '-'), ])
  neg_profile <- neg_profile/max(neg_profile)

  if(strand == '+'){
    plot(pos_profile[1:(length(pos_profile)/2)], type = 'l', col = 'darkblue', xaxt='n',
         xlab = 'Relative position', ylab = 'Normalized cuts', lwd = 2, ylim = c(0, 1), main = title)
    lines(pos_profile[(length(pos_profile)/2+1):length(pos_profile)], type = 'l', col = 'darkred', lwd = 2)
    axis(1, at = c(1, round(length(pos_profile)/4-mlen/2), round(length(pos_profile)/4+mlen/2), length(pos_profile)/2), labels = c('-100bp', '', '','100bp'))
    legend('topright',
           legend = c('Fwd profile of + strand sites', 'Rev profile of + strand sites'),
           lty = 1, col = c('darkblue', 'darkred'), bty = 'n', lwd = 2, cex = 0.7)
  }else if (strand == '-'){
    plot(neg_profile[1:(length(neg_profile)/2)], type = 'l', col = 'cyan', xaxt='n',
         xlab = 'position', ylab = 'Normalized cuts', lwd = 2, ylim = c(0, 1), main = title)
    lines(neg_profile[(length(neg_profile)/2+1):length(neg_profile)], type = 'l', col = 'purple', lwd = 2)
    axis(1, at = c(1, round(length(pos_profile)/4-mlen/2), round(length(pos_profile)/4+mlen/2), length(pos_profile)/2), labels = c('-100bp', '', '','100bp'))
    legend('topright',
           legend = c('Fwd profile of - strand sites', 'Rev profile of - strand sites'),
           lty = 1, col = c( 'cyan', 'purple'), bty = 'n', lwd = 2, cex = 0.7)
  }else{
    plot(pos_profile[1:(length(pos_profile)/2)], type = 'l', col = 'darkblue', xaxt='n',
         xlab = 'position', ylab = 'Normalized cuts', lwd = 2, ylim = c(0, 1), main = title)
    lines(pos_profile[(length(pos_profile)/2+1):length(pos_profile)], type = 'l', col = 'darkred', lwd = 2)
    lines(neg_profile[1:(length(neg_profile)/2)], type = 'l', col = 'cyan', lwd = 1)
    lines(neg_profile[(length(neg_profile)/2+1):length(neg_profile)], type = 'l', col = 'purple', lwd = 2)
    axis(1, at = c(1, round(length(pos_profile)/4-mlen/2), round(length(pos_profile)/4+mlen/2), length(pos_profile)/2), labels = c('-100bp', '', '','100bp'))
    legend('topright',
           legend = c('Fwd profile of + strand sites', 'Rev profile of + strand sites', 'Fwd profile of - strand sites', 'Rev profile of - strand sites'),
           lty = 1, col = c('darkblue', 'darkred', 'cyan', 'purple'), bty = 'n', lwd = 2, cex = 0.7)
  }

}
