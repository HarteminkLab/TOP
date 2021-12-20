#' @title Sort and merge the overlapping sites,
#' remove the overlapping sites with lower occupancy
#'
#' @param predicted_bedgraph.df predicted TF occupancy in bedgraph format
#'
sort_merge_overlap_bedgraph <- function(predicted_bedgraph.df){

  predicted_bedgraph_sorted.df <- predicted_bedgraph.df[with(predicted_bedgraph.df, order(chr, start)), ]
  starts.df <- predicted_bedgraph_sorted.df[-1,]
  ends.df <- predicted_bedgraph_sorted.df[-nrow(predicted_bedgraph_sorted.df),]

  idx_overlap <- which(starts.df$start - ends.df$end <= 0 & starts.df$chr == ends.df$chr) + 1

  if(length(idx_overlap) == 0){
    return(predicted_bedgraph_sorted.df)
  }else{
    # remove the overlapping sites with lower occupancy
    idx_min_occ <- c()
    for(i in idx_overlap){
      if(predicted_bedgraph_sorted.df[i-1, 4] < predicted_bedgraph_sorted.df[i, 4]){
        idx_min_occ <- c(idx_min_occ, i-1)
      }else{
        idx_min_occ <- c(idx_min_occ, i)
      }
    }
    return(sort_merge_overlap_bedgraph(predicted_bedgraph_sorted.df[-idx_min_occ, ]))

  }

}

#' @title Set genome browser track parameters for predictions in bedGraph format
#'
#' @param tf_name TF name
#' @param pwm_id PWM ID
#' @param cell_type cell type
#' @param rep_name replicate name
#' @param type_model model type
#' @param viewMax upper limit to view TF occupancy in genome browser
#' @param mycolor color of the track
#'
track_def_line <- function(tf_name, pwm_id, cell_type, rep_name, type_model,
                           viewMax = 100, mycolor = '0,0,0'){
  track_name <- paste(tf_name, cell_type, rep_name, 'predicted occupancy')
  track_description <- paste(tf_name, pwm_id, cell_type, rep_name, type_model, 'predicted occupancy')
  track_options <- paste0(' visibility=full color=', mycolor,
                          ' altColor=0,0,0 priority=100 autoScale=off alwaysZero=on graphType=bar maxHeightPixels=100:32:16 viewLimits=0:', viewMax)

  track_def <- paste0('track type=bedGraph ', 'name="', track_name, '" description="', track_description,'"', track_options)
  return(track_def)
}
