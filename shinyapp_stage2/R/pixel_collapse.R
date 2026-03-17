# ##
# # @title find_adj_neighbors
# # @param coords_df the coordinates data frame (pixel_id, x_coord, y_coord)
# #
# find_adj_neighbors_OLD = function(coords_df) {
#   neighbors_list = list()
#
#   x = coords_df[['x_coord']][1]
#   y = coords_df[['y_coord']][1]
#   while(nrow(coords_df) >= 1) {
#     if(any(coords_df[['x_coord']] == x+3 & coords_df[['y_coord']] == y)){
#       x = x + 3
#       y = y
#     } else if(any(coords_df[['x_coord']] == x & coords_df[['y_coord']] == y+3)){
#       x = x
#       y = y + 3
#     } else if(any(coords_df[['x_coord']] == x+3 & coords_df[['y_coord']] == y+1)){
#       x = x + 3
#       y = y + 1
#     } else if(any(coords_df[['x_coord']] == x+1 & coords_df[['y_coord']] == y+3)){
#       x = x + 1
#       y = y + 3
#     } else if(any(coords_df[['x_coord']] == x+3 & coords_df[['y_coord']] == y+2)){
#       x = x + 3
#       y = y + 2
#     } else if(any(coords_df[['x_coord']] == x+2 & coords_df[['y_coord']] == y+3)){
#       x = x + 2
#       y = y + 3
#     } else if(any(coords_df[['x_coord']] == x+3 & coords_df[['y_coord']] == y+3)){
#       x = x + 3
#       y = y + 3
#     } else{
#       x = coords_df[['x_coord']][1]
#       y = coords_df[['y_coord']][1]
#     }
#
#     # Define the 8 possible neighbors (include ninth pixel - center)
#     neighbors = data.frame(
#       x_coord = c(x, x+1, x+2, x  , x+1, x+2, x  , x+1, x+2),
#       y_coord = c(y, y  , y  , y+1, y+1, y+1, y+2, y+2, y+2)
#     )
#
#     # Filter out neighbors that are not in the data
#     valid_neighbors <- neighbors %>%
#       inner_join(coords_df, by = c("x_coord", "y_coord"))
#
#     # Add the valid neighbors to the list
#     name_pix = valid_neighbors[[3]][ valid_neighbors[["x_coord"]] == x & valid_neighbors[["y_coord"]] == y ]
#     neighbors_list[[name_pix]] = valid_neighbors[[3]]
#
#     coords_df = coords_df[ !(coords_df[[1]] %in% valid_neighbors[[3]]), , drop=F]
#   }
#
#   return(neighbors_list)
# }

##
#
#' @title find_adj_neighbors_nxn
#' @param coords_df Data frame with columns: pixel_id, x_coord, y_coord
#' @param n Block size (default 3 for 3x3 blocks)
#' @return A named list where each element contains pixel_ids in the same block
#'
find_adj_neighbors = function(coords_df=NULL, n=3){

  # Validate inputs
  if (!all(c("pixel_id", "x_coord", "y_coord") %in% names(coords_df))) {
    stop("coords_df must contain columns: pixel_id, x_coord, y_coord")
  }
  if (n < 1 || n != as.integer(n)) {
    stop("n must be a positive integer")
  }

  # Create working copy
  df = coords_df

  # Assign each pixel to an nxn block using integer division
  df$block_x = floor(df$x_coord / n)
  df$block_y = floor(df$y_coord / n)

  # Create unique block identifier
  df$block_id = paste(df$block_x, df$block_y, sep = "_")

  # Group pixels by block
  neighbors_list = split(df$pixel_id, df$block_id)

  # Name each list element by the first pixel_id in that block
  names(neighbors_list) <- sapply(neighbors_list, function(x) x[1])

  return(neighbors_list)
}


##
# Function to find the most common category (mode) for a vector of pixel_ids
#
#
get_most_common_category = function(pixel_ids=NULL, cat_df=NULL){
  cat_df_tmp = cat_df[cat_df$pixel_id %in% pixel_ids, ]

  # Get categories for the given pixel_ids
  categories = cat_df_tmp[['hc_orig']]

  # Count occurrences of each category
  category_counts = table(categories)

  # Check if all categories are different (all counts equal to 1)
  if (all(category_counts == 1)) {
    # Return the category of the first pixel_id
    return(cat_df_tmp[['hc_orig']][1])
  } else {
    # Return the most common category
    return(names(which.max(category_counts)))
  }
}

