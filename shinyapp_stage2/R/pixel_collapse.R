##
# @title find_adj_neighbors
# @param coords_df the coordinates data frame (pixel_id, x_coord, y_coord)
#
find_adj_neighbors = function(coords_df) {
  neighbors_list = list()

  x = coords_df[['x_coord']][1]
  y = coords_df[['y_coord']][1]
  while(nrow(coords_df) >= 1) {
    if(any(coords_df[['x_coord']] == x+3 & coords_df[['y_coord']] == y)){
      x = x + 3
      y = y
    } else if(any(coords_df[['x_coord']] == x & coords_df[['y_coord']] == y+3)){
      x = x
      y = y + 3
    } else if(any(coords_df[['x_coord']] == x+3 & coords_df[['y_coord']] == y+1)){
      x = x + 3
      y = y + 1
    } else if(any(coords_df[['x_coord']] == x+1 & coords_df[['y_coord']] == y+3)){
      x = x + 1
      y = y + 3
    } else if(any(coords_df[['x_coord']] == x+3 & coords_df[['y_coord']] == y+2)){
      x = x + 3
      y = y + 2
    } else if(any(coords_df[['x_coord']] == x+2 & coords_df[['y_coord']] == y+3)){
      x = x + 2
      y = y + 3
    } else if(any(coords_df[['x_coord']] == x+3 & coords_df[['y_coord']] == y+3)){
      x = x + 3
      y = y + 3
    } else{
      x = coords_df[['x_coord']][1]
      y = coords_df[['y_coord']][1]
    }

    # Define the 8 possible neighbors (include ninth pixel - center)
    neighbors = data.frame(
      x_coord = c(x, x+1, x+2, x  , x+1, x+2, x  , x+1, x+2),
      y_coord = c(y, y  , y  , y+1, y+1, y+1, y+2, y+2, y+2)
    )

    # Filter out neighbors that are not in the data
    valid_neighbors <- neighbors %>%
      inner_join(coords_df, by = c("x_coord", "y_coord"))

    # Add the valid neighbors to the list
    name_pix = valid_neighbors[[3]][ valid_neighbors[["x_coord"]] == x & valid_neighbors[["y_coord"]] == y ]
    neighbors_list[[name_pix]] = valid_neighbors[[3]]

    coords_df = coords_df[ !(coords_df[[1]] %in% valid_neighbors[[3]]), , drop=F]
  }

  return(neighbors_list)
}
