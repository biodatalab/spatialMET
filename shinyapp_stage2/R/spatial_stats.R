##
# @title spatial_stats: Computes global spatial autocorrelation statistics (Moran's I / Geary's C) on intensity data
#
# @param x a matrix of intensity values
# @param spw a list of spatial weights among neighbors
#
spatial_stats = function(x=NULL, spw=NULL){
  # Record time
  zero_t = Sys.time()
  cat('\tSpatial statistics calculation started...\n')

  morans = gene_moran_i_notest(x=x, listw=spw)
  geary = gene_geary_c_notest(x=x, listw=spw)

  df_tmp = morans %>% dplyr::left_join(., geary, by='feature')

  # Print time
  end_t = difftime(Sys.time(), zero_t, units='min')
  cat(paste0('\tSpatial statistics calculation in ', round(end_t, 2), ' min.\n'))

  return(df_tmp)
}


# Helpers ----------------------------------------------------------------------

##
# @ calculate_spatial_weights
# @param coords a data frame with pixel IDs and coordinates (x, y)
# @param k the number of neighbors to estimate weights. By default NULL, meaning that
# spatial weights will be estimated from Euclidean distances. If an positive integer is
# entered, then the faster k nearest-neighbors approach is used. Please keep in mind
# that estimates are not as accurate as when using the default distance-based method
#
calculate_spatial_weights = function(coords=NULL, k=NULL){
  # Check whether or not a list of weights have been created
  cat("\tCalculating spatial weights...\n") ## Mostly added to make sure calculation is happening only when needed.
  if(!is.null(k)){
    k = as.integer(k)
    cat('\t\tUsing KNN method...\n')
    spw = create_listw_from_knn(coords, ks=k)
  } else{
    cat('\t\tUsing distance-based method...\n')
    spw = create_listw_from_dist(coords)
  }

  return(spw)
}

##
# @title gene_moran_i_notest
# @description Calculates Moran's I from intensity data
#
# @param x intensity data matrix
#
gene_moran_i_notest = function(x=NULL, listw=NULL){
  # Use method to compute autocorrelation described in tutorial of https://rspatial.org/
  df_tmp = data.frame(feature=rownames(x), moran_i=NA)
  for(j in rownames(x)){
    # Extract expression data for a given gene.
    itx_tmp = as.vector(x[j, ])

    # Estimate statistic.
    stat_est = spdep::moran(x=itx_tmp,
                            listw=listw,
                            n=length(listw$neighbours),
                            S0=spdep::Szero(listw))

    df_tmp[['moran_i']][ df_tmp[['feature']] == j ] = stat_est[['I']]
  }

  return(df_tmp)
}

##
# @title gene_geary_c_notest
# @description Calculates Geary's C from intensity data
#
# @param x intensity data matrix
#
gene_geary_c_notest = function(x=NULL, listw=NULL){
  # Use method to compute autocorrelation described in tutorial of https://rspatial.org/
  df_tmp = data.frame(feature=rownames(x), geary_c=NA)
  for(j in rownames(x)){
    # Extract expression data for a given gene.
    itx_tmp = as.vector(x[j, ])

    # Estimate statistic.
    stat_est = spdep::geary(x=itx_tmp,
                            listw=listw,
                            n=length(listw$neighbours),
                            n1=length(listw$neighbours)-1,
                            S0=spdep::Szero(listw))

    df_tmp[['geary_c']][ df_tmp[['feature']] == j ] = stat_est[['C']]
  }

  return(df_tmp)
}

##
# @title create_listw_from_knn
# @param coords coordinates data frame
# @param ks
#
create_listw_from_knn = function(coords=NULL, ks=NULL){
  # Create distance matrix based on the coordinates of each sampled location.
  coords_mtx = coords
  rownames(coords_mtx) = coords_mtx[[1]]
  coords_mtx = coords_mtx[, -1]
  coords_mtx = as.matrix(coords_mtx)
  subj_listw = spdep::nb2listw(spdep::knn2nb(spdep::knearneigh(coords_mtx, k=ks, longlat=F)), style='B')

  return(subj_listw)
}

##
# @title create_listw_from_dist
# @param coords coordinates data frame
#
create_listw_from_dist = function(coords=NULL){
  # Create distance matrix based on the coordinates of each sampled location.
  coords_mtx = coords
  rownames(coords_mtx) = coords_mtx[[1]]
  coords_mtx = coords_mtx[, -1]
  coords_mtx = as.matrix(coords_mtx)
#assign('test_coords_mtx', coords_mtx, envir=.GlobalEnv)
  adj = as.matrix(stats::dist(coords_mtx))
#assign('test_adj', adj, envir=.GlobalEnv)
  adj = adj/max(adj)
  diag(adj) = NA
  colnames(adj) = coords[[1]]
  rownames(adj) = coords[[1]]
#assign('test_adj', adj, envir=.GlobalEnv)
  neighbours = list()
  weights = list()
  for(id in 1:nrow(adj)){
    idx = 1:nrow(adj)
    neighbours[[id]] = idx[!(rownames(adj) %in% rownames(adj)[id])]
    names(neighbours[[id]]) = idx[!(rownames(adj) %in% rownames(adj)[id])]
    weights[[id]] = as.vector(adj[id, ])
    weights[[id]] = 1/(weights[[id]][!(rownames(adj) %in% rownames(adj)[id])])
  }
  class(neighbours) = 'nb'

  subj_listw = list(style='B',
                    neighbours=neighbours,
                    weights=weights)
  class(subj_listw) = c("listw", "nb")

  return(subj_listw)
}

