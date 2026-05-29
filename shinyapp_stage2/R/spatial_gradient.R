##
# @title spatial_gradient
# @description Calculates Spearman's coefficients to detect metabolites showing
# spatial gradients
# @details
# The `spatial_gradient` function fits regression models and calculates Spearman
# coefficients between the metabolite abundances and the minimum or average distance of
# pixels or to a reference tissue domain. In other words the `spatial_gradient` function
# can be used to investigate if a gene is more abundant in pixels closer to
# a specific reference tissue domain, compared to pixels farther from the
# reference domain (or viceversa as indicated by the Spearman's cofficient).
#
# @param x intensity matrix (rows=m/z, columns=pixels)
# @param sp_df coordinates data frame with pixel names in first column matching column
# names in `x`
# @param ref a vector with pixel names assigned to the reference spatial domain
# corresponding to the "reference" cluster or domain. Spearman's correlations will
# be calculated using pixels not included in this list (or domains specified in `exclude`).
# @param exclude optional, a vector with pixel names assigned to domains to be excluded
# from the analysis
# @param out_rm logical (optional), remove outliers defined by the interquartile
# method. This option is only valid when `robust=F`
# @param limit limits the analysis to pixels with distances to `ref` shorter
# than the value specified here. Useful when gradients might occur at smaller scales
# or when the domain in `ref` is scattered through the tissue. Caution must be used
# due to difficult interpretation of imposed limits. It is suggested to run analysis
# without restricted distances in addition for comparison.
# @param distsumm the distance summary metric to use in correlations. One of `min` or `avg`
# @param min_nb the minimum number of immediate neighbors a pixel has to
# have in order to be included in the analysis. This parameter seeks to reduce the
# effect of isolated `ref` spots on the correlation
# @param robust logical, whether to use robust regression (`MASS` and `sfsmisc` packages)
# @param nbs_ls a list with pixel neighborghoods for collapsing
# @param log_dist whethet to log-transform the pixel-pixel distances
# @return a dataframe with the results of the tests
#
#
spatial_gradient = function(x=NULL, sp_df=NULL, ref=NULL, exclude=NULL,
                            out_rm=F, limit=NULL, distsumm='avg', min_nb=9,
                            robust=T, nbs_ls=NULL, log_dist=F){
  # Record time
  zero_t = Sys.time()

  # Make sure the reference cluster is character
  ref = as.character(ref)

  # Define neighborhood tolerance
  nb_dist_thr = c(0.25, 3)

  # Calculate euclidean distances
  coords_tmp = sp_df

  # If list of neighbors is passed, subset coordinates
  if(!is.null(nbs_ls)){
    coords_tmp = coords_tmp[ coords_tmp[[1]] %in% names(nbs_ls), ]
  }
  rm(sp_df) # Clean env

  # Calculate distance matrix
  rownames(coords_tmp) = coords_tmp[[1]]
  coords_tmp = coords_tmp[, -1]
  dist_tmp = as.matrix(stats::dist(coords_tmp, method='euclidean'))

  # Save spots in the different categories (ref, nonref, excl)
  ref_tmp = ref
  nonref_tmp = rownames(dist_tmp)[ !(rownames(dist_tmp) %in% c(ref_tmp, exclude)) ]
  # Subset ref and non_ref spots if collapsing was requested
  if(!is.null(nbs_ls)){
    ref_tmp = ref_tmp[ref_tmp %in% rownames(dist_tmp)]
    nonref_tmp = nonref_tmp[nonref_tmp %in% rownames(dist_tmp)]
  }

  # Identify spots to be removed from reference if not enough neighbors
  # Get minimum distance among all spots within a sample (for Visium would be approximately the same for any sample)
  min_sample = min(as.data.frame(dist_tmp[lower.tri(dist_tmp)]))
#assign('min_sample', min_sample, envir = .GlobalEnv)
  # Get distances among reference spots
  dists_ref_tmp = dist_tmp[ref_tmp, ref_tmp, drop=FALSE]
#assign('dists_ref_tmp', dists_ref_tmp, envir = .GlobalEnv)
  # Get number of neighbors within minimum distance
  # NOTE: When dealing with other technologies like SMI, will need to be more flexible with
  # minimum distances as not an array of equally distant spots. In this case, allowed a "buffer"
  # of a quarter of the minimum distance
  #nbs = colSums(dists_ref_tmp >= (min_sample * nb_dist_thr[1]) & dists_ref_tmp <= (min_sample * nb_dist_thr[2]) )
  nbs = colSums(dists_ref_tmp >= (min_sample * nb_dist_thr[1]) & dists_ref_tmp <= quantile(dists_ref_tmp, 0.2) )
#assign('nbs', nbs, envir = .GlobalEnv)
  if(sum(nbs >= min_nb) < 1){ # At least 1 cluster of spots to continue with analysis
    nbs_keep = c()
  } else{
    nbs_keep = names(nbs)[nbs >= min_nb] # Save spots to be kept (enough neighbors)
  }
  rm(nbs, dists_ref_tmp) # Clean environment

  # Get summarized distances from the reference for each spot in the non-reference
  # Select spots in analysis (non reference in rows, reference in columns)
  dists_nonref_tmp = as.data.frame(dist_tmp[nonref_tmp, ref_tmp, drop=F])
  # Remove columns corresponding to spots without enough neighbors
  dists_nonref_tmp = dists_nonref_tmp[, colnames(dists_nonref_tmp) %in% nbs_keep, drop=F]

  # Check that distances are available for the comparison
  # Number of rows larger than 1, because cannot compute variable genes with a single non-reference spot
  if(nrow(dists_nonref_tmp) > 1 & ncol(dists_nonref_tmp) > 0){
    if(distsumm == 'avg'){ # Summarize non-reference spots using minimun or mean distance
      dists_summ_tmp = tibble::tibble(pixel_id=rownames(dists_nonref_tmp),
                                      dist2ref=apply(dists_nonref_tmp, 1, mean))
    } else{
      dists_summ_tmp = tibble::tibble(pixel_id=rownames(dists_nonref_tmp),
                                      dist2ref=apply(dists_nonref_tmp, 1, min))
    }
  } else{
    dists_summ_tmp = tibble::tibble()
  }
  rm(dists_nonref_tmp) # Clean environment

  # Remove distances if outside user-specified limit
  if(!is.null(limit) & nrow(dists_summ_tmp) > 1){
    # Get lower and upper distance limits
    # If lower limit is higher than user limit, then set lower limit as upper limit
    if(!all(is.na(dists_summ_tmp[['dist2ref']]))){
      dist2reflower = min(dists_summ_tmp[['dist2ref']], na.rm=T)
      if(dist2reflower > limit){
        dist2refupper = dist2reflower
      } else{
        dist2refupper = limit
      }
      # Make NA the distances outside range
      dists_summ_tmp = dists_summ_tmp %>%
        dplyr::mutate(dist2ref=dplyr::case_when(dist2ref <= dist2refupper ~ as.numeric(dist2ref)))

      rm(dist2reflower, dist2refupper) # Clean environment
    }
  }

  # Get abundance data from metabolites
  dist_cor = tibble::tibble() # Initialize result data frame in case no computations can be done (e.g., sample without non-ref pixels)
  if(nrow(dists_summ_tmp) > 1){
    # Collapse abundance if list of neighbors is provided
    if(!is.null(nbs_ls)){
      nbs_ls_tmp = nbs_ls[names(nbs_ls) %in% nonref_tmp] # Subset neighors to those in non-reference
      raw_cts = x[, unlist(nbs_ls_tmp), drop=F]
      raw_cts = lapply(1:length(nbs_ls_tmp), function(i){
        tmp = raw_cts[, colnames(raw_cts) %in% nbs_ls_tmp[[i]], drop=F]
        tmp = as.data.frame(apply(tmp, 1, mean))
        colnames(tmp) = names(nbs_ls_tmp[i])
        return(tmp)
      })
      raw_cts = do.call(cbind, raw_cts)
    } else{
      raw_cts = x[, nonref_tmp, drop=FALSE]
    }

    # Get spots that have at least 1 distance value
    raw_cts = raw_cts[, dists_summ_tmp[['pixel_id']][ !is.na(dists_summ_tmp[['dist2ref']]) ], drop=F]

    # Number of rows larger than 1
    if(ncol(raw_cts) > 1){
      # Get transformed gene expression data (will be used for the correlations with distance)
      # Matrices will contain only the non-reference spots (as defined by non-NA value in distance)
      abund_df = t(raw_cts) %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var='pixel_id') %>%
        dplyr::left_join(dists_summ_tmp, ., by='pixel_id') %>%
        dplyr::filter(pixel_id %in% colnames(raw_cts)) %>%
        dplyr::inner_join(coords_tmp %>% tibble::rownames_to_column('pixel_id'), ., by='pixel_id') %>%
        tibble::column_to_rownames('pixel_id')

      # Detect gene expression outlier spots for each sample and gene
      if(out_rm & !robust){
        outs_dist2ref = list()
        dfdist2ref = abund_df[!is.na(abund_df[['dist2ref']]), ] %>%
          dplyr::select(-c('y_coord', 'x_coord', 'dist2ref'))

        for(gene in colnames(dfdist2ref)){
          # Calculate gene expression quartiles
          quarts = quantile(dfdist2ref[[gene]], probs=c(0.25, 0.75))
          # Calculate inter-quartile range
          iqr_dist2ref = IQR(dfdist2ref[[gene]])
          # Calculate distribution lower and upper limits
          low_up_limits = c((quarts[1]-1.5*iqr_dist2ref),
                            (quarts[2]+1.5*iqr_dist2ref))
          # Save outliers (barcodes)
          outs_dist2ref[[gene]] = rownames(dfdist2ref)[ dfdist2ref[[gene]] < low_up_limits[1] | dfdist2ref[[gene]] > low_up_limits[2] ]
        }
        rm(list=grep("iqr|quarts|low_up|dfdist2ref", ls(), value=T)) # Clean environment
      }

      # Calculate Spearman's correlations
      # Initialize data frame to store results
      dist_cor = tibble::tibble(molecule=character(),
                                lm_coef=numeric(),
                                lm_pval=numeric(),
                                spearman_r=numeric(),
                                spearman_r_pval=numeric(),
                                pval_comment=character())

      # CORRELATIONS DISTANCE TO REFERENCE CLUSTER
      mols_sample = colnames(abund_df %>% dplyr::select(-c('y_coord', 'x_coord', 'dist2ref')))
      for(m in mols_sample){
        tibble_tmp = tibble::tibble(molecule=character(),
                                    lm_coef=numeric(),
                                    lm_pval=numeric(),
                                    spearman_r=numeric(),
                                    spearman_r_pval=numeric(),
                                    pval_comment=character())

        df_mol = abund_df %>% dplyr::select(dist2ref, !!!m)

        lm_res = list(estimate=NA, estimate_p=NA)
        cor_res = list(estimate=NA, p.value=NA)
        if(out_rm & !robust){ # Regular linear models after removal of outliers
          # Remove outliers
          if(length(outs_dist2ref[[m]]) > 0){
            df_mol_outrm = df_mol %>%
              dplyr::filter( !(rownames(.) %in% outs_dist2ref[[m]]) )
          } else{
            df_mol_outrm = df_mol
          }

          if(nrow(df_mol_outrm) > 1){
            # log-transform distances if selected by user
            if(log_dist){
              df_mol_outrm[['dist2ref']] = log(df_mol_outrm[['dist2ref']] + 1e-200)
            }

            # Run linear model and get summary
            lm_tmp = lm(df_mol_outrm[[m]] ~ df_mol_outrm[['dist2ref']])
            lm_summ_tmp = summary(lm_tmp)[['coefficients']]
            if(nrow(lm_summ_tmp) > 1){ # Test a linear model could be run
              lm_res = list(estimate=lm_summ_tmp[2,1],
                            estimate_p=lm_summ_tmp[2,4])
            }
            # Calculate Spearman correlation
            cor_res = tryCatch({cor.test(df_mol_outrm[['dist2ref']], df_mol_outrm[[m]], method='spearman')}, warning=function(w){return(w)})
            pval_warn = NA_character_
            if(any(class(cor_res) == 'simpleWarning')){ # Let known user if p-value could not be exactly calculated
              if(grepl('standard deviation is zero', cor_res$message)){
                pval_warn = 'zero_st_deviation'
              }
              cor_res = cor.test(df_mol_outrm[['dist2ref']], df_mol_outrm[[m]], method='spearman', exact=F)
            }
          }
        } else {
          if(robust){ # Robust linear models?
            df_mol_range = df_mol
            if(nrow(df_mol_range) > 1){
              pval_warn = NA_character_

              # log-transform distances if selected by user
              if(log_dist){
                df_mol_range[['dist2ref']] = log(df_mol_range[['dist2ref']] + 1e-200)
              }

#if(m == 'mz_279.2223'){assign('df_mol_range', df_mol_range, envir= .GlobalEnv)}

              # Run robust linear model and get summary
              lm_tmp = MASS::rlm(df_mol_range[[m]] ~ df_mol_range[['dist2ref']], maxit=100)
              if(lm_tmp[['converged']] & lm_tmp[['coefficients']][2] != 0){ # Check the model converged and an effect was estimated
                # Run Wald test (MASS::rlm does not provide a p-value)
                lm_test_tmp = sfsmisc::f.robftest(lm_tmp)
                lm_res = list(estimate=summary(lm_tmp)[['coefficients']][2,1],
                              estimate_p=lm_test_tmp[['p.value']])
                # Calculate Spearman correlation
                cor_res = tryCatch({cor.test(df_mol_range[['dist2ref']], df_mol_range[[m]], method='spearman')}, warning=function(w){return(w)})
                if(any(class(cor_res) == 'simpleWarning')){ # Let known user if p-value could not be exactly calculated
                  if(grepl('standard deviation is zero', cor_res$message)){
                    pval_warn = 'zero_st_deviation'
                  } #else if(grepl('Cannot compute exact p-value with ties', cor_res$message)){ ## WARNING REMOVED AS MOST GENES WILL HAVE TIES = NON EXACT P-VAL
                  #pval_warn = 'non_exact_pvalue'
                  #}
                  cor_res = cor.test(df_mol_range[['dist2ref']], df_mol_range[[m]], method='spearman', exact=F)
                }
              } else{
                pval_warn = 'rob_regr_no_convergence'
              }
            }
          } else{ # Regular linear models without outlier removal
            df_mol_range = df_mol
            if(nrow(df_mol_range) > 1){

              # log-transform distances if selected by user
              if(log_dist){
                df_mol_range[['dist2ref']] = log(df_mol_range[['dist2ref']] + 1e-200)
              }

              lm_tmp = lm(df_mol_range[[gene]] ~ df_mol_range[['dist2ref']])
              lm_summ_tmp = summary(lm_tmp)[['coefficients']]
              if(nrow(lm_summ_tmp) > 1){ # Test a linear model could be run
                lm_res = list(estimate=lm_summ_tmp[2,1],
                              estimate_p=lm_summ_tmp[2,4])
              }
              # Calculate Spearman correlation
              cor_res = tryCatch({cor.test(df_mol_range[['dist2ref']], df_mol_range[[gene]], method='spearman')}, warning=function(w){return(w)})
              pval_warn = NA_character_
              if(any(class(cor_res) == 'simpleWarning')){ # Let known user if p-value could not be exactly calculated
                if(grepl('standard deviation is zero', cor_res$message)){
                  pval_warn = 'zero_st_deviation'
                }
                cor_res = cor.test(df_mol_range[['dist2ref']], df_mol_range[[m]], method='spearman', exact=F)
              }
            }
          }
        }

        # Create row with results
        tibble_tmp = tibble::tibble(molecule=m,
                                    lm_coef=lm_res[['estimate']],
                                    lm_pval=lm_res[['estimate_p']],
                                    spearman_r=as.vector(cor_res[['estimate']]),
                                    spearman_r_pval=cor_res[['p.value']],
                                    pval_comment=pval_warn)

        rm(list=grep("lm_|_res|_test|cor_|df_mol|exact_p", ls(), value=T)) # Clean environment

        # Add row to result table if there is one row with results
        if(nrow(tibble_tmp) == 1){
          dist_cor = dplyr::bind_rows(dist_cor, tibble_tmp)
          rm(tibble_tmp) # Clean environment
        }
      }
      rm(mols_sample) # Clean environment
    }
  }

    if(nrow(dist_cor) > 0){
      # Adjust p-values for multiple comparison
      dist_cor[['spearman_r_pval_adj']] = p.adjust(dist_cor[['spearman_r_pval']], method='BH')
      dist_cor = dist_cor %>%
        dplyr::relocate(spearman_r_pval_adj, .after=spearman_r_pval) %>%
        dplyr::arrange(spearman_r_pval_adj)
    }

  # Print time
  verbose = 1L
  end_t = difftime(Sys.time(), zero_t, units='min')
  if(verbose > 0L){
    cat(paste0('STgradient completed in ', round(end_t, 2), ' min.\n'))
  }

  return(dist_cor)
}

