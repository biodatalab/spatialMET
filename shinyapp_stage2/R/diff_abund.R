##
# @title da_tests: Differential abundance of spatial metabolomics data
#
# @param x a matrix of transformed intensity values per pixel (rows=m/z values and
# columns=pixels)
# @param annot_df a data frame with cluster assignments for each pixel (for example, the
# output of `spatial_clust`). The first column contains the pixel names matching the
# columns in `x`
# @param pval_adj Method to adjust p-values. Defaults to `FDR`. Other options as
# available from `p.adjust`
# @param method
# @param cores Number of cores to use in parallelization. If `NULL`, the number of
# cores to use is detected automatically
#
# @return a data frame with results
#
#
da_tests = function(x=NULL, group1=NULL, group2=NULL, group1_lbl=NULL, group2_lbl=NULL, method=NULL, pval_adj=NULL){
  
  #annot_df[[3]] = as.character(annot_df[[3]])
  #clusters_tmp = unique(annot_df[[3]])

  cluster_test_res = data.frame()
  #Loop through clusters
  #for(i in clusters_tmp){
    #annot_df_tmp = annot_df
    #annot_df_tmp[[3]] = ifelse(annot_df_tmp[[3]] == i, i, 'other')
    # Loop through molecules
    test_ls = lapply(1:nrow(x), function(j){
      group1_dat = log1p(as.vector(x[j, group1]))
      group2_dat = log1p(as.vector(x[j, group2]))
        
      # Calculate log-fold change
      avgabund_test1 = mean(group1_dat)
      avgabund_test2 = mean(group2_dat)
      avglogfold_tmp = (avgabund_test1 - avgabund_test2)
      
      # Perform test for differences in mean
      if(method == 'hellinger'){
        res_test = statip::hellinger(group1_dat, group2_dat)
      } else if(method == 'wilcoxon'){
        res_test = wilcox.test(group1_dat, group2_dat)
        statistic_tmp = as.vector(res_test[['statistic']])
        pvalue_tmp = as.vector(res_test[['p.value']])
      } else if(method == 'ttest'){
        res_test = t.test(group1_dat, group2_dat)
        statistic_tmp = as.vector(res_test[['statistic']])
        pvalue_tmp = as.vector(res_test[['p.value']])
      }
    
      res_row = c(
        metabolite = rownames(x)[j],
        cluster = group1_lbl,
        avg_log_itx1 = as.numeric(avgabund_test1),
        avg_log_itx2 = as.numeric(avgabund_test2),
        log_fc = as.numeric(avglogfold_tmp),
        wilcox_est = statistic_tmp,
        p_value = as.numeric(pvalue_tmp)
      )
      
      return(res_row)
    })
    
    # Compile cluster results in a data frame
    cluster_test_res = as.data.frame(do.call(rbind, test_ls))
    cluster_test_res[['adj_p_value']] = p.adjust(cluster_test_res[['p_value']], method=pval_adj)
    cluster_test_res = cluster_test_res %>% 
      dplyr::mutate(dplyr::across(c("avg_log_itx1", "avg_log_itx2", "log_fc", "p_value", "adj_p_value"), ~ as.numeric(.x))) %>% 
      #dplyr::mutate_at(c("avg_log_itx1", "avg_log_itx2", "log_fc", "p_value", "adj_p_value"), as.numeric) %>% 
      dplyr::arrange(adj_p_value, dplyr::desc(log_fc))
    
    # Add cluster results to overall results
    #cluster_test_res = rbind(cluster_test_res, cluster_test_res_tmp)
    
    #rm(cluster_test_res_tmp, test_ls, annot_df_tmp) # Clean env
  #}
  
  return(cluster_test_res)
}

