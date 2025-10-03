#define UPGMA   1
#define WARD2   2
#define CLINK   3
#define WPGMA   4
#define WPGMC2  5
#define UPGMC2  6
#define WARDU   7
#define WPGMCU  8
#define UPGMCU  9
#define SLINK   10

/* merge unclustered sub-clusters together in cluster output (not in tree) */
#define MERGE_UNCLUSTERED 1

/* 0: THG 1994 (CLUSTALW), 1: GSC 1994, 2: experimental, DO NOT USE */
#define WEIGHT_TYPE   1

/* modify GSC weights to floor leaf edges to DBL_MIN,
 * so that multiple duplicate leaves are better down-weighted;
 * yields results similar to if the zeroes were very small values instead
 */
#define GSC_ZERO_HACK 1

/* attempt to counteract deterimental effects of asymmetry on GSC weights */
#define GSC_ASYMMETRY_HACK	1

#define MAX_THREADS 128      /* maximum allowable number of threads */


struct tree_node
{
    struct tree_node  *parent_ptr;
    struct tree_node  *left_ptr;
    struct tree_node  *right_ptr;
    struct tree_node **child_ptrs;
    char   *name;

    double  edge_length;
    double  avg_length;        /* avg distance to leaves */
    double  weight;
    int32_t clevel;
    int32_t num_members;
    int32_t num_children;
#if GSC_ZERO_HACK_currently_not_used
    int32_t num_zero_leaves;
#endif

    int8_t  walk_flag;
    int32_t walk_child_idx;
    
    int32_t cluster_num;       /* within cluster # or unassigned region # */
    char    cluster_good_flag;
    char    cluster_flag;      /* 1: in cluster   0: in unassigned region */
};


struct name_index_pair
{
    char     *name;      /* pointer to name string allocated elsewhere */
    int64_t   index;     /* original index in unsorted array */
};


struct node_data_map
{
    /* map to the leaf pointer in the tree */
    struct  tree_node *node_ptr;
    
    /* map to the row or col index within the data file */
    int32_t data_index;

    /* map to the leaf index within the tree leaf array */
    int32_t leaf_index;

    /* store x,y coordinates parsed from pixel names */
    int32_t x;
    int32_t y;
};


struct options
{
    char linkage_method;
    char ties_fewest_flag;
    char ties_order_flag;
    char ties_random_flag;

    char calc_euclidian_flag;    /* set if we need Euclidian/RMSD */
    char pearson_flag;
    char cosine_flag;
    char ignore_weak_flag;
    char mean_center_flag;
    char unit_variance_flag;
    char maxmag_flag;
    char dump_data_matrix_flag;
    char heuristic_flag;
    char geomean_flag;
    char similarity_flag;
    char absolute_flag;
    char scale_ten_thousand_flag;
    char rmsd_vs_euclidian_flag;
    char cityblock_flag;        /* overrides RMSD / Euclidian */
    char mostly_trig_flag;
    
    char set_any_distance_flag;
    char build_tree_flag;
    char dists_from_file_flag;
    char tree_from_file_flag;
    
    char tree_flip_size_flag;
    char tree_flip_edge_flag;
    char tree_flip_avg_flag;
    
    char tree_weights_flag;
    char tree_defaults_flag;
    char clusters_flag;
    char one_line_flag;
    char clusters_no_merge_flag;    /* don't merge too-small clusters */
    
    char spatial_flag;
    char norm_median_flag;
    char impute_col_min_half_flag;
    char impute_col_global_min_half_flag;
    char impute_global_min_half_flag;
    char impute_row_col_global_min_half_flag;
    char floor_to_one_flag;
    char floor_to_lod_flag;
    char impute_later_flag;
    char log2_prior_flag;
    char unlog2_last_flag;
    
    char transpose_first_flag;
    char transpose_last_flag;
    
    int32_t target_n_clusters;
    int32_t threads;

    double  floor_to_value;
    double  filter_unlog_max_signal_cutoff;
    double  filter_unlog_mean_cutoff;
    double  filter_unlog_sd_cutoff;
    double  filter_log2_mean_cutoff;
    double  filter_log2_sd_cutoff;
    double  filter_present_cutoff;
    double  filter_present_mag;    /* optional, non-zero magnitude */
    double  minkowski;
    double  dist_pow;              /* power to raise distances to */
    double  unit_stretch;          /* 0.99 would be 99th percentile */
    int32_t filter_log2_flag;
    int32_t filter_unlog_flag;
    char    depower_flag;
    
    char   *weight_filename;
};


struct dist_thread_data
{
    pthread_mutex_t mutex_flags;
    pthread_mutex_t mutex_vars;
    pthread_mutex_t mutex_wait;
    pthread_cond_t  cond_wait;
    int     *t_flags_array;
    int      t_index;

    double  **data_matrix;
    double  **dist_matrix;
    double   *weights;
    double   *flip_row1;
    int32_t  *row_good_counts;
    int32_t  *x_coord_array;
    int32_t  *y_coord_array;
    struct    options *opt;
    int32_t   num_cols;
    int32_t   row2_inc;

    double   *row_ptr1;
    int32_t   n1;
    int32_t   row1;
    int32_t   row2_start;

    double    max_dist;         /* return maxdist */
    uint64_t  n_scanned;        /* return n_scanned */
};


struct linkage_thread_data
{
    double  **dist_ptrs;
    double   *new_dists;
    int32_t  *order_array;
    int32_t  *orig_i_array;
    int32_t  *nmemb_array;

    double    best_dist;
    int32_t   method_flag;
    int32_t   best_ci;
    int32_t   best_cj;
    int32_t   best_i;
    int32_t   best_j;

    int32_t   i_start;
    int32_t   i_end;
    int32_t   i_inc;
};


/* tree wrapper function (modified main() from hctree.c) */
extern int32_t launch_tree(double **dist_matrix, char **name_array,
                           int32_t num_leaves,
                           int32_t argc, char **argv,
                           struct options *opt);
extern char * read_in_tree_string(char *tree_file_name);
extern struct tree_node * create_tree_from_string(char *tree_string,
                                  struct tree_node ***return_node_ptr_array,
                                  int32_t *return_num_nodes,
                                  struct tree_node ***return_leaf_ptr_array,
                                  int32_t *return_num_leaves);
extern void bless_tree(struct tree_node *tree_root);
extern void calc_node_clevels(struct tree_node **leaf_ptr_array,
                              int32_t num_leaves);
extern void calc_node_avg_lengths(struct tree_node **node_ptr_array,
                                int32_t num_nodes,
                                struct tree_node **leaf_ptr_array,
                                int32_t num_leaves);
extern int32_t cmp_node_ptr_clevel(const void *keyval, const void *datum);
extern double calc_tree_weights_gsc(struct tree_node **node_ptr_array,
                                    int32_t num_nodes,
                                    struct tree_node **leaf_ptr_array,
                                    int32_t num_leaves, int32_t print_flag);
extern double calc_tree_weights_clustalw(struct tree_node **node_ptr_array,
                                    int32_t num_nodes,
                                    struct tree_node **leaf_ptr_array,
                                    int32_t num_leaves, int32_t print_flag);
extern double calc_tree_weights_experimental(struct tree_node **node_ptr_array,
                                    int32_t num_nodes,
                                    struct tree_node **leaf_ptr_array,
                                    int32_t num_leaves, int32_t print_flag);
extern char * create_newick_string(struct tree_node *root_node,
                            int32_t one_line_per_tree,
                            char overwrite_names_flag);
extern char * create_newick_string_multi(struct tree_node *root_node,
                                  int32_t one_line_per_tree,
                                  char overwrite_names_flag);
extern void free_tree_stuff(char *tree_string,
                     struct tree_node **node_ptr_array, int32_t num_nodes,
                     struct tree_node **leaf_ptr_array, int32_t num_leaves,
                     int32_t free_node_flag);
extern double scan_clusters(struct  tree_node *root_node,
                            struct  tree_node **node_ptr_array,
                            int32_t num_nodes, double min_tiny_fraction,
                            int32_t tiny_size, double max_size_fraction,
                            double  max_dist_fraction,
                            int32_t target_n_clusters,
                            int32_t *return_n_clusters,
                            struct  options *opt,
                            int32_t print_flag);
extern void sum_node_num_members(struct tree_node **node_ptr_array,
                          int32_t num_nodes);
extern void clean_tree(struct tree_node *root_node,
                       struct tree_node ***return_node_ptr_array,
                       int32_t *return_num_nodes,
                       struct tree_node ***return_leaf_ptr_array,
                       int32_t *return_num_leaves);
extern void flip_nodes_multi(struct tree_node **node_ptr_array,
                             int32_t num_nodes, struct options *opt);
extern void depower_tree(struct tree_node **node_ptr_array,
                         int32_t num_nodes, struct options *opt);
extern void output_clusters_multi(struct  tree_node *root_node,
                                  int32_t print_flag, char *filename);

extern double * read_weights_for_data(char     *filename,
                                      double   *weight_array,
                                      char    **data_names,
                                      int32_t   num_points);



/* tree.c, may want to move it to a separate file eventually? */
extern int32_t cmp_double(const void *keyval, const void *datum);
extern int32_t cmp_int32(const void *keyval, const void *datum);
