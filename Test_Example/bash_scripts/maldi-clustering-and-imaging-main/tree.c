#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <float.h>
#include <pthread.h>
#include "xoshiro256.h"
#include "text.h"
#include "hcdist.h"

/* TODO -- singleton clusters are left unflagged by '*' or '**' in tree output,
 *         decide if we want to flag them or not
 */

/* 2025-05-14:  correct typo in GSC_ASYMMETRY_HACK documentation comment
 * 2025-04-28:  enable breaking edge length ties by cluster size via options
 * 2025-04-17:  move WARDU to top of linkage methods tests
 * 2025-04-17:  special case --depower for p = [1, 1.5, 2, 3]
 * 2025-04-16:  add --depower option
 * 2025-04-15:  improved leaf ordering via improved edge-length node flipping
 * 2024-07-10:  support flipping nodes on average length
 * 2024-07-07:  disable some cluster score scaling when target n specified
 * 2024-07-02:  cleaned up some comments
 * 2024-06-24:  reinstate small cluster skipping when initializing scan
 * 2024-06-24:  fix some edge cases for largeish --nclusters
 * 2024-06-20:  changed many output printf %.8g to %.14g
 * 2024-06-20:  support --cluster-merge and --cluster-no-merge
 * 2024-06-19:  extend newer tree flips/names features to jumbles as well
 * 2024-06-19:  default tree building to flip nodes on edge length
 * 2024-06-18:  print initial number of nodes during tree building
 * 2024-06-17:  moved WARD2 to top of the if statement chain
 * 2024-06-07:  added (not so accurate) tree building time estimates
 * 2024-06-05:  tweak auto-cluster scoring to better handle too-small clusters
 * 2024-06-04:  disable auto flipping of nodes by size on tree import
 * 2024-05-30:  change row weight output format for better re-input
 * 2024-05-28:  add --tree-flip-size --tree-flip-dist
 * 2024-05-09:  added basic pthreads multithreading to update new distances
 * 2024-05-03:  pass various flags for MALDI program
 * 2024-04-30:  added *return_n_clusters to cluster functions
 * 2024-04-29:  print EFF_COUNT to stderr as well with --leaf-weights, etc.
 * 2024-04-26:  merge unclustered sub-clusters together in cluster output table
 * 2023-09-07:  better handle weights for highly assymetric branches
 * 2023-08-18:  fixed bug in too-small cluster optimization from 2023-08-03
 * 2023-08-06:  add support for CLUSTALW style weights
 * 2023-08-03:  sub-cluster too-small unclustered regions
 * 2023-08-02:  rename height to clevel (cluster level above leaves)
 * 2023-07-26:  more early back-out during auto-clustering
 * 2023-07-25:  name clusters using letters if <= 26 clusters
 * 2023-07-24:  further adjustments to auto-clustering, all-zeroes clustering
 * 2023-07-19:  added random tie breaking
 * 2023-07-18:  fix for auto-clustering on zero-distance trees
 * 2023-07-14:  change all int to int32_t
 * 2023-07-14:  better integrate hcdist.c and tree.c
 * 2023-07-13:  new auto-clustering and associated improvements
 * 2023-07-13:  fix rare uninitialized read during tree building
 * 2020-12-10:  fix jumble code to work with lower triangular matrix
 * 2020-12-08:  begin converting everything to lower trianglular to save mem
 */


/* #define DEBUG_TREE_INPUT */

#define MEM_OVERHEAD 1.01    /* speed hack -- overallocate to avoid reallocs */

#define STALE_DISTANCE       DBL_MAX
#define STALE_LOWER_BOUNDS  -DBL_MAX

#define DEBUG_COUNTS		0
#define SWAP_ALL		1    /* required for lower trianglular mem */
#define NO_COPY_COLS            1    /* required for lower trianglular mem */

/* CLINK is *incredibly* sensitive to missing data ties
 * SLINK is moderately sensitive to missing data ties
 * everything else is only slightly sensitive to missing data ties
 */

/* --ward2 and --wardu test_rand1.txt currently result in the most updates
 *  per iteration
 */


#define BEST_GLOBAL_DIST_MACRO(tmp_dist) \
if (tmp_dist <= best_dist) \
{ \
    max_order = order_array[i]; \
    delta = max_order - best_order_array[i]; \
\
    nmemb = nmemb_array[i] + \
            nmemb_array[orig_to_cur_array[best_j_array[i]]]; \
\
    /* break ties */ \
    if (tmp_dist == best_dist) \
    { \
        /* break ties randomly */ \
        if (ties_random_flag) \
        { \
            ties_ci_array[num_ties_ci++] = i; \
            continue; \
        } \
\
        /* take merged node with fewer children */ \
        if (nmemb > best_nmemb && \
            ties_fewest_flag) \
        { \
            continue; \
        } \
\
        /* break ties using input order */ \
        if (nmemb == best_nmemb || \
            ties_fewest_flag == 0) \
        { \
            /* take the nearest node when tied */ \
            if (delta > best_delta) \
                continue; \
\
            /* take the earlier node when still tied */ \
            if (delta == best_delta && \
                max_order > best_max_order) \
            { \
                continue; \
            } \
        } \
    } \
    else if (ties_random_flag)\
    { \
        num_ties_ci = 1; \
        ties_ci_array[0] = i; \
    } \
\
    best_nmemb     = nmemb; \
    best_dist      = tmp_dist; \
    best_ci        = i; \
\
    best_max_order = max_order; \
    best_delta     = delta; \
}


int32_t cmp_double(const void *keyval, const void *datum)
{
    double value1, value2;
    
    value1 = * (double *) keyval;
    value2 = * (double *) datum;
    
    if (value1 < value2) return -1;
    if (value1 > value2) return  1;
    
    return 0;
}


int32_t cmp_int32(const void *keyval, const void *datum)
{
    int32_t value1, value2;
    
    value1 = * (int32_t *) keyval;
    value2 = * (int32_t *) datum;
    
    if (value1 < value2) return -1;
    if (value1 > value2) return  1;
    
    return 0;
}


int32_t cmp_node_ptr_clevel(const void *keyval, const void *datum)
{
    struct tree_node *ptr1, *ptr2;
    
    ptr1 = * (struct tree_node **) keyval;
    ptr2 = * (struct tree_node **) datum;
    
    if (ptr1->clevel < ptr2->clevel) return -1;
    if (ptr1->clevel > ptr2->clevel) return  1;

    /* break ties on distance */
    if (ptr1->avg_length < ptr2->avg_length) return -1;
    if (ptr1->avg_length > ptr2->avg_length) return  1;

    /* break ties on number of members */
    if (ptr1->num_members < ptr2->num_members) return -1;
    if (ptr1->num_members > ptr2->num_members) return  1;

    /* shorter edges first */
    if (ptr1->edge_length < ptr2->edge_length) return -1;
    if (ptr1->edge_length > ptr2->edge_length) return  1;
    
    return 0;
}


int32_t cmp_node_ptr_avg_length(const void *keyval, const void *datum)
{
    struct tree_node *ptr1, *ptr2;
    
    ptr1 = * (struct tree_node **) keyval;
    ptr2 = * (struct tree_node **) datum;
    
    if (ptr1->avg_length < ptr2->avg_length) return -1;
    if (ptr1->avg_length > ptr2->avg_length) return 1;
    
    /* break ties on number of members */
    if (ptr1->num_members < ptr2->num_members) return -1;
    if (ptr1->num_members > ptr2->num_members) return 1;

    /* shorter edges first */
    if (ptr1->edge_length < ptr2->edge_length) return -1;
    if (ptr1->edge_length > ptr2->edge_length) return  1;
    
    return 0;
}


int32_t cmp_node_ptr(const void *keyval, const void *datum)
{
    struct tree_node *ptr1, *ptr2;
    int32_t temp_int32;
    
    ptr1 = * (struct tree_node **) keyval;
    ptr2 = * (struct tree_node **) datum;

    /* larger number of members first */
    if (ptr1->num_members > ptr2->num_members) return -1;
    if (ptr1->num_members < ptr2->num_members) return  1;
    
    /* then shorter distances */
    if (ptr1->avg_length < ptr2->avg_length)   return -1;
    if (ptr1->avg_length > ptr2->avg_length)   return  1;

    /* then shorter edges */
    if (ptr1->edge_length < ptr2->edge_length) return -1;
    if (ptr1->edge_length > ptr2->edge_length) return  1;
    
    /* names: as numeric integers (strings of floats are treated as text) */
    if (is_all_digits(ptr1->name) &&
        is_all_digits(ptr2->name))
    {
        if (atol(ptr1->name) < atol(ptr2->name)) return -1;
        if (atol(ptr1->name) > atol(ptr2->name)) return  1;
    }
    /* names: as non-integer strings */
    else if (ptr1->name && ptr2->name)
    {
        temp_int32 = strcmp(ptr1->name, ptr2->name);
        
        if (temp_int32 != 0) return temp_int32;
    }
    
    return 0;
}


/* leave ties as-is */
int32_t cmp_node_ptr_size(const void *keyval, const void *datum)
{
    struct tree_node *ptr1, *ptr2;
    
    ptr1 = * (struct tree_node **) keyval;
    ptr2 = * (struct tree_node **) datum;

    /* larger number of members first */
    if (ptr1->num_members > ptr2->num_members) return -1;
    if (ptr1->num_members < ptr2->num_members) return  1;
    
    return 0;
}


/* leave ties as-is */
int32_t cmp_node_ptr_edge(const void *keyval, const void *datum)
{
    struct tree_node *ptr1, *ptr2;
    
    ptr1 = * (struct tree_node **) keyval;
    ptr2 = * (struct tree_node **) datum;

    /* shorter edges first */
    if (ptr1->edge_length < ptr2->edge_length) return -1;
    if (ptr1->edge_length > ptr2->edge_length) return  1;
    
    return 0;
}


/* leave ties as-is */
int32_t cmp_node_ptr_edge_longer(const void *keyval, const void *datum)
{
    struct tree_node *ptr1, *ptr2;
    
    ptr1 = * (struct tree_node **) keyval;
    ptr2 = * (struct tree_node **) datum;

    /* shorter edges first */
    if (ptr1->edge_length > ptr2->edge_length) return -1;
    if (ptr1->edge_length < ptr2->edge_length) return  1;
    
    return 0;
}


/* free_node_flag: 0 for tree built from scratch, 1 for from file */
void free_tree_stuff(char *tree_string,
                     struct tree_node **node_ptr_array, int32_t num_nodes,
                     struct tree_node **leaf_ptr_array, int32_t num_leaves,
                     int32_t free_node_flag)
{
    struct tree_node *node_ptr;
    int32_t i;

    if (tree_string)
        free(tree_string);

    if (node_ptr_array)
    {
        for (i = 0; i < num_nodes; i++)
        {
            node_ptr = node_ptr_array[i];
        
            if (node_ptr)
            {
                if (node_ptr->name)
                    free(node_ptr->name);

                if (node_ptr->child_ptrs)
                    free(node_ptr->child_ptrs);
        
                /* nodes allocated separately, rather than in array */
                if (free_node_flag)
                    free(node_ptr);
            }
        }

        free(node_ptr_array);
    }
    
    if (leaf_ptr_array)
    {
        /* names/names allocated separately, rather than in array */
        if (free_node_flag)
        {
            for (i = 0; i < num_leaves; i++)
            {
                node_ptr = leaf_ptr_array[i];

                if (node_ptr)
                {
                    if (node_ptr->name)
                        free(node_ptr->name);

                    /* nodes allocated separately, rather than in array */
                    free(node_ptr);
                }
            }
        }

        free(leaf_ptr_array);
    }
}


/*
 * Input matrix format is tab-delimited full distance matrix, but with
 *  labels added as the first column.
 * There can be no blank/extra rows or columns
 *
 * return pointer to distance matrix
 *
 */
double ** read_distance_matrix(char *filename, char ***return_name_array,
                               int32_t *return_num_nodes)
{
    FILE     *infile;
    char     *string         = NULL;
    char    **fields         = NULL;
    int32_t   max_string_len = 0;
    int32_t   num_fields     = 0;
    int32_t   max_num_fields = 0;
    char     *buffer         = NULL;
    int32_t   line_num       = 0;

    double  **dist_matrix    = NULL;
    double   *dptr           = NULL;
    char    **name_array     = NULL;
    int32_t   num_nodes      = 0;
    
    int32_t   i;

    /* allocate my own i/o buffer, since we can't trust the system and/or
     * compiler to allocate a decently large one...
     */
    
    buffer = (char *) malloc(1000000 * sizeof(char));

    /* standard input */
    if (filename == NULL || strcmp(filename, "-") == 0)
    {
        infile = stdin;
    }
    else
    {
        infile = fopen(filename, "rb");
        if (!infile)
        {
            fprintf(stderr, "ERROR -- can't open input file %s\n", filename);

            if (buffer) free(buffer);
            if (string) free(string);
            if (fields) free(fields);
            
            return NULL;
        }
        setvbuf(infile, buffer, _IOFBF, 1000000);
    }

    while(fgets_strip_realloc(&string, &max_string_len, infile))
    {
        num_fields = split_tabs(string, &fields, &max_num_fields);
        
        num_nodes = line_num + 1;

        /* reallocate arrays */
        name_array  =
            (char   **) realloc(name_array,  num_nodes * sizeof(char *));
        dist_matrix =
            (double **) realloc(dist_matrix, num_nodes * sizeof(double *));
            
        name_array[line_num] = malloc((strlen(fields[0]) + 1) * sizeof(char));
        strcpy(name_array[line_num], fields[0]);

        /* need +1 to store the diagonal as well, since we do not
         * skip reading it during the new dist update loop
         * (the extra check would slow the loop down)
         */
        dist_matrix[line_num] =
            (double *) malloc((line_num+1) * sizeof(double));
        
        dptr = dist_matrix[line_num];
        for (i = 1; i <= line_num && i < num_fields; i++)
            *dptr++ = atof(fields[i]);

        /* set diagonal to zero */
        *dptr = 0.0;
        
        line_num++;
    }
    

    fclose(infile);

    if (buffer) free(buffer);
    if (string) free(string);
    if (fields) free(fields);

    *return_name_array = name_array;
    *return_num_nodes  = num_nodes;
    
    return dist_matrix;
}


void *linkage_new_dist_thread(void *passed_ptr)
{
    struct linkage_thread_data *tdata_ptr =
        (struct linkage_thread_data *) passed_ptr;

    double  **dist_ptrs    = tdata_ptr->dist_ptrs;
    double   *new_dists    = tdata_ptr->new_dists;
    int32_t  *order_array  = tdata_ptr->order_array;
    int32_t  *orig_i_array = tdata_ptr->orig_i_array;
    int32_t  *nmemb_array  = tdata_ptr->nmemb_array;

    double    best_dist    = tdata_ptr->best_dist;
    int32_t   method_flag  = tdata_ptr->method_flag;
    int32_t   best_ci      = tdata_ptr->best_ci;
    int32_t   best_cj      = tdata_ptr->best_cj;
    int32_t   best_i       = tdata_ptr->best_i;
    int32_t   best_j       = tdata_ptr->best_j;

    int32_t   i_start      = tdata_ptr->i_start;
    int32_t   i_end        = tdata_ptr->i_end;
    int32_t   i_inc        = tdata_ptr->i_inc;


    /* for new dist calculations */
    double  num_j;
    double  num_i;
    double  num_k;
    double  sum_ij;
    double  dist_ci, dist_cj;
    double  ki2, kj2, ij2;
    double  last_term;
    
    int32_t i;


    /* See Data Clustering in C++: An Object-Oriented Approach
     * by Guojun Gan, p187 for discussion of squared euclidian distances
     * for ward, median, centroid.
     *
     * It is also mentioned, albiet confusingly, in Lance-Williams 1967
     *
     * The Lance-Williams equations for Ward/Centroid/Median are
     * commonly mis-implemented, due to various mistakes throughout
     * the literature over time.  This is due to confusion between using
     * squared vs. unsquared distances.  Lance-Williams 1967 uses
     * unsquared distances for most linkage methods, but squared distances
     * for Ward/Centroid/Median.  This is commonly overlooked, resulting
     * in incorrect equations for unsquared distance inputs, and/or
     * incorrect usage of functions that assume squared inputs.
     *
     * Wikipedia (https://en.wikipedia.org/wiki/Ward%27s_method) is
     * a good example, as they discuss distance d (unsquared), give the
     * Lance-Williams equation for already-squared distances, then
     * never mention anything about requiring input distances to be
     * squared.  R hclust <= v 3.0.3 was similarly incorrect.  I'll
     * want to take a look at the hclust source code to see if they
     * have median and centroid implemented correctly....
     *
     * Many publications give the Lance-Williams formula for squared
     * distances.  We need one for unsquared distances, so that it
     * can be used on the same input distance matrix as the unsquared
     * linkage methods without having to modify the input distance matrix.
     *
     * Correct:   Murtaugh 2011-12-13; Ward2 equation with discussion
     * Correct:   Mullner 2013; has both sqrt and d^2
     * Incorrect: Mullner 2011; close, but forgot to square d
     */
    
    /* sensitive to round-off error asymmetries in input matrix
     *
     * Ward's minimum variance, assumes pre-exponentiated input distances
     *
     * I know, it is confusing that WARDU is squared inputs and WARD2
     * is unsquared inputs, but the U and 2 have to do with whether
     * the distance update calculation is squared or unsquared, not whether
     * the input is squared or unsquared.
     *
     * --wardu --distpow=2 --depower is faster than --ward2, with
     * near-identical distances due to roundoff error at machine precision
     */
    if (method_flag == WARDU)
    {
        num_j  = nmemb_array[best_cj];
        num_i  = nmemb_array[best_ci];
        sum_ij = num_j + num_i;

        for (i = i_start; i <= i_end; i += i_inc)
        {
            num_k   = nmemb_array[i];

#if NO_COPY_COLS
            if (order_array[i] >= order_array[best_cj])
                dist_cj = dist_ptrs[i][best_j];
            else
                dist_cj = dist_ptrs[best_cj][orig_i_array[i]];

            if (order_array[i] >= order_array[best_ci])
                dist_ci = dist_ptrs[i][best_i];
            else
                dist_ci = dist_ptrs[best_ci][orig_i_array[i]];
#else
            dist_cj = dist_ptrs[best_cj][orig_i_array[i]];
            dist_ci = dist_ptrs[best_ci][orig_i_array[i]];
#endif
            new_dists[i] = ((num_j + num_k) * dist_cj +
                            (num_i + num_k) * dist_ci -
                            (num_k        ) * best_dist) /
                           (num_k + sum_ij);
        }
    }
    /* sensitive to round-off error asymmetries in input matrix
     * only truly valid for Euclidian distances, but actually works reasonably
     *  for other distances as well, especially metric distances
     * this is the *ONLY* good linkage method for clustering very large N
     *
     * --wardu --distpow=2 --depower is faster, and will differ only due to
     * machine precision round off error
     */
    /* Ward's minimum variance, unsquared input distances */
    else if (method_flag == WARD2)
    {
        num_j  = nmemb_array[best_cj];
        num_i  = nmemb_array[best_ci];
        sum_ij = num_j + num_i;
        ij2    = best_dist * best_dist;

        for (i = i_start; i <= i_end; i += i_inc)
        {
            num_k = nmemb_array[i];
            
#if NO_COPY_COLS
            if (order_array[i] >= order_array[best_cj])
                kj2 = dist_ptrs[i][best_j];
            else
                kj2 = dist_ptrs[best_cj][orig_i_array[i]];
            kj2 *= kj2;

            if (order_array[i] >= order_array[best_ci])
                ki2 = dist_ptrs[i][best_i];
            else
                ki2 = dist_ptrs[best_ci][orig_i_array[i]];
            ki2 *= ki2;
#else
            kj2  = dist_ptrs[best_cj][orig_i_array[i]];
            kj2 *= kj2;

            ki2  = dist_ptrs[best_ci][orig_i_array[i]];
            ki2 *= ki2;
#endif
            new_dists[i] = sqrt(((num_j + num_k) * kj2 +
                                 (num_i + num_k) * ki2 -
                                 (num_k        ) * ij2) /
                                (num_k + sum_ij));
        }
    }
    /* robust, but tends to not work well for clustering */
    /* Average linkage */
    else if (method_flag == UPGMA)
    {
        num_j  = nmemb_array[best_cj];
        num_i  = nmemb_array[best_ci];
        sum_ij = num_j + num_i;

        for (i = i_start; i <= i_end; i += i_inc)
        {
#if NO_COPY_COLS
            if (order_array[i] >= order_array[best_cj])
                dist_cj = dist_ptrs[i][best_j];
            else
                dist_cj = dist_ptrs[best_cj][orig_i_array[i]];

            if (order_array[i] >= order_array[best_ci])
                dist_ci = dist_ptrs[i][best_i];
            else
                dist_ci = dist_ptrs[best_ci][orig_i_array[i]];
#else
            dist_cj = dist_ptrs[best_cj][orig_i_array[i]];
            dist_ci = dist_ptrs[best_ci][orig_i_array[i]];
#endif

            new_dists[i] = (num_j * dist_cj + num_i * dist_ci) / sum_ij;
        }
    }
    /* sensitive to round-off error asymmetries in input matrix */
    /* Centroid linkage, incorrect update formula */
    else if (method_flag == UPGMCU)
    {
        num_j     = nmemb_array[best_cj];
        num_i     = nmemb_array[best_ci];
        sum_ij    = num_j + num_i;
        last_term = (num_j * num_i * best_dist) /
                    (sum_ij * sum_ij);

        for (i = i_start; i <= i_end; i += i_inc)
        {
#if NO_COPY_COLS
            if (order_array[i] >= order_array[best_cj])
                dist_cj = dist_ptrs[i][best_j];
            else
                dist_cj = dist_ptrs[best_cj][orig_i_array[i]];

            if (order_array[i] >= order_array[best_ci])
                dist_ci = dist_ptrs[i][best_i];
            else
                dist_ci = dist_ptrs[best_ci][orig_i_array[i]];
#else
            dist_cj = dist_ptrs[best_cj][orig_i_array[i]];
            dist_ci = dist_ptrs[best_ci][orig_i_array[i]];
#endif
            new_dists[i] = (num_j * dist_cj + num_i * dist_ci) /
                           sum_ij - last_term;
        }
    }
    /* sensitive to round-off error asymmetries in input matrix */
    /* Centroid linkage, correct update formula */
    else if (method_flag == UPGMC2)
    {
        num_j     = nmemb_array[best_cj];
        num_i     = nmemb_array[best_ci];
        sum_ij    = num_j + num_i;
        ij2       = best_dist * best_dist;
        last_term = (num_j * num_i * ij2) /
                    (sum_ij * sum_ij);

        for (i = i_start; i <= i_end; i += i_inc)
        {
#if NO_COPY_COLS
            if (order_array[i] >= order_array[best_cj])
                kj2 = dist_ptrs[i][best_j];
            else
                kj2 = dist_ptrs[best_cj][orig_i_array[i]];
            kj2 *= kj2;

            if (order_array[i] >= order_array[best_ci])
                ki2 = dist_ptrs[i][best_i];
            else
                ki2 = dist_ptrs[best_ci][orig_i_array[i]];
            ki2 *= ki2;
#else
            kj2  = dist_ptrs[best_cj][orig_i_array[i]];
            kj2 *= kj2;

            ki2  = dist_ptrs[best_ci][orig_i_array[i]];
            ki2 *= ki2;
#endif
            new_dists[i] = sqrt((num_j * kj2 + num_i * ki2) / sum_ij -
                                last_term);
        }
    }
    /* sensitive to round-off error asymmetries in input matrix */
    /* McQuitty's weighted average linkage */
    else if (method_flag == WPGMA)
    {
        for (i = i_start; i <= i_end; i += i_inc)
        {
#if NO_COPY_COLS
            if (order_array[i] >= order_array[best_cj])
                dist_cj = dist_ptrs[i][best_j];
            else
                dist_cj = dist_ptrs[best_cj][orig_i_array[i]];

            if (order_array[i] >= order_array[best_ci])
                dist_ci = dist_ptrs[i][best_i];
            else
                dist_ci = dist_ptrs[best_ci][orig_i_array[i]];
            new_dists[i] = 0.5 * (dist_cj + dist_ci);
#else
            new_dists[i] = 0.5 * (dist_ptrs[best_cj][orig_i_array[i]] +
                                  dist_ptrs[best_ci][orig_i_array[i]]);
#endif
        }
    }
    /* sensitive to round-off error asymmetries in input matrix */
    /* Gower's median linkage, incorrect update formula */
    else if (method_flag == WPGMCU)
    {
        last_term = 0.25 * best_dist;
    
        for (i = i_start; i <= i_end; i += i_inc)
        {
#if NO_COPY_COLS
            if (order_array[i] >= order_array[best_cj])
                dist_cj = dist_ptrs[i][best_j];
            else
                dist_cj = dist_ptrs[best_cj][orig_i_array[i]];

            if (order_array[i] >= order_array[best_ci])
                dist_ci = dist_ptrs[i][best_i];
            else
                dist_ci = dist_ptrs[best_ci][orig_i_array[i]];
            new_dists[i] = 0.5 * (dist_cj + dist_ci) - last_term;
#else
            new_dists[i] = 0.5 * (dist_ptrs[best_cj][orig_i_array[i]] +
                                  dist_ptrs[best_ci][orig_i_array[i]]) -
                            last_term;
#endif
        }
    }
    /* sensitive to round-off error asymmetries in input matrix */
    /* Gower's median linkage, correct update formula */
    else if (method_flag == WPGMC2)
    {
        ij2       = best_dist * best_dist;
        last_term = 0.25 * ij2;
    
        for (i = i_start; i <= i_end; i += i_inc)
        {
#if NO_COPY_COLS
            if (order_array[i] >= order_array[best_cj])
                kj2 = dist_ptrs[i][best_j];
            else
                kj2 = dist_ptrs[best_cj][orig_i_array[i]];
            kj2 *= kj2;

            if (order_array[i] >= order_array[best_ci])
                ki2 = dist_ptrs[i][best_i];
            else
                ki2 = dist_ptrs[best_ci][orig_i_array[i]];
            ki2 *= ki2;
#else
            kj2  = dist_ptrs[best_cj][orig_i_array[i]];
            kj2 *= kj2;

            ki2  = dist_ptrs[best_ci][orig_i_array[i]];
            ki2 *= ki2;
#endif
            new_dists[i] = sqrt(0.5 * (kj2 + ki2) - last_term);
        }
    }
    /* sensitive to round-off error asymmetries in input matrix */
    else if (method_flag == SLINK)
    {
        for (i = i_start; i <= i_end; i += i_inc)
        {
#if NO_COPY_COLS
            if (order_array[i] >= order_array[best_cj])
                dist_cj = dist_ptrs[i][best_j];
            else
                dist_cj = dist_ptrs[best_cj][orig_i_array[i]];

            if (order_array[i] >= order_array[best_ci])
                dist_ci = dist_ptrs[i][best_i];
            else
                dist_ci = dist_ptrs[best_ci][orig_i_array[i]];

            new_dists[i] = dist_cj;
            if (dist_ci < new_dists[i])
                new_dists[i] = dist_ci;
#else
            new_dists[i] = dist_ptrs[best_cj][orig_i_array[i]];
            if (dist_ptrs[best_ci][orig_i_array[i]] < new_dists[i])
                new_dists[i] = dist_ptrs[best_ci][orig_i_array[i]];
#endif
        }
    }
    /* sensitive to round-off error asymmetries in input matrix */
    /* lower triangle matrix is generally slower than full matrix for
     *  CLINK, so use full matrix for CLINK.  The time saved by
     *  skipping updating the matrix columns is less than the extra time
     *  from (probably) lower triangle induced cache misses.
     */
    else if (method_flag == CLINK)
    {
        for (i = i_start; i <= i_end; i += i_inc)
        {
#if NO_COPY_COLS
            if (order_array[i] >= order_array[best_cj])
                dist_cj = dist_ptrs[i][best_j];
            else
                dist_cj = dist_ptrs[best_cj][orig_i_array[i]];

            if (order_array[i] >= order_array[best_ci])
                dist_ci = dist_ptrs[i][best_i];
            else
                dist_ci = dist_ptrs[best_ci][orig_i_array[i]];

            new_dists[i] = dist_cj;
            if (dist_ci > new_dists[i])
                new_dists[i] = dist_ci;
#else
            new_dists[i] = dist_ptrs[best_cj][orig_i_array[i]];
            if (dist_ptrs[best_ci][orig_i_array[i]] > new_dists[i])
                new_dists[i] = dist_ptrs[best_ci][orig_i_array[i]];
#endif
        }
    }

    /* return as normal function, since it wasn't launched as thread */
    if (i_start == 0)
        return(NULL);

    pthread_exit(NULL);
}


/* WARNING -- clobbers distance matrix */
/*
 * When tied distances are encountered, the node with the smallest number
 * of child leaves (recurse through the sub-tree and count all terminal
 * leaves) is chosen.  When still tied, take the nearest input order,
 * and if still tied, take the earlier input order.
 *
 * Selecting for the smallest number of child leaves is desirable for several
 * reasons:
 *
 *    A) For non-ideal distances, the more we merge nodes, the more error we
 *       incur due to violations of assumptions.  Therefore, fewer children
 *       generally implies less error in the distance.
 *
 *    B) Encourages a more evenly balanced tree (left/right branches should
 *       generally be more similar in their number of children).
 *
 *    C) The final tree is less dependent on input order
 *
 * Breaking ties with the nearest order, rather than simply the first order,
 * results in much better worst-case performance on ties with, for example,
 * SLINK.  Taking the first order can result in all or most ties having the
 * same best choice, pushing it towards O(n3).  Taking the nearest makes it
 * less likely that most ties share the same best choice, so that ties wind
 * up closer to the optimal best-case O(n2) run time.  All-zeroes distance
 * matrix makes for a good pathological worst-case test.
 *
 * We also do some book keeping to keep track of the global lower bound
 * distance, as well as the lower bound distance for each row.  We can then
 * postpone updating the best distance for that row once it is invalidated
 * until its best distance could affect the global best distance.
 */
struct tree_node * build_hctree(double  **dist_matrix,
                                char    **name_array,
                                int32_t   num_nodes_orig,
                                int32_t   method_flag,
                                int32_t   njumb,
                                int32_t   ties_fewest_flag,
                                int32_t   ties_random_flag,
                                int32_t   max_threads)
{
    struct tree_node **clust_ptrs = NULL;
    struct tree_node  *clust_pool = NULL;
    struct tree_node  *node_ptr;
    struct tree_node  *left_ptr, *right_ptr;
    double **dist_ptrs = NULL;
    double  *best_dists = NULL;
    double  *lower_bounds = NULL;
    double  *new_dists = NULL;
    double  *dptr;
#if NO_COPY_COLS == 0
    double **dptr2;
#endif
    double   best_dist, tmp_best_dist, tmp_second;
    double   tmp_double;
    int32_t *orig_i_array = NULL;	/* old indices */
    int32_t *best_j_array = NULL;	/* old indices */
    int32_t *orig_to_cur_array = NULL;  /* map orig index to cur index */
    int32_t *order_array = NULL;	/* order we've added nodes in */
    int32_t *best_order_array = NULL;	/* order of best hit */
    int32_t *nmemb_array = NULL;	/* total descendant leaves */
    int32_t *ties_k_array = NULL;       /* store tied indicies */
    int32_t *ties_ci_array = NULL;      /* store tied indicies */
    int32_t  node_order = 0;		/* number of node about to be added */
    int32_t  max_order, best_max_order;
    int32_t  num_nodes = num_nodes_orig;
    int32_t  clust_pool_idx = num_nodes_orig;
    int32_t  nm1, nm2;
    int32_t  best_i, best_j;	/* old indices */
    int32_t  best_ci, best_cj;	/* current indicies */
    int32_t  swap_ci, swap_cj;	/* current indicies */
    int32_t  tmp_best_k;
    int32_t  delta, best_delta;
    int32_t  num_ties_k, num_ties_ci;
    int32_t  nmemb, best_nmemb;
    int32_t  i, j, k;

    /* start and end times for system clock
     * we're using 1 second resolution due to ANSI-compatibility
     * we could use clock_gettime() for high-precision, but that is POSIX
     */
    time_t start_time, end_time;
    double delta_time, time_left;
    uint64_t n_scanned, n_left_to_go;
    char time_str1[20], time_str2[20];    /* printf can't re-eval one str */
    
    pthread_t threads[MAX_THREADS];
    int       pt_ret[MAX_THREADS];
    struct    linkage_thread_data *tdata_array = NULL;
    int32_t   t, cur_threads;

#if DEBUG_COUNTS
    uint64_t counts_nest1 = 0;
    uint64_t counts_nest2 = 0;
    uint64_t counts_nest3 = 0;
    uint64_t counts_nest4 = 0;
#endif

#ifdef GPERF
ProfilerStart("gperf.out");
#endif


    /* too many threads causes too many penalties
     * ~8 threads is optimal on my test system
     */
    if (max_threads > 8)
        max_threads = 8;

    /* sanity check on number of threads */
    if (max_threads > MAX_THREADS)
        max_threads = MAX_THREADS;


    /* NOTE -- orig_i_array[], best_ci, best_cj, best_i, best_j
     *
     * best_j_array[] is relative to the ORIGINAL matrix cooridinates.
     * Thus, we need to keep track of the original locations as we swap
     * rows around (since columns are not swapped)
     *
     * best_ci, best_cj are the current reduced/swapped around row indices
     * best_i, best_j are the indicies in the original un-swapped matrix
     */
    
    best_dists        = (double  *) malloc(num_nodes * sizeof(double));
    lower_bounds      = (double  *) malloc(num_nodes * sizeof(double));
    new_dists         = (double  *) malloc(num_nodes * sizeof(double));
    orig_i_array      = (int32_t *) malloc(num_nodes * sizeof(int32_t));
    best_j_array      = (int32_t *) malloc(num_nodes * sizeof(int32_t));
    orig_to_cur_array = (int32_t *) malloc(num_nodes * sizeof(int32_t));
    order_array       = (int32_t *) malloc(num_nodes * sizeof(int32_t));
    best_order_array  = (int32_t *) malloc(num_nodes * sizeof(int32_t));
    nmemb_array       = (int32_t *) malloc(num_nodes * sizeof(int32_t));
    
    if (ties_random_flag)
    {
        ties_k_array  = (int32_t *) malloc(num_nodes * sizeof(int32_t));
        ties_ci_array = (int32_t *) malloc(num_nodes * sizeof(int32_t));
    }
    
    dist_ptrs = (double **) calloc(num_nodes, sizeof(double *));
    for (i = 0; i < num_nodes; i++)
        dist_ptrs[i] = dist_matrix[i];

    /* allocate initial + future nodes */
    /* must cast 2 as uint32_t, otherwise gcc -O1 gives a false warning */
    clust_pool = (struct tree_node  *) calloc(2U * num_nodes,
                                              sizeof(struct tree_node));
    /* initial leaves */
    clust_ptrs = (struct tree_node **) calloc(num_nodes,
                                              sizeof(struct tree_node *));

    /* initialize leaves */
    for (i = 0; i < num_nodes; i++)
    {
        clust_ptrs[i]         = &clust_pool[i];
        orig_i_array[i]       = i;
        order_array[i]        = i;
        orig_to_cur_array[i]  = i;
        nmemb_array[i]        = 1;

        node_ptr              = clust_ptrs[i];
        node_ptr->name        = name_array[i];
        node_ptr->clevel      = 0;
        node_ptr->num_members = 1;
    }
    node_order = num_nodes - 1;
    
    
    /* initialize best j */
    /* first node is perma-stale due to lower order triangle */
    best_j_array[0] = best_order_array[0] = 0;
    best_dists[0]   = lower_bounds[0]     = STALE_DISTANCE;
    for (i = 1; i < num_nodes; i++)
    {
        /* If I initialize num_ties_k inside the ties_random_flag check
         * where it belongs, gcc thinks that num_ties_k may be unintialized
         * when I use it below, so I have to put it outside the
         * ties_random_flag check.
         */
        num_ties_k = 1;
        if (ties_random_flag)
        {
            ties_k_array[0] = 0;
        }
        best_j        = 0;
        dptr          = dist_ptrs[i];
        tmp_best_dist = dptr[0];
        tmp_second    = STALE_DISTANCE;

        /* only need to scan lower triangle; same half we scan later on */
        for (j = 1; j < i; j++)
        {
            tmp_double = dptr[j];

            /* lower triangle enforces j < i, so no need to check for j!=i */
            if (tmp_double <= tmp_best_dist)
            {
                tmp_second = tmp_best_dist;

                if (ties_random_flag)
                {
                    if (tmp_double < tmp_best_dist)
                    {
                        num_ties_k      = 1;
                        ties_k_array[0] = j;
                    }
                    /* break ties randomly */
                    else
                    {
                        ties_k_array[num_ties_k++] = j;
                        continue;
                    }
                }
                
                /* break ties */
                /* lower triangle means j can't ever be > i,
                 *  and increasing scan means equal dists = better delta,
                 *  so ties are already broken
                 */
                
                tmp_best_dist = tmp_double;
                best_j        = j;
            }
            else if (tmp_double < tmp_second)
                tmp_second = tmp_double;
        }
        
        /* break ties randomly */
        if (ties_random_flag && num_ties_k > 1)
        {
            best_j = ties_k_array[(int32_t) (num_ties_k * dxoshiro256p())];
        }

        best_dists[i]       = tmp_best_dist;
        lower_bounds[i]     = tmp_second;
        best_j_array[i]     = best_j;
        best_order_array[i] = best_j;
    }


    /* initialize best dist */
    /* break ties randomly */
    if (ties_random_flag)
    {
        num_ties_ci      = 1;
        ties_ci_array[0] = 1;
    }
    best_ci = 1;
    best_max_order = order_array[1];
    if (best_order_array[1] > best_max_order)
        best_max_order = best_order_array[1];
    best_delta = order_array[1] - best_order_array[1];
    best_nmemb = nmemb_array[1] +
                 nmemb_array[orig_to_cur_array[best_j_array[1]]];
    best_dist = best_dists[1];

    /* find best dist */
    for (i = 2; i < num_nodes; i++)
    {
        tmp_double = best_dists[i];

        BEST_GLOBAL_DIST_MACRO(tmp_double)
    }
    /* break ties randomly */
    if (ties_random_flag && num_ties_ci > 1)
    {
        best_ci = ties_ci_array[(int32_t) (num_ties_ci * dxoshiro256p())];
    }


    /* allocate and initialize thread data */
    tdata_array = (struct linkage_thread_data *)
        calloc(max_threads, sizeof(struct linkage_thread_data));
    for (t = 0; t < max_threads; t++)
    {
        tdata_array[t].dist_ptrs    = dist_ptrs;
        tdata_array[t].new_dists    = new_dists;
        tdata_array[t].order_array  = order_array;
        tdata_array[t].orig_i_array = orig_i_array;
        tdata_array[t].nmemb_array  = nmemb_array;
        tdata_array[t].method_flag  = method_flag;
        tdata_array[t].i_inc        = max_threads;
    }
    

    /* find nearest pair of nodes */
    fprintf(stderr, "NumNodesLeft: %d %d\n", njumb, num_nodes);
    n_scanned  = 0;
    start_time = time(NULL);
    while (num_nodes > 1)
    {
#if DEBUG_COUNTS
        counts_nest1++;
#endif

#if 1
        if (n_scanned && num_nodes % (max_threads * 100) == 0)
        {
            end_time     = time(NULL);
            delta_time   = difftime(end_time, start_time);
            n_left_to_go = 0.5 * num_nodes * (num_nodes + 1);
            
            /* we know it slows down as it goes, so fudge it */
            n_left_to_go *= 2.0;

            /* accumulate at least 15 seconds before estimating */
            if (delta_time >= 15.0)
            {
                time_left = (double) n_left_to_go * delta_time /
                            (double) n_scanned;

                fprintf(stderr, "NumNodesLeft: %d %d   TimeElapsed: %s   TimeLeft: %s\n",
                    njumb, num_nodes,
                    seconds_to_str(time_str1, 20, delta_time),
                    seconds_to_str(time_str2, 20, time_left));
            }
            else
            {
                fprintf(stderr, "NumNodesLeft: %d %d   TimeElapsed: %s   TimeLeft: %s\n",
                    njumb, num_nodes,
                    seconds_to_str(time_str1, 20, delta_time),
                    "estimating...");
            }
        }
        n_scanned += num_nodes - 1;
#endif

        best_j  = best_j_array[best_ci];
        best_cj = orig_to_cur_array[best_j];

        /* make sure that best_cj is < best_ci
         * various optimizations in the code require this assumption to hold
         */
        if (best_cj > best_ci)
        {
            k       = best_ci;
            best_ci = best_cj;
            best_cj = k;
        }
        best_j = orig_i_array[best_cj];
        best_i = orig_i_array[best_ci];


        /* calculate new distances */

        /* initialize each thread
         *
         * NOTE -- Due to the severe performance penaly of mutexes in this
         *         context, it is faster to create/terminate all threads
         *         each iteration, rather than permanent worker threads that
         *         process new data when ready.  The mutex performance penalty
         *         is worse than the thread creation/termination penalty.
         *
         *         pthread_join() should result in new_dists[] being
         *         fully written by the time it is used afterwards?
         */
        for (t = 1; t < max_threads && t <= num_nodes - 1; t++)
        {
            tdata_array[t].best_dist = best_dist;
            tdata_array[t].best_ci   = best_ci;
            tdata_array[t].best_cj   = best_cj;
            tdata_array[t].best_i    = best_i;
            tdata_array[t].best_j    = best_j;
            
            tdata_array[t].i_start   = t;
            tdata_array[t].i_end     = num_nodes - 1;
        }
        cur_threads = t;

        
        /* launch each thread */
        for (t = 1; t < cur_threads; t++)
        {
            pt_ret[t] = pthread_create(&threads[t],
                            NULL,
                            linkage_new_dist_thread,
                            (void *) &tdata_array[t]);
        }

        
        /* use current thread as a worker thread while others finish */
        t = 0;
        tdata_array[t].best_dist = best_dist;
        tdata_array[t].best_ci   = best_ci;
        tdata_array[t].best_cj   = best_cj;
        tdata_array[t].best_i    = best_i;
        tdata_array[t].best_j    = best_j;
        tdata_array[t].i_start   = t;
        tdata_array[t].i_end     = num_nodes - 1;
        linkage_new_dist_thread((void *) &tdata_array[t]);

        
        /* wait for all to finish before proceeding */
        for (t = 1; t < cur_threads; t++)
            pthread_join(threads[t], NULL);


        nm1 = num_nodes - 1;
        nm2 = num_nodes - 2;

        /* standard swap, preserve both nm1, nm2 */
        /* if we skip swapping swap_ci later on in a few places (!SWAP_FULL),
         *  we must assign swap_cj = nm2 in order to not break,
         *  since, otherwise, swap_ci will sometimes be the one re-used
         *  instead of swap_cj, and partial swapping is then invalid
         */
        if (best_cj < nm2 && best_ci < nm2)
        {
            swap_cj = nm2;
            swap_ci = nm1;
        }
        /* no swap, preserve neither */
        else if (best_cj == nm2)
        {
            swap_cj = best_cj;
            swap_ci = best_ci;
        }
        /* preserve nm1 */
        else if (best_ci == nm2)
        {
            swap_cj = nm1;
            swap_ci = best_ci;
        }
        /* preserve nm2 */
        else
        {
            swap_cj = nm2;
            swap_ci = best_ci;
        }
        
        /* swap out old orig_i */
        k                     = orig_i_array[best_cj];
        orig_i_array[best_cj] = orig_i_array[swap_cj];
        orig_i_array[swap_cj] = k;
#if SWAP_ALL
        k                     = orig_i_array[best_ci];
#endif
        orig_i_array[best_ci] = orig_i_array[swap_ci];
#if SWAP_ALL
        orig_i_array[swap_ci] = k;
#endif

        /* map orig index to cur index */
        orig_to_cur_array[orig_i_array[best_cj]] = best_cj;
        orig_to_cur_array[orig_i_array[best_ci]] = best_ci;
        orig_to_cur_array[orig_i_array[nm2]]     = nm2;

        /* swap out old rows */
        dptr               = dist_ptrs[best_cj];
        dist_ptrs[best_cj] = dist_ptrs[swap_cj];
        dist_ptrs[swap_cj] = dptr;
#if SWAP_ALL
        dptr               = dist_ptrs[best_ci];
#endif
        dist_ptrs[best_ci] = dist_ptrs[swap_ci];
#if SWAP_ALL
        dist_ptrs[swap_ci] = dptr;
#endif
        
        /* Distance columns are not swapped, to save an N loop.
         * Thus, all distance columns are indexed by their original order!
         */

        /* swap out new dists */
        /* no need to actually swap, both swap_cj and swap_ci are clobbered */
        new_dists[best_cj] = new_dists[swap_cj];
        new_dists[best_ci] = new_dists[swap_ci];

        /* zero out new self-self */
        new_dists[nm2] = 0.0;

        /* swap out best_dists */
        /* no need to actually swap, both swap_cj and swap_ci are clobbered */
        best_dists[best_cj] = best_dists[swap_cj];
        best_dists[best_ci] = best_dists[swap_ci];

        /* initialize best_dist of new merged node */
        best_dists[nm2] = STALE_DISTANCE;

        /* swap out lower_bounds */
        /* no need to actually swap, both swap_cj and swap_ci are clobbered */
        lower_bounds[best_cj] = lower_bounds[swap_cj];
        lower_bounds[best_ci] = lower_bounds[swap_ci];

        /* initialize lower_bound of new merged node */
        lower_bounds[nm2] = STALE_LOWER_BOUNDS;

        /* swap out order */
        /* no need to actually swap, both swap_cj and swap_ci are clobbered */
        order_array[best_cj] = order_array[swap_cj];
        order_array[best_ci] = order_array[swap_ci];

        /* assign order of new merged node */
        order_array[nm2] = ++node_order;

        /* swap out best_order */
        /* no need to actually swap, both swap_cj and swap_ci are clobbered */
        best_order_array[best_cj] = best_order_array[swap_cj];
        best_order_array[best_ci] = best_order_array[swap_ci];

        /* swap out nmemb */
        /* no need to actually swap, both swap_cj and swap_ci are clobbered */
        k                    = nmemb_array[best_cj] + nmemb_array[best_ci];
        nmemb_array[best_cj] = nmemb_array[swap_cj];
        nmemb_array[best_ci] = nmemb_array[swap_ci];
        nmemb_array[nm2]     = k;

        /* swap out best_j_array */
        /* no need to actually swap, both swap_cj and swap_ci are clobbered */
        best_j_array[best_cj] = best_j_array[swap_cj];
        best_j_array[best_ci] = best_j_array[swap_ci];

        /* initialize nm2 for updating later */
        best_j_array[nm2]     = best_j;



        /* store new node */
        node_ptr = &clust_pool[clust_pool_idx++];
        node_ptr->left_ptr    = clust_ptrs[best_cj];
        node_ptr->right_ptr   = clust_ptrs[best_ci];
        node_ptr->num_members = nmemb_array[nm2];
        
        /* take lower-sorted name string */

        if (strcmp(node_ptr->left_ptr->name, node_ptr->right_ptr->name) <= 0)
        {
            node_ptr->name = malloc((strlen(node_ptr->left_ptr->name) + 1) *
                                    sizeof(char));
            strcpy(node_ptr->name, node_ptr->left_ptr->name);
        }
        else
        {
            node_ptr->name = malloc((strlen(node_ptr->right_ptr->name) + 1) *
                                    sizeof(char));
            strcpy(node_ptr->name, node_ptr->right_ptr->name);
        }

        left_ptr  = node_ptr->left_ptr;
        right_ptr = node_ptr->right_ptr;

        if (left_ptr->clevel > right_ptr->clevel)
            node_ptr->clevel = left_ptr->clevel + 1;
        else
            node_ptr->clevel = right_ptr->clevel + 1;
        

        /* store edges, etc. */
        
        /* leaves */
        node_ptr->avg_length = 0.5 * best_dist;
        if (left_ptr->num_members == 1)
        {
            left_ptr->avg_length  = 0.0;
            left_ptr->edge_length = node_ptr->avg_length;
        }
        if (right_ptr->num_members == 1)
        {
            right_ptr->avg_length  = 0.0;
            right_ptr->edge_length = node_ptr->avg_length;
        }

        /* nodes */
        if (left_ptr->num_members > 1)
        {
            left_ptr->edge_length =
                node_ptr->avg_length - left_ptr->avg_length;
        }
        if (right_ptr->num_members > 1)
        {
            right_ptr->edge_length =
                node_ptr->avg_length - right_ptr->avg_length;
        }
        

        /* swap nodes */
        /* no need to actually swap, both swap_cj and swap_ci are clobbered */
        clust_ptrs[best_cj] = clust_ptrs[swap_cj];
        clust_ptrs[best_ci] = clust_ptrs[swap_ci];
        
        clust_ptrs[nm2]     = node_ptr;


        /* free dead row */
        if (dist_ptrs[nm1])
        {
            free(dist_ptrs[nm1]);
            dist_matrix[orig_i_array[nm1]] = NULL;
        }

        /* free and allocate new row */
        /* must allocate full row, due to book keeping */
        if (dist_ptrs[nm2])
            free(dist_ptrs[nm2]);
        dist_ptrs[nm2] = (double *) malloc(num_nodes_orig * sizeof(double));

        /* update original matrix pointer */
        dist_matrix[orig_i_array[nm2]] = dist_ptrs[nm2];


#if NO_COPY_COLS == 0
        /* fill in new dist col */
        dptr  = new_dists;
        dptr2 = dist_ptrs;
        do
            (*dptr2++)[orig_i_array[nm2]] = *dptr++;
        while (dptr < new_dists + nm2);
#endif

        /* update newly merged node nm2 using new_dists[],
         *  since that is faster than reading the dist matrix
         */

        /* If I initialize num_ties_k inside the ties_random_flag check
         * where it belongs, gcc thinks that num_ties_k may be unintialized
         * when I use it below in the dxoshiro256p() calculation, so I have to
         * put it outside the ties_random_flag check.
         */
        num_ties_k = 1;
        if (ties_random_flag)
        {
            ties_k_array[0] = 0;
        }
        tmp_best_k    = 0;
        tmp_best_dist = new_dists[0];
        tmp_second    = STALE_DISTANCE;

        for (k = 1; k < nm2; k++)
        {
            if (new_dists[k] <= tmp_best_dist)
            {
                tmp_second = tmp_best_dist;

                /* break ties */
                if (new_dists[k] == tmp_best_dist)
                {
                    /* break ties randomly */
                    if (ties_random_flag)
                    {
                        ties_k_array[num_ties_k++] = k;
                        continue;
                    }
                
                    /* take merged node with fewer children */
                    if (nmemb_array[k] > nmemb_array[tmp_best_k] &&
                        ties_fewest_flag)
                    {
                        continue;
                    }

                    /* break ties using input order */
                    if (nmemb_array[k] == nmemb_array[tmp_best_k] ||
                        ties_fewest_flag == 0)
                    {
                        /* take the nearest node when tied */
                        if (order_array[k] < order_array[tmp_best_k])
                            continue;
                    
                        /* take the earlier node when still tied */
                        /* lower triangle guarantees this already */
                    }
                }
                else if (ties_random_flag)
                {
                    num_ties_k = 1;
                    ties_k_array[0] = k;
                }
                
                tmp_best_dist = new_dists[k];
                tmp_best_k    = k;
            }
            else if (new_dists[k] < tmp_second)
                tmp_second = new_dists[k];
        }

        /* break ties randomly */
        if (ties_random_flag && num_ties_k > 1)
        {
            tmp_best_k =
                ties_k_array[(int32_t) (num_ties_k * dxoshiro256p())];
        }

        best_order_array[nm2] = order_array[tmp_best_k];
        best_j_array[nm2]     = orig_i_array[tmp_best_k];
        lower_bounds[nm2]     = tmp_second;
        best_dists[nm2]       = tmp_best_dist;


        /* initialize new best dist */
        if (ties_random_flag)
        {
            num_ties_ci      = 1;
            ties_ci_array[0] = nm2;
        }
        best_ci        = nm2;
        best_max_order = node_order;
        best_delta     = node_order - order_array[tmp_best_k];
        best_dist      = tmp_best_dist;
        best_nmemb     = nmemb_array[tmp_best_k] + nmemb_array[nm2];


        /* Update distance matrix with newly merged distances
         * Flag newly stale nodes
         * Refine best dist using old best distances
         */
        dptr = dist_ptrs[nm2];
        dptr[orig_i_array[nm2]] = 0.0;   /* initialize self-self dist */
        for (i = 0; i < nm2; i++)
        {
            /* update merged node dists */
            /* it is faster to put it here than do it in other loops */
            /* note that the diagonal will NOT be zero, but doesn't matter */
            dptr[orig_i_array[i]] = new_dists[i];

            /* only use non-stale distances */
            if (best_j_array[i] != best_j &&
                best_j_array[i] != best_i)
            {
                tmp_double = best_dists[i];

                BEST_GLOBAL_DIST_MACRO(tmp_double)
            }
            /* flag as stale */
            else
            {
                best_dists[i] = STALE_DISTANCE;
            }
        }


        /* update any nodes that need updating */
        for (i = 0; i < nm2; i++)
        {
#if DEBUG_COUNTS
          counts_nest2++;
#endif

          /* since we only scan the lower triangle, there is no need
           *  to check or update against the newly merged node,
           *  not even to check for inversions
           */

          /* only recalculate it now if it might matter now */
          if (lower_bounds[i] <= best_dist &&
              best_dists[i] == STALE_DISTANCE)
          {
#if DEBUG_COUNTS
              counts_nest3++;
#endif

              /* recycle j to hold order_array[i] */
              j = order_array[i];

              /* it is possible that no lower triangle points are left,
               *  which will result in invalid best k's
               * this is ok, since it's all flagged as stale
               */
              tmp_best_dist = STALE_DISTANCE;
              tmp_second    = STALE_DISTANCE;

              dptr = dist_ptrs[i];
              
              /* we aren't initializing to the first element this time */
              num_ties_k = 0;

              /* nm2 is scanned earlier, so we can skip it here */
              for (k = 0; k < nm2; k++)
              {
                  /* skip higher order triangle; REQUIRED */
                  if (order_array[k] >= j)
                      continue;

#if DEBUG_COUNTS
                  counts_nest4++;
#endif

                  tmp_double = dptr[orig_i_array[k]];
                  if (tmp_double <= tmp_best_dist)
                  {
                      tmp_second = tmp_best_dist;

                      /* break ties */
                      if (tmp_double == tmp_best_dist)
                      {
                          if (ties_random_flag)
                          {
                              ties_k_array[num_ties_k++] = k;
                              continue;
                          }
                      
                          /* take merged node with fewer children */
                          if (nmemb_array[k] > nmemb_array[tmp_best_k] &&
                              ties_fewest_flag)
                          {
                              continue;
                          }

                          /* break ties using input order */
                          if (nmemb_array[k] == nmemb_array[tmp_best_k] ||
                              ties_fewest_flag == 0)
                          {
                              /* take the nearest node when tied */
                              if (order_array[k] < order_array[tmp_best_k])
                                  continue;
                          
                              /* take the earlier node when still tied */
                              /* lower triangle guarantees this already */
                          }
                      }
                      else if (ties_random_flag)
                      {
                          num_ties_k      = 1;
                          ties_k_array[0] = k;
                      }
                      
                      tmp_best_dist = tmp_double;
                      tmp_best_k    = k;
                  }
                  else if (tmp_double < tmp_second)
                  {
                      tmp_second = tmp_double;
                  }
              }

              lower_bounds[i] = tmp_second;
              best_dists[i]   = tmp_best_dist;

              /* we must check tmp_best_dist, rather than tmp_second */
              if (tmp_best_dist != STALE_DISTANCE)
              {
                  /* break ties randomly */
                  if (ties_random_flag && num_ties_k > 1)
                  {
                      tmp_best_k =
                          ties_k_array[(int32_t) (num_ties_k * dxoshiro256p())];
                  }

                  best_j_array[i]     = orig_i_array[tmp_best_k];
                  best_order_array[i] = order_array[tmp_best_k];

                  /* update best dist */
                  BEST_GLOBAL_DIST_MACRO(tmp_best_dist)
              }
          }
        }

        /* break ties randomly */
        if (ties_random_flag && num_ties_ci > 1)
        {
            best_ci = ties_ci_array[(int32_t) (num_ties_ci * dxoshiro256p())];
        }

        --num_nodes;
    }

#if 1
    fprintf(stderr, "NumMembers: %d\n", clust_ptrs[0]->num_members);
#endif

#if DEBUG_COUNTS
    fprintf(stderr, "Nested loop fractions:   %ld   %lf   %lf   %lf   %ld\n",
               counts_nest1,
      (double) counts_nest2 / ((double) counts_nest1 * (double) counts_nest1),
      (double) counts_nest3 /  (double) counts_nest1,
      (double) counts_nest4 /  (double) counts_nest2,
               counts_nest4);
#endif

    if (best_dists)
        free(best_dists);
    if (lower_bounds)
        free(lower_bounds);
    if (new_dists)
        free(new_dists);
    if (dist_ptrs)
        free(dist_ptrs);
    
    if (clust_ptrs)
        free(clust_ptrs);
    if (orig_i_array)
        free(orig_i_array);
    if (orig_to_cur_array)
        free(orig_to_cur_array);
    if (best_j_array)
        free(best_j_array);
    if (order_array)
        free(order_array);
    if (best_order_array)
        free(best_order_array);
    if (nmemb_array)
        free(nmemb_array);

    if (ties_k_array)
        free(ties_k_array);
    if (ties_ci_array)
        free(ties_ci_array);
    
    if (tdata_array)
        free(tdata_array);

#ifdef GPERF
ProfilerStop();
#endif
        
    return clust_pool;
}


/* reorder nodes in a multifurcating tree
 * leave ties as-is
 */
void flip_nodes_multi(struct tree_node **node_ptr_array,
                      int32_t num_nodes, struct options *opt)
{
    struct tree_node **node_ptr_array2 = NULL;    /* sorted on clevel */
    struct tree_node  *node_ptr;
    int32_t i;
    

    for (i = 0; i < num_nodes; i++)
    {
        node_ptr = node_ptr_array[i];

        /* trees can come from two routes at the moment:
         *   hierarchical clustering function using only left/right ptrs
         *   read from file, which uses child_ptrs
         *
         * detect which one and fill in the other
         */
        /* binary tree, left/right ptrs */
        if (node_ptr->left_ptr && node_ptr->right_ptr &&
            node_ptr->child_ptrs == NULL)
        {
            node_ptr->num_children  = 2;
            node_ptr->child_ptrs    = calloc(2, sizeof(struct tree_node *));
            node_ptr->child_ptrs[0] = node_ptr->left_ptr;
            node_ptr->child_ptrs[1] = node_ptr->right_ptr;
        }
        /* tree read in from file, fill in left/right with first two */
        else if (node_ptr->child_ptrs && node_ptr->num_children >= 2 &&
                 node_ptr->left_ptr == NULL && node_ptr->right_ptr == NULL)
        {
            node_ptr->left_ptr  = node_ptr->child_ptrs[0];
            node_ptr->right_ptr = node_ptr->child_ptrs[1];
        }
        /* must be a leaf node */
        else if (node_ptr->child_ptrs == NULL &&
                 node_ptr->left_ptr   == NULL &&
                 node_ptr->right_ptr  == NULL)
        {
            node_ptr->num_children = 0;
            node_ptr->num_members  = 1;
            node_ptr->clevel       = 0;
        }

        /* reorder branches */
        if (node_ptr->num_children >= 2 && node_ptr->child_ptrs)
        {
            if (opt->tree_flip_size_flag)
            {
                qsort(node_ptr->child_ptrs, node_ptr->num_children,
                      sizeof(struct tree_node *), cmp_node_ptr_size);
            }
            if (opt->tree_flip_edge_flag)
            {
                qsort(node_ptr->child_ptrs, node_ptr->num_children,
                      sizeof(struct tree_node *), cmp_node_ptr_edge);
            }
            if (opt->tree_flip_avg_flag)
            {
                qsort(node_ptr->child_ptrs, node_ptr->num_children,
                      sizeof(struct tree_node *), cmp_node_ptr_avg_length);
            }

            /* sync up left/right child pointers */
            node_ptr->left_ptr  = node_ptr->child_ptrs[0];
            node_ptr->right_ptr = node_ptr->child_ptrs[1];
        }
    }
    
    /* more flipping for edges */
    if (opt->tree_flip_edge_flag)
    {
        node_ptr_array2 = (struct tree_node **) malloc(num_nodes *
                          sizeof(struct tree_node *));
        memcpy(node_ptr_array2, node_ptr_array,
               num_nodes * sizeof(struct tree_node *));

        /* sort nodes by clevel, so that children are before parents */
        qsort(node_ptr_array2, num_nodes, sizeof(struct tree_node *),
              cmp_node_ptr_clevel);

        /* flip nodes if node is a left node of parent node */
        for (i = num_nodes - 1; i >= 0; --i)
        {
            node_ptr = node_ptr_array2[i];

            /* reorder branches */
            if (node_ptr->num_children >= 2 && node_ptr->child_ptrs &&
                node_ptr->parent_ptr &&
                node_ptr == node_ptr->parent_ptr->left_ptr)
            {
                /* reverse sort */
                qsort(node_ptr->child_ptrs, node_ptr->num_children,
                      sizeof(struct tree_node *), cmp_node_ptr_edge_longer);

                /* sync up left/right child pointers */
                node_ptr->left_ptr  = node_ptr->child_ptrs[0];
                node_ptr->right_ptr = node_ptr->child_ptrs[1];
            }
        }

        if (node_ptr_array2)
        {
            free(node_ptr_array2);
            node_ptr_array2 = NULL;
        }
    }
}


/* depower the distances from wardu method
 * output tree distances we would have gotten with ward2 equivalent
 */
void depower_tree(struct tree_node **node_ptr_array,
                  int32_t num_nodes, struct options *opt)
{
    struct tree_node **node_ptr_array2 = NULL;    /* sorted on clevel */
    struct tree_node  *node_ptr, *child_ptr;
    double power, depower, dist;
    int32_t i, j;
    

    /* --distpow option */
    power   = opt->dist_pow;
    
    /* sanity check, set power to 1 if invalid power was given */
    if (power <= 0)
        power = 1.0;

    /* nothing to do, return */
    if (power == 1.0)
    {
        return;
    }
    
    depower = 1.0 / power;
    

    /* allocate temporary node pointer array for sorting */
    node_ptr_array2 = (struct tree_node **) malloc(num_nodes *
                      sizeof(struct tree_node *));
    memcpy(node_ptr_array2, node_ptr_array,
           num_nodes * sizeof(struct tree_node *));


    /* first, make sure we have all our tree structures set up correctly */
    for (i = 0; i < num_nodes; i++)
    {
        node_ptr = node_ptr_array[i];

        /* trees can come from two routes at the moment:
         *   hierarchical clustering function using only left/right ptrs
         *   read from file, which uses child_ptrs
         *
         * detect which one and fill in the other
         */
        /* binary tree, left/right ptrs */
        if (node_ptr->left_ptr && node_ptr->right_ptr &&
            node_ptr->child_ptrs == NULL)
        {
            node_ptr->num_children  = 2;
            node_ptr->child_ptrs    = calloc(2, sizeof(struct tree_node *));
            node_ptr->child_ptrs[0] = node_ptr->left_ptr;
            node_ptr->child_ptrs[1] = node_ptr->right_ptr;
        }
        /* tree read in from file, fill in left/right with first two */
        else if (node_ptr->child_ptrs && node_ptr->num_children >= 2 &&
                 node_ptr->left_ptr == NULL && node_ptr->right_ptr == NULL)
        {
            node_ptr->left_ptr  = node_ptr->child_ptrs[0];
            node_ptr->right_ptr = node_ptr->child_ptrs[1];
        }
        /* must be a leaf node */
        else if (node_ptr->child_ptrs == NULL &&
                 node_ptr->left_ptr   == NULL &&
                 node_ptr->right_ptr  == NULL)
        {
            node_ptr->num_children = 0;
            node_ptr->num_members  = 1;
            node_ptr->clevel       = 0;
        }
    }


    /* sort nodes by clevel, so that children are before parents */
    qsort(node_ptr_array2, num_nodes, sizeof(struct tree_node *),
          cmp_node_ptr_clevel);

    /* first, recalculate all the average distances in case we haven't yet */
    for (i = 0; i < num_nodes; i++)
    {
        node_ptr = node_ptr_array2[i];

        /* sum distances of child nodes */
        dist = 0;
        for (j = 0; j < node_ptr->num_children; j++)
        {
            child_ptr = node_ptr->child_ptrs[j];
            
            /* leaves */
            if (child_ptr->num_members == 1)
            {
                child_ptr->avg_length = 0.0;
            }

            dist += child_ptr->edge_length + child_ptr->avg_length;
        }
        
        if (node_ptr->num_children)
        {
            dist /= node_ptr->num_children;
        }
        
        node_ptr->avg_length = dist;
    }

    for (i = 0; i < num_nodes; i++)
    {
        node_ptr = node_ptr_array2[i];

        dist = node_ptr->avg_length;

        if (node_ptr->num_children)
        {
            dist *= node_ptr->num_children;
        }

        /* use more accurate evaluations where easily available */
        if (power == 2.0)
        {
            dist  = sqrt(dist);
        }
        else if (power == 1.5)
        {
            dist  = cbrt(dist);
            dist *= dist;
        }
        else if (power == 3.0)
        {
            dist  = cbrt(dist);
        }
        else
        {
            dist  = pow(dist, depower);
        }

        if (node_ptr->num_children)
        {
            dist /= node_ptr->num_children;
        }
        
        node_ptr->avg_length = dist;
    }

    /* recalculate edges after we've depowered average lengths */
    for (i = 0; i < num_nodes; i++)
    {
        node_ptr = node_ptr_array2[i];

        for (j = 0; j < node_ptr->num_children; j++)
        {
            child_ptr = node_ptr->child_ptrs[j];

            child_ptr->edge_length = node_ptr->avg_length -
                                     child_ptr->avg_length;
        }
    }

    if (node_ptr_array2)
    {
        free(node_ptr_array2);
        node_ptr_array2 = NULL;
    }
}


/* only to be run on generated bifurcated trees */
/* *NOT* to be run on trees read in from a file */
/* should be run *BEFORE* creating the newick string */
void clean_tree(struct tree_node *root_node,
                struct tree_node ***return_node_ptr_array,
                int32_t *return_num_nodes,
                struct tree_node ***return_leaf_ptr_array,
                int32_t *return_num_leaves)
{
    struct tree_node **node_stack = NULL;
    int32_t max_stack_size = 0;
    int32_t stack_idx = 0;
    
    struct tree_node *node_ptr;
    struct tree_node *left_ptr;
    struct tree_node *right_ptr;

    int32_t num_nodes = 1;
    int32_t num_leaves = 0;
    
    struct tree_node **node_ptr_array = *return_node_ptr_array;
    struct tree_node **leaf_ptr_array = *return_leaf_ptr_array;

    node_ptr_array = (struct tree_node **) realloc(node_ptr_array,
                        num_nodes * sizeof(struct tree_node *));
    node_ptr_array[num_nodes - 1] = root_node;
    
    /* push root node onto stack */
    node_stack = (struct tree_node **) malloc(sizeof(struct tree_node *));
    node_stack[0] = root_node;
    root_node->walk_flag = 0;
    max_stack_size = 1;
    
    while (stack_idx >= 0)
    {
        node_ptr = node_stack[stack_idx];

        /* this is going to be confusing, since it is a hack
         *
         * swap left_ptr and right_ptr so that left_ptr always has the
         *  the larger number of children, so that it is output first
         */
        if (node_ptr->left_ptr && node_ptr->right_ptr)
        {
            struct tree_node *temp_ptr;
        
            if (node_ptr->left_ptr->num_members <
                node_ptr->right_ptr->num_members)
            {
                temp_ptr            = node_ptr->left_ptr;
                node_ptr->left_ptr  = node_ptr->right_ptr;
                node_ptr->right_ptr = temp_ptr;
            }
            else if (node_ptr->left_ptr->num_members ==
                     node_ptr->right_ptr->num_members)
            {
              /* break ties on distance */
              if (node_ptr->left_ptr->avg_length !=
                  node_ptr->right_ptr->avg_length)
              {
                  if (node_ptr->left_ptr->avg_length >
                      node_ptr->right_ptr->avg_length)
                  {
                      temp_ptr            = node_ptr->left_ptr;
                      node_ptr->left_ptr  = node_ptr->right_ptr;
                      node_ptr->right_ptr = temp_ptr;
                  }
              }
              /* break ties on name */
              else if (node_ptr->left_ptr->name && node_ptr->right_ptr->name)
              {
                /* numeric integers (strings of floats are treated as text) */
                if (is_all_digits(node_ptr->left_ptr->name) &&
                    is_all_digits(node_ptr->right_ptr->name))
                {
                    if (atol(node_ptr->left_ptr->name) >
                        atol(node_ptr->right_ptr->name))
                    {
                        temp_ptr            = node_ptr->left_ptr;
                        node_ptr->left_ptr  = node_ptr->right_ptr;
                        node_ptr->right_ptr = temp_ptr;
                    }
                }
                /* non-integer strings */
                else if (strcmp(node_ptr->left_ptr->name,
                                node_ptr->right_ptr->name) > 0)
                {
                    temp_ptr            = node_ptr->left_ptr;
                    node_ptr->left_ptr  = node_ptr->right_ptr;
                    node_ptr->right_ptr = temp_ptr;
                }
              }
            }
        }

        left_ptr = node_ptr->left_ptr;
        right_ptr = node_ptr->right_ptr;
        
        if (left_ptr && node_ptr->walk_flag == 0)
        {
            /* initialize walk flags */
            left_ptr->walk_flag = 0;
            if (right_ptr)
            {
                right_ptr->walk_flag = 0;
            }
            
            left_ptr->parent_ptr = node_ptr;
            
            if (node_ptr->left_ptr->left_ptr ||
                node_ptr->left_ptr->right_ptr)
            {
                node_ptr_array = (struct tree_node **) realloc(node_ptr_array,
                            (num_nodes+1) * sizeof(struct tree_node *));
                node_ptr_array[num_nodes++] = node_ptr->left_ptr;
            }

            node_ptr->walk_flag = 1;

            /* push new node */
            stack_idx++;
            if (stack_idx >= max_stack_size)
            {
                max_stack_size = stack_idx + 1;
                node_stack = realloc(node_stack,
                            max_stack_size * sizeof(struct tree_node *));
            }
            node_stack[stack_idx] = left_ptr;
                
            continue;
        }

        /* walk down right_ptr */
        else if (right_ptr && node_ptr->walk_flag == 1)
        {
            right_ptr->parent_ptr = node_ptr;

            if (node_ptr->right_ptr->left_ptr ||
                node_ptr->right_ptr->right_ptr)
            {
                node_ptr_array = (struct tree_node **) realloc(node_ptr_array,
                            (num_nodes+1) * sizeof(struct tree_node *));
                node_ptr_array[num_nodes++] = node_ptr->right_ptr;
            }

            node_ptr->walk_flag = 2;

            /* push new node */
            stack_idx++;
            if (stack_idx >= max_stack_size)
            {
                max_stack_size = stack_idx + 1;
                node_stack = realloc(node_stack,
                            max_stack_size * sizeof(struct tree_node *));
            }
            node_stack[stack_idx] = right_ptr;
            
            continue;
        }

        /* back up a node on the stack */
        else
        {
            /* print leaf node */
            if (left_ptr == NULL && right_ptr == NULL)
            {
                leaf_ptr_array = (struct tree_node **) realloc(leaf_ptr_array,
                 (num_leaves+1) * sizeof(struct tree_node *));
                leaf_ptr_array[num_leaves++] = node_ptr;
            }
        
            stack_idx--;
        }
    }
    
    if (node_stack) free(node_stack);

    *return_node_ptr_array = node_ptr_array;
    *return_leaf_ptr_array = leaf_ptr_array;
    *return_num_nodes = num_nodes;
    *return_num_leaves = num_leaves;
}


/* clobbers original string with blessed string */
void bless_tree_names(struct tree_node *root_node,
                      struct tree_node **node_ptr_array,
                      int32_t num_nodes,
                      struct tree_node **leaf_ptr_array,
                      int32_t num_leaves)
{
    int32_t i;
    char *sptr;

    for (i = 0; i < num_nodes; i++)
    {
        for (sptr = node_ptr_array[i]->name; *sptr; sptr++)
        {
            switch (*sptr)
            {
                case '(' :
                case ')' :
                case ';' :
                case ':' :
                case ',' :
                    *sptr = '-';
                    break;
            }
        }
    }

    for (i = 0; i < num_leaves; i++)
    {
        for (sptr = leaf_ptr_array[i]->name; *sptr; sptr++)
        {
            switch (*sptr)
            {
                case '(' :
                case ')' :
                case ';' :
                case ':' :
                case ',' :
                    *sptr = '-';
                    break;
            }
        }
    }
}


/* can handle multifating trees */
char * create_newick_string_multi(struct tree_node *root_node,
                                  int32_t one_line_per_tree,
                                  char overwrite_names_flag)
{
    struct tree_node **node_stack = NULL;
    struct tree_node *node_ptr;
    int32_t max_stack_size = 0;
    int32_t stack_idx = 0;
    int32_t node_index = 0;

    char *string         = NULL;
    char *str_clust_good = "*";
    char *str_clust_bad  = "**";    /* ** adds clutter, so use for bad */
    char *str_blank      = "";
    char *str_ptr        = NULL;
    uint64_t string_len  = 0, buf_len = 1, new_len;

    uint32_t len_str_clust_good = strlen(str_clust_good);
    uint32_t len_str_clust_bad  = strlen(str_clust_bad);


    /* allocate initial null-terminated buffer
     * estimate initial size to save on future reallocs
     * assume root string is representative of label lengths
     *     ':' + ~8 characters per 2n-1 leaves
     *     1 ',' per n-1 leaf
     *     1 '(' per n-1 leaf
     *     1 ')' per n-1 leaf
     *     1 '\n' per n-1 leaf
     *
     *     = name_len + 18 + 4
     *
     * ;\n at the end
     * terminal \0
     */
    
    new_len = root_node->num_members *
              strlen(root_node->name) +
              (root_node->num_members - 1) * 22 +
              9 + 2 + 1;
    if (one_line_per_tree)	/* subtract newlines if not printed */
        new_len -= (root_node->num_members - 1);
    buf_len = new_len;
    
    string = (char *) malloc(new_len * sizeof(char));
    string[0] = '\0';


    /* push root node onto stack */
    node_stack = (struct tree_node **) malloc(sizeof(struct tree_node *));
    node_stack[0] = root_node;
    root_node->walk_child_idx = -1;
    max_stack_size = 1;
    
    while (stack_idx >= 0)
    {
        node_ptr = node_stack[stack_idx];
        
        /* trees can come from two routes at the moment:
         *   hierarchical clustering function using only left/right ptrs
         *   read from file, which uses child_ptrs
         *
         * detect which one and fill in the other
         */
        /* binary tree, left/right ptrs */
        if (node_ptr->left_ptr && node_ptr->right_ptr &&
            node_ptr->child_ptrs == NULL)
        {
            node_ptr->num_children  = 2;
            node_ptr->child_ptrs    = calloc(2, sizeof(struct tree_node *));
            node_ptr->child_ptrs[0] = node_ptr->left_ptr;
            node_ptr->child_ptrs[1] = node_ptr->right_ptr;
        }
        /* tree read in from file, fill in left/right with first two */
        else if (node_ptr->child_ptrs && node_ptr->num_children >= 2 &&
                 node_ptr->left_ptr == NULL && node_ptr->right_ptr == NULL)
        {
            node_ptr->left_ptr  = node_ptr->child_ptrs[0];
            node_ptr->right_ptr = node_ptr->child_ptrs[1];
        }
        /* must be a leaf node */
        else if (node_ptr->child_ptrs == NULL &&
                 node_ptr->left_ptr   == NULL &&
                 node_ptr->right_ptr  == NULL)
        {
            node_ptr->num_children = 0;
            node_ptr->num_members  = 1;
            node_ptr->clevel       = 0;
        }

        
        /* advance to next child */
        if (node_ptr->num_children &&
            node_ptr->walk_child_idx < node_ptr->num_children - 1)
        {
            /* store which child we're about to walk down */
            node_ptr->walk_child_idx++;


            /* open new branch */
            if (node_ptr->walk_child_idx == 0)
            {
                new_len = string_len + 1;
                if (new_len + 1 > buf_len)
                {
                    buf_len = (new_len + 1) * MEM_OVERHEAD;
                    string = (char *) realloc(string, buf_len * sizeof(char));
                }
                string[new_len - 1] = '(';
                string_len          = new_len;
            }

            /* sibling node */
            else if (node_ptr->walk_child_idx)
            {
                if (one_line_per_tree)
                {
                    new_len = string_len + 1;
                    if (new_len + 1 > buf_len)
                    {
                        buf_len = (new_len + 1) * MEM_OVERHEAD;
                        string = (char *) realloc(string, buf_len *
                                                          sizeof(char));
                    }
                    string[new_len - 1] = ',';
                    string_len          = new_len;
                }
                else
                {
                    new_len = string_len + 2;
                    if (new_len + 1 > buf_len)
                    {
                        buf_len = (new_len + 1) * MEM_OVERHEAD;
                        string = (char *) realloc(string, buf_len * sizeof(char));
                    }
                    string[new_len - 2] = ',';
                    string[new_len - 1] = '\n';
                    string_len          = new_len;
                }
            }


            /* push new node */
            stack_idx++;
            if (stack_idx >= max_stack_size)
            {
                max_stack_size = stack_idx + 1;
                node_stack = realloc(node_stack,
                            max_stack_size * sizeof(struct tree_node *));
            }
            node_stack[stack_idx] =
                node_ptr->child_ptrs[node_ptr->walk_child_idx];
            node_ptr->child_ptrs[node_ptr->walk_child_idx]->walk_child_idx=-1;
        }

        /* back up a node on the stack */
        else
        {
            /* print leaf node */
            if (node_ptr->num_children == 0)
            {
                new_len = string_len +
                          snprintf(NULL, 0, "%s:%.14g",
                                   node_ptr->name, node_ptr->edge_length);

                if (new_len + 1 > buf_len)
                {
                    buf_len = (new_len + 1) * MEM_OVERHEAD;
                    string = (char *) realloc(string, buf_len * sizeof(char));
                }

                snprintf(string + string_len,
                         (new_len - string_len + 1) * sizeof(char),
                         "%s:%.14g",
                         node_ptr->name, node_ptr->edge_length);

                string_len = new_len;
            }
            
            /* close branch and back up */
            if (node_ptr->parent_ptr &&
                node_ptr->parent_ptr->walk_child_idx ==
                 node_ptr->parent_ptr->num_children - 1)
            {
                node_index++;

                /* replace name with node number for debugging */
                if (overwrite_names_flag && node_ptr->parent_ptr->name)
                {
                    free(node_ptr->parent_ptr->name);

                    /* allocate enough room for the name string */
                    new_len = snprintf(NULL, 0, "Node%05d", node_index);

                    node_ptr->parent_ptr->name =
                        malloc((new_len+1) * sizeof(char));
                    node_ptr->parent_ptr->name[new_len] = '\0';

                    sprintf(node_ptr->parent_ptr->name,
                        "Node%05d", node_index);
                }
            
                new_len = string_len +
                          snprintf(NULL, 0,
                                   ")Node%05d:%.14g",
                                   node_index,
                                   node_ptr->parent_ptr->edge_length);
                
                /* increase length if it is in a cluster */
                str_ptr = str_blank;
                if (node_ptr->parent_ptr->cluster_flag)
                {
                    if (node_ptr->parent_ptr->cluster_good_flag)
                    {
                        str_ptr  = str_clust_good;
                        new_len += len_str_clust_good;
                    }
                    else
                    {
                        str_ptr  = str_clust_bad;
                        new_len += len_str_clust_bad;
                    }
                }

                if (new_len + 1 > buf_len)
                {
                    buf_len = (new_len + 1) * MEM_OVERHEAD;
                    string = (char *) realloc(string, buf_len * sizeof(char));
                }

                snprintf(string + string_len,
                         (new_len - string_len + 1) * sizeof(char),
                         ")%sNode%05d:%.14g",
                         str_ptr,
                         node_index,
                         node_ptr->parent_ptr->edge_length);

                string_len = new_len;
            }

            stack_idx--;
        }
    }

    new_len = string_len + 2;
    if (new_len + 1 > buf_len)
    {
        buf_len = (new_len + 1) * MEM_OVERHEAD;
        string = (char *) realloc(string, buf_len * sizeof(char));
    }
    string[new_len - 2] = ';';
    string[new_len - 1] = '\n';
    string_len          = new_len;
    
    string[string_len] = '\0';
    
    if (node_stack) free(node_stack);
    
    return string;
}


/* can handle multifating trees */
double flag_clusters_multi(struct   tree_node *root_node,
                           double   min_tiny_fraction, int32_t tiny_size,
                           double   max_size_fraction,
                           double   max_dist,
                           int32_t  target_n_clusters,
                           int32_t *return_n_clusters,
                           struct   options *opt,
                           int32_t  store_flag,
                           int32_t  print_flag)
{
    struct  tree_node **node_stack        = NULL;
    struct  tree_node *cluster_ptr        = NULL;
    struct  tree_node *node_ptr, *child_ptr;
    int32_t max_stack_size = 0;
    int32_t stack_idx = 0;

    int32_t max_size;
    int32_t min_size;
    
    int32_t no_merge_flag = opt->clusters_no_merge_flag;
    
    int32_t num_clusters;
    int32_t num_good_clusters        = 0;
    int32_t num_good_leaves          = 0;
    int32_t num_unclustered_regions  = 0;
    int32_t num_unclustered_leaves   = 0;
    int32_t in_cluster_flag          = 0;
    int32_t good_cluster_flag        = 0;
    int32_t out_of_good_cluster_flag = 1;
    
    int32_t num_best_clusters        = 0;
    int32_t num_best_leaves          = 0;
    int32_t num_worst_clusters       = 0;
    int32_t num_worst_leaves         = 0;
    int32_t node_ok_flag, child_ok_count, all_ok_flag;
    int32_t i;
    
    double dist, dist_fraction, dist_child;
    double score;
    double avg_dist_score   = 0.0;
    double avg_size_score   = 0.0;
    double avg_clevel_score = 0.0;
    double n_modified;

    /* cap target number of clusters at number of leaves */
    if (target_n_clusters > root_node->num_members)
        target_n_clusters = root_node->num_members;
    
    /* reducing the multiplier from 0.025 to 0.02 causes bad behavior */
    min_size = (int32_t) (root_node->num_members * min_tiny_fraction + 0.5);
    if (min_size < tiny_size)
    {
        min_size = tiny_size;
    }
    
    /* disable minimum sizes, and therefore small cluster merging/penalties */
    if (no_merge_flag)
    {
        min_size  = 0;
        tiny_size = 0;
    }

    max_size = (int32_t) (max_size_fraction * root_node->num_members + 0.5);


    /* push root node onto stack */
    node_stack = (struct tree_node **) malloc(sizeof(struct tree_node *));
    node_stack[0] = root_node;
    root_node->walk_child_idx = -1;
    max_stack_size = 1;
    
    num_clusters = 0;
    while (stack_idx >= 0)
    {
        node_ptr = node_stack[stack_idx];
        
#if 0
        /* trees can come from two routes at the moment:
         *   hierarchical clustering function using only left/right ptrs
         *   read from file, which uses child_ptrs
         *
         * detect which one and fill in the other
         */
        /* binary tree, left/right ptrs */
        if (node_ptr->left_ptr && node_ptr->right_ptr &&
            node_ptr->child_ptrs == NULL)
        {
            node_ptr->num_children  = 2;
            node_ptr->child_ptrs    = calloc(2, sizeof(struct tree_node *));
            node_ptr->child_ptrs[0] = node_ptr->left_ptr;
            node_ptr->child_ptrs[1] = node_ptr->right_ptr;
        }
        /* tree read in from file, fill in left/right with first two */
        else if (node_ptr->child_ptrs && node_ptr->num_children >= 2 &&
                 node_ptr->left_ptr == NULL && node_ptr->right_ptr == NULL)
        {
            node_ptr->left_ptr  = node_ptr->child_ptrs[0];
            node_ptr->right_ptr = node_ptr->child_ptrs[1];
        }
        /* must be a leaf node */
        else if (node_ptr->child_ptrs == NULL &&
                 node_ptr->left_ptr   == NULL &&
                 node_ptr->right_ptr  == NULL)
        {
            node_ptr->num_children = 0;
            node_ptr->num_members  = 1;
            node_ptr->clevel       = 0;
        }
#endif


        dist = node_ptr->avg_length;
        if (root_node->avg_length == 0.0)
            dist = node_ptr->num_members;

        node_ok_flag   = 0;
        child_ok_count = 0;
        all_ok_flag    = 0;

        if (node_ptr != root_node &&
            in_cluster_flag == 0 &&
            node_ptr->walk_child_idx == -1)
        {
            if (dist <= max_dist)
            {
                node_ok_flag = 1;
            }

            /* allow less-good clusters if all children are OK */
            if (min_size == 0 ||
                (node_ptr->num_children >= 2 &&
                 (dist > max_dist || node_ptr->num_members < min_size)))
            {
                for (i = 0; i < node_ptr->num_children; i++)
                {
                    child_ptr = node_ptr->child_ptrs[i];

                    dist_child = child_ptr->avg_length;
                    if (root_node->avg_length == 0.0)
                        dist_child = child_ptr->num_members;

                    /* good child, must all be smaller than min size */
                    if (dist_child             <= max_dist &&
                        child_ptr->num_members <= max_size &&
                        child_ptr->num_members <  min_size)
                    {
                        child_ok_count++;
                    }
                }
                
                if (child_ok_count == node_ptr->num_children)
                    all_ok_flag = 1;
            }
        }


        /* found a cluster */
        /* allow larger distance if all child nodes are OK */
        if (node_ok_flag || all_ok_flag)
        {
            cluster_ptr              = node_ptr;
            in_cluster_flag          = 1;

            good_cluster_flag = 0;
            if (node_ptr->num_members <= max_size)
            {
                if (node_ptr->num_members >= min_size || all_ok_flag)
                {
                    good_cluster_flag        = 1;
                    out_of_good_cluster_flag = 0;
                    
                    if (all_ok_flag && min_size > 0)
                    {
                        num_worst_clusters++;
                        num_worst_leaves += node_ptr->num_members;
                    }
                    else
                    {
                        num_best_clusters++;
                        num_best_leaves += node_ptr->num_members;
                    }
                }
            }

            if (good_cluster_flag)
            {
                if (node_ptr->num_members)
                {
                    /* even X^0.95 isn't sufficient
                     * for ties-fewest/orderall-ties, needs X^1,
                     * but sqrt(X) seems to work well for not-all-ties
                     *
                     * sqrt(X) is too loose for all-ties ties-fewest/order
                     */
                    avg_size_score += sqrt(node_ptr->num_members);
                }
                if (node_ptr->clevel)
                {
                    avg_clevel_score += node_ptr->clevel;
                }

                /* down-scale large distances, in case highly dissimilar leaves
                 * have MUCH larger distances than similar leaves,
                 * such as when too-dissimilar nodes are set to a fixed
                 * very large distance to indicate no relation whatsoever
                 */
                score = 0;
                if (node_ptr->avg_length)
                {
                    score           = log10(root_node->num_members *
                                            node_ptr->avg_length);
                    avg_dist_score += score;
                }

                num_good_leaves += node_ptr->num_members;

                num_good_clusters++;
            }


            num_clusters++;

            if (print_flag && good_cluster_flag)
            {
                dist_fraction = 0.0;
                if (root_node->avg_length)
                {
                    dist_fraction = node_ptr->avg_length /
                                    root_node->avg_length;
                }
                
                fprintf(stderr, "%s\t%d\t%0.6f\t%0.6f\n",
                    node_ptr->name,
                    node_ptr->num_members,
                    dist_fraction,
                    (double) node_ptr->num_members / root_node->num_members);
            }
        }
        /* we're backing out of the current cluster */
        else if (node_ptr == cluster_ptr &&
                 node_ptr->walk_child_idx == node_ptr->num_children - 1)
        {
            cluster_ptr       = NULL;
            in_cluster_flag   = 0;
            good_cluster_flag = 0;
        }
        
        
#if 0
        /* currently bugged, doesn't take merged too-small clusters under
         * the parent cluster into account, parent can be flagged UNCLUSTERED
         */
        if ((node_ptr->num_members < min_size ||
             node_ptr->num_children == 0) &&
            node_ptr->walk_child_idx == -1 &&
            good_cluster_flag == 0 &&
            out_of_good_cluster_flag == 0)
        {
            out_of_good_cluster_flag = 1;
            num_unclustered_regions++;
            
            if (print_flag)
                fprintf(stderr, "UNCLUSTERED\t%d\t%s\n",
                    num_unclustered_regions, node_ptr->name);
        }
#endif


        /* back out early
         * no need to recurse further, we're not flagging leaves
         * must back out after unclustered bookkeeping
         */
        if (store_flag == 0 &&
            (good_cluster_flag || node_ptr->num_members < min_size))
        {
            cluster_ptr       = NULL;
            in_cluster_flag   = 0;
            good_cluster_flag = 0;
        
            node_ptr->walk_child_idx = node_ptr->num_children - 1;

            stack_idx--;
            continue;
        }


        /* store cluster information */
        if (store_flag && node_ptr->walk_child_idx == -1)
        {
            /* only store for leaves and parent cluster node */
            if (node_ptr == cluster_ptr || node_ptr->num_children == 0)
            {
                if (in_cluster_flag)
                {
                     node_ptr->cluster_num       = num_clusters;
                     node_ptr->cluster_flag      = 1;
                     node_ptr->cluster_good_flag = 0;
                     
                     if (good_cluster_flag)
                         node_ptr->cluster_good_flag = 1;
                }
                /* should no longer ever happen */
                else
                {
                     node_ptr->cluster_num       = num_unclustered_regions;
                     node_ptr->cluster_flag      = 0;
                     node_ptr->cluster_good_flag = 0;
                }
            }
        }

        
        /* advance to next child */
        if (node_ptr->num_children &&
            node_ptr->walk_child_idx < node_ptr->num_children - 1)
        {
            /* store which child we're about to walk down */
            node_ptr->walk_child_idx++;

            /* push new node */
            stack_idx++;
            if (stack_idx >= max_stack_size)
            {
                max_stack_size = stack_idx + 1;
                node_stack = realloc(node_stack,
                               max_stack_size * sizeof(struct tree_node *));
            }
            node_stack[stack_idx] =
                node_ptr->child_ptrs[node_ptr->walk_child_idx];
            node_ptr->child_ptrs[node_ptr->walk_child_idx]->walk_child_idx=-1;
        }
        /* back up from singleton leaf cluster */
        else if (node_ptr == cluster_ptr && node_ptr->num_children == 0)
        {
            cluster_ptr       = NULL;
            in_cluster_flag   = 0;
            good_cluster_flag = 0;

            stack_idx--;
        }
        /* back up a node on the stack */
        else
        {
            stack_idx--;
        }
    }
    
    num_unclustered_leaves = root_node->num_members - num_good_leaves;
    if (num_good_clusters)
    {
        avg_dist_score   /= num_good_clusters;
        avg_size_score   /= num_good_clusters;
        avg_clevel_score /= num_good_clusters;
        
        avg_clevel_score /= root_node->clevel;
    }
    else
    {
        avg_dist_score   = 0.0;
        avg_size_score   = 0.0;
        avg_clevel_score = 0.0;
    }


#if 1
    n_modified = (num_good_clusters - (1.0 / 3.0) *
                  (num_unclustered_regions + num_worst_clusters)) *
                 (num_best_leaves / (double) root_node->num_members);
#else
    /* 0.5x is too risky, generates too many small clusters on too many
     * random all-zeros trees.  2/3rds is much safer.
     */
    n_modified = (num_good_clusters - (2.0 / 3.0) * num_unclustered_regions) *
                 (num_good_leaves / (double) root_node->num_members);
#endif

    score = n_modified;
    
#if 1
    if (target_n_clusters == 0)
        score *= avg_clevel_score;
#endif

#if 1
    /* improves clustering on some trees, even with log10 distances */
    if (target_n_clusters == 0 && root_node->avg_length)
    {
        score *= avg_dist_score;
        score /= log10(root_node->num_members * root_node->avg_length);
    }
#endif

#if 1
    if (target_n_clusters == 0 && root_node->avg_length == 0.0)
    {
        score *= avg_size_score;


        /* scale it last, otherwise unexpected things happen */
        /* score /= root_node->num_members; */
    }
#endif

    /* penalize the further away we are from the target number of clusters */
    dist = num_good_clusters + num_unclustered_regions;
    if (target_n_clusters && dist)
    {
        dist   = labs(dist - target_n_clusters);
        dist   = 1 + dist * dist * dist;

        score /= dist;
    }
    
    /* flag as bad cluster */
    if (n_modified < 0)
    {
        score = -999.99;
    }

    if (print_flag)
    {
        fprintf(stderr, "Cluster statistics:\t%d\t%d\t%d\t%d\t%0.5f\n",
                num_best_clusters,
                num_unclustered_regions + num_worst_clusters,
                num_best_leaves,
                num_unclustered_leaves + num_worst_leaves, score);
    }

    *return_n_clusters = num_good_clusters + num_unclustered_regions;

    if (node_stack)
        free(node_stack);

    return score;
}


double scan_clusters(struct   tree_node *root_node,
                     struct   tree_node **node_ptr_array,
                     int32_t  num_nodes, double min_tiny_fraction,
                     int32_t  tiny_size, double max_size_fraction,
                     double   max_dist_fraction,
                     int32_t  target_n_clusters,
                     int32_t *return_n_clusters,
                     struct   options *opt,
                     int32_t  print_flag)
{
    struct  tree_node *node_ptr;
    double  score, old_score, old_best_score;
    double  best_score = -9999.9999;
    double  dist_cutoff;
    
    int32_t i, best_i, old_best_i;
    int32_t min_size;


    min_size = (int32_t) (root_node->num_members * min_tiny_fraction + 0.5);
    if (min_size < tiny_size)
    {
        min_size = tiny_size;
    }

    /* we need to start with a reasonably sized cluster,
     * otherwise, gigantic trees take too long
     */
    if (opt->clusters_no_merge_flag)
    {
        /* rather than disable, set to reasonably small percentage */
        if (min_size > (int32_t) (0.001 * root_node->num_members + 0.5))
            min_size = (int32_t) (0.001 * root_node->num_members + 0.5);
    }

    /* sort node by increasing average distance to leaves */
    qsort(node_ptr_array, num_nodes, sizeof(struct tree_node *),
          cmp_node_ptr_avg_length);

    
    if (print_flag)
    {
        fprintf(stderr, "Max cluster size cutoff:\t%f\t%d\n",
                max_size_fraction,
                (int32_t) (max_size_fraction * root_node->num_members + 0.5));
    }


    old_score        = -9999.9999;
    old_best_score   = -9999.9999;
    old_best_i       = -9999;

    best_i = 0;

    for (i = 0; i < num_nodes; i++)
    {
        node_ptr = node_ptr_array[i];

#if 1
        /* unlikely to find many clusters of sufficient size */
        /* generally results in a significant speedup */
        if (node_ptr->num_members < min_size)
            continue;
#endif
    
        if (root_node->avg_length &&
            node_ptr->avg_length / root_node->avg_length > max_dist_fraction)
        {
            break;
        }

        dist_cutoff = node_ptr->avg_length;
        if (root_node->avg_length == 0.0)
            dist_cutoff = node_ptr->num_members;

        score = flag_clusters_multi(root_node,
            min_tiny_fraction, 5,
            max_size_fraction, dist_cutoff,
            target_n_clusters, return_n_clusters, opt, 0, 0);

        if (score >= best_score)
        {
            old_best_i = best_i;
            best_i     = i;
            best_score = score;
        }

        old_score = score;
    }

    /* For both cutoffs, higher cutoffs = bigger clusters.
     * The new scoring system, which takes number of unclustered
     * leaves into account, should take care of mots bad clustering
     * behaveio, but we'll still take advantage of the increasing nature of
     * the scan * and just use the last-best values, just in case.
     */
    
    dist_cutoff = node_ptr_array[best_i]->avg_length;
    if (root_node->avg_length == 0.0)
        dist_cutoff = node_ptr_array[best_i]->num_members;

    score = flag_clusters_multi(root_node,
                                min_tiny_fraction, 5,
                                max_size_fraction, dist_cutoff,
                                target_n_clusters, return_n_clusters,
                                opt, 1, print_flag);

    return score;
}


void output_clusters_multi(struct  tree_node *root_node,
                           int32_t print_flag, char *outfile_name)
{
    FILE              *outfile                    = NULL;
    char              *buffer_out                 = NULL;

    struct tree_node **node_stack                 = NULL;
    struct tree_node  *node_ptr;
    int32_t            max_stack_size             = 0;
    int32_t            stack_idx                  = 0;

    struct tree_node **leaf_ptr_array             = NULL;
    int32_t           *cluster_size_array         = NULL;
    int32_t            num_leaves                 = 0;
    int32_t            num_clusters               = 0;
    int32_t            cur_leaf                   = 0;
    int32_t            cur_cluster                = -1;
    
    /* how the cluster was stored in the leaf, not the final cluster number */
    int32_t            leaf_cluster_num           = -42;
    int32_t            leaf_cluster_flag          = -42;
    int32_t            leaf_cluster_good_flag     = -42;
    int32_t            old_leaf_cluster_num       = -42;
    int32_t            old_leaf_cluster_flag      = -42;
    int32_t            old_leaf_cluster_good_flag = -42;


    /* open file for writing */    
    if (outfile_name)
    {
        /* open as text, so it will translate EOL automatically */
        outfile = fopen(outfile_name, "wt");

        if (!outfile)
        {
            fprintf(stderr, "ERROR -- can't open output cluster file %s\n",
                    outfile_name);

            return;
        }
    }
    /* write to stdout */
    else
    {
        outfile = stdout;
    }

    buffer_out = (char *) malloc(1048576 * sizeof(char));
    setvbuf(outfile, buffer_out, _IOFBF, 1048576);


    /* allocate leaf pointer array */
    if (root_node && root_node->num_members)
    {
        leaf_ptr_array =
            (struct tree_node **) calloc(root_node->num_members,
                                         sizeof(struct tree_node *));
        num_leaves = root_node->num_members;
    }
    /* no tree, exit */
    else
    {
        return;
    }

    /* push root node onto stack */
    node_stack = (struct tree_node **) malloc(sizeof(struct tree_node *));
    node_stack[0] = root_node;
    root_node->walk_child_idx = -1;
    max_stack_size = 1;
    
    while (stack_idx >= 0)
    {
        node_ptr = node_stack[stack_idx];
        
        /* trees can come from two routes at the moment:
         *   hierarchical clustering function using only left/right ptrs
         *   read from file, which uses child_ptrs
         *
         * detect which one and fill in the other
         */
        /* binary tree, left/right ptrs */
        if (node_ptr->left_ptr && node_ptr->right_ptr &&
            node_ptr->child_ptrs == NULL)
        {
            node_ptr->num_children  = 2;
            node_ptr->child_ptrs    = calloc(2, sizeof(struct tree_node *));
            node_ptr->child_ptrs[0] = node_ptr->left_ptr;
            node_ptr->child_ptrs[1] = node_ptr->right_ptr;
        }
        /* tree read in from file, fill in left/right with first two */
        else if (node_ptr->child_ptrs && node_ptr->num_children >= 2 &&
                 node_ptr->left_ptr == NULL && node_ptr->right_ptr == NULL)
        {
            node_ptr->left_ptr  = node_ptr->child_ptrs[0];
            node_ptr->right_ptr = node_ptr->child_ptrs[1];
        }
        /* must be a leaf node */
        else if (node_ptr->child_ptrs == NULL &&
                 node_ptr->left_ptr   == NULL &&
                 node_ptr->right_ptr  == NULL)
        {
            node_ptr->num_children = 0;
            node_ptr->num_members  = 1;
            node_ptr->clevel       = 0;
        }
        
        
        /* leaf node information */
        if (node_ptr->num_children == 0)
        {
            leaf_cluster_num       = node_ptr->cluster_num;
            leaf_cluster_flag      = node_ptr->cluster_flag;
            leaf_cluster_good_flag = node_ptr->cluster_good_flag;

            
            /* new cluster */
#if MERGE_UNCLUSTERED
            if (leaf_cluster_num != old_leaf_cluster_num &&
                (leaf_cluster_good_flag || old_leaf_cluster_good_flag))
#else
            if (leaf_cluster_num  != old_leaf_cluster_num ||
                leaf_cluster_flag != old_leaf_cluster_flag)
#endif
            {
                num_clusters++;
                cur_cluster++;
                
                cluster_size_array = (int32_t *) realloc(cluster_size_array,
                                     num_clusters * sizeof(int32_t));
                cluster_size_array[cur_cluster] = 0;
            }
            

            cluster_size_array[cur_cluster]++;
            leaf_ptr_array[cur_leaf++] = node_ptr;

            /* overwrite original cluster num with new cluster num */
            node_ptr->cluster_num = cur_cluster;

            
            old_leaf_cluster_num  = leaf_cluster_num;
            old_leaf_cluster_flag = leaf_cluster_flag;
            old_leaf_cluster_good_flag = leaf_cluster_good_flag;
        }

        /* advance to next child */
        if (node_ptr->num_children &&
            node_ptr->walk_child_idx < node_ptr->num_children - 1)
        {
            /* store which child we're about to walk down */
            node_ptr->walk_child_idx++;

            /* push new node */
            stack_idx++;
            if (stack_idx >= max_stack_size)
            {
                max_stack_size = stack_idx + 1;
                node_stack = realloc(node_stack,
                            max_stack_size * sizeof(struct tree_node *));
            }
            node_stack[stack_idx] =
                node_ptr->child_ptrs[node_ptr->walk_child_idx];
            node_ptr->child_ptrs[node_ptr->walk_child_idx]->walk_child_idx=-1;
        }

        /* back up a node on the stack */
        else
        {
            stack_idx--;
        }
    }


    if (print_flag)
    {
        fprintf(outfile, "%s\t%s\t%s\n",
            "Leaf", "Cluster", "ClusterSize");
    }

    for (cur_leaf = 0; cur_leaf < num_leaves; cur_leaf++)
    {
        node_ptr    = leaf_ptr_array[cur_leaf];
        cur_cluster = node_ptr->cluster_num;

        if (print_flag)
        {
            /* name clusters using letters */
            if (0 && num_clusters <= 26)
            {
                fprintf(outfile, "%s\t%c\t%d\n",
                       node_ptr->name,
                       'A' + cur_cluster,
                       cluster_size_array[cur_cluster]);
            }
            /* name clusters using numbers, since we have too many */
            else
            {
                fprintf(outfile, "%s\tC%0*d\t%d\n",
                       node_ptr->name,
                       (int32_t) log10(num_clusters) + 1,
                       cur_cluster + 1,
                       cluster_size_array[cur_cluster]);
            }
        }
    }
    

    if (outfile_name)       fclose(outfile);

    if (node_stack)         free(node_stack);
    if (leaf_ptr_array)     free(leaf_ptr_array);
    if (cluster_size_array) free(cluster_size_array);
}


void jumble_matrix(double ***ret_jumb_matrix,
                   char   ***ret_jumb_names,
                   double  **dist_matrix,
                   char    **name_array,
                   int32_t   n)
{
    double   **jumb_matrix = *ret_jumb_matrix;
    char     **jumb_names  = *ret_jumb_names;
    uint32_t  *jumb_order  = NULL;
    int32_t    row, col;
    int32_t    i, j;
    
    jumb_order = (uint32_t *) calloc(n, sizeof(uint32_t));
    for (i = 0; i < n; i++)
        jumb_order[i] = i;

    
    jumb_matrix = *ret_jumb_matrix;

    /* allocate matrix */
    if (jumb_matrix == NULL)
    {
        jumb_matrix = (double **) calloc(n, sizeof(double *));
    }
    /* free old matrix rows */
    else
    {
        for (i = 0; i < n; i++)
        {
            if (jumb_matrix[i])
                free(jumb_matrix[i]);
        }
    }
    /* allocate new matrix rows */
    for (i = 0; i < n; i++)
        jumb_matrix[i] = (double *) malloc((i + 1) * sizeof(double));


    jumb_names = *ret_jumb_names;
    if (jumb_names == NULL)
        jumb_names = (char **) calloc(n, sizeof(char *));

    /* Fisher-Yates-Durstenfeld shuffle */
    for (i = n-1; i; --i)
    {
        uint32_t temp;
    
        /* j = [0,i]: 0 <= j <= i */
        j = (int32_t) ((i+1) * dxoshiro256p());

        temp          = jumb_order[j];
        jumb_order[j] = jumb_order[i];
        jumb_order[i] = temp;
    }

    /* populate jumbled matrix */
    for (i = 0; i < n; i++)
    {
        row = jumb_order[i];

        for (j = 0; j < i; j++)
        {
            col = jumb_order[j];
            
            if (col <= row)
                jumb_matrix[i][j] = dist_matrix[row][col];
            else
                jumb_matrix[i][j] = dist_matrix[col][row];
        }
        jumb_matrix[i][i] = 0.0;
    }

    /* populate jumbled names */
    for (i = 0; i < n; i++)
        jumb_names[i] = name_array[jumb_order[i]];

    *ret_jumb_matrix = jumb_matrix;
    *ret_jumb_names  = jumb_names;

    if (jumb_order)
        free(jumb_order);
}


void copy_matrix(double ***ret_jumb_matrix,
                 char   ***ret_jumb_names,
                 double  **dist_matrix,
                 char    **name_array,
                 uint32_t   n)
{
    double   **jumb_matrix = *ret_jumb_matrix;
    char     **jumb_names  = *ret_jumb_names;
    int32_t    i, j;
    
    jumb_matrix = *ret_jumb_matrix;

    /* allocate matrix */
    if (jumb_matrix == NULL)
    {
        jumb_matrix = (double **) calloc(n, sizeof(double *));
    }
    /* free old matrix rows */
    else
    {
        for (i = 0; i < n; i++)
        {
            if (jumb_matrix[i])
                free(jumb_matrix[i]);
        }
    }
    /* allocate new matrix rows */
    for (i = 0; i < n; i++)
        jumb_matrix[i] = (double *) malloc((i + 1) * sizeof(double));


    jumb_names = *ret_jumb_names;
    if (jumb_names == NULL)
        jumb_names = (char **) calloc(n, sizeof(char *));

    /* populate jumbled matrix */
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < i; j++)
        {
            jumb_matrix[i][j] = dist_matrix[i][j];
        }
        jumb_matrix[i][i] = 0.0;
    }

    /* populate jumbled names */
    for (i = 0; i < n; i++)
        jumb_names[i] = name_array[i];

    *ret_jumb_matrix = jumb_matrix;
    *ret_jumb_names  = jumb_names;
}


#if 0
/* not needed for this program, as the tree building process handles this */
void sum_node_num_members(struct tree_node **node_ptr_array,
                          int32_t num_nodes)
{
    struct tree_node *node_ptr;
    int32_t i, j;
    
    for (i = 0; i < num_nodes; i++)
    {
        node_ptr              = node_ptr_array[i];
        node_ptr->num_members = 0;
        
        for (j = 0; j < node_ptr->num_children; j++)
        {
            if (node_ptr->child_ptrs[j]->clevel == 0)
            {
                node_ptr->num_members++;
            }
            else
            {
                node_ptr->num_members += node_ptr->child_ptrs[j]->num_members;
            }
        }
    }
}
#endif


/* Return values:
 *
 *   -1  failure, error
 *    0  no tree-related options specified
 *    1  success, tree generated
 */
int32_t launch_tree(double **dist_matrix, char **name_array,
                    int32_t num_leaves,
                    int32_t argc, char *argv[],
                    struct options *opt)
{
    double          **jumb_matrix = NULL;
    char            **jumb_names  = NULL;
    char             *filename    = NULL;
    char             *newick_str  = NULL;
    struct tree_node *root_node   = NULL;
    struct tree_node *clust_pool  = NULL;
    int32_t           n_clusters  = 0;
    int32_t           i, j;

    struct tree_node **node_ptr_array = NULL;
    struct tree_node **leaf_ptr_array = NULL;
    int32_t            num_nodes;
    
    int32_t method_flag            = UPGMA;
    int32_t arg_jumble_flag        = 0;
    int32_t arg_seed_flag          = 0;

    uint64_t random_seed = 0;
    int32_t  num_jumbles = 0;

    /* start and end times for system clock
     * we're using 1 second resolution due to ANSI-compatibility
     * we could use clock_gettime() for high-precision, but that is POSIX
     */
    time_t start_time, end_time;
    double delta_time;


    /* default to UPGMA if no linkage method specified */
    if (opt->linkage_method)
        method_flag = opt->linkage_method;


    /* TODO -- handle this in hcdist.c where it belongs,
     *         rather than this left-over hack from
     *         merging two separate programs together...
     */
    if (argc > 1)
    {
        for (i = 1; i < argc; i++)
        {
            if (strncmp(argv[i], "--", 2) == 0)
            {
            }
            /* take only the first filename */
            else if (filename == NULL)
            {
                filename = argv[i];
            }
            /* already have the filename, assume it is number of jumbles */
            else if (arg_jumble_flag == 0 && is_all_digits(argv[i]))
            {
                num_jumbles = atol(argv[i]);
                arg_jumble_flag = 1;
            }
            /* already specified number of jumbles, treat as RNG seed */
            else if (arg_seed_flag == 0   && is_all_digits(argv[i]))
            {
                random_seed = atol(argv[i]);
                arg_seed_flag = 1;
            }
        }
    }
    

    /* if seed is zero, use current system time */
    /* if --ties-random and no seed given, use current system time */
    if (random_seed == 0.0)
    {
        random_seed = time(NULL);
    }
    
    /* print PRNG seed to STDERR */
    if ((arg_jumble_flag && random_seed) ||
        opt->ties_random_flag)
    {
        fprintf(stderr, "PRNG seed:\t%ld\n", random_seed);
    }


    initialize_xoshiro256(random_seed);

    if (opt->dists_from_file_flag)
      dist_matrix = read_distance_matrix(filename, &name_array, &num_leaves);


    if (num_jumbles == 0)
    {
        start_time = time(NULL);

        clust_pool = build_hctree(dist_matrix, name_array,
                                  num_leaves, method_flag,
                                  0, opt->ties_fewest_flag,
                                  opt->ties_random_flag,
                                  opt->threads);

        end_time   = time(NULL);
        delta_time = difftime(end_time, start_time);

        fprintf(stderr, "Time spent building trees:\t%.0lfs\n", delta_time);

        /* root_node is at the end of the cluster nodes array */
        root_node = &clust_pool[2 * (num_leaves - 1)];

        
        /* clean up after tree building, fill final arrays properly */
        clean_tree(root_node, &node_ptr_array, &num_nodes,
                   &leaf_ptr_array, &num_leaves);
        bless_tree(root_node);

#if 0
        /* handle number of descendants properly if we didn't earlier */
        /* count the number of leaves below each node */
        qsort(node_ptr_array, num_nodes, sizeof(struct tree_node *),
              cmp_node_ptr_clevel);
        sum_node_num_members(node_ptr_array, num_nodes);

        clean_tree(root_node, &node_ptr_array, &num_nodes,
                   &leaf_ptr_array, &num_leaves);
#endif

        /* default to flip tree on edge */
        if (opt->tree_flip_size_flag == 0 && opt->tree_flip_edge_flag == 0)
        {
            /* HACK -- set the unset command line option */
            opt->tree_flip_edge_flag = 1;
        }
        
        /* HACK -- sort of size first, before sorting on edge length */
        if (opt->tree_flip_edge_flag)
        {
            /* flip_nodes_multi() is coded to sort size before sorting edges;
             * enabling size effectively breaks edge length ties using size
             */
            opt->tree_flip_size_flag = 1;
        }

        /* de-exponentiate tree distances */
        if (opt->dist_pow && opt->depower_flag)
            depower_tree(node_ptr_array, num_nodes, opt);

        /* flip nodes to reorder on size or edge length */
        if (opt->tree_flip_size_flag || opt->tree_flip_edge_flag)
            flip_nodes_multi(node_ptr_array, num_nodes, opt);


        /* clobber node/leaf names with safe strings for output */
        bless_tree_names(root_node, node_ptr_array, num_nodes,
                         leaf_ptr_array, num_leaves);


        /* only print tree if we aren't printing tree weights */
        if (opt->tree_weights_flag == 0)
        {
            /* MUST run before clustering, since it creates node labels
             * TODO -- fix this so that it happens elsewhere
             */
            newick_str = create_newick_string_multi(root_node,
                                                    opt->one_line_flag, 1);
            scan_clusters(root_node, node_ptr_array, num_nodes,
                          0.025, 5, 1.0, 1.0,
                          opt->target_n_clusters, &n_clusters, opt, 1);

            /* output clusters */
            if (opt->clusters_flag)
            {
                output_clusters_multi(root_node, 1, NULL);
            }
            /* output tree */
            else
            {
                /* re-create output string, with cluster nodes flagged */
                if (newick_str) free(newick_str);
                newick_str = create_newick_string_multi(root_node,
                                                        opt->one_line_flag, 0);
                if (newick_str)
                    printf("%s", newick_str);
            }
        }


        /* print effective row count */
#if WEIGHT_TYPE == 1
        fprintf(stderr, "Effective count:\t%0.8g\n",
                calc_tree_weights_gsc(node_ptr_array, num_nodes,
                leaf_ptr_array, num_leaves,
                opt->tree_weights_flag));
#elif WEIGHT_TYPE == 2
        fprintf(stderr, "Effective count:\t%0.8g\n",
                calc_tree_weights_experimental(node_ptr_array, num_nodes,
                leaf_ptr_array, num_leaves,
                opt->tree_weights_flag));
#else
        fprintf(stderr, "Effective count:\t%0.8g\n",
                calc_tree_weights_clustalw(node_ptr_array, num_nodes,
                leaf_ptr_array, num_leaves,
                opt->tree_weights_flag));
#endif


#if 1
        /* free everything in the tree */
        if (newick_str)
            free(newick_str);

        if (leaf_ptr_array)
            free(leaf_ptr_array);

        if (node_ptr_array)
        {
            for (i = 0; i < num_nodes; i++)
            {
                if (node_ptr_array[i])
                {
                    if (node_ptr_array[i]->name)
                        free(node_ptr_array[i]->name);
                    if (node_ptr_array[i]->child_ptrs)
                        free(node_ptr_array[i]->child_ptrs);
                }
            }

            free(node_ptr_array);
        }
#else
        free_tree_stuff(newick_str, node_ptr_array, num_nodes,
                        leaf_ptr_array, num_leaves, 0);
#endif

        if (clust_pool)
            free(clust_pool);

        newick_str     = NULL;
        leaf_ptr_array = NULL;
        node_ptr_array = NULL;
        clust_pool     = NULL;
    }
    else
    {
      for (j = 0; j < num_jumbles; j++)
      {
        if (opt->ties_random_flag == 0)
        {
            jumble_matrix(&jumb_matrix, &jumb_names,
                          dist_matrix, name_array,
                          num_leaves);
        }
        /* break ties randomly, no need to actually jumble */
        else
        {
            copy_matrix(&jumb_matrix, &jumb_names,
                        dist_matrix, name_array,
                        num_leaves);
        }

        start_time = time(NULL);

        clust_pool = build_hctree(jumb_matrix, jumb_names,
                                  num_leaves, method_flag,
                                  j+1, opt->ties_fewest_flag,
                                  opt->ties_random_flag,
                                  opt->threads);

        end_time   = time(NULL);
        delta_time = end_time - start_time;
        
        fprintf(stderr, "Time spent building trees:\t%.0lfs\n", delta_time);

        /* root_node is at the end of the cluster nodes array */
        root_node = &clust_pool[2 * (num_leaves - 1)];

        /* clean up after tree building, fill final arrays properly */
        clean_tree(root_node, &node_ptr_array, &num_nodes,
                   &leaf_ptr_array, &num_leaves);
        bless_tree(root_node);

#if 0
        /* handle number of descendants properly if we didn't earlier */
        /* count the number of leaves below each node */
        qsort(node_ptr_array, num_nodes, sizeof(struct tree_node *),
              cmp_node_ptr_clevel);
        sum_node_num_members(node_ptr_array, num_nodes);

        clean_tree(root_node, &node_ptr_array, &num_nodes,
                   &leaf_ptr_array, &num_leaves);
#endif

        /* default to flip tree on edge */
        if (opt->tree_flip_size_flag == 0 && opt->tree_flip_edge_flag == 0)
        {
            /* HACK -- set the unset command line option */
            opt->tree_flip_edge_flag = 1;
        }

        /* de-exponentiate tree distances */
        if (opt->dist_pow && opt->depower_flag)
            depower_tree(node_ptr_array, num_nodes, opt);

        /* flip nodes to reorder on size or edge length */
        if (opt->tree_flip_size_flag || opt->tree_flip_edge_flag)
            flip_nodes_multi(node_ptr_array, num_nodes, opt);


        /* clobber node/leaf names with safe strings for output */
        bless_tree_names(root_node, node_ptr_array, num_nodes,
                         leaf_ptr_array, num_leaves);


        /* MUST run before clustering, since it creates node labels
         * TODO -- fix this so that it happens elsewhere
         */
        newick_str = create_newick_string_multi(root_node,
                                                opt->one_line_flag, 1);

        /* only one tree, so we'll allow all usual output */
        if (num_jumbles == 1)
        {
            scan_clusters(root_node, node_ptr_array, num_nodes,
                          0.025, 5, 1.0, 1.0,
                          opt->target_n_clusters, &n_clusters, opt, 1);

            /* output clusters */
            if (opt->clusters_flag)
            {
                output_clusters_multi(root_node, 1, NULL);
            }
            /* output tree */
            else
            {
                /* re-create output string, with cluster nodes flagged */
                if (newick_str) free(newick_str);
                newick_str = create_newick_string_multi(root_node,
                                                        opt->one_line_flag, 0);
                if (newick_str)
                    printf("%s", newick_str);
            }
        }
        else
        {
            if (newick_str)
                printf("%s", newick_str);
        }

        /* print effective row count */
#if WEIGHT_TYPE == 1
        fprintf(stderr, "Effective count:\t%0.8g\n",
                calc_tree_weights_gsc(node_ptr_array, num_nodes,
                                      leaf_ptr_array, num_leaves,
                                      opt->tree_weights_flag));
#elif WEIGHT_TYPE == 2
        fprintf(stderr, "Effective count:\t%0.8g\n",
                calc_tree_weights_experimental(node_ptr_array, num_nodes,
                                      leaf_ptr_array, num_leaves,
                                      opt->tree_weights_flag));
#else
        fprintf(stderr, "Effective count:\t%0.8g\n",
                calc_tree_weights_clustalw(node_ptr_array, num_nodes,
                                      leaf_ptr_array, num_leaves,
                                      opt->tree_weights_flag));
#endif

        /* free everything in the tree */
        if (newick_str)
            free(newick_str);

        if (leaf_ptr_array)
            free(leaf_ptr_array);

        if (node_ptr_array)
        {
            for (i = 0; i < num_nodes; i++)
            {
                if (node_ptr_array[i])
                {
                    if (node_ptr_array[i]->name)
                        free(node_ptr_array[i]->name);
                    if (node_ptr_array[i]->child_ptrs)
                        free(node_ptr_array[i]->child_ptrs);
                }
            }

            free(node_ptr_array);
        }


        if (clust_pool)
            free(clust_pool);

        newick_str     = NULL;
        node_ptr_array = NULL;
        leaf_ptr_array = NULL;
        clust_pool     = NULL;
      }
    }


    /* free distance matrix */
    if (dist_matrix)
    {
        for (i = 0; i < num_leaves; i++)
            if (dist_matrix[i])
                free(dist_matrix[i]);
        
        free(dist_matrix);
    }

    /* free row names */
    if (name_array)
    {
        for (i = 0; i < num_leaves; i++)
            if (name_array[i])
                free(name_array[i]);
        
        free(name_array);
    }


    /* free jumbled distance matrix */
    if (jumb_matrix)
    {
        for (i = 0; i < num_leaves; i++)
            if (jumb_matrix[i])
                free(jumb_matrix[i]);
        
        free(jumb_matrix);
    }

    /* free jumbled name array */
    if (jumb_names)
        free(jumb_names);


    return 0;
}


char * read_in_tree_string(char *tree_file_name)
{
    FILE *in_tree;
    char *string = NULL;
    char *tree_string = NULL;
    char *sptr;
    char *dptr;
    int32_t max_string_len = 0;

    int32_t max_tree_string_len = 0;
    int32_t temp_tree_string_len = 0;
    int32_t tree_string_len = 0;


    if (tree_file_name == NULL || strcmp(tree_file_name, "-") == 0)
    {
        in_tree = stdin;
    }
    else
    {
        in_tree = fopen(tree_file_name, "rb");

        if (!in_tree)
        {
            printf("Can't open treefile %s\n", tree_file_name);
            return NULL;
        }
    }


    tree_string = calloc(1, sizeof(char));

    /* read in the tree string */
    while (fgets_strip_realloc(&string, &max_string_len, in_tree) != NULL)
    {
        /* max_string_len alredy takes into account extra \0 at the end */
        temp_tree_string_len = tree_string_len +
                               max_string_len;
        if (temp_tree_string_len > max_tree_string_len)
        {
            /* overallocate by 1% for a *MASSIVE* speedup in valgrind
             *  due to saving on the number of reallocs
             */
            max_tree_string_len = temp_tree_string_len * 1.01;
            tree_string = (char *) realloc(tree_string,
                                           max_tree_string_len * sizeof(char));
        }
        
        sptr = string;
        dptr = tree_string + tree_string_len;

        /* concatenate new line, skipping whitespace */
        while(*sptr)
        {
            if (!isspace(*sptr))
                *dptr++ = *sptr;
            sptr++;
        }
        tree_string_len = dptr - tree_string;
        *dptr++ = '\0';
    }

    /* trim it down to size, since we overallocated */
    tree_string = (char*) realloc(tree_string,
                                  (tree_string_len + 1) * sizeof(char));

    fclose(in_tree);
    if (string)
        free(string);
    
    return tree_string;
}


struct tree_node * create_tree_from_string(char *tree_string,
                                  struct tree_node ***return_node_ptr_array,
                                  int32_t *return_num_nodes,
                                  struct tree_node ***return_leaf_ptr_array,
                                  int32_t *return_num_leaves)
{
    struct tree_node **node_ptr_array = NULL;
    struct tree_node **leaf_ptr_array = NULL;

    struct tree_node *tree_root;
    struct tree_node *node_ptr;
    struct tree_node *parent_ptr;
    char *sptr, *sptr2;
    char *temp_string = NULL;
    char *str_clust   = "*";    /* I use this to denote clusters */
    char *sptr3, *sptr4;
    int32_t temp_length;
    int32_t num_children;
    char paren_state = ' ';

    int32_t num_nodes  = 1;
    int32_t num_leaves = 0;
    int32_t i, j;

    tree_root = calloc(1, sizeof(struct tree_node));
    /* tree_root->name = strdup("root"); */

    node_ptr_array = (struct tree_node **) realloc(node_ptr_array,
                         num_nodes * sizeof(struct tree_node *));
    node_ptr_array[num_nodes - 1] = tree_root;
    
    node_ptr = tree_root;
    parent_ptr = tree_root;
    sptr = tree_string;
    while (*sptr && *sptr != ';')
    {
        /* skip whitespace */
        if (isspace(*sptr))
        {
            sptr++;
            continue;
        }

        paren_state = *sptr;
    
        /* add a new child and sibling nodes */
        if (*sptr == '(' || *sptr == ',')
        {
            /* drop down a level, going to add a new level of children */
            if (*sptr == '(')
            {
                parent_ptr = node_ptr;
            }
            else
            {
                parent_ptr = node_ptr->parent_ptr;
            }

            if (parent_ptr->child_ptrs == NULL)
            {
                parent_ptr->num_children = 0;
            }
            parent_ptr->num_children++;
            
            num_children = parent_ptr->num_children;

            parent_ptr->child_ptrs = (struct tree_node **)
                                      realloc(parent_ptr->child_ptrs,
                                              num_children *
                                              sizeof(struct tree_node *));
            node_ptr = parent_ptr->child_ptrs[num_children - 1] =
                (struct tree_node *) calloc(1, sizeof(struct tree_node));

            node_ptr->parent_ptr   = parent_ptr;
            node_ptr->num_children = 0;
        }

        /* back up a level */
        if (*sptr == ')')
        {
            node_ptr = node_ptr->parent_ptr;

            /* add new nodes to node list */
            if (node_ptr != tree_root)
            {
                num_nodes++;
                node_ptr_array = (struct tree_node **) realloc(node_ptr_array,
                                    num_nodes * sizeof(struct tree_node *));
                node_ptr_array[num_nodes - 1] = node_ptr;
                
            }
        }

        sptr++;
        while(isspace(*sptr))
            sptr++;
        sptr2 = sptr;

#ifdef DEBUG_TREE_INPUT
        printf ("PAREN_STATE %c\n", paren_state);
#endif

        /* read names, add siblings, edge length, etc. */
        while (*sptr2 && *sptr2 != ';' &&
               (*sptr2 != ',' && *sptr2 != '(' && *sptr2 != ')' &&
                *sptr2 != ':'))
        {
            sptr2++;
        }
        if (sptr2 - sptr > 0)
        {
            if (temp_string)
            {
                free(temp_string);
            }
            temp_length = sptr2 - sptr;
            temp_string = (char *) malloc((temp_length + 1) * sizeof(char));
            strncpy(temp_string, sptr, temp_length);
            temp_string[temp_length] = '\0';
            
            /* edge length */
            if (paren_state == ':')
            {
                node_ptr->edge_length = atof(temp_string);

#ifdef DEBUG_TREE_INPUT
                printf("EDGE_LENGTH %s %.14g\n",
                    temp_string,
                    node_ptr->edge_length);
#endif
            }
            /* names */
            else
            {
#ifdef DEBUG_TREE_INPUT
                /* first leaf name */
                if (paren_state == '(')
                {
                    printf("FIRST_LEAF_NAME %s\n", temp_string);
                }
                else if (paren_state == ',')
                {
                    printf("MORE_LEAF_NAMES %s\n", temp_string);
                }
#endif
                
                if (paren_state == '(' || paren_state == ',')
                {
                    if (node_ptr->name)
                        free(node_ptr->name);
                    node_ptr->name = strdup(temp_string);
                    
                    num_leaves++;
                    leaf_ptr_array = (struct tree_node **)
                                      realloc(leaf_ptr_array,
                                      num_leaves * sizeof(struct tree_node *));
                    leaf_ptr_array[num_leaves - 1] = node_ptr;
                }

                /* node name, given after a ) */
                else if (paren_state == ')')
                {
#ifdef DEBUG_TREE_INPUT
                    printf("NODE_NAME %s\n", temp_string);
#endif
                    if (node_ptr->name)
                        free(node_ptr->name);
                    
                    /* strip leading cluster characters */
                    sptr3 = temp_string;
                    sptr4 = strstr(sptr3, str_clust);
                    while (sptr4)
                    {
                       sptr3 = sptr4 + 1;
                       sptr4 = strstr(sptr3, str_clust);
                    }
                    
                    node_ptr->name = strdup(sptr3);
                }
            }
        }

        sptr = sptr2;
    }


    /* calculate various tree properties */
    /* *MUST* do this before sorting nodes */
    calc_node_clevels(leaf_ptr_array, num_leaves);

    /* sum number of members */
    /* must first sort nodes by increasing clevel */
    qsort(node_ptr_array, num_nodes, sizeof(struct tree_node *),
          cmp_node_ptr_clevel);
    for (i = 0; i < num_nodes; i++)
    {
        node_ptr              = node_ptr_array[i];
        node_ptr->num_members = 0;
        
        for (j = 0; j < node_ptr->num_children; j++)
        {
            if (node_ptr->child_ptrs[j]->num_children == 0)
            {
                node_ptr->child_ptrs[j]->num_members = 1;
            }

            node_ptr->num_members += node_ptr->child_ptrs[j]->num_members;
        }
    }

    calc_node_avg_lengths(node_ptr_array, num_nodes,
                          leaf_ptr_array, num_leaves);
    
#if 0
    /* sort children */
    for (i = 0; i < num_nodes; i++)
    {
        node_ptr = node_ptr_array[i];
        
        if (node_ptr->num_children)
        {
            qsort(node_ptr->child_ptrs, node_ptr->num_children,
                  sizeof(struct tree_node *), cmp_node_ptr);
        }
    }
#endif

fprintf(stderr, "ROOT\t%d\n", tree_root->num_members);

    fprintf(stderr, "NumLeaves\t%d\n", num_leaves);
    fprintf(stderr, "NumNodes\t%d\n", num_nodes);

    if (temp_string)
        free(temp_string);

    *return_node_ptr_array = node_ptr_array;
    *return_leaf_ptr_array = leaf_ptr_array;
    *return_num_nodes = num_nodes;
    *return_num_leaves = num_leaves;

    return(tree_root);
}


void bless_tree(struct tree_node *tree_root)
{
    struct tree_node **node_stack = NULL;
    struct tree_node *node_ptr;
    int32_t max_stack_size = 0;
    int32_t stack_idx = 0;
    int32_t num_not_leaves = 0;
    int32_t new_len;
    int32_t node_index = 0;

    /* push root node onto stack */
    node_stack = (struct tree_node **) malloc(sizeof(struct tree_node *));
    node_stack[0] = tree_root;
    tree_root->walk_child_idx = -1;
    max_stack_size = 1;
    
    while (stack_idx >= 0)
    {
        node_ptr = node_stack[stack_idx];
        
        /* trees can come from two routes at the moment:
         *   hierarchical clustering function using only left/right ptrs
         *   read from file, which uses child_ptrs
         *
         * detect which one and fill in the other
         */
        /* binary tree, left/right ptrs */
        if (node_ptr->left_ptr && node_ptr->right_ptr &&
            node_ptr->child_ptrs == NULL)
        {
            node_ptr->num_children  = 2;
            node_ptr->child_ptrs    = calloc(2, sizeof(struct tree_node *));
            node_ptr->child_ptrs[0] = node_ptr->left_ptr;
            node_ptr->child_ptrs[1] = node_ptr->right_ptr;
        }
        /* tree read in from file, fill in left/right with first two */
        else if (node_ptr->child_ptrs && node_ptr->num_children >= 2 &&
                 node_ptr->left_ptr == NULL && node_ptr->right_ptr == NULL)
        {
            node_ptr->left_ptr  = node_ptr->child_ptrs[0];
            node_ptr->right_ptr = node_ptr->child_ptrs[1];
        }
        /* must be a leaf node */
        else if (node_ptr->child_ptrs == NULL &&
                 node_ptr->left_ptr   == NULL &&
                 node_ptr->right_ptr  == NULL)
        {
            node_ptr->num_children = 0;
            node_ptr->num_members  = 1;
            node_ptr->clevel       = 0;
        }
        
        /* deal with new node the first time we come to it */
        if (node_ptr->walk_child_idx == -1)
        {
            /* not a leaf */
            if (node_ptr->num_children)
            {
                num_not_leaves++;
            }
        }
        
        /* advance to next child */
        if (node_ptr->num_children &&
            node_ptr->walk_child_idx < node_ptr->num_children - 1)
        {
            /* store which child we're about to walk down */
            node_ptr->walk_child_idx++;

            /* push new node */
            stack_idx++;
            if (stack_idx >= max_stack_size)
            {
                max_stack_size = stack_idx + 1;
                node_stack = realloc(node_stack,
                            max_stack_size * sizeof(struct tree_node *));
            }
            node_stack[stack_idx] =
                node_ptr->child_ptrs[node_ptr->walk_child_idx];
            node_ptr->child_ptrs[node_ptr->walk_child_idx]->walk_child_idx=-1;
        }

        /* back up a node on the stack */
        else
        {
            /* close branch and back up */
            if (node_ptr->parent_ptr &&
                node_ptr->parent_ptr->walk_child_idx ==
                 node_ptr->parent_ptr->num_children - 1)
            {
                node_index++;

                /* fill empty names with node number */
                if (node_ptr->parent_ptr->name == NULL)
                {
                    /* allocate enough room for the name string */
                    new_len = snprintf(NULL, 0, "Node%05d", node_index);

                    node_ptr->parent_ptr->name =
                        malloc((new_len+1) * sizeof(char));
                    node_ptr->parent_ptr->name[new_len] = '\0';

                    sprintf(node_ptr->parent_ptr->name,
                        "Node%05d", node_index);
                }
            }

            stack_idx--;
        }
    }
    
    if (node_stack) free(node_stack);
}


void calc_node_clevels(struct tree_node **leaf_ptr_array, int32_t num_leaves)
{
    struct tree_node *node_ptr;
    int32_t i;
    int32_t clevel;
    
    for (i = 0; i < num_leaves; i++)
    {
        node_ptr = leaf_ptr_array[i];
        
        /* walk up tree */
        while (node_ptr->parent_ptr && node_ptr->parent_ptr != node_ptr)
        {
            /* new trial clevel of parent node */
            clevel = node_ptr->clevel + 1;
            
            node_ptr = node_ptr->parent_ptr;

            /* keep new clevel if higher */
            if (clevel > node_ptr->clevel)
            {
                node_ptr->clevel = clevel;
            }
            /* node is already higher, abort walk */
            else
            {
                break;
            }
        }
    }
}


/* fill in average distance to leaves */
/* call this after blessing a newly read in tree and calculating clevels */
void calc_node_avg_lengths(struct tree_node **node_ptr_array,
                           int32_t num_nodes,
                           struct tree_node **leaf_ptr_array,
                           int32_t num_leaves)
{
    struct tree_node *node_ptr, *node_ptr2;
    int32_t i, j;
    int32_t count;

    /* sort nodes by clevel, so that children are before parents */
    qsort(node_ptr_array, num_nodes, sizeof(struct tree_node *),
          cmp_node_ptr_clevel);

    /* initialize leaves */
    for (i = 0; i < num_leaves; i++)
    {
        node_ptr = leaf_ptr_array[i];
        node_ptr->avg_length = 0.0;
    }

    for (i = 0; i < num_nodes; i++)
    {
        node_ptr = node_ptr_array[i];
        
        /* average child distance to leaves */
        node_ptr->avg_length = 0.0;
        count                = 0;
        for (j = 0; j < node_ptr->num_children; j++)
        {
            node_ptr2 = node_ptr->child_ptrs[j];

#if 1
            /* give each child branch equal weight */
            node_ptr->avg_length += node_ptr2->edge_length +
                                    node_ptr2->avg_length;
            count++;
#else
            /* weight by number of leaves within each child branch */
            node_ptr->avg_length += node_ptr2->num_members *
                                    (node_ptr2->edge_length +
                                     node_ptr2->avg_length);
            count += node_ptr2->num_members;
#endif
        }
        if (count)
            node_ptr->avg_length /= count;
    }
}


/* Gerstein-Sonnhammer-Chothia tree weights
 *
 * Gerstein M, Sonnhammer ELL, Chothia C;
 * Volume Changes in Protein Evolution;
 * J. Mol. Biol. (1994) 236. 1067-1078
 *
 * Ward's/Centroid/Median tend to not work so well for effective counts
 * and/or weights.  This can be demonstrated by N random vectors resulting in
 * effective counts noticably smaller or larger than N.  They also have
 * trouble on the dupe/variants effective counts example trees.
 * Additionally, Centroid/Median can produce inversions, which may cause
 * further strangeness.  So, don't use Ward's/Centroid/Median for tree
 * weights.
 *
 * UPGMA (default) is recommended for tree weights.
 *
 * Nicola De Maio et al. 2021 BMC Bioinformatics 2021 compare the effects
 * of several different weighting schemes on the accuracy of their
 * particular use case.  For ultrametric trees, GSC generally performs well.
 * For very non-ultrametric trees, GSC isn't so great.
 *
 *      Leaf      CLUSTALW       GSC       GSC-floored
 *    --------   ----------   ----------   -----------
 *    B_dupe1    0.24573131   0.24594504   0.16396334    
 *    B_dupe2    0.24573131   0.24594504   0.16396334
 *    B_dupe3    0.24573131   0.24594504   0.16396334
 *    B_var1     0.24573135   0.24594504   0.49189002
 *    A_dupe1    0.24573135   0.30948831   0.22576900
 *    A_dupe2    0.24573135   0.30948831   0.22576900
 *    A_var1     0.24573136   0.18240185   0.26612109
 *    A_var2     0.24573136   0.18240185   0.26612109
 *    D_dupe1    0.16494015   0.16465524   0.098793134
 *    D_dupe2    0.16494015   0.16465524   0.098793134
 *    D_dupe5    0.16494015   0.16465524   0.098793134
 *    D_dupe3    0.16494015   0.16465524   0.098793134
 *    D_dupe4    0.16494015   0.16465524   0.098793134
 *    D_var1     0.16494016   0.16465524   0.49396567
 *    C_var2     0.20000003   0.12723290   0.17723717
 *    C_var3     0.20000003   0.12723290   0.17723717
 *    C_var1     0.20000004   0.18127020   0.25251188
 *    C_dupe1    0.20000003   0.28213213   0.19650697
 *    C_dupe2    0.20000003   0.28213213   0.19650697
 *    --------   ----------   ----------   -----------
 *    EFFCOUNT   3.9554917    3.9554922    3.9554917
 *
 *    GSC can overemphasise the effects of VERY tiny differences between
 *    leaves.  In the above example, in a perfect weighting scheme,
 *    all the dupes for each letter would be down-weighted to sum to roughly
 *    the same weight as each individual _var#, which are ever so slightly
 *    different from the dupes (but slightly more different from the dupes
 *    than from each other).  GSC can exacerbate the problem, resuting in
 *    the dupes being even more heavily overweighted vs. the variants.
 *    CLUSTALW doesn't condense the dupes relative to the variants, so it
 *    does not achieve our ideal weighting, but a least it doesn't potentially
 *    make things worse.  Modifying GSC to floor weights to DBL_MIN, but only
 *    where needed to deal with dupes, may be a reasonable compromise?
 *
 *    GSC, as originally published, adds the full current edge weight whenever
 *    the current leaf weight is zero.  This can result in unexpected relative
 *    weights between sets of duplicate nodes and their close neighbors (see
 *    example table above), since the duplicate nodes aren't down-weighted
 *    relative to close neighbors -- they all receive equal weight.
 *    Additionally, highly non-ultrametric trees can suffer from a similar,
 *    somewhat related problem, due to overly-low weighting of branches with
 *    shorter average node-to-leaf distances relative to neighboring branches.
 *
 *    I have implemented two hacks to attempt to deal with these issues,
 *    which aim to result in both more intuitively correct weights, as well
 *    as more consistent weights between highly similar trees differing only
 *    in zero vs. short edge lengths:
 *
 *    The first hack floors all edge lengths to DBL_MIN, the smallest non-
 *    zero value representable by double precision floating point, and
 *    modifies the calculations to appropriately take this into account.  This
 *    fixes the incorrect relative weighting of duplicate nodes vs. their
 *    close neighbors (see example table above).  However, it actually
 *    exacerbates the problem for highly non-ultrametric trees where the
 *    leaf edges are very small (say, 1E-10), but not exactly zero.  This is
 *    dealt with by the 2nd hack, which up-weights leaves that appear to have
 *    been overly down-weighted due to asymmetry in the non-ultrametric tree.
 *    This fixes the zero / near-zero issues with non-ultrametric trees, and
 *    appears to benefit non-ultrametric trees in general -- at least in my toy
 *    examples I tested it on.  It could very well break some non-ultrametric
 *    trees worse than they were already, but I don't work much with non-
 *    ultrametric trees, and this software generates only ultrametric trees
 *    anyways, so improving the non-ultrametric hack isn't high on my list of
 *    priorities.
 */
double calc_tree_weights_gsc(struct tree_node **node_ptr_array,
                             int32_t num_nodes,
                             struct tree_node **leaf_ptr_array,
                             int32_t num_leaves,
                             int32_t print_flag)
{
    struct tree_node **node_stack     = NULL;
    double           *leaf_root_dists = NULL;
    struct tree_node *node_ptr, *node_ptr2, *node_ptr_wsum;
    int32_t max_stack_size = 0;
    int32_t stack_idx = 0;

    double  weight_sum, weight_inc;
    double  max_weight, max_dist, weight_scale, eff_count;
    double  edge_length, weight_frac;
#if GSC_ASYMMETRY_HACK
    double  avg_length, avg_length_child, max_avg;
    double  asymmetry, adjust;
    int32_t nmemb_asym;
#endif
    int32_t i;

    /* initialize leaf-to-root distances */
    leaf_root_dists = (double *) calloc(num_leaves, sizeof(double));

    /* initialize stack */
    node_stack = (struct tree_node **) malloc(sizeof(struct tree_node *));
    max_stack_size = 1;

    /* sort nodes by clevel, so that children are before parents */
    qsort(node_ptr_array, num_nodes, sizeof(struct tree_node *),
          cmp_node_ptr_clevel);

    /* initialize leaf weights */
    for (i = 0; i < num_leaves; i++)
    {
        node_ptr         = leaf_ptr_array[i];
        node_ptr->weight = 0.0;
    }
    
    /* walk up the nodes, walking back down to the leaves each time */
    for (i = 0; i < num_nodes; i++)
    {
        node_ptr = node_ptr_array[i];

        /* push parent node onto stack */
        node_stack[0] = node_ptr;
        stack_idx     = 0;
        node_ptr->walk_child_idx = -1;

        /* GCC warns it might be uninitialized unless we set it here,
         * but, due to how the tree is traversed, it should be fine if we
         * didn't.
         */
        node_ptr_wsum = node_ptr;
        
        while (stack_idx >= 0)
        {
            node_ptr2 = node_stack[stack_idx];

            /* first child node */
            if (node_ptr2->parent_ptr == node_ptr)
            {
                node_ptr_wsum = node_ptr2;

                /* initialize sum of subtree weights */
                if (node_ptr2->walk_child_idx == -1)
                {
                    node_ptr_wsum->weight = 0.0;

#if GSC_ZERO_HACK_currently_not_used
                    node_ptr_wsum->num_zero_leaves = 0;
#endif
                }
            }


            /* add leaf weight to subtree sum */
            if (node_ptr2->num_children == 0)
            {
                node_ptr_wsum->weight += node_ptr2->weight;

#if GSC_ZERO_HACK_currently_not_used
                if (node_ptr2->edge_length < DBL_MIN)
                    node_ptr_wsum->num_zero_leaves++;
#endif
            }
            
            /* advance to next child */
            if (node_ptr2->num_children &&
                node_ptr2->walk_child_idx < node_ptr2->num_children - 1)
            {
                /* store which child we're about to walk down */
                node_ptr2->walk_child_idx++;

                /* push new node */
                stack_idx++;
                if (stack_idx >= max_stack_size)
                {
                    max_stack_size = stack_idx + 1;
                    node_stack = realloc(node_stack,
                                max_stack_size * sizeof(struct tree_node *));
                }
                node_stack[stack_idx] =
                    node_ptr2->child_ptrs[node_ptr2->walk_child_idx];
                node_ptr2->child_ptrs[node_ptr2->walk_child_idx]->walk_child_idx=-1;
            }
            /* back up a node on the stack */
            else
            {
                stack_idx--;
            }
        }


        /* increment leaf weights */

        /* push parent node onto stack */
        node_stack[0] = node_ptr;
        stack_idx     = 0;
        node_ptr->walk_child_idx = -1;

        edge_length   = node_ptr->edge_length;
#if GSC_ASYMMETRY_HACK
        asymmetry     = 0.0;
        nmemb_asym    = node_ptr->num_members;
#endif

        while (stack_idx >= 0)
        {
            node_ptr2 = node_stack[stack_idx];

            /* first child node */
            if (node_ptr2->parent_ptr == node_ptr)
            {
                edge_length   = node_ptr2->edge_length;
                node_ptr_wsum = node_ptr2;
            }

#if GSC_ASYMMETRY_HACK
            /* calculate a measure of how asymmetric the node is */
            if (node_ptr2->parent_ptr == node_ptr_wsum)
            {
                asymmetry  = 0.0;
                nmemb_asym = node_ptr2->num_members;

                if (node_ptr_wsum)
                {
                    /* average length of this branch */
                    avg_length_child  = node_ptr2->edge_length +
                                        node_ptr2->avg_length;

#if 1
                    /* average length of sibling branches */
                    /* works better when adjust further scaled by asymmetry */
                    avg_length        = (node_ptr_wsum->avg_length *
                                         node_ptr_wsum->num_members -
                                         avg_length_child *
                                         node_ptr2->num_members) /
                                        (node_ptr_wsum->num_members -
                                         node_ptr2->num_members);
#else
                    /* works better on test_eff4_upgma_near-zero.tree ?? */
                    avg_length        = node_ptr_wsum->avg_length;
#endif
                    
                    /* absolute value, in case of negative edges */
                    /* negative edges may still cause weird behavior... */
                    max_avg = fabs(avg_length);
                    if (fabs(avg_length_child) > fabs(max_avg))
                        max_avg = fabs(avg_length_child);

                    if (max_avg)
                    {
                        asymmetry = (avg_length - avg_length_child) /
                                    max_avg;
                    }
                }
            }
#endif

            /* leaf node */
            if (node_ptr2->num_children == 0)
            {
                weight_inc  = edge_length;
                
                /* add full edge length for 0 weights, per GSC 1994 */
                weight_frac = 1.0;

                if (node_ptr2->weight)
                {
                    weight_frac = fabs(node_ptr2->weight /
                                       node_ptr_wsum->weight);
                    weight_inc  = edge_length * weight_frac;

#if GSC_ZERO_HACK
                    /* add full edge length for ~0 weights */
                    if (weight_frac <
                        DBL_MIN * node_ptr_wsum->num_members)
                    {
                        weight_inc  = edge_length;
                        /* weight_frac = 1.0; */
                    }
#endif
                }

#if GSC_ZERO_HACK
                /* initialize weights to DBL_MIN */
                if (fabs(weight_inc) < DBL_MIN &&
                    fabs(node_ptr2->weight) < DBL_MIN)
                {
                     weight_inc = DBL_MIN;
                }
#endif

#if GSC_ASYMMETRY_HACK
                /* Adjust weight_inc based on amount of asymmetry.
                 *
                 * Helps with cases where one edge is near, but not
                 * near enough, to zero, which causes major failures
                 * in the weighting algorithm.
                 *
                 * Still need to tweak the conditional so that it might
                 * benefit other asymmetric branches without making them
                 * worse.
                 */
                if (weight_inc != edge_length && asymmetry > 0)
                {
#if 0
                    /* adjust shorter leaves more than longer ones */
                    adjust = asymmetry *
                             (asymmetry / nmemb_asym - weight_frac);
#else
                    /* distribute asymmetry amongst shorter branch leaves */
                    /* unsure whether ^3 or ^2 is better */
                    adjust = asymmetry * asymmetry * asymmetry / nmemb_asym;
#endif

                    /* do not allow negative adjustments */
                    if (adjust < 0)
                        adjust = 0.0;
                    
                    /* do not sum to more than 1.0 */
                    if (weight_frac + adjust > 1.0)
                        adjust = 1.0 - weight_frac;

                    /* undo standard zero-weight GSC rule */
                    weight_inc  = edge_length * weight_frac;

                    /* add in new adjustment amount */
                    weight_inc += edge_length * adjust;


                    if (adjust && adjust > 1E-10)
                        fprintf(stderr, "ASYM  %s   %s\t%g   %g   %g\n",
                            node_ptr_wsum->name, node_ptr2->name,
                            asymmetry,
                            weight_frac, adjust);
                }
#endif

#if 0
                if (print_flag)
                {
                    printf("WeightIter   %s   %.8g   %s   %.8g   %.8g   %.8g\n",
                        node_ptr->name,
                        edge_length,
                        node_ptr2->name,
                        node_ptr2->weight,
                        node_ptr_wsum->weight,
                        weight_inc);
                }
#endif

                node_ptr2->weight += weight_inc;
            }
            
            /* advance to next child */
            if (node_ptr2->num_children &&
                node_ptr2->walk_child_idx < node_ptr2->num_children - 1)
            {
                /* store which child we're about to walk down */
                node_ptr2->walk_child_idx++;

                /* push new node */
                stack_idx++;
                if (stack_idx >= max_stack_size)
                {
                    max_stack_size = stack_idx + 1;
                    node_stack = realloc(node_stack,
                                max_stack_size * sizeof(struct tree_node *));
                }
                node_stack[stack_idx] =
                    node_ptr2->child_ptrs[node_ptr2->walk_child_idx];
                node_ptr2->child_ptrs[node_ptr2->walk_child_idx]->walk_child_idx=-1;
            }
            /* back up a node on the stack */
            else
            {
#if GSC_ASYMMETRY_HACK
                /* reset asymmetry in case we're too low in the tree */
                /* this is necessary, and does matter sometimes */
                if (node_ptr2->parent_ptr == node_ptr_wsum &&
                    node_ptr2->num_children &&
                    node_ptr2->walk_child_idx >= node_ptr2->num_children - 1)
                {
                    asymmetry = 0.0;
                }
#endif

                stack_idx--;
            }
        }
    }
    
    /* walk up tree from leaves to sum leaf-to-root distances */
    max_dist = 0.0;
    for (i = 0; i < num_leaves; i++)
    {
        node_ptr2 = leaf_ptr_array[i];
        node_ptr  = node_ptr2->parent_ptr;
        
        leaf_root_dists[i] += node_ptr2->edge_length;

        node_ptr2 = node_ptr;

        while (node_ptr2)
        {
            leaf_root_dists[i] += node_ptr2->edge_length;
            node_ptr2           = node_ptr2->parent_ptr;
        }
    }
    

    max_weight = 0.0;
    weight_sum = 0.0;
    max_dist   = 0.0;
    for (i = 0; i < num_leaves; i++)
    {
        /* HACK -- deal with all-zero trees */
        if (leaf_ptr_array[i]->weight < DBL_MIN)
            leaf_ptr_array[i]->weight = DBL_MIN;

        weight_inc  = leaf_ptr_array[i]->weight;
        weight_sum += weight_inc;

        if (weight_inc > max_weight)
            max_weight = weight_inc;

        if (leaf_root_dists[i] > max_dist)
            max_dist   = leaf_root_dists[i];
    }

    
    /* NOTE -- we have a problem with nearly identical vectors
     *
     * For example: (A,B,C) (X,Y,Z)
     *
     * ABC are nearly identical, XYZ are nearly identical, ABC and XYZ differ
     *
     * The effective count should be ~2, not ~6
     *
     * (((A2:2.74971e-05,
     *    A3:2.74971e-05):1.37488e-05,
     *    A1:4.12459e-05):0.367125,
     * ((B2:2.75335e-05,
     *   B3:2.75335e-05):1.37669e-05,
     *   B1:4.13003e-05):0.367125):0;
     *
     * normalizing by average root-to-leaf length appears to work well
     */

    /* If we have a non-ultrametric tree, such as from neighbor-joining,
     * the distance between nodes averaged between child nodes as we walk
     * up the tree isn't going to be the same as the maximum root-to-leaf
     * distance.  For ultrametric trees, they should be the same.
     */
#if 0
    weight_scale = node_ptr_array[num_nodes-1]->avg_length;
#else
    weight_scale = max_dist;
#endif

    /* edge case for 100% identical */
    if (weight_scale < DBL_MIN)
        weight_scale = weight_sum;

    /* normalize weights */
    for (i = 0; i < num_leaves; i++)
    {
        node_ptr          = leaf_ptr_array[i];
        node_ptr->weight /= weight_scale;
    }
    eff_count = weight_sum / weight_scale;

    if (print_flag)
    {
        printf("%s\t%s\n", "Leaf", "Weight");
    
        for (i = 0; i < num_leaves; i++)
        {
            node_ptr = leaf_ptr_array[i];
        
            printf("%s\t%.14g\n",
                   node_ptr->name, node_ptr->weight);
        }
    }

    if (leaf_root_dists)
        free(leaf_root_dists);
    
    if (node_stack) free(node_stack);
    
    return eff_count;
}


/* basically GSC, flooring leaf weights to zero, optimized for no recursion
 *
 * __NOT__ suitable for non-ultrametric trees!!
 */
/* return effective number of leaves */
double calc_tree_weights_experimental(struct tree_node **node_ptr_array,
                               int32_t num_nodes,
                               struct tree_node **leaf_ptr_array,
                               int32_t num_leaves, int32_t print_flag)
{
    struct tree_node *node_ptr, *node_ptr2, *node_ptr_wsum;
    double           *leaf_weights    = NULL;
    double           *leaf_root_dists = NULL;

    double  weight_scale, weight_sum, eff_count;
    double  max_weight, max_dist;
    double  weight_inc, weight_frac;
    double  edge_length;
#if GSC_ASYMMETRY_HACK
    double  avg_length, avg_length_child, max_avg;
    double  asymmetry, adjust;
    int32_t nmemb_asym;
#endif
    int32_t i, j;

    /* initialize leaf weights and root-to-leaf distances */
    leaf_weights    = (double *) calloc(num_leaves, sizeof(double));
    leaf_root_dists = (double *) calloc(num_leaves, sizeof(double));

    /* sort nodes by clevel, so that children are before parents */
    qsort(node_ptr_array, num_nodes, sizeof(struct tree_node *),
          cmp_node_ptr_clevel);
    
    /* initialize leaf edge sums */
    for (i = 0; i < num_leaves; i++)
    {
        node_ptr = leaf_ptr_array[i];
        node_ptr->weight = 0.0;
    }

    /* initialize edge sums, store them in the node weights */
    for (i = 0; i < num_nodes; i++)
    {
        node_ptr = node_ptr_array[i];
        node_ptr->weight = 0.0;

        for (j = 0; j < node_ptr->num_children; j++)
        {
            node_ptr_wsum = node_ptr->child_ptrs[j];
            edge_length   = node_ptr_wsum->edge_length;

            /* avoid divide by zero later, floor leaves to DBL_MIN */
            /* since these are leaves, there cannot be any inversions yet */
            if (node_ptr_wsum->num_children == 0 && edge_length < DBL_MIN)
                edge_length = DBL_MIN;

            node_ptr->weight += node_ptr_wsum->weight + edge_length;
        }
    }

    /* walk up tree to sum weights */
    for (i = 0; i < num_leaves; i++)
    {
        node_ptr2   = leaf_ptr_array[i];
        edge_length = node_ptr2->edge_length;
        node_ptr    = node_ptr2->parent_ptr;

        /* must be >= DBL_MIN in order for all-zero trees to weight properly */
        if (edge_length < DBL_MIN)
            edge_length = DBL_MIN;

        /* initialize leaf weights */
        leaf_weights[i]    = edge_length;
        leaf_root_dists[i] = edge_length;

        
        node_ptr_wsum = node_ptr;
        
        /* walk up tree */
        while (node_ptr_wsum)
        {
            edge_length = node_ptr_wsum->edge_length;
            node_ptr    = node_ptr_wsum->parent_ptr;

#if GSC_ASYMMETRY_HACK
            asymmetry  = 0.0;
            nmemb_asym = node_ptr2->num_members;

            if (node_ptr)
            {
                avg_length        = node_ptr_wsum->avg_length;
                avg_length_child  = node_ptr2->edge_length +
                                    node_ptr2->avg_length;
                
                /* in case negative edges cause unexpected behavior */
                max_avg = fabs(avg_length);
                if (fabs(avg_length_child) > fabs(max_avg))
                    max_avg = fabs(avg_length_child);

                if (max_avg)
                {
                    asymmetry = (avg_length - avg_length_child) /
                                max_avg;
                }
            }
#endif

#if 0
            /* set inversions to zero */
            if (edge_length < 0)
                edge_length = 0;
#endif

            weight_inc = edge_length;

            /* add full edge length for 0 weights, per GSC 1994 */
            weight_frac = 1.0;
            
            if (leaf_weights[i])
            {
                weight_frac = fabs(leaf_weights[i] / node_ptr_wsum->weight);
                weight_inc  = edge_length * weight_frac;

                /* add full edge length for ~0 weights */
                if (weight_frac < DBL_MIN * node_ptr_wsum->num_members)
                    weight_inc  = edge_length;
            }

#if GSC_ASYMMETRY_HACK
            /* Adjust weight_inc based on amount of asymmetry.
             *
             * Helps with cases where one edge is near, but not
             * near enough, to zero, which causes major failures
             * in the weighting algorithm.
             *
             * Still need to tweak the conditional so that it might
             * benefit other asymmetric branches without making them
             * worse.
             */

            /* if (asymmetry) */
            /* if (asymmetry > 0 && weight_frac <= asymmetry) */
            if (weight_inc != edge_length &&
                asymmetry > 0 &&
                weight_frac < asymmetry)
            {
#if 0
                adjust = fabs(asymmetry) *
                         (asymmetry / node_ptr_wsum->num_members -
                          weight_frac);

#else
                /* adjust=fabs(asymmetry) * (asymmetry - weight_frac); */

                adjust = asymmetry - weight_frac;
#endif

                /* undo standard zero-weight GSC rule */
                weight_inc  = edge_length * weight_frac;

                /* add in new adjustment amount */
                weight_inc += edge_length * adjust;

                fprintf(stderr, "ASYMMETRY  %s  %s\t%f\t%g\t%f\n",
                    node_ptr_wsum->name, node_ptr2->name,
                    asymmetry,
                    weight_frac, weight_frac + adjust);
            }
#endif


            leaf_weights[i]    += weight_inc;
            leaf_root_dists[i] += edge_length;

            node_ptr2     = node_ptr_wsum;
            node_ptr_wsum = node_ptr_wsum->parent_ptr;
        }
    }

    
    /* store unnormalized weights, find max and sum*/
    max_weight = 0.0;
    weight_sum = 0.0;
    max_dist   = 0.0;
    for (i = 0; i < num_leaves; i++)
    {
        weight_sum += leaf_weights[i];

        if (leaf_weights[i] > max_weight)
            max_weight = leaf_weights[i];
        if (leaf_root_dists[i] > max_dist)
            max_dist   = leaf_root_dists[i];

        node_ptr = leaf_ptr_array[i];
        node_ptr->weight = leaf_weights[i];
    }
    
    /* NOTE -- we have a problem with nearly identical vectors
     *
     * For example: (A,B,C) (X,Y,Z)
     *
     * ABC are nearly identical, XYZ are nearly identical, ABC and XYZ differ
     *
     * The effective count should be ~2, not ~6
     *
     * (((A2:2.74971e-05,
     *    A3:2.74971e-05):1.37488e-05,
     *    A1:4.12459e-05):0.367125,
     * ((B2:2.75335e-05,
     *   B3:2.75335e-05):1.37669e-05,
     *   B1:4.13003e-05):0.367125):0;
     *
     * normalizing by average root-to-leaf length appears to work well
     */

    /* If we have a non-ultrametric tree, such as from neighbor-joining,
     * the distance between nodes averaged between child nodes as we walk
     * up the tree isn't going to be the same as the maximum root-to-leaf
     * distance.  For ultrametric trees, they should be the same.
     */
#if 0
    weight_scale = node_ptr_array[num_nodes-1]->avg_length;
#else
    weight_scale = max_dist;
#endif
    
    /* edge case for 100% identical */
    if (weight_scale < DBL_MIN)
        weight_scale = weight_sum;

    /* normalize weights */
    for (i = 0; i < num_leaves; i++)
    {
        node_ptr          = leaf_ptr_array[i];
        node_ptr->weight /= weight_scale;
    }
    eff_count = weight_sum / weight_scale;

    if (print_flag)
    {
        printf("%s\t%s\n", "Leaf", "Weight");

        for (i = 0; i < num_leaves; i++)
        {
            node_ptr = leaf_ptr_array[i];
        
            printf("%s\t%.14g\n",
                   node_ptr->name, node_ptr->weight);
        }
    }

    if (leaf_weights)
        free(leaf_weights);
    if (leaf_root_dists)
        free(leaf_root_dists);
    
    return eff_count;
}


/* Thompson-Higgins-Gibson (CLUSTALW) tree weights
 *
 * Thompson JD, Higgins DG, Gibson TJ
 * Improved sensitivity of profile searches through the use of sequence
 *  weights and gap excision
 * CABIOS (1994) 10(1) 19-29
 *
 * Analogy: apply voltage to root, zero potential to leaves, weights
 *          are the currents flowing from the leaves.
 */
double calc_tree_weights_clustalw(struct tree_node **node_ptr_array,
                         int32_t num_nodes,
                         struct tree_node **leaf_ptr_array,
                         int32_t num_leaves, int32_t print_flag)
{
    struct tree_node *node_ptr;
    double           *leaf_weights    = NULL;
    double           *leaf_root_dists = NULL;

    double  weight_scale, weight_sum, eff_count;
    double  max_weight, max_dist;
    double  temp_edge;
    int32_t i;

    /* initialize leaf weights and root-to-leaf distances */
    leaf_weights    = (double *) calloc(num_leaves, sizeof(double));
    leaf_root_dists = (double *) calloc(num_leaves, sizeof(double));

    /* sort nodes by clevel, so that children are before parents */
    qsort(node_ptr_array, num_nodes, sizeof(struct tree_node *),
          cmp_node_ptr_clevel);
    
    /* walk up tree to sum weights */
    for (i = 0; i < num_leaves; i++)
    {
        node_ptr = leaf_ptr_array[i];

        /* initialize leaf weights */
        leaf_weights[i]    = node_ptr->edge_length;
        leaf_root_dists[i] = node_ptr->edge_length;
        
#if 0
        /* leaves need to be non-zero, in case entire tree is zeroes */
        /* since these are leaves, there cannot be any inversions yet */
        if (leaf_weights[i] < DBL_MIN)
            leaf_weights[i] = DBL_MIN;
#endif
        
        /* walk up tree */
        node_ptr = node_ptr->parent_ptr;
        while (node_ptr)
        {
            temp_edge           = node_ptr->edge_length;

#if 0
            /* set inversions to zero */
            if (temp_edge < 0)
                temp_edge = 0;
#endif

            leaf_weights[i]    += temp_edge / node_ptr->num_members;
            leaf_root_dists[i] += temp_edge;

            node_ptr = node_ptr->parent_ptr;
        }
    }

    
    /* store unnormalized weights, find max and sum*/
    max_weight = 0.0;
    weight_sum = 0.0;
    max_dist   = 0.0;
    for (i = 0; i < num_leaves; i++)
    {
        /* HACK -- deal with all-zero trees */
        if (leaf_weights[i] < DBL_MIN)
            leaf_weights[i] = DBL_MIN;

        weight_sum += leaf_weights[i];

        if (leaf_weights[i] > max_weight)
            max_weight = leaf_weights[i];
        if (leaf_root_dists[i] > max_dist)
            max_dist   = leaf_root_dists[i];

        node_ptr = leaf_ptr_array[i];
        node_ptr->weight = leaf_weights[i];
    }
    
    /* NOTE -- we have a problem with nearly identical vectors
     *
     * For example: (A,B,C) (X,Y,Z)
     *
     * ABC are nearly identical, XYZ are nearly identical, ABC and XYZ differ
     *
     * The effective count should be ~2, not ~6
     *
     * (((A2:2.74971e-05,
     *    A3:2.74971e-05):1.37488e-05,
     *    A1:4.12459e-05):0.367125,
     * ((B2:2.75335e-05,
     *   B3:2.75335e-05):1.37669e-05,
     *   B1:4.13003e-05):0.367125):0;
     *
     * normalizing by average root-to-leaf length appears to work well
     */

    /* If we have a non-ultrametric tree, such as from neighbor-joining,
     * the distance between nodes averaged between child nodes as we walk
     * up the tree isn't going to be the same as the maximum root-to-leaf
     * distance.  For ultrametric trees, they should be the same.
     */
#if 0
    weight_scale = node_ptr_array[num_nodes-1]->avg_length;
#else
    weight_scale = max_dist;
#endif
    
    /* edge case for 100% identical */
    if (weight_scale < DBL_MIN)
        weight_scale = weight_sum;

    /* normalize weights */
    for (i = 0; i < num_leaves; i++)
    {
        node_ptr          = leaf_ptr_array[i];
        node_ptr->weight /= weight_scale;
    }
    eff_count = weight_sum / weight_scale;

    if (print_flag)
    {
        printf("%s\t%s\n", "Leaf", "Weight");
    
        for (i = 0; i < num_leaves; i++)
        {
            node_ptr = leaf_ptr_array[i];
        
            printf("%s\t%.14g\n",
                   node_ptr->name, node_ptr->weight);
        }
    }

    if (leaf_weights)
        free(leaf_weights);
    if (leaf_root_dists)
        free(leaf_root_dists);

    return eff_count;
}
