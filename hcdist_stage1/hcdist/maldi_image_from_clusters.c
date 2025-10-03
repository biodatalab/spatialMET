#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <stdint.h>
#ifdef __GLIBC__
#include <malloc.h>
#endif
#include "cet_colors.h"
#include "FreeImage.h"
#include "text.h"
#include "hcdist.h"

#define MALDI_MISSING DBL_MAX

/* put these here to make valgrind happy with tree.c,
 * even though I'm not calling any pthreaded functions in this program
 */
pthread_mutex_t mutex_wait_main;
pthread_cond_t  cond_wait_main;
int             n_threads_left;


/* Ported from original F90 code written by John Burkardt in 1998.
 *
 * Reference:
 *    Foley, van Dam, Feiner, and Hughes,
 *    Computer Graphics, Principles and Practice,
 *    Addison Wesley, Second Edition, 1990.
 */
double hls_value(double n1, double n2, double hue)
{
    if (hue < 0.0)
        hue += 360.0;
    else if (hue >= 360.0)
        hue -= 360.0;

    if (hue < 60.0)
        return n1 + (n2 - n1) * hue / 60.0;
    else if(hue < 180.0)
        return n2;
    else if (hue < 240.0)
        return n1 + (n2 - n1) * (240.0 - hue) / 60.0;

    return n1;
}


/* Ported from original F90 code written by John Burkardt in 1998.
 *
 * Reference:
 *    Foley, van Dam, Feiner, and Hughes,
 *    Computer Graphics, Principles and Practice,
 *    Addison Wesley, Second Edition, 1990.
 *
 * R,G,B range from 0.0 - 1.0
 * H:0 - 360, L/S: 0.0 - 1.0
 */
void hls_to_rgb(double h, double l, double s,
                double *r, double *g, double *b)
{
    double m1, m2;

    if (l <= 0.5)
        m2 = l + l * s;
    else
        m2 = l + s - l * s;

    m1 = l + l - m2;

    if (s == 0.0)
        *r = *g = *b = l;
    else
    {
        *r = hls_value(m1, m2, h + 120.0);
        *g = hls_value(m1, m2, h);
        *b = hls_value(m1, m2, h - 120.0);
    }
}


/* Ported from original F90 code written by John Burkardt in 1998.
 *
 * Reference:
 *    Foley, van Dam, Feiner, and Hughes,
 *    Computer Graphics, Principles and Practice,
 *    Addison Wesley, Second Edition, 1990.
 *
 * R,G,B range from 0.0 - 1.0
 * H:0 - 360, L/S: 0.0 - 1.0
 */
void rgb_to_hls(double r, double g, double b,
	        double *h, double *l, double *s)
{
    double bc, gc, rc;
    double rgbmax, rgbmin;

    /* Lightness */
    rgbmax = r;
    if (g > rgbmax) rgbmax = g;
    if (b > rgbmax) rgbmax = b;
    rgbmin = r;
    if (g < rgbmin) rgbmin = g;
    if (b < rgbmin) rgbmin = b;
    *l = (rgbmax + rgbmin) * 0.5;

    /* Saturation */
    if (rgbmax == rgbmin)
        *s = 0.0;
    else
    {
        if (*l <= 0.5)
            *s = (rgbmax - rgbmin) / (rgbmax + rgbmin);
        else
            *s = (rgbmax - rgbmin) / (2.0 - rgbmax - rgbmin);
    }

    /* Hue */
    if (rgbmax == rgbmin)
        *h = 0.0;
    else
    {
        rc = (rgbmax - r) / (rgbmax - rgbmin);
        gc = (rgbmax - g) / (rgbmax - rgbmin);
        bc = (rgbmax - b) / (rgbmax - rgbmin);

        if (r == rgbmax)
            *h = bc - gc;
        else if (g == rgbmax)
            *h = 2.0 + rc - bc;
        else
            *h = 4.0 + gc - rc;

        *h *= 60.0;

        if (*h < 0.0)
            *h += 360.0;
        else if (*h >= 360.0)
            *h -= 360.0;
    }
}


int cmp_sort_double_ptr(const void *vptr1, const void *vptr2)
{
    double val1 = **(double **) vptr1;
    double val2 = **(double **) vptr2;
    
    if (val1 < val2) return -1;
    if (val1 > val2) return 1;
    return 0;
}


int cmp_sort_int32(const void *vptr1, const void *vptr2)
{
    int32_t value1 = *(int32_t *) vptr1;
    int32_t value2 = *(int32_t *) vptr2;
    
    if (value1 < value2) return -1;
    if (value1 > value2) return 1;
    return 0;
}

int cmp_sort_float(const void *vptr1, const void *vptr2)
{
    float value1 = *(float *) vptr1;
    float value2 = *(float *) vptr2;
    
    if (value1 < value2) return -1;
    if (value1 > value2) return 1;
    return 0;
}


int cmp_sort_string(const void *vptr1, const void *vptr2)
{
    char *data1_ptr = *(char **) vptr1;
    char *data2_ptr = *(char **) vptr2;
    
    return strcmp(data1_ptr, data2_ptr);
}


/* sort node_data_map structures (array of structures, not pointers) */
int cmp_sort_mapping_by_node(const void *vptr1, const void *vptr2)
{
    struct node_data_map *map1_ptr = (struct node_data_map *) vptr1;
    struct node_data_map *map2_ptr = (struct node_data_map *) vptr2;

    /* sort by node index */
    if (map1_ptr->leaf_index < map2_ptr->leaf_index) return -1;
    if (map1_ptr->leaf_index > map2_ptr->leaf_index) return  1;

    /* sort by original allocated order */
    if (vptr1 < vptr2) return -1;
    if (vptr1 > vptr2) return  1;
    
    return 0;
}


/* sort node_data_map structures (array of structures, not pointers) */
int cmp_sort_mapping_by_data(const void *vptr1, const void *vptr2)
{
    struct node_data_map *map1_ptr = (struct node_data_map *) vptr1;
    struct node_data_map *map2_ptr = (struct node_data_map *) vptr2;

    /* sort by node index */
    if (map1_ptr->data_index < map2_ptr->data_index) return -1;
    if (map1_ptr->data_index > map2_ptr->data_index) return  1;

    /* sort by original allocated order */
    if (vptr1 < vptr2) return -1;
    if (vptr1 > vptr2) return  1;
    
    return 0;
}


void log2_data(double **data_matrix, int32_t num_rows, int32_t num_cols)
{
    double *row_ptr;
    double  value;
    int32_t row, col;
    
    for (row = 0; row < num_rows; row++)
    {
        row_ptr = data_matrix[row];

        for (col = 0; col < num_cols; col++)
        {
            value = row_ptr[col];

            /* leave MISSING and bad spots as-is */
            if (value == MALDI_MISSING || value >= DBL_MAX)
                continue;
            
            /* uh oh, negative value or zero, set to missing */
            if (value <= 0.0)
            {
                row_ptr[col] = MALDI_MISSING;
            }
            /* take the log2 */
            else
            {
                row_ptr[col] = log2(value);
            }
        }
    }
}


void fill_log2_row_stats(double **data_matrix,
                         int32_t num_rows, int32_t num_cols,
                         double *row_min_array,
                         double *row_max_array,
                         double *row_mean_array,
                         double *row_sd_array)
{
    double *value_array = NULL;
    double *row_ptr;
    double  value;
    int32_t row, col;
    int32_t n, i;
    
    double  min, max, mean, sd;
    
    /* allocate array for the values, so we don't need to re-log them, etc. */
    value_array = (double *) malloc(num_cols * sizeof(double));
    
    for (row = 0; row < num_rows; row++)
    {
        row_ptr = data_matrix[row];
        
        min  = DBL_MAX;
        max  = DBL_MIN;
        mean = 0.0;
        sd   = 0.0;
        n    = 0;

        for (col = 0; col < num_cols; col++)
        {
            value = row_ptr[col];

            /* leave MISSING and bad spots as-is */
            if (value == MALDI_MISSING || value >= DBL_MAX)
                continue;
            
            /* uh oh, negative value or zero, set to missing */
            if (value <= 0.0)
            {
                row_ptr[col] = MALDI_MISSING;
                continue;
            }

            /* store the log2 value */
            value_array[n++] = log2(value);
        }
        
        for (i = 0; i < n; i++)
        {
            value = value_array[i];
            
            if (value < min) min = value;
            if (value > max) max = value;
            
            mean += value;
        }
        
        if (n)
            mean /= n;
        
        for (i = 0; i < n; i++)
        {
            value  = value_array[i];
            value  = mean - value;
            sd    += value * value;
        }
        
        if (n)
            sd = sqrt(sd / n);
        
        row_min_array[row]  = min;
        row_max_array[row]  = max;
        row_mean_array[row] = mean;
        row_sd_array[row]   = sd;
    }
    
    if (value_array)
        free(value_array);
}


/* This program uses the population stdev ("n"),
 * not the sample stdev ("n-1"), for standardizing to
 * unit variance, so that the resulting vector lengths = n
 * instead of n-1.  This is nice, in that cosine calculations
 * can then use RMSD^2 / 2, rather than the odd-looking
 * Euclidian^2 / 2(n-1), as a speed optimization.
 */
void mean_center_uv(double **data_matrix, int32_t num_rows, int32_t num_cols)
{
    double *row_ptr;
    double  value;
    double  avg, sd;
    int32_t row, col;
    int32_t n;
    
    for (row = 0; row < num_rows; row++)
    {
        row_ptr = data_matrix[row];
        
        avg = 0.0;
        n = 0;
        for (col = 0; col < num_cols; col++)
        {
            value = row_ptr[col];

            if (value == MALDI_MISSING)
                continue;

            /* skip bad spots */
            if (value >= DBL_MAX)
                continue;
            
            avg += value;
            n++;
        }
        
        if (n)
            avg /= n;

        sd = 0;
        n  = 0;
        for (col = 0; col < num_cols; col++)
        {
            value = row_ptr[col];

            if (value == MALDI_MISSING)
                continue;

            /* skip bad spots */
            if (value >= DBL_MAX)
                continue;
            
            value -= avg;
            sd += value * value;
            n++;
        }
        
        if (n)
            sd = sqrt(sd / n);

        /* subtract the mean */
        if (avg)
        {
            for (col = 0; col < num_cols; col++)
            {
                if (row_ptr[col] == MALDI_MISSING)
                    continue;

                if (row_ptr[col] < DBL_MAX)
                    row_ptr[col] -= avg;
            }
        }

        /* scale to unit variance */
        if (sd)
        {
            for (col = 0; col < num_cols; col++)
            {
                if (row_ptr[col] == MALDI_MISSING)
                    continue;

                if (row_ptr[col] < DBL_MAX)
                    row_ptr[col] /= sd;
            }
        }
    }
}


void map_leaves_to_data(double  **data_matrix,
                        char    **data_names,
                        struct  tree_node **leaf_ptr_array,
                        struct  node_data_map *node_map_array,
                        int32_t num_points, int32_t num_leaves)
{
    struct name_index_pair *name_index_pairs = NULL;
    struct tree_node       *node_ptr;
    struct node_data_map   *node_map_ptr;
    struct name_index_pair *pair_ptr;
    struct name_index_pair  query_pair;
    int32_t i;


    /* allocate array to use with qsort/bsearch */
    /* initialize relative to data, not relative to tree */
    name_index_pairs = calloc(num_points, sizeof(struct name_index_pair));
    
    /* initialize name:index pairs */
    for (i = 0; i < num_points; i++)
    {
        name_index_pairs[i].name  = data_names[i];
        name_index_pairs[i].index = i;
    }
    
    /* sort name:index pairs on name */
    qsort(name_index_pairs, num_points, sizeof(struct name_index_pair),
          cmp_name_index_by_name);


    /* store row or col index */
    for (i = 0; i < num_points; i++)
    {
        node_map_ptr = &node_map_array[i];

        /* store data row or col index */
        node_map_ptr->data_index = i;

        /* initialize to NULL (unmapped) */
        node_map_ptr->node_ptr   = NULL;
    }
    

    /* map to tree */
    query_pair.index = -1;
    for (i = 0; i < num_leaves; i++)
    {
        /* initialize the query */
        node_ptr        = leaf_ptr_array[i];
        query_pair.name = node_ptr->name;

        /* search for tree pointer by name */
        pair_ptr = bsearch(&query_pair, name_index_pairs, num_points,
                           sizeof(struct name_index_pair),
                           cmp_name_index_by_name);

        /* found it */
        if (pair_ptr)
        {
            node_map_ptr             = &node_map_array[pair_ptr->index];
            node_map_ptr->node_ptr   = node_ptr;
            node_map_ptr->leaf_index = i;
        }
    }


    if (name_index_pairs)
        free(name_index_pairs);
}


/* x#y#
 *
 * set x or y to -1 if unable to parse it
 */
void parse_maldi_xy(char **pixel_name_array,
                    struct node_data_map *node_map_pixel_array,
                    int32_t num_pixels)
{
    struct node_data_map *node_map_ptr;
    char                 *pixel_name;
    char                 *str_ptr;
    int32_t x, y;
    int32_t i;
    
    for (i = 0; i < num_pixels; i++)
    {
        node_map_ptr = &node_map_pixel_array[i];
        pixel_name   =  pixel_name_array[i];
        
        x = y = -1;
        
        str_ptr = strstr(pixel_name, "x");
        if (str_ptr && strlen(str_ptr) > 1)
        {
            x = atol(str_ptr + 1);
        }

        str_ptr = strstr(pixel_name, "y");
        if (str_ptr && strlen(str_ptr) > 1)
        {
            y = atol(str_ptr + 1);
        }
        
        node_map_ptr->x = x;
        node_map_ptr->y = y;
    }
}


#if 1
/* Required minimal error handling function for FreeImage */
void FreeImageErrorHandler(FREE_IMAGE_FORMAT fif, const char *message)
{
    printf("\n*** ");
    printf("%s Format\n", FreeImage_GetFormatFromFIF(fif));
    printf(message);
    printf(" ***\n");
}
#endif



/* needed by qsort() comparison function, so it must be global */
RGBQUAD *global_pal       = NULL;
int32_t  global_pal_count = 256;    /* # colors in original palette */

int cmp_sort_pal_index(const void *vptr1, const void *vptr2)
{
    BYTE index1 = *(BYTE *) vptr1;
    BYTE index2 = *(BYTE *) vptr2;
    double r_1, g_1, b_1, h_1, l_1, s_1;
    double r_2, g_2, b_2, h_2, l_2, s_2;

    /* sort unused colors last */
    if (index1 <  global_pal_count && index2 >= global_pal_count) return -1;
    if (index1 >= global_pal_count && index2 <  global_pal_count) return  1;
    
    r_1 = global_pal[index1].rgbRed   / 255.0;
    g_1 = global_pal[index1].rgbGreen / 255.0;
    b_1 = global_pal[index1].rgbBlue  / 255.0;

    r_2 = global_pal[index2].rgbRed   / 255.0;
    g_2 = global_pal[index2].rgbGreen / 255.0;
    b_2 = global_pal[index2].rgbBlue  / 255.0;

    rgb_to_hls(r_1, g_1, b_1, &h_1, &l_1, &s_1);
    rgb_to_hls(r_2, g_2, b_2, &h_2, &l_2, &s_2);

    /* sort pure grey to last */
    if (s_1 > 1E-5 && s_2 < 1E-5) return -1;
    if (s_1 < 1E-5 && s_2 > 1E-5) return  1;

    /* then in reverse hue order */
    if (h_1 > h_2)         return -1;
    if (h_1 < h_2)         return  1;
    
    /* then by lightness */
    if (l_1 > l_2)         return -1;
    if (l_1 < l_2)         return  1;
    
    /* then saturation */
    if (s_1 > s_2)         return -1;
    if (s_1 < s_2)         return  1;
    
    /* then index */
    if (index1 < index2)   return -1;
    if (index1 > index2)   return  1;

    return 0;
}


void sort_image_palette_by_hue(FIBITMAP *image, int32_t count)
{
    BYTE    orig_indices[256];
    BYTE    new_indices[256];
    int32_t i;
    RGBQUAD pal_new[256];

    global_pal = FreeImage_GetPalette(image);
    
    /* assume 0 means entire palette */
    if (count == 0)
        count = 256;

    global_pal_count = count;

    for (i = 0; i < 256; i++)
    {
        orig_indices[i] = i;
        new_indices[i]  = i;
    }

    /* sort palette to new order */
    qsort(new_indices, 256, sizeof(BYTE), cmp_sort_pal_index);

    /* reassign indices to new order */
    FreeImage_ApplyPaletteIndexMapping(image, new_indices, orig_indices,
                                       count, FALSE);

    /* reorder the palette to new order */
    for (i = 0; i < 256; i++)
    {
        memcpy(&pal_new[i], &global_pal[new_indices[i]], sizeof(RGBQUAD));
    }
    memcpy(global_pal, pal_new, 256 * sizeof(RGBQUAD));
}


void render_maldi_images(double  **data_matrix,
                         struct  node_data_map *node_map_array_pixels,
                         int32_t   num_pixels,
                         struct  node_data_map *node_map_array_mz,
                         int32_t   num_mz,
                         int32_t   n_pixel_clusters,
                         int32_t   n_mz_clusters,
                         uint32_t *pixel_cluster_colors,
                         char     *base_str,
                         int       opt_reverse_palette,
                         int       opt_mz_images,
                         int      *selected_clusters_array,
                         int       n_selected_clusters,
                         int       pixel_size)
{
    FIBITMAP *cluster_image = NULL;
    FIBITMAP *quantized    = NULL;
    struct    tree_node     *node_ptr;
    struct    node_data_map *node_map_ptr_col, *node_map_ptr_row;
    int32_t   x_max, x_min, y_max, y_min;
    uint32_t  col, col_data;
    int32_t   cluster_pixel;
    int32_t   cluster_pixel_bin;
    int32_t   missing_brightness = 0x7f;
    
    FIBITMAP *row_image       = NULL;
    double   *double_2d_array = NULL;
    double   *value_array     = NULL;
    double   *row_ptr, min_global;
    double    value, value_median, value_geomean;
    double    value_high, value_low, value_range, value_scale;
    int32_t   n_values;
    int32_t   cluster_mz;
    int32_t   n_mz_selected;
    uint32_t  row, row_data;
    uint32_t  index;
    uint64_t  size;
    
    /* arrays for storing start/end indices for m/z clusters */
    uint64_t  *mz_cluster_start_array = NULL;
    uint64_t  *mz_cluster_end_array   = NULL;
    uint64_t  start, end;

    BYTE     *src, *sptr, *colored_src;
    WORD      r_word, g_word, b_word;

    int32_t   w, h, pitch;
    int32_t   x, y;	/* x = col, y = row */
    int32_t   i;
    int32_t   x2, y2;

    FIBITMAP *contact_image = NULL;
    BYTE     *contact_src, *contact_sptr;
    int32_t   contact_w, contact_h, contact_pitch;
    int32_t   contact_x, contact_y;
    int32_t   grid_w, grid_h, grid_max;
    int32_t   grid_x, grid_y, grid_n, contact_n;
    int32_t   grid_y_old;
    double    target_aspect = 16.0 / 9.0;
    double    aspect, score, best_grid_w, best_grid_h, best_score;
    int       pixel_size_half = pixel_size >> 1;

    FILE     *outfile_mz = NULL;


    /* sprintf() returns (int), and %*d must be a short to prevent
     * bogus GCC -Wformat-truncation warnings
     */
    unsigned short print_number_width;
    int       str_len;
    char     *output_file_name = NULL;


    /* store pixel cluster colors for later */
    for (cluster_pixel = 0; cluster_pixel < n_pixel_clusters; cluster_pixel++)
    {
#if 1
        /* CET-R1 or CET-R2 */
        cluster_pixel_bin = 0;
        if (n_pixel_clusters > 1)
            cluster_pixel_bin = (int32_t) (0.5 + 255.0 *
                                           (cluster_pixel /
                                           (n_pixel_clusters - 1.0L)));
        if (opt_reverse_palette)
           cluster_pixel_bin = 255.0 - cluster_pixel_bin;

        r_word = cet_r2_palette_255[cluster_pixel_bin][0];
        g_word = cet_r2_palette_255[cluster_pixel_bin][1];
        b_word = cet_r2_palette_255[cluster_pixel_bin][2];
#else
       /* CET-C6, ~300 degrees, stop before wrap-around
        * while it looks nice it is *NOT* perceptually isoluminant,
        * since the colors between the primaries and secondaries are darker
        */
        cluster_pixel_bin = 0;
        if (n_pixel_clusters > 1)
            cluster_pixel_bin = 210.0 *
                                (cluster_pixel / (n_pixel_clusters - 1.0L));
        if (opt_reverse_palette)
            cluster_pixel_bin = 210.0 - cluster_pixel_bin;

        r_word = cet_c6_palette_255[cluster_pixel_bin][0];
        g_word = cet_c6_palette_255[cluster_pixel_bin][1];
        b_word = cet_c6_palette_255[cluster_pixel_bin][2];
#endif

        pixel_cluster_colors[cluster_pixel] =
            b_word + (g_word << 8) + (r_word << 16);
    }

    
    x_max = y_max = -42;
    x_min = y_min = 999999999;
    
    /* find min/max x,y coordinates */
    for (col = 0; col < num_pixels; col++)
    {
        node_map_ptr_col = &node_map_array_pixels[col];
        
        x = node_map_ptr_col->x;
        y = node_map_ptr_col->y;
        
        if (x > x_max)
            x_max = x;

        if (x < x_min)
            x_min = x;

        if (y > y_max)
            y_max = y;

        if (y < y_min)
            y_min = y;
    }
    

    /* sanity check make sure we have some pixels to plot */
    if (x_max < 0 || y_max < 0)
        return;

    /* crop blank image borders */
    w = x_max - x_min + 1;
    h = y_max - y_min + 1;


    cluster_image = FreeImage_Allocate(w, h, 24,
                    FI_RGBA_RED_MASK, FI_RGBA_GREEN_MASK, FI_RGBA_BLUE_MASK);
    src           = FreeImage_GetBits(cluster_image);
    pitch         = FreeImage_GetPitch(cluster_image);


    /* fill image with grey to indicate missing values */
    for (y = 0; y < h; y++)
    {
        for (x = 0; x < w; x++)
        {
            /* calculate the address of the pixel to write */
            sptr = src + y * pitch + 3*x;

            sptr[FI_RGBA_RED]   = missing_brightness;
            sptr[FI_RGBA_GREEN] = missing_brightness;
            sptr[FI_RGBA_BLUE]  = missing_brightness;
        }
    }

    /* color pixel clusters */
    for (col = 0; col < num_pixels; col++)
    {
        node_map_ptr_col = &node_map_array_pixels[col];
        node_ptr         = node_map_ptr_col->node_ptr;
        
        /* uh oh, this pixel wasn't mapped */
        if (node_ptr == NULL)
        {
            fprintf(stderr, "WARNING -- unmapped col %d\n", col);
            continue;
        }
        
        x = node_map_ptr_col->x - x_min;
        y = node_map_ptr_col->y - y_min;
        
        cluster_pixel = node_ptr->cluster_num;

        r_word = (pixel_cluster_colors[cluster_pixel] & 0xFF0000) >> 16;
        g_word = (pixel_cluster_colors[cluster_pixel] & 0x00FF00) >> 8;
        b_word = (pixel_cluster_colors[cluster_pixel] & 0x0000FF);

        if (pixel_size == 1)
        {
            /* calculate the address of the pixel to write */
            sptr = src + y * pitch + 3*x;

            sptr[FI_RGBA_RED]   = r_word;
            sptr[FI_RGBA_GREEN] = g_word;
            sptr[FI_RGBA_BLUE]  = b_word;
        }
        else
        {
            for (y2  = y - pixel_size_half;
                 y2 <= y - pixel_size_half + pixel_size; y2++)
            {
                if (y2 >= 0 && y2 < h)
                {
                    for (x2  = x - pixel_size_half;
                         x2 <= x - pixel_size_half + pixel_size; x2++)
                    {
                        if (x2 >= 0 && x2 < w)
                        {
                            /* calculate the address of the pixel to write */
                            sptr = src + y2 * pitch + 3*x2;

                            sptr[FI_RGBA_RED]   = r_word;
                            sptr[FI_RGBA_GREEN] = g_word;
                            sptr[FI_RGBA_BLUE]  = b_word;
                        }
                    }
                }
            }
        }

#if 0
        fprintf(stderr, "%s\t%ld\t%ld\n",
            node_ptr->name,
            (long) node_map_ptr_col->x, (long) node_map_ptr_col->y);
#endif
    }
    
    
    /* save src pointer for later */
    colored_src = src;


    print_number_width = (int) log10(n_mz_clusters) + 1;
    if (print_number_width < 5)
        print_number_width = 5;


    /* write pixel cluster image */
    if (cluster_image)
    {
        /* allocate enough room for the name string */
        str_len                   = snprintf(NULL, 0, "%s_%0*d.png", base_str,
                                      print_number_width, 0);
        output_file_name          = realloc(output_file_name,
                                            (str_len+1) * sizeof(char));
        output_file_name[str_len] = '\0';

        sprintf(output_file_name, "%s_%0*d.png",
            base_str, print_number_width, 0);

        /* palletize the output image
         *
         * FIQ_WUQUANT throws out too many colors
         * FIQ_NNQUANT is much better, but still loses some, and
         *  produces a medium grey that isn't quite grey...
         * FIQ_LFPQUANT is perfect, but we must do a sanity check
         *  on the number of pixel clusters earlier to limit it to 255
         *
         * The most robust solution would be to specify our exact palette to
         * use, so that FIQ_WUQUANT doesn't auto-generate poor colors,
         * but FIQ_LFPQUANT should work fine in the mean time.
         */
        quantized = FreeImage_ColorQuantizeEx(cluster_image, FIQ_LFPQUANT,
                        n_pixel_clusters + 1, 0, NULL);

        sort_image_palette_by_hue(quantized, n_pixel_clusters + 1);

        FreeImage_FlipVertical(quantized);
        FreeImage_SetDotsPerMeterX(quantized, 0);
        FreeImage_SetDotsPerMeterY(quantized, 0);
        FreeImage_Save(FIF_PNG, quantized, output_file_name,
                       PNG_Z_BEST_COMPRESSION);
        FreeImage_Unload(quantized);
    }


    /* allocate start/end arrays for m/z clusters */
    mz_cluster_start_array = (uint64_t *) malloc(num_mz * sizeof(uint64_t));
    mz_cluster_end_array   = (uint64_t *) calloc(num_mz, sizeof(uint64_t));
    memset(mz_cluster_start_array, 0xFF, num_mz * sizeof(uint64_t));

    /* !!! assume m/z mappings were sorted in tree order previously !!! */
    for (row = 0; row < num_mz; row++)
    {
        node_map_ptr_row = &node_map_array_mz[row];
        node_ptr         = node_map_ptr_row->node_ptr;

        if (node_ptr == NULL)
        {
            fprintf(stderr, "WARNING -- unmapped row %d\n", row);
            continue;
        }

        cluster_mz = node_ptr->cluster_num;
        
        /* we could be clever and set start/end as we change clusters,
         * but I'm opting for the slightly slower, but more robust to
         * coding errors method of min/max checking for now
         */
        if (row < mz_cluster_start_array[cluster_mz])
            mz_cluster_start_array[cluster_mz] = row;
        if (row > mz_cluster_end_array[cluster_mz])
            mz_cluster_end_array[cluster_mz] = row;
    }

    /* find global minimum */
    min_global = DBL_MAX;
    for (row = 0; row < num_mz; row++)
    {
        row_ptr = data_matrix[row];
    
        for (col = 0; col < num_pixels; col++)
        {
            value = data_matrix[row][col];
            
            if (value < min_global)
                min_global = value;
        }
    }


    /* allocate maximum size of temp array */
    size = (uint64_t) w * h * sizeof(double);
    if ((uint64_t) num_mz * sizeof(double) > size)
    {
        size = (uint64_t) num_mz * sizeof(double);
    }
    value_array = (double *) malloc(size);

    /* allocate array for storing values to go into individual images */
    double_2d_array = (double *) malloc((uint64_t) w * h * sizeof(double));


    /* regular cluster contact sheets */
    if (n_selected_clusters == 0)
    {
        /* set up the contact sheet */
        contact_n   = n_mz_clusters + 1;
        best_grid_w = contact_n;
        best_grid_h = 1;
        best_score  = DBL_MAX;
        for (grid_w = 1; grid_w <= contact_n; grid_w++)
        {
            for (grid_h = 1; grid_h <= contact_n; grid_h++)
            {
                grid_n = grid_w * grid_h;
                
                /* grid too small, skip it */
                if (grid_n < contact_n)
                    continue;
                
                /* grid too big, skip it */
                if (grid_n - contact_n > ceil(0.5 * grid_w) &&
                    grid_n - contact_n > ceil(0.5 * grid_h))
                {
                    continue;
                }
                
                grid_max = grid_w;
                if (grid_h > grid_max)
                    grid_max = grid_h;
                    
                aspect = (double) (grid_w * w) / (double) (grid_h * h);

#if 0
                score  = (1.0 + (grid_n - contact_n)) *
                         exp(fabs(log(aspect / target_aspect)));
#else
                score  = (double)(grid_max+(grid_n-contact_n))/grid_max *
                         exp(fabs(log(aspect / target_aspect)));
#endif
                
                /* smaller score is better score */
                if (score < best_score)
                {
                    best_grid_w = grid_w;
                    best_grid_h = grid_h;
                    best_score  = score;
                }
            }
        }
        

        grid_w = best_grid_w;
        grid_h = best_grid_h;
        
        fprintf(stderr, "Contact grid:\t%d\t%d\n",
                (int) grid_w, (int) grid_h);
        
        contact_w     = grid_w * w;
        contact_h     = grid_h * h;
        contact_image = FreeImage_Allocate(contact_w, contact_h, 24,
                        FI_RGBA_RED_MASK, FI_RGBA_GREEN_MASK, FI_RGBA_BLUE_MASK);
        contact_src   = FreeImage_GetBits(contact_image);
        contact_pitch = FreeImage_GetPitch(contact_image);


        /* copy cluster image into contact sheet */
        grid_x = grid_y = 0;
        for (y = 0; y < h; y++)
        {
            contact_y = y + grid_y * h;
        
            for (x = 0; x < w; x++)
            {
                contact_x = x + grid_x * w;
            
                /* calculate the address of the pixel to write */
                sptr         = src + y * pitch + 3*x;
                contact_sptr = contact_src + contact_y * contact_pitch +
                               3*contact_x;

                contact_sptr[FI_RGBA_RED]   = sptr[FI_RGBA_RED];
                contact_sptr[FI_RGBA_GREEN] = sptr[FI_RGBA_GREEN];
                contact_sptr[FI_RGBA_BLUE]  = sptr[FI_RGBA_BLUE];
            }
        }


        /* generate row cluster images */
        for (cluster_mz = 0; cluster_mz < n_mz_clusters; cluster_mz++)
        {
            row_image = FreeImage_Clone(cluster_image);
            src       = FreeImage_GetBits(row_image);
            pitch     = FreeImage_GetPitch(row_image);

            start     = mz_cluster_start_array[cluster_mz];
            end       = mz_cluster_end_array[cluster_mz];


            /* initialize 2D array to missing data */
            for (y = 0; y < h; y++)
            {
                for (x = 0; x < w; x++)
                {
                    double_2d_array[y * w + x] = MALDI_MISSING;
                }
            }
        
            /* each pixel has its own spectra */
            for (col = 0; col < num_pixels; col++)
            {
                node_map_ptr_col = &node_map_array_pixels[col];
                col_data         = node_map_ptr_col->data_index;

                x = node_map_ptr_col->x - x_min;
                y = node_map_ptr_col->y - y_min;


                value_geomean = 0.0;
                n_values      = 0;
                for (row = start; row <= end; row++)
                {
                    node_map_ptr_row = &node_map_array_mz[row];
                    node_ptr         = node_map_ptr_row->node_ptr;
                    row_data         = node_map_ptr_row->data_index;

                    if (node_ptr == NULL)
                    {
                        fprintf(stderr, "WARNING -- unmapped row %d\n", row);
                        continue;
                    }
                    
                    value                    = data_matrix[row_data][col_data];
                    value_array[n_values++]  = value;

                    if (value > 0)
                        value_geomean       += log(value);
                }
                
                /* sort the values from low to high */
                qsort(value_array, n_values, sizeof(double), cmp_double);
                
                if (n_values)
                {
                    index = n_values >> 1;
                
                    /* odd */
                    if (n_values % 2)
                    {
                        value_median = value_array[index];
                    }
                    /* even */
                    else
                    {
                        value_median = 0.5 * (value_array[index-1] +
                                              value_array[index]);
                    }
                    
                    value_geomean = exp(value_geomean / n_values);
                    
                    /* store the value to render for this pixel */
                    double_2d_array[y * w + x] = value_geomean;
                }
            }
            
            /* set high/low thresholds for entire image */
            n_values = 0;
            for (y = 0; y < h; y++)
            {
                for (x = 0; x < w; x++)
                {
                    value = double_2d_array[y * w + x];
                    
                    if (value != MALDI_MISSING)
                        value_array[n_values++] = value;
                }
            }

            /* sort the values from low to high */
            qsort(value_array, n_values, sizeof(double), cmp_double);
            

            value_low   = value_high = MALDI_MISSING;
            value_range = 1.0;

            if (n_values)
            {
                /*
                index       = 0.00 * (n_values - 1);
                value_low   = value_array[index];
                */
                value_low   = min_global;
                
                /* there can be a few outlier pixels, even after normalization */
                /* 0.999 is too high, can swamp out too much */
                index       = 0.99 * (n_values - 1);
                value_high  = value_array[index];

                value_range = value_high - value_low;
            }

            /* sanity check checks */
            if (value_range == 0.0)
                value_range = value_high;
            if (value_range == 0.0)
                value_range = 1.0;
            
            /* generate image */
            for (col = 0; col < num_pixels; col++)
            {
                node_map_ptr_col = &node_map_array_pixels[col];
                col_data         = node_map_ptr_col->data_index;

                x = node_map_ptr_col->x - x_min;
                y = node_map_ptr_col->y - y_min;

                /* scale values to display range */
                value = double_2d_array[y * w + x];
                if (value != MALDI_MISSING)
                {
                    if (value < value_low)  value = value_low;
                    if (value > value_high) value = value_high;

                    value -= value_low;

                    /* calculate the address of the pixel to write */
                    sptr = src + y * pitch + 3*x;

                    r_word = sptr[FI_RGBA_RED];
                    g_word = sptr[FI_RGBA_GREEN];
                    b_word = sptr[FI_RGBA_BLUE];

                    value_scale = 1.0;
                    if (value_range)
                        value_scale = value / value_range;
                    
                    /* sanity check */
                    if (value_scale > 1.0) value_scale = 1.0;

                    r_word = (WORD) (value_scale * r_word + 0.5);
                    g_word = (WORD) (value_scale * g_word + 0.5);
                    b_word = (WORD) (value_scale * b_word + 0.5);

                    if (pixel_size == 1)
                    {
                        sptr[FI_RGBA_RED]   = r_word;
                        sptr[FI_RGBA_GREEN] = g_word;
                        sptr[FI_RGBA_BLUE]  = b_word;
                    }
                    else
                    {
                        for (y2  = y - pixel_size_half;
                             y2 <= y - pixel_size_half + pixel_size; y2++)
                        {
                            if (y2 >= 0 && y2 < h)
                            {
                                for (x2  = x - pixel_size_half;
                                     x2 <= x - pixel_size_half + pixel_size; x2++)
                                {
                                    if (x2 >= 0 && x2 < w)
                                    {
                                        /* calculate the address of the pixel to write */
                                        sptr = src + y2 * pitch + 3*x2;

                                        sptr[FI_RGBA_RED]   = r_word;
                                        sptr[FI_RGBA_GREEN] = g_word;
                                        sptr[FI_RGBA_BLUE]  = b_word;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            /* copy row image into contact sheet */
            if (row_image)
            {
                grid_x = (cluster_mz + 1) % grid_w;
                grid_y = (cluster_mz + 1) / grid_w;

                for (y = 0; y < h; y++)
                {
                    contact_y = y + grid_y * h;

                    for (x = 0; x < w; x++)
                    {
                        contact_x = x + grid_x * w;

                        /* calculate the address of the pixel to write */
                        sptr         = src + y * pitch + 3*x;
                        contact_sptr = contact_src + contact_y * contact_pitch +
                                       3*contact_x;

                        contact_sptr[FI_RGBA_RED]   = sptr[FI_RGBA_RED];
                        contact_sptr[FI_RGBA_GREEN] = sptr[FI_RGBA_GREEN];
                        contact_sptr[FI_RGBA_BLUE]  = sptr[FI_RGBA_BLUE];
                    }
                }
            }

            /* write row image */
            if (opt_mz_images && row_image)
            {
                /* allocate enough room for the name string */
                str_len                   = snprintf(NULL, 0, "%s_%0*ld.png",
                                              base_str, print_number_width,
                                              (long) (cluster_mz + 1));
                output_file_name          = realloc(output_file_name,
                                                    (str_len+1) * sizeof(char));
                output_file_name[str_len] = '\0';

                /* use snprintf() to shut up the bogus
                 * GCC v7 -Wformat-overflow warnings
                 *
                 * but then it still issues bogus -Wformat-truncation warnings...
                 * turns out that %* needs to be passed a short instead of int
                 */
                snprintf(output_file_name, str_len+1, "%s_%0*ld.png",
                         base_str, print_number_width,
                         (long) (cluster_mz + 1));

                FreeImage_FlipVertical(row_image);
                FreeImage_SetDotsPerMeterX(row_image, 0);
                FreeImage_SetDotsPerMeterY(row_image, 0);
                FreeImage_Save(FIF_PNG, row_image, output_file_name,
                               PNG_Z_BEST_COMPRESSION);

            }

            if (row_image)
            {
                FreeImage_Unload(row_image);
                row_image = NULL;
            }
        }

        /* write contact image */
        if (contact_image)
        {
            /* allocate enough room for the name string */
            str_len                   = snprintf(NULL, 0, "%s_contact.png",
                                                 base_str);
            output_file_name          = realloc(output_file_name,
                                                (str_len+1) * sizeof(char));
            output_file_name[str_len] = '\0';

            /* use snprintf() to shut up the bogus
             * GCC v7 -Wformat-overflow warnings
             *
             * but then it still issues bogus -Wformat-truncation warnings...
             * turns out that %* needs to be passed a short instead of int
             */
            snprintf(output_file_name, str_len+1, "%s_contact.png", base_str);

            FreeImage_FlipVertical(contact_image);
            FreeImage_SetDotsPerMeterX(contact_image, 0);
            FreeImage_SetDotsPerMeterY(contact_image, 0);
            FreeImage_Save(FIF_PNG, contact_image, output_file_name,
                           PNG_Z_BEST_COMPRESSION);

            FreeImage_Unload(contact_image);
            contact_image = NULL;
        }
    }
    /* output contact sheets of all m/z within selected clusters */
    else
    {
        for (i = 0; i < n_selected_clusters; i++)
        {
            /* command line is base-1, arrays are base-0 */
            cluster_mz = selected_clusters_array[i] - 1;
            
            /* we don't have that many clusters, skip it */
            if (cluster_mz >= n_pixel_clusters)
                continue;
            
            start = mz_cluster_start_array[cluster_mz];
            end   = mz_cluster_end_array[cluster_mz];
            
            n_mz_selected = end - start + 1;

            /* set up the contact sheet */
            contact_n   = n_mz_selected + 1;
            best_grid_w = contact_n;
            best_grid_h = 1;
            best_score  = DBL_MAX;
            for (grid_w = 1; grid_w <= contact_n; grid_w++)
            {
                for (grid_h = 1; grid_h <= contact_n; grid_h++)
                {
                    grid_n = grid_w * grid_h;
                    
                    /* grid too small, skip it */
                    if (grid_n < contact_n)
                        continue;
                    
                    /* grid too big, skip it */
                    if (grid_n - contact_n > ceil(0.5 * grid_w) &&
                        grid_n - contact_n > ceil(0.5 * grid_h))
                    {
                        continue;
                    }
                    
                    grid_max = grid_w;
                    if (grid_h > grid_max)
                        grid_max = grid_h;
                    
                    aspect = (double) (grid_w * w) / (double) (grid_h * h);

#if 0
                    score  = (1.0 + (grid_n - contact_n)) *
                             exp(fabs(log(aspect / target_aspect)));
#else
                    score  = (double)(grid_max+(grid_n-contact_n))/grid_max *
                             exp(fabs(log(aspect / target_aspect)));
#endif
                    
                    /* smaller score is better score */
                    if (score < best_score)
                    {
                        best_grid_w = grid_w;
                        best_grid_h = grid_h;
                        best_score  = score;
                    }
                }
            }
            

            grid_w = best_grid_w;
            grid_h = best_grid_h;
            
            fprintf(stderr, "Contact grid %d:\t%d\t%d\n",
                    (int) cluster_mz + 1,
                    (int) grid_w, (int) grid_h);
            
            contact_w     = grid_w * w;
            contact_h     = grid_h * h;
            contact_image = FreeImage_Allocate(contact_w, contact_h, 24,
                            FI_RGBA_RED_MASK, FI_RGBA_GREEN_MASK, FI_RGBA_BLUE_MASK);
            contact_src   = FreeImage_GetBits(contact_image);
            contact_pitch = FreeImage_GetPitch(contact_image);


            /* copy cluster image into contact sheet */
            grid_x = grid_y = 0;
            for (y = 0; y < h; y++)
            {
                contact_y = y + grid_y * h;
            
                for (x = 0; x < w; x++)
                {
                    contact_x = x + grid_x * w;
                
                    /* calculate the address of the pixel to write */
                    sptr         = colored_src + y * pitch + 3*x;
                    contact_sptr = contact_src + contact_y * contact_pitch +
                                   3*contact_x;

                    contact_sptr[FI_RGBA_RED]   = sptr[FI_RGBA_RED];
                    contact_sptr[FI_RGBA_GREEN] = sptr[FI_RGBA_GREEN];
                    contact_sptr[FI_RGBA_BLUE]  = sptr[FI_RGBA_BLUE];
                }
            }


            /* generate m/z images */
            for (row = start; row <= end; row++)
            {
                row_image = FreeImage_Clone(cluster_image);
                src       = FreeImage_GetBits(row_image);
                pitch     = FreeImage_GetPitch(row_image);

                /* initialize 2D array to missing data */
                for (y = 0; y < h; y++)
                {
                    for (x = 0; x < w; x++)
                    {
                        double_2d_array[y * w + x] = MALDI_MISSING;
                    }
                }


                node_map_ptr_row = &node_map_array_mz[row];
                node_ptr         = node_map_ptr_row->node_ptr;
                row_data         = node_map_ptr_row->data_index;


                /* store intensity data for each pixel */
                for (col = 0; col < num_pixels; col++)
                {
                    node_map_ptr_col = &node_map_array_pixels[col];
                    col_data         = node_map_ptr_col->data_index;

                    x = node_map_ptr_col->x - x_min;
                    y = node_map_ptr_col->y - y_min;

                    if (node_ptr == NULL)
                    {
                        fprintf(stderr, "WARNING -- unmapped row %d\n", row);
                        continue;
                    }
                        
                    value = data_matrix[row_data][col_data];

                    /* store the value to render for this pixel */
                    double_2d_array[y * w + x] = value;
                }
                
                /* set high/low thresholds for entire image */
                n_values = 0;
                for (y = 0; y < h; y++)
                {
                    for (x = 0; x < w; x++)
                    {
                        value = double_2d_array[y * w + x];
                        
                        if (value != MALDI_MISSING)
                            value_array[n_values++] = value;
                    }
                }

                /* sort the values from low to high */
                qsort(value_array, n_values, sizeof(double), cmp_double);
                

                value_low   = value_high = MALDI_MISSING;
                value_range = 1.0;

                if (n_values)
                {
                    /*
                    index       = 0.00 * (n_values - 1);
                    value_low   = value_array[index];
                    */
                    value_low   = min_global;
                    
                    /* there can be a few outlier pixels, even after normalization */
                    /* 0.999 is too high, can swamp out too much */
                    index       = 0.99 * (n_values - 1);
                    value_high  = value_array[index];

                    value_range = value_high - value_low;
                }

                /* sanity check checks */
                if (value_range == 0.0)
                    value_range = value_high;
                if (value_range == 0.0)
                    value_range = 1.0;
                
                /* generate image */
                for (col = 0; col < num_pixels; col++)
                {
                    node_map_ptr_col = &node_map_array_pixels[col];
                    col_data         = node_map_ptr_col->data_index;

                    x = node_map_ptr_col->x - x_min;
                    y = node_map_ptr_col->y - y_min;

                    /* calculate the address of the pixel to write */
                    sptr = src + y * pitch + 3*x;
                    
                    /* scale values to display range */
                    value = double_2d_array[y * w + x];
                    if (value != MALDI_MISSING)
                    {
                        if (value < value_low)  value = value_low;
                        if (value > value_high) value = value_high;

                        value -= value_low;

                        r_word = sptr[FI_RGBA_RED];
                        g_word = sptr[FI_RGBA_GREEN];
                        b_word = sptr[FI_RGBA_BLUE];

                        value_scale = 1.0;
                        if (value_range)
                            value_scale = value / value_range;
                        
                        /* sanity check */
                        if (value_scale > 1.0) value_scale = 1.0;

                        r_word = (WORD) (value_scale * r_word + 0.5);
                        g_word = (WORD) (value_scale * g_word + 0.5);
                        b_word = (WORD) (value_scale * b_word + 0.5);

                        sptr[FI_RGBA_RED]   = r_word;
                        sptr[FI_RGBA_GREEN] = g_word;
                        sptr[FI_RGBA_BLUE]  = b_word;
                    }
                }

                /* copy row image into contact sheet */
                if (row_image)
                {
                    grid_x = (row - start + 1) % grid_w;
                    grid_y = (row - start + 1) / grid_w;

                    for (y = 0; y < h; y++)
                    {
                        contact_y = y + grid_y * h;

                        for (x = 0; x < w; x++)
                        {
                            contact_x = x + grid_x * w;

                            /* calculate the address of the pixel to write */
                            sptr         = src + y * pitch + 3*x;
                            contact_sptr = contact_src + contact_y * contact_pitch +
                                           3*contact_x;

                            contact_sptr[FI_RGBA_RED]   = sptr[FI_RGBA_RED];
                            contact_sptr[FI_RGBA_GREEN] = sptr[FI_RGBA_GREEN];
                            contact_sptr[FI_RGBA_BLUE]  = sptr[FI_RGBA_BLUE];
                        }
                    }
                }

                if (row_image)
                {
                    FreeImage_Unload(row_image);
                    row_image = NULL;
                }
            }


            print_number_width = (int) log10(n_mz_selected) + 1;
            if (print_number_width < 5)
                print_number_width = 5;


            /* write contact image */
            if (contact_image)
            {
                /* allocate enough room for the name string */
                str_len = snprintf(NULL, 0, "%s_%0*ld_contact.png",
                         base_str, print_number_width,
                         (long) (cluster_mz + 1));

                output_file_name          = realloc(output_file_name,
                                                    (str_len+1) * sizeof(char));
                output_file_name[str_len] = '\0';

                /* use snprintf() to shut up the bogus
                 * GCC v7 -Wformat-overflow warnings
                 *
                 * but then it still issues bogus -Wformat-truncation warnings,
                 * turns out that %* needs to be passed a short instead of int
                 */
                snprintf(output_file_name, str_len+1, "%s_%0*ld_contact.png",
                         base_str, print_number_width,
                         (long) (cluster_mz + 1));

                FreeImage_FlipVertical(contact_image);
                FreeImage_SetDotsPerMeterX(contact_image, 0);
                FreeImage_SetDotsPerMeterY(contact_image, 0);
                FreeImage_Save(FIF_PNG, contact_image, output_file_name,
                               PNG_Z_BEST_COMPRESSION);

                FreeImage_Unload(contact_image);
                contact_image = NULL;
            }
            
            
            /* write contact sheet grid legend */
            
            /* allocate enough room for the name string */
            str_len = snprintf(NULL, 0, "%s_%0*ld_contact_legend.txt",
                     base_str, print_number_width,
                     (long) (cluster_mz + 1));

            output_file_name          = realloc(output_file_name,
                                                (str_len+1) * sizeof(char));
            output_file_name[str_len] = '\0';

            /* use snprintf() to shut up the bogus
             * GCC v7 -Wformat-overflow warnings
             *
             * but then it still issues bogus -Wformat-truncation warnings...
             * turns out that %* needs to be passed a short instead of int
             */
            snprintf(output_file_name, str_len+1, "%s_%0*ld_contact_legend.txt",
                     base_str, print_number_width,
                     (long) (cluster_mz + 1));

            /* open as text, so it will translate EOL automatically */
            /* this should generally be a small file, so we shouldn't need
             * to muck around with a custom buffer size
             */
            outfile_mz = fopen(output_file_name, "wt");
            if (!outfile_mz)
            {
                fprintf(stderr, "ERROR -- can't open output file %s\n",
                        output_file_name);
            }
            
            /* output header row */
            fprintf(outfile_mz, "Grid");
            for (grid_x = 0; grid_x < grid_w; grid_x++)
                fprintf(outfile_mz, "\t%d", (int) grid_x + 1);
            fprintf(outfile_mz, "\n");
            
            fprintf(outfile_mz, "1");
            fprintf(outfile_mz, "\tCluster %d", (int) cluster_mz + 1);

            grid_y_old = 0;
            for (row = start; row <= end; row++)
            {
                node_map_ptr_row = &node_map_array_mz[row];
                node_ptr         = node_map_ptr_row->node_ptr;

                grid_x = (row - start + 1) % grid_w;
                grid_y = (row - start + 1) / grid_w;

                if (grid_y != grid_y_old)
                {
                    fprintf(outfile_mz, "\n");
                    fprintf(outfile_mz, "%d", (int) grid_y + 1);
                }
                
                fprintf(outfile_mz, "\t%s", node_ptr->name);

                grid_y_old = grid_y;
            }
            fprintf(outfile_mz, "\n");
            
            fclose(outfile_mz);
        }
    }


    if (value_array)
        free(value_array);

    if (double_2d_array)
        free(double_2d_array);

    if (output_file_name)
        free(output_file_name);

    if (mz_cluster_start_array)
        free(mz_cluster_start_array);
    if (mz_cluster_end_array)
        free(mz_cluster_end_array);

    /* deallocate images */
    if (cluster_image)
        FreeImage_Unload(cluster_image);
    if (row_image)
        FreeImage_Unload(row_image);
}


int guess_reverse_palette(double **data_matrix,
                          struct  node_data_map *node_map_array_pixels,
                          int32_t num_pixels,
                          int32_t n_pixel_clusters)
{
    struct     tree_node     *node_ptr;
    struct     node_data_map *node_map_ptr_col;
    int32_t    col;
    int32_t    cluster_pixel;
    
    double    *cluster_mean_dists = NULL;    /* mean distance from center */
    int32_t   *cluster_counts     = NULL;    /* number of pixels per cluster */
    double     x_mean, y_mean, dist;
    double     x_min, x_max, y_min, y_max;
    double     best_dist;
    double     sum_half1, sum_half2;
    int32_t    x, y;
    int32_t    n_mapped;
    int32_t    i;
    int32_t    n_clusters_half;
    int32_t    best_i;

    x_min    = y_min  =  DBL_MAX;
    x_max    = y_max  = -DBL_MAX;
    x_mean   = y_mean =  0.0;
    n_mapped = 0;
    
    for (col = 0; col < num_pixels; col++)
    {
        node_map_ptr_col = &node_map_array_pixels[col];
        node_ptr         = node_map_ptr_col->node_ptr;
        
        /* uh oh, this pixel wasn't mapped */
        if (node_ptr == NULL)
        {
            fprintf(stderr, "WARNING -- unmapped col %d\n", col);
            continue;
        }
        
        x = node_map_ptr_col->x;
        y = node_map_ptr_col->y;
        
        if (x < x_min) x_min = x;
        if (x > x_max) x_max = x;
        if (y < y_min) y_min = y;
        if (y > y_max) y_max = y;
        
        x_mean += x;
        y_mean += y;

        n_mapped++;
    }
    
    /* uh oh, we're probably going to crash anyways... */
    if (n_mapped == 0)
        return -1;


    cluster_mean_dists = (double *)  calloc(n_mapped, sizeof(double));
    cluster_counts     = (int32_t *) calloc(n_mapped, sizeof(int32_t));


    /* x/y means are the coordinates of the center of the shape */
    x_mean /= n_mapped;
    y_mean /= n_mapped;


    for (col = 0; col < num_pixels; col++)
    {
        node_map_ptr_col = &node_map_array_pixels[col];
        node_ptr         = node_map_ptr_col->node_ptr;
        
        /* uh oh, this pixel wasn't mapped */
        if (node_ptr == NULL)
        {
            fprintf(stderr, "WARNING -- unmapped col %d\n", col);
            continue;
        }
        
        x = node_map_ptr_col->x;
        y = node_map_ptr_col->y;
        
        dist = sqrt((x - x_mean) * (x - x_mean) +
                    (y - y_mean) * (y - y_mean));

        cluster_pixel = node_ptr->cluster_num;
        
        cluster_mean_dists[cluster_pixel] += dist;
        cluster_counts[cluster_pixel]++;
    }
    

    best_dist = -DBL_MAX;
    best_i = 0;
    n_clusters_half = n_pixel_clusters >> 1;

    /* HACK - we may want to weight these by num pixels? */
    sum_half1 = sum_half2 = 0.0;

    for (i = 0; i < n_pixel_clusters; i++)
    {
        if (cluster_counts[i])
            cluster_mean_dists[i] /= cluster_counts[i];

        dist = cluster_mean_dists[i];
        
        if (dist > best_dist)
        {
            best_dist = dist;
            best_i    = i;
        }
        
        if (i < n_clusters_half)
            sum_half1 += dist;
        
        if (i >= n_pixel_clusters - n_clusters_half)
            sum_half2 += dist;
    }
    
    if (cluster_mean_dists)
        free(cluster_mean_dists);
    
    if (cluster_counts)
        free(cluster_counts);


    if (sum_half2 > sum_half1)
        return 1;

    return 0;
}


/* red = 0; green = 120; blue = 240 */
void generate_heatmap_blue_white_red(FIBITMAP *image,
                         double **data_matrix,
                         struct  node_data_map *node_map_array_pixels,
                         int32_t num_pixels,
                         struct  node_data_map *node_map_array_mz,
                         int32_t num_mz,
                         double forced_upper_bound)
{
    BYTE    *src = FreeImage_GetBits(image);
    BYTE    *sptr;
    struct   node_data_map *node_map_ptr_col, *node_map_ptr_row;
    uint64_t w, h, pitch;
    uint64_t n;

    uint64_t x, y, row, col, row_data, col_data;    /* x = col, y = row */
    double data_value, data_value_orig;
    double *dptr;
    
    double sd;
    double dist;
    double max_value = -9E99;
    double min_orig_value = 9E99;
    double upper_bound;

    double hue, l, s;
    double r, g, b;
    WORD r_word, g_word, b_word;
    
    w = FreeImage_GetWidth(image);
    h = FreeImage_GetHeight(image);
    pitch = FreeImage_GetPitch(image);


    fprintf(stderr, "w:%ld h:%ld p:%ld\n", w, h, pitch);
    
    sd = 0.0;
    n  = 0;
    for (row = 0; row < h; row++)
    {
        node_map_ptr_row = &node_map_array_mz[row];
        row_data = node_map_ptr_row->data_index;

        dptr = data_matrix[row_data];
    
        for (col = 0; col < w; col++)
        {
            node_map_ptr_col = &node_map_array_pixels[col];
            col_data = node_map_ptr_col->data_index;
        
            if (dptr[col_data] != MALDI_MISSING)
            {
                dist = fabs(dptr[col_data]);
                sd += dist * dist;

                if (dist > max_value)
                    max_value = dist;
                if (dptr[col_data] < min_orig_value)
                    min_orig_value = dptr[col_data];

                n++;
            }
        }
    }
    if (n)
        sd = sqrt(sd / n);


    if (forced_upper_bound)
    {
        upper_bound = forced_upper_bound;
    }
    else
    {
        upper_bound = sqrt(2.0) * sd;
    
        if (upper_bound > max_value)
            upper_bound = max_value;
    }

    fprintf(stderr, "|CappedValue|: %lf\n", upper_bound);
    fprintf(stderr, "MinValue: %lf\n", min_orig_value);
    


    for (row = 0; row < h; row++)
    {
        node_map_ptr_row = &node_map_array_mz[row];
        row_data         = node_map_ptr_row->data_index;
        y                = node_map_ptr_row->leaf_index;

        dptr             = data_matrix[row_data];

        for (col = 0; col < w; col++)
        {
            node_map_ptr_col = &node_map_array_pixels[col];
            col_data         = node_map_ptr_col->data_index;
            x                = node_map_ptr_col->leaf_index;

            /* calculate the address of the pixel to write */
            sptr = src + y * pitch + 3*x;

            data_value = dptr[col_data];
            data_value_orig = data_value;
            
            s = 1.0;
            l = 0.5;
            
            /* scale value to [0, 1] */
            if (data_value >  upper_bound) data_value =  upper_bound;
            if (data_value < -upper_bound) data_value = -upper_bound;
            
            /* special case for all positive */
            if (min_orig_value >= -1.0E-5)
            {
		/* red / white / blue */
                data_value = 2.0 * (data_value - 0.5 * upper_bound);
                l = 1.0 - 0.5 * fabs(data_value) / upper_bound;

                hue = 120;
                if (data_value < 0)
                    hue = 240;
                else if (data_value > 0)
                    hue = 0;
            }
            else
            {
                l = 1.0 - 0.5 * fabs(data_value) / upper_bound;

                hue = 120;
                if (data_value < 0)
                    hue = 240;
                else if (data_value > 0)
                    hue = 0;
            }


            /* regular non-missing data */
            if (data_value_orig != MALDI_MISSING)
            {
                hls_to_rgb(hue, l, s, &r, &g, &b);
            
                r_word = (WORD) (256.0 * r);
                g_word = (WORD) (256.0 * g);
                b_word = (WORD) (256.0 * b);
            
                if (r_word > 255) r_word = 255;
                if (g_word > 255) g_word = 255;
                if (b_word > 255) b_word = 255;
            }
            /* set missing/bad data to grey */
            else
            {
               r_word = 0x7F;
               g_word = 0x7F;
               b_word = 0x7F;
            }

            
            sptr[FI_RGBA_RED]   = r_word;
            sptr[FI_RGBA_GREEN] = g_word;
            sptr[FI_RGBA_BLUE]  = b_word;
        }
    }
}


void write_heatmap(char *base_str, double **data_matrix,
                   struct  node_data_map *node_map_array_pixels,
                   int32_t num_pixels,
                   struct  node_data_map *node_map_array_mz,
                   int32_t num_mz,
                   double forced_upper_bound)
{
    FIBITMAP *heatmap_image = NULL;
    FIBITMAP *resized_image;
    char *new_filename = NULL;
    char *sptr;

    /* 11" x 8.5", 600 DPI */
    WORD w_target = 6600;
    WORD h_target = 5100;
    WORD w = num_pixels, h = num_mz;
    
    /* shrink the width if there are very few samples */
    if (num_pixels * 300 < w_target)
        w_target = 300 * num_pixels;
    
    heatmap_image = FreeImage_Allocate(num_pixels, num_mz, 24,
                    FI_RGBA_RED_MASK, FI_RGBA_GREEN_MASK, FI_RGBA_BLUE_MASK);
    
    generate_heatmap_blue_white_red(heatmap_image, data_matrix,
                                    node_map_array_pixels, num_pixels,
                                    node_map_array_mz, num_mz,
                                    forced_upper_bound);

    /* output heatmap image */
    new_filename = (char *) realloc(new_filename,
            (strlen(base_str) + 1 + 12) * sizeof(char));
    strcpy(new_filename, base_str);
    
    /* null out the last . */
    for (sptr = new_filename + strlen(new_filename) - 1;
         sptr >= new_filename; sptr--)
    {
        if (*sptr == '.')
        {
            *sptr = '\0';
            break;
        }
    }
    
    strcat(new_filename, "_heatmap.jpg");
    
    fprintf(stderr, "Saving file: %s\n", new_filename);

    /* resize dimensions separately, for finer control of resize method
     * Use Box filter if up-sizing, Catmull-Rom if down-sizing)
     */
    if (w <= w_target && h <= h_target)
    {
        resized_image = FreeImage_Rescale(heatmap_image, w_target, h_target,
                          FILTER_BOX);
    }
    else if (w <= w_target && h > h_target)
    {
        resized_image = FreeImage_Rescale(heatmap_image, w_target, h,
                          FILTER_BOX);

        FreeImage_Unload(heatmap_image);
        heatmap_image = FreeImage_Clone(resized_image);

        resized_image = FreeImage_Rescale(heatmap_image, w_target, h_target,
                          FILTER_CATMULLROM);
    }
    else if (h <= h_target && w > w_target)
    {
        resized_image = FreeImage_Rescale(heatmap_image, w, h_target,
                          FILTER_BOX);

        FreeImage_Unload(heatmap_image);
        heatmap_image = FreeImage_Clone(resized_image);

        resized_image = FreeImage_Rescale(heatmap_image, w_target, h_target,
                          FILTER_CATMULLROM);
    }
    else
    {
        resized_image = FreeImage_Rescale(heatmap_image, w_target, h_target,
                          FILTER_CATMULLROM);
    }

    FreeImage_FlipVertical(resized_image);
    FreeImage_Save(FIF_JPEG, resized_image, new_filename,
                   (85 | JPEG_OPTIMIZE | JPEG_BASELINE));
    
    FreeImage_Unload(resized_image);
    FreeImage_Unload(heatmap_image);
}


void fill_pixel_clusters_means(double **data_matrix,
                               double **col_cluster_spectra,
                               struct  node_data_map *node_map_array_pixels,
                               int32_t num_pixels,
                               struct  node_data_map *node_map_array_mz,
                               int32_t num_mz,
                               int32_t n_pixel_clusters)
{
    struct    tree_node     *node_ptr;
    struct    node_data_map *node_map_ptr_col;
    uint32_t *cluster_pixel_counts = NULL;
    double    value;
    uint32_t  col, col_data, row;
    int32_t   cluster_pixel;


    cluster_pixel_counts = calloc(n_pixel_clusters, sizeof(uint32_t));


    for (col = 0; col < num_pixels; col++)
    {
        node_map_ptr_col = &node_map_array_pixels[col];
        node_ptr         = node_map_ptr_col->node_ptr;
        col_data         = node_map_ptr_col->data_index;
        
        /* uh oh, this pixel wasn't mapped */
        if (node_ptr == NULL)
        {
            fprintf(stderr, "WARNING -- unmapped col %d\n", col);
            continue;
        }
        
        cluster_pixel = node_ptr->cluster_num;
        cluster_pixel_counts[cluster_pixel]++;
        
        /* sum m/z in original data file order */
        for (row = 0; row < num_mz; row++)
        {
            value = data_matrix[row][col_data];
            
            /* for these purposes, convert any missing data to zero */
            if (value == MALDI_MISSING)
                value = 0;
            
            /* log2(1) and log2(0) will be the same;
             * ideally, we've already imputed all the 0's previously
             */
            if (value > 0)
            {
                col_cluster_spectra[cluster_pixel][row] += log2(value);
            }
        }
    }

    /* divide by counts to get the average */
    for (cluster_pixel = 0; cluster_pixel < n_pixel_clusters; cluster_pixel++)
    {
        if (cluster_pixel_counts[cluster_pixel])
        {
            for (row = 0; row < num_mz; row++)
            {
                col_cluster_spectra[cluster_pixel][row] /=
                    cluster_pixel_counts[cluster_pixel];
            }
        }
    }
    
    /* leave the geometric means in log2 space for now */
    

    if (cluster_pixel_counts)
        free(cluster_pixel_counts);
}


void write_clusters_to_file(double  **col_cluster_spectra,
                            char    **col_name_array,
                            char    **row_name_array,
                            struct node_data_map *node_map_array_pixels,
                            int32_t   num_pixels,
                            struct    node_data_map *node_map_array_mz,
                            int32_t   num_mz,
                            int32_t   n_pixel_clusters, int32_t n_mz_clusters,
                            uint32_t *pixel_cluster_colors,
                            double   *row_min_array,  double *row_max_array,
                            double   *row_mean_array, double *row_sd_array,
                            char     *base_str)
{
    FILE     *outfile_pixels  = NULL;
    FILE     *outfile_mz      = NULL;
    char     *buffer_pixels   = NULL;
    char     *buffer_mz       = NULL;
    char     *filename_pixels = NULL;
    char     *filename_mz     = NULL;

    struct    tree_node     *node_ptr;
    struct    node_data_map *node_map_ptr_col, *node_map_ptr_row;

    uint32_t *cluster_counts_pixels = NULL;
    uint32_t *cluster_counts_mz     = NULL;
    uint32_t  row, col, row_data, col_data;
    uint32_t  i;
    
    double    min, max, mean, sd, value;


    /* allocate my own i/o buffers, since we can't trust the system and/or
     * compiler to allocate a decently large one...
     */
    buffer_pixels = (char *) malloc(1000000 * sizeof(char));
    buffer_mz     = (char *) malloc(1000000 * sizeof(char));

    cluster_counts_pixels = calloc(n_pixel_clusters, sizeof(uint32_t));
    cluster_counts_mz     = calloc(n_mz_clusters,    sizeof(uint32_t));

    filename_pixels = (char *) calloc(strlen(base_str) + 1 +
                                      strlen("_pixel_clusters.txt"),
                                      sizeof(char));
    strcpy(filename_pixels, base_str);
    strcpy(filename_pixels + strlen(base_str), "_pixel_clusters.txt");

    /* open as text, so it will translate EOL automatically */
    outfile_pixels = fopen(filename_pixels, "wt");
    if (!outfile_pixels)
    {
        fprintf(stderr, "ERROR -- can't open output file %s\n",
                filename_pixels);
    }
    if (outfile_pixels)
        setvbuf(outfile_pixels, buffer_pixels, _IOFBF, 1000000);


    filename_mz     = (char *) calloc(strlen(base_str) + 1 +
                                      strlen("_mz_clusters.txt"),
                                      sizeof(char));
    strcpy(filename_mz, base_str);
    strcpy(filename_mz + strlen(base_str), "_mz_clusters.txt");

    /* open as text, so it will translate EOL automatically */
    outfile_mz = fopen(filename_mz, "wt");
    if (!outfile_mz)
    {
        fprintf(stderr, "ERROR -- can't open output file %s\n",
                filename_mz);
    }
    if (outfile_mz)
        setvbuf(outfile_mz, buffer_mz, _IOFBF, 1000000);


    /* output the m/z clusters */
    if (outfile_mz)
    {
        /* count number of members in each cluster */
        for (row = 0; row < num_mz; row++)
        {
            node_map_ptr_row = &node_map_array_mz[row];
            node_ptr         = node_map_ptr_row->node_ptr;
            i                = node_ptr->cluster_num;
            
            cluster_counts_mz[i]++;
        }

        /* output the header line */
        fprintf(outfile_mz, "%s\t%s\t%s\t%s",
                "m/z", "tree order", "cluster", "cluster_size");
        fprintf(outfile_mz, "\t%s\t%s\t%s\t%s",
                "RowMin", "RowMean", "RowMax", "RowStDev");
        fprintf(outfile_mz, "\t%s\t%s\t%s\t%s",
                "CMin", "CMean", "CMax", "CStDev");

        /* append hexadecimal color to the end */
        for (i = 0; i < n_pixel_clusters; i++)
        {
            fprintf(outfile_mz, "\tC%03d_x%06X",
                i+1, pixel_cluster_colors[i]);
        }
        fprintf(outfile_mz, "\n");
        
        
        /* output the m/z rows */
        for (row = 0; row < num_mz; row++)
        {
            node_map_ptr_row = &node_map_array_mz[row];
            node_ptr         = node_map_ptr_row->node_ptr;
            row_data         = node_map_ptr_row->data_index;
            
            
            /* calculate pixel cluster stats */
            mean =  sd = 0.0;
            min  =  DBL_MAX;
            max  = -DBL_MAX;

            for (i = 0; i < n_pixel_clusters; i++)
            {
                value = col_cluster_spectra[i][row_data];
                
                if (value < min) min = value;
                if (value > max) max = value;
                
                mean += value;
            }
            if (n_pixel_clusters)
                mean /= n_pixel_clusters;

            for (i = 0; i < n_pixel_clusters; i++)
            {
                value  = col_cluster_spectra[i][row_data];
                value  = value - mean;
                sd    += value * value;
            }
            if (n_pixel_clusters)
                sd = sqrt(sd / n_pixel_clusters);


            fprintf(outfile_mz, "%s\t%d\t%03d\t%d",
                    row_name_array[row_data],
                    (int) row + 1,
                    (int) node_ptr->cluster_num + 1,
                    cluster_counts_mz[node_ptr->cluster_num]);

            fprintf(outfile_mz, "\t%lf\t%lf\t%lf\t%lf",
                    row_min_array[row], row_mean_array[row],
                    row_max_array[row], row_sd_array[row]);

            fprintf(outfile_mz, "\t%lf\t%lf\t%lf\t%lf",
                    min, mean, max, sd);
            
            for (i = 0; i < n_pixel_clusters; i++)
            {
                fprintf(outfile_mz, "\t%lf",
                        col_cluster_spectra[i][row_data]);
            }
            fprintf(outfile_mz, "\n");
        }
        
        
        fclose(outfile_mz);
    }
    
    
    if (outfile_pixels)
    {
        int32_t x, y;
    
        /* count number of members in each cluster */
        for (col = 0; col < num_pixels; col++)
        {
            node_map_ptr_col = &node_map_array_pixels[col];
            node_ptr         = node_map_ptr_col->node_ptr;
            i                = node_ptr->cluster_num;
            
            cluster_counts_pixels[i]++;
        }


        /* output the header line */
        fprintf(outfile_pixels, "%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                "Pixel", "X", "Y", "tree order", "cluster", "cluster_size",
                "color");


        /* output the pixels lines */
        for (col = 0; col < num_pixels; col++)
        {
            node_map_ptr_col = &node_map_array_pixels[col];
            node_ptr         = node_map_ptr_col->node_ptr;
            col_data         = node_map_ptr_col->data_index;
            x                = node_map_ptr_col->x;
            y                = node_map_ptr_col->y;
            
            fprintf(outfile_pixels, "%s\t%d\t%d\t%d\t%03d\t%d\tx%06X\n",
                    col_name_array[col_data],
                    (int) x, (int) y,
                    (int) col + 1,
                    (int) node_ptr->cluster_num + 1,
                    cluster_counts_pixels[node_ptr->cluster_num],
                    pixel_cluster_colors[node_ptr->cluster_num]);
        }

        fclose(outfile_pixels);
    }


    if (cluster_counts_pixels)
        free(cluster_counts_pixels);
    
    if (cluster_counts_mz)
        free(cluster_counts_mz);

    if (filename_pixels)
        free(filename_pixels);

    if (filename_mz)
        free(filename_mz);
    
    if (buffer_pixels)
        free(buffer_pixels);

    if (buffer_mz)
        free(buffer_mz);
}


/* Very aggressive in parsing of clusters and range...
 * Returns NULL if no clusters were parsed from string
 */
int * parse_selected_clusters_string(char *selected_clusters_str,
                                     int  *return_n_clusters)
{
    char *sptr, *sptr_hyphen, *sptr_delim;
    int  *selected_clusters_array = NULL;
    int   n = 0, n_old;
    int   i, j, start, end;
    int   found_flag;
    
    /* initialize return n clusters, in case we abort early */
    *return_n_clusters = 0;

    sptr = selected_clusters_str;
    
    /* scan forward until we find the first digit */
    while (!isdigit(*sptr))
        sptr++;
     
    while (isdigit(*sptr))
    {
        selected_clusters_array = (int *) realloc(selected_clusters_array,
                                                  (n+1) * sizeof(int));
        start = selected_clusters_array[n++] = atoi(sptr);

        
        /* advance to next non-digit */
        sptr++;

        while (isdigit(*sptr))
            sptr++;

        /* check for range(s) */
        while (*sptr == '-')
        {
            sptr_hyphen = sptr;
        
            /* read through multiple hyhpens in a row */
            sptr = sptr_hyphen + 1;
            while (*sptr++ == '-')
                sptr_hyphen = sptr;

            sptr = sptr_hyphen + 1;

            /* non-digit encountered, exit early */
            if (!isdigit(sptr))
            {
                break;
            }

            end = atoi(sptr);

            /* swap start/end */
            if (end < start)
            {
                start ^= end;
                end   ^= start;
                start ^= end;
            }
            
            n_old = n;
            for (i = start; i <= end; i++)
            {
                /* skip any clusters we've already stored */
                found_flag = 0;
                for (j = 0; j < n_old; j++)
                {
                    if (selected_clusters_array[j] == i)
                    {
                        found_flag = 1;
                        break;
                    }
                }
                if (found_flag)
                    continue;

                selected_clusters_array =
                    (int *) realloc(selected_clusters_array,
                                    (n+1) * sizeof(int));
                selected_clusters_array[n++] = i;
            }
            
            /* set start to original end, in case we have more hyphens */
            start = atoi(sptr);
            
            /* advance to next non-digit */
            sptr++;
            while (isdigit(*sptr))
                sptr++;
        }

        /* advance to next delimiter
         * accept comma, semicolon, vertical bar, space
         */
        sptr_delim = strchr(sptr, ',');
        if (!sptr_delim)
            sptr_delim = strchr(sptr, ';');
        if (!sptr_delim)
            sptr_delim = strchr(sptr, '|');
        if (!sptr_delim)
            sptr_delim = strchr(sptr, ' ');

        if (sptr_delim)
        {
            /* read through multiple delimiters in a row */
            sptr = sptr_delim + 1;
            while (*sptr == ',' || *sptr == ';' ||
                   *sptr == '|' || *sptr == ' ')
            {
                sptr_delim = sptr++;
            }
            
            sptr = sptr_delim + 1;

            /* scan past any garbage that might be there */
            while (*sptr != '\0' && !isdigit(*sptr))
                sptr++;
        }
        /* we're either done, or there is some garbage to scan past */
        else
        {
            while (*sptr != '\0' && !isdigit(*sptr))
                sptr++;
        }
    }

    /* sort the values from low to high */
    if (n)
        qsort(selected_clusters_array, n, sizeof(int32_t), cmp_int32);


    *return_n_clusters = n;

    return selected_clusters_array;
}


int main(int argc, char *argv[])
{
    char *data_file_name      = NULL;
    char *tree_file_name_rows = NULL;
    char *tree_file_name_cols = NULL;
    char *output_file_prefix  = NULL;
    char *tree_string_rows    = NULL;
    char *tree_string_cols    = NULL;
    struct tree_node *tree_root_rows;
    struct tree_node *tree_root_cols;
    
    char **row_name_array = NULL;
    char **col_name_array = NULL;

    double            **data_matrix          = NULL;
    double            **col_cluster_spectra  = NULL;
    struct tree_node  **node_ptr_array_rows  = NULL;
    struct tree_node  **leaf_ptr_array_rows  = NULL;
    struct tree_node  **node_ptr_array_cols  = NULL;
    struct tree_node  **leaf_ptr_array_cols  = NULL;
    uint32_t           *pixel_cluster_colors = NULL;   /* RGB byte order */
    int32_t num_cols_chunk = 0;    /* #cols in file last read in */
    int32_t num_nodes_rows;
    int32_t num_leaves_rows;
    int32_t num_nodes_cols;
    int32_t num_leaves_cols;

    struct node_data_map *node_map_array_rows = NULL;
    struct node_data_map *node_map_array_cols = NULL;
    int32_t num_rows = 0, num_cols = 0;
    int32_t i;
    
    /* struct tree_node *node_ptr; */
    int32_t target_n_pixel_clusters = 0;
    int32_t n_pixel_clusters        = 0;
    int32_t target_n_mz_clusters    = 0;
    int32_t n_mz_clusters           = 0;
    int32_t guess_reverse           = 0;

    /* stats for row output */
    double *row_min_array  = NULL;
    double *row_max_array  = NULL;
    double *row_mean_array = NULL;
    double *row_sd_array   = NULL;

    /* command line parsing */
    int     error_flag              = 0;
    int     num_files               = 0;
    char   *sptr;
    struct  options opt;

    int     opt_reverse_palette     = 0;
    int     opt_mz_images           = 0;  /* output individual m/z clusters */

    /* for selecting clusters to generate m/z contact sheets */
    char   *opt_selected_clusters_str = NULL;    /* point within argv[] */
    int    *selected_clusters_array   = NULL;    /* will allocate */
    int     n_selected_clusters       = 0;
    int     pixel_size                = 1;


    /* set various non-default malloc options */
#ifdef __GLIBC__

#if (__GLIBC__ >= 2 && __GLIBC_MINOR__ >= 15)
    mallopt(M_ARENA_MAX, 1);   /* we're not worried about thread contention */
#endif

#if (__GLIBC__ >= 2 && __GLIBC_MINOR__ >=  3)
    mallopt(M_MXFAST, 0);      /* reduce memory footprint and fragmentation */
#endif

#if 1
    /* size of pointer or double = 8 bytes */
    mallopt(M_TRIM_THRESHOLD, 8);
#endif

#if 1
    /* default of 65536 may be too small
     *
     * 2147483648 acts like it may be disabling MMAP allocations,
     * potentially represented by signed int32 that goes negative at 2^31
     */
    mallopt(M_MMAP_MAX, 2147483647);
#endif

#endif


    /* initialize FreeImage library stuff */
    FreeImage_Initialise(0);
    FreeImage_SetOutputMessage(FreeImageErrorHandler);


    target_n_pixel_clusters    = 22;    /* (7*n + 1); n=3    */
    target_n_mz_clusters       = 23;    /* 6x4 = 24 = 1 + 23 */
    opt.clusters_no_merge_flag = 1;


    if (argc > 1)
    {
        for (i = 1; i < argc; i++)
        {
            if (strncmp(argv[i], "--", 2) == 0)
            {
                if (strcmp(argv[i], "--reverse-palette") == 0)
                {
                    opt_reverse_palette = 1;
                }
                else if (strncmp(argv[i], "--nclusters-pixel=",
                         strlen("--nclusters-pixel=")) == 0 &&
                         strlen(argv[i]) > strlen("--nclusters-pixel="))
                {
                    target_n_pixel_clusters =
                        atol(argv[i] + strlen("--nclusters-pixel="));
                }
                else if (strncmp(argv[i], "--nclusters-mz=",
                         strlen("--nclusters-mz=")) == 0 &&
                         strlen(argv[i]) > strlen("--nclusters-mz="))
                {
                    target_n_mz_clusters =
                        atol(argv[i] + strlen("--nclusters-mz="));
                }
                else if (strncmp(argv[i], "--mz-contact=",
                         strlen("--mz-contact=")) == 0 &&
                         strlen(argv[i]) > strlen("--mz-contact="))
                {
                    opt_selected_clusters_str =
                        argv[i] + strlen("--mz-contact=");
                }
                else if (strncmp(argv[i], "--pixel-size=",
                         strlen("--pixel-size=")) == 0 &&
                         strlen(argv[i]) > strlen("--pixel-size="))
                {
                    pixel_size = atol(argv[i] + strlen("--pixel-size="));
                    
                    if (pixel_size < 1)
                        pixel_size = 1;
                }
                else if (strcmp(argv[i], "--cluster-no-merge") == 0)
                {
                    opt.clusters_no_merge_flag = 1;
                }
                /* merge small clusters and penalize them */
                else if (strcmp(argv[i], "--cluster-merge") == 0)
                {
                    opt.clusters_no_merge_flag = 0;
                }
                /* merge small clusters and penalize them */
                else if (strcmp(argv[i], "--mz-images") == 0)
                {
                    opt_mz_images = 1;
                }
                else
                {
                    fprintf(stderr, "ABORT -- unknown option %s\n",
                            argv[i]);
                
                    error_flag = 1;
                }
            }
            else
            {
                if (num_files == 0)
                {
                    data_file_name      = strdup(argv[i]);
                }
                else if (num_files == 1)
                {
                    tree_file_name_rows = strdup(argv[i]);
                }
                else if (num_files == 2)
                {
                    tree_file_name_cols = strdup(argv[i]);
                }
                else if (num_files == 3)
                {
                    output_file_prefix  = strdup(argv[i]);
                }
                
                num_files++;
            }
        }
    }
    else
    {
        error_flag = 1;
    }


    /* parse any selected clusters for m/z contact sheets */
    if (opt_selected_clusters_str)
    {
        selected_clusters_array =
            parse_selected_clusters_string(opt_selected_clusters_str,
                                           &n_selected_clusters);
    }

    
    if (error_flag)
    {
        printf("Usage: maldi_program [options] unlogged_data.txt mz_tree.tree pixel_tree.tree [output_file_prefix]\n");
        printf("  Options:\n");
        printf("    --cluster-merge     merge and penalize overly small clusters\n");
        printf("    --cluster-no-merge  do not merge/penalize overly small clusters (default)\n");
        printf("    --mz-contact=\"\"     comma-delimited cluster list for m/z contact sheets\n");
        printf("    --mz-images         output individual m/z cluster images\n");
        printf("    --nclusters-mz=N    target number of m/z clusters\n");
        printf("    --nclusters-pixel=N target number of pixel clusters\n");
        printf("    --pixel-size=N      NxN box to draw around each pixel\n");
        printf("    --reverse-palette   reverse color palette in output images\n");

        exit(1);
    }
    
    
    /* sanity check, don't allow more than 255 pixel colors */
    if (target_n_pixel_clusters > 255)
        target_n_pixel_clusters = 255;
    
    
    /* create default output file prefix if none was provided */
    if (output_file_prefix == NULL)
    {
        output_file_prefix = strdup(data_file_name);
        
        /* strip common file extensions */
        if ((sptr = strstr(output_file_prefix, ".txt")))
        {
            *sptr = '\0';
        }
        else if ((sptr = strstr(output_file_prefix, ".tsv")))
        {
            *sptr = '\0';
        }
        else if ((sptr = strstr(output_file_prefix, ".csv")))
        {
            *sptr = '\0';
        }
        if ((sptr = strstr(output_file_prefix, ".TXT")))
        {
            *sptr = '\0';
        }
        else if ((sptr = strstr(output_file_prefix, ".TSV")))
        {
            *sptr = '\0';
        }
        else if ((sptr = strstr(output_file_prefix, ".CSV")))
        {
            *sptr = '\0';
        }
    }
    

    data_matrix = read_data_matrix(data_file_name, data_matrix,
                                   &row_name_array, &col_name_array,
                                   &num_rows, &num_cols, &num_cols_chunk);

    tree_string_rows = read_in_tree_string(tree_file_name_rows);
    tree_string_cols = read_in_tree_string(tree_file_name_cols);

    tree_root_rows = create_tree_from_string(tree_string_rows,
                                    &node_ptr_array_rows, &num_nodes_rows,
                                    &leaf_ptr_array_rows, &num_leaves_rows);
    tree_root_cols = create_tree_from_string(tree_string_cols,
                                    &node_ptr_array_cols, &num_nodes_cols,
                                    &leaf_ptr_array_cols, &num_leaves_cols);

    bless_tree(tree_root_rows);
    bless_tree(tree_root_cols);

    calc_node_clevels(leaf_ptr_array_rows, num_leaves_rows);
    calc_node_clevels(leaf_ptr_array_cols, num_leaves_cols);


    /* allocate mapping stuff, now that everything is read in */
    node_map_array_rows = (struct node_data_map *) calloc(num_rows,
                              sizeof(struct node_data_map));
    node_map_array_cols = (struct node_data_map *) calloc(num_cols,
                              sizeof(struct node_data_map));

    fprintf(stderr, "Mapping data to trees...\n");
    map_leaves_to_data(data_matrix, row_name_array, leaf_ptr_array_rows,
                       node_map_array_rows,
                       num_rows, num_leaves_rows);
    map_leaves_to_data(data_matrix, col_name_array, leaf_ptr_array_cols,
                       node_map_array_cols,
                       num_cols, num_leaves_cols);

    /* parse x#y# image coordinates from pixel names
     * must be done BEFORE sorting the node_map arrays
     */
    parse_maldi_xy(col_name_array, node_map_array_cols, num_cols);

    /* sort tree mappings by tree order
     * must be done AFTER parsing the X,Y coordinates
     *
     * even with 222k pixels, sorting structures is virtually instantaneous,
     * so there is no need sort a pointer array for any minor speed increase
     */
    qsort(node_map_array_rows, num_rows, sizeof(struct node_data_map),
          cmp_sort_mapping_by_node);

    fprintf(stderr, "Finished mapping data to trees\n");




    /* assign clusters to pixel nodes */
    /* unclustered regions still need to be dealt with later */
    scan_clusters(tree_root_cols, node_ptr_array_cols, num_nodes_cols,
                  0.025, 5, 1.0, 1.0,
                  target_n_pixel_clusters, &n_pixel_clusters, &opt, 0);

#if 1
    /* shrink the target number of clusters if we go over */
    while (1)
    {
        scan_clusters(tree_root_rows, node_ptr_array_rows, num_nodes_rows,
                      0.025, 5, 1.0, 1.0,
                      target_n_mz_clusters, &n_mz_clusters, &opt, 0);
        
        /* we've reached or undershot our target, so everything is fine */
        if (n_mz_clusters <= target_n_mz_clusters)
            break;

        /* we've overshot, lower the target number and try again */
        target_n_mz_clusters--;
    }
#else
    scan_clusters(tree_root_rows, node_ptr_array_rows, num_nodes_rows,
                  0.025, 5, 1.0, 1.0,
                  target_n_mz_clusters, &n_mz_clusters, &opt, 0);
#endif


    /* merge good clusters and unassigned regions into final clusters */
    output_clusters_multi(tree_root_cols, 0, NULL);
    output_clusters_multi(tree_root_rows, 0, NULL);

    fprintf(stderr, "Number of col clusters:\t%ld\n",
        (long) n_pixel_clusters);
    fprintf(stderr, "Number of row clusters:\t%ld\n",
        (long) n_mz_clusters);


    /* allocate data for pixel geometric mean cluster spectra */
    col_cluster_spectra = (double **) calloc(n_pixel_clusters,
                                             sizeof(double *));
    for (i = 0; i < n_pixel_clusters; i++)
        col_cluster_spectra[i] = (double *) calloc(num_leaves_rows,
                                                   sizeof(double));


    /* generate pixel cluster spectra geometric means */
    fill_pixel_clusters_means(data_matrix,
                              col_cluster_spectra,
                              node_map_array_cols, num_leaves_cols,
                              node_map_array_rows, num_leaves_rows,
                              n_pixel_clusters);


    /* guess whether the cluster colors should be reversed or not */
    guess_reverse = guess_reverse_palette(data_matrix, node_map_array_cols,
                                          num_leaves_cols, n_pixel_clusters);
    fprintf(stderr, "Guess reverse palette:\t%d\n", (int) guess_reverse);

    /* flip the reverse palette setting */
    if (guess_reverse)
    {
        if (opt_reverse_palette)
            opt_reverse_palette = 0;
        else
            opt_reverse_palette = 1;
    }
    
    
    /* allocate pixel colors array */
    pixel_cluster_colors =
        (uint32_t *) calloc(n_pixel_clusters, sizeof(uint32_t));


    /* render the images */
    render_maldi_images(data_matrix,
                        node_map_array_cols, num_leaves_cols,
                        node_map_array_rows, num_leaves_rows,
                        n_pixel_clusters,    n_mz_clusters,
                        pixel_cluster_colors,
                        output_file_prefix,  opt_reverse_palette,
                        opt_mz_images,
                        selected_clusters_array,
                        n_selected_clusters, pixel_size);


    /* allocate row stats arrays */
    row_min_array  = (double *) malloc(num_rows * sizeof(double));
    row_max_array  = (double *) malloc(num_rows * sizeof(double));
    row_mean_array = (double *) malloc(num_rows * sizeof(double));
    row_sd_array   = (double *) malloc(num_rows * sizeof(double));

    /* calculate stats for row cluster output */
    fill_log2_row_stats(data_matrix, num_rows, num_cols,
                        row_min_array,  row_max_array,
                        row_mean_array, row_sd_array);


    /* sort columns after image generation, for improved cache hits earlier */
    qsort(node_map_array_cols, num_cols, sizeof(struct node_data_map),
          cmp_sort_mapping_by_node);


    /* output the clusters, after sorting the mapping arrays */
    write_clusters_to_file(col_cluster_spectra,
                           col_name_array,
                           row_name_array,
                           node_map_array_cols, num_leaves_cols,
                           node_map_array_rows, num_leaves_rows,
                           n_pixel_clusters, n_mz_clusters,
                           pixel_cluster_colors,
                           row_min_array,  row_max_array,
                           row_mean_array, row_sd_array,
                           output_file_prefix);


    /* only output the heatmap if we aren't doing individual cluster sheets */
    if (n_selected_clusters == 0)
    {
        /* transform data matrix to unit variance, now that imaging is finished */
        /* then write the heatmap */
        log2_data(data_matrix, num_rows, num_cols);
        mean_center_uv(data_matrix, num_rows, num_cols);
        write_heatmap(output_file_prefix, data_matrix,
                      node_map_array_cols, num_cols,
                      node_map_array_rows, num_rows, 0.0);
    }


    /* don't free newick string before this without setting it to NULL,
     * otherwise, free_tree_stuff() frees it again
     */
    free_tree_stuff(tree_string_rows, node_ptr_array_rows, num_nodes_rows,
                    leaf_ptr_array_rows, num_leaves_rows, 1);
    free_tree_stuff(tree_string_cols, node_ptr_array_cols, num_nodes_cols,
                    leaf_ptr_array_cols, num_leaves_cols, 1);


    if (data_file_name)
        free(data_file_name);
    if (tree_file_name_rows)
        free(tree_file_name_rows);
    if (tree_file_name_cols)
        free(tree_file_name_cols);
    if (output_file_prefix)
        free(output_file_prefix);


    if (data_matrix)
    {
        for (i = 0; i < num_rows; i++)
            if (data_matrix[i])
                free(data_matrix[i]);
        free(data_matrix);
    }

    if (node_map_array_rows)
        free(node_map_array_rows);
    if (node_map_array_cols)
        free(node_map_array_cols);

    if (pixel_cluster_colors)
        free(pixel_cluster_colors);

    /* free row/col name arrays */
    if (row_name_array)
    {
        for (i = 0; i < num_rows; i++)
            if (row_name_array[i])
                free(row_name_array[i]);

        free(row_name_array);
    }
    if (col_name_array)
    {
        for (i = 0; i < num_cols; i++)
            if (col_name_array[i])
                free(col_name_array[i]);

        free(col_name_array);
    }

    /* free row stats array */
    if (row_min_array)  free(row_min_array);
    if (row_max_array)  free(row_max_array);
    if (row_mean_array) free(row_mean_array);
    if (row_sd_array)   free(row_sd_array);

    /* free memory allocated by FreeImage */
    FreeImage_DeInitialise();
    
    return 0;
}
