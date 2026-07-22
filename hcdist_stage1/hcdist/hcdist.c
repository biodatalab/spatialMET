#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdint.h>
#include <float.h>
#include <time.h>
#include <inttypes.h>
#include <pthread.h>
#ifdef __GLIBC__
#include <malloc.h>
#endif
#include <sys/mman.h>
#include "text.h"
#include "hcdist.h"

/* TODO -- fix lazy tree jumble/seed parsing */


/* 2025-05-14:  print effective count to STDERR when outputting tree weights
 * 2025-05-06:  more --unit-stretch=N bug fixes
 * 2025-05-05:  add --unit-stretch=N stretch [min, fraction Nth precentile]
 * 2025-05-03:  --filter-present=N uses global min cutoff if detects imputed
 * 2025-05-02:  change all-zero peak pixel warning to col/row terminology
 * 2025-05-02:  skip worst distance loop entirely if no missing data
 * 2025-05-01:  minor optimization to worst distance loop after calculations
 * 2025-04-28:  enable breaking edge length ties by cluster size via options
 * 2024-04-24:  update --filter-present argument description
 * 2024-04-07:  use C99 log2() function instead of log() / log2
 * 2025-04-16:  add --depower option
 * 2025-04-11:  several bug fixes to median normaliztion when >= 50% zeroes
 * 2025-04-11:  change wording of all-zero peaks pixel warning
 * 2024-12-19:  print time spent calculating distance matrix
 * 2024-12-03:  changed --minowski to --minkowski to correct typo,
 *              still accept --minowski for backwards compatability
 * 2024-10-01;  --floor-lod-ub now warns on terminal empty columns
 * 2024-07-15:  optimize --minowski=# for p=1 and p=2 (City Block, Euclidian)
 * 2024-07-10:  edit usage statement, RMSD --> RM#D to reflect Minowski
 * 2024-07-10:  add some code for testing leaf reordering methods
 * 2024-07-08:  change flag description for *2 and *u linkage methods
 * 2024-07-08:  optimize p=1.25 and p=1.75 for --distpow and --minowski
 * 2024-07-07:  add --distpow=#  and optimize --distpow=1.5
 * 2024-07-07:  add --minowski=# and optimize --minowski=1.5
 * 2024-07-02:  add --floor-value=N
 * 2024-07-01:  add --floor-lod-ub
 * 2024-06-28:  print estimated lod noise floor statistics
 * 2024-06-20:  don't print currently used non-default mallopt() options
 * 2024-06-20:  add madvise(row_ptr, ...) to pre-processing functions
 * 2024-06-20:  move madvise(row_ptr1, ...) inside loop
 * 2024-06-19:  hopefully fix perma-sleep/deadlocked threads this time
 * 2024-06-18:  posix_madvise(,,POSIX_MADV_WILLNEED) in distance calculations
 * 2024-06-18:  non-default malloc() settings to save footprint/fragmentation
 * 2024-06-17:  sleep main thread to hopefully avoid mutex deadlock
 * 2024-06-13:  minor optimization to waiting on threads loop
 * 2024-06-12:  fix memory leak in filtering, mutex dereference/ptr mistakes
 * 2024-06-11:  add other (disabled) options for mean/median of medians
 * 2024-06-10:  add --floor-to-lod
 * 2024-06-10:  change a few 0 to 0.0 to encourage fewer recasts to (double)
 * 2024-06-10:  begin cleaning up various MISSING checks
 * 2024-06-10:  convert DBL_MAX to MISSING, remove most DBL_MAX code
 * 2024-06-08:  add experimental --impute-row-col (not in usage statement)
 * 2024-06-07:  rename --floor-later to --impute-later, applies to both
 * 2024-06-07:  add --impute-global, rename --impute to --impute-col-min
 * 2024-06-06:  output data: choose shorter of printf %.6f or %.14g
 * 2024-06-06:  add --floor-to-one --floor-later
 * 2024-06-05:  experimental median norm code for >= half zeroes, off for now
 * 2024-06-03:  add ,M option to --filter-present=N,M
 * 2024-06-01:  disable cosine->Euclidian optimization when --transpose-last
 * 2024-06-01:  add --cityblock Manhattan distances divided by n
 * 2024-05-31:  don't calculate unused Euclidian dist with weighted correlation
 * 2024-05-31:  use correct weighted cosine/pearson functions for missing data
 * 2024-05-30:  --row-weights=file
 * 2024-05-28:  add --tree-flip-size --tree-flip-edge
 * 2024-05-23:  check for and skip some processing of empty rows
 * 2024-05-23:  include empty rows in counts of scanned data
 * 2024-05-22:  add --filter-present, change <= N filters to < N
 * 2024-05-21:  various cleanup and features for release
 * 2024-05-16:  add --transpose-first --transpose-last flags
 * 2024-05-16:  add --log2 --impute --norm-median flags
 * 2024-05-10:  add --spatial flag
 * 2024-05-10:  added pthreads multithreading to distance calculations
 * 2024-05-10:  various new flags and minor changes for MALDI clustering
 * 2024-04-30:  move read_data_matrix() to text.c
 * 2024-04-29:  add --geomean flag to usage statement
 * 2023-09-07:  tweak --similarity distance to auto-cluster better
 * 2023-08-02:  rename height to clevel (cluster level above leaves)
 * 2023-07-26:  only calculate Euclidian/RMSD when we need it
 * 2023-07-26:  add no_missing versions of cosine/pearson calculations
 * 2023-07-14:  change all int to int32_t
 * 2023-07-14:  better integrate hcdist.c and tree.c
 * 2023-07-14:  auto-clustering, memory freeing cleanup
 * 2023-07-03:  add --mostly-trig; 999:1 trig:euclidian to avoid ties
 * 2023-02-23:  add --unit-max-mag to scale to unit 1 maximum magnitude
 * 2021-04-09:  fix file input error on blank initial field (ie. R output)
 * 2021-02-15:  additional memory free'ing clean up after printing output
 * 2021-01-29:  remove some checks for n == 0, since we already return early
 * 2021-01-28:  remove now-invalid write of transposed worst_dist values (bug)
 * 2021-01-27:  rewrite similarity a,b,c,d calculations to be more readable
 * 2021-01-27:  only apply various correlation speed-ups when no missing data
 * 2021-01-27:  calculate uv cosine similarity using Euclidian equation
 * 2020-12-23:  create options structure, add --geomean
 * 2020-12-23:  print ddd:hh:mm:ss format time elapsed/remaining strings
 * 2020-12-22:  print TimeLeft estimate during distance calculations
 * 2020-12-16:  improve --help descriptions
 * 2020-12-10:  free input matrix prior to building tree
 *              free data matrix as distance matrix is calculated
 * 2020-12-09:  use lower half triangular matrix to save memory
 *              extend precision to printf %.13g so the difference between
 *              direct tree building and intermediate distance matrix output
 *              file is lessened (identical results for 1000 row test case)
 * 2020-11-10:  dealt with blank input vector rows, blank trailing columns;
 *              fixed pearson/cosine functions to have less round-off error
 *              1E-16 error turns into 1E-8 with sqrt transform, so avoiding
 *              round-off error is quite important
 */

#define MEM_OVERHEAD 1.01    /* speed hack -- overallocate to avoid reallocs */

#define MISSING -DBL_MIN


pthread_mutex_t mutex_wait_main;
pthread_cond_t  cond_wait_main;
int             n_threads_left;


void free_filenames(char **filename_array, int32_t num_files)
{
    int32_t i;

    if (filename_array)
    {
        for (i = 0; i < num_files; i++)
            if (filename_array[i])
                free(filename_array[i]);

        free(filename_array);
    }

    return;
}


double calc_pearson_r(double *array1, double *array2, int32_t n)
{
  double x_avg = 0.0, y_avg = 0.0;
  double sum_xy_diff = 0.0, sum_x2_diff = 0.0, sum_y2_diff = 0.0;
  double x_diff, y_diff;
  double temp;
  int32_t i;
  int32_t count = 0;

  for (i = 0; i < n; i++)
  {
    /* skip missing */
    if (array1[i] != MISSING && array2[i] != MISSING)
    {
      x_avg += array1[i];
      y_avg += array2[i];

      count++;
    }
  }
  if (count)
  {
    x_avg /= count;
    y_avg /= count;
  }

  for (i = 0; i < n; i++)
  {
    /* skip missing */
    if (array1[i] != MISSING && array2[i] != MISSING)
    {
      x_diff = array1[i] - x_avg;
      sum_x2_diff += x_diff * x_diff;

      y_diff = array2[i] - y_avg;
      sum_y2_diff += y_diff * y_diff;

      sum_xy_diff += x_diff * y_diff;
    }
  }


  if (sum_x2_diff && sum_y2_diff)
  {
    /* avoid further round off error entirely, return +/- 1 */
    if (sum_x2_diff == sum_y2_diff && fabs(sum_x2_diff) == sum_xy_diff)
    {
        if (sum_xy_diff >= 0)
            return 1.0;
        return -1.0;
    }
  
    /* sqrt(A*B), instead of sqrt(A)*sqrt*(B), avoids round off errors
     * However, it can also lead to overflows...
     *
     * So, we'll sqrt sum_xy_diff, then square it again, to put the errors
     * on an even footing
     */
#if 0
    temp         = sqrt(sum_x2_diff * sum_y2_diff);
#else
    sum_x2_diff  = sqrt(sum_x2_diff);
    sum_y2_diff  = sqrt(sum_y2_diff);
    temp         = sum_x2_diff * sum_y2_diff;

    if (sum_xy_diff >= 0)
    {
        sum_xy_diff = sqrt(sum_xy_diff);
        sum_xy_diff = sum_xy_diff * sum_xy_diff;
    }
    else
    {
        sum_xy_diff = sqrt(-sum_xy_diff);
        sum_xy_diff = -sum_xy_diff * sum_xy_diff;
    }
#endif

    if (temp)
    {
      temp = sum_xy_diff / temp;
      
      /* round off errors can lead to slightly greater than 1 */
      if (temp >  1.0)
        return    1.0;
      if (temp < -1.0)
        return   -1.0;

      return temp;
    }
  }

  return 0;
}


/* only use when no missing/weak data */
double calc_pearson_r_no_missing(double *array1, double *array2, int32_t n)
{
  double x_avg = 0.0, y_avg = 0.0;
  double sum_xy_diff = 0.0, sum_x2_diff = 0.0, sum_y2_diff = 0.0;
  double x_diff, y_diff;
  double temp;
  int32_t i;
  int32_t count = 0;

  for (i = 0; i < n; i++)
  {
    x_avg += array1[i];
    y_avg += array2[i];
    
    count++;
  }
  if (count)
  {
    x_avg /= count;
    y_avg /= count;
  }

  for (i = 0; i < n; i++)
  {
    x_diff = array1[i] - x_avg;
    sum_x2_diff += x_diff * x_diff;

    y_diff = array2[i] - y_avg;
    sum_y2_diff += y_diff * y_diff;

    sum_xy_diff += x_diff * y_diff;
  }


  if (sum_x2_diff && sum_y2_diff)
  {
    /* avoid further round off error entirely, return +/- 1 */
    if (sum_x2_diff == sum_y2_diff && fabs(sum_x2_diff) == sum_xy_diff)
    {
        if (sum_xy_diff >= 0)
            return 1.0;
        return -1.0;
    }
  
    /* sqrt(A*B), instead of sqrt(A)*sqrt*(B), avoids round off errors
     * However, it can also lead to overflows...
     *
     * So, we'll sqrt sum_xy_diff, then square it again, to put the errors
     * on an even footing
     */
#if 0
    temp         = sqrt(sum_x2_diff * sum_y2_diff);
#else
    sum_x2_diff  = sqrt(sum_x2_diff);
    sum_y2_diff  = sqrt(sum_y2_diff);
    temp         = sum_x2_diff * sum_y2_diff;

    if (sum_xy_diff >= 0)
    {
        sum_xy_diff = sqrt(sum_xy_diff);
        sum_xy_diff = sum_xy_diff * sum_xy_diff;
    }
    else
    {
        sum_xy_diff = sqrt(-sum_xy_diff);
        sum_xy_diff = -sum_xy_diff * sum_xy_diff;
    }
#endif

    if (temp)
    {
      temp = sum_xy_diff / temp;
      
      /* round off errors can lead to slightly greater than 1 */
      if (temp >  1.0)
        return    1.0;
      if (temp < -1.0)
        return   -1.0;

      return temp;
    }
  }

  return 0;
}


double calc_pearson_r_weighted(double *array1, double *array2,
                               double *weights, int32_t n)
{
  double x_avg = 0.0, y_avg = 0.0;
  double sum_xy_diff = 0.0, sum_x2_diff = 0.0, sum_y2_diff = 0.0;
  double x_diff, y_diff;
  double temp;
  int32_t i;
  double w_count = 0.0;     /* weighted count, sum of weights */

  for (i = 0; i < n; i++)
  {
    /* skip missing */
    if (array1[i] != MISSING && array2[i] != MISSING)
    {
      x_avg   += weights[i] * array1[i];
      y_avg   += weights[i] * array2[i];

      w_count += weights[i];
    }
  }
  if (w_count)
  {
    x_avg /= w_count;
    y_avg /= w_count;
  }

  for (i = 0; i < n; i++)
  {
    /* skip missing */
    if (array1[i] != MISSING && array2[i] != MISSING)
    {
      x_diff = array1[i] - x_avg;
      sum_x2_diff += weights[i] * x_diff * x_diff;

      y_diff = array2[i] - y_avg;
      sum_y2_diff += weights[i] * y_diff * y_diff;

      sum_xy_diff += weights[i] * x_diff * y_diff;
    }
  }

  if (sum_x2_diff && sum_y2_diff)
  {
    /* avoid further round off error entirely, return +/- 1 */
    if (sum_x2_diff == sum_y2_diff && fabs(sum_x2_diff) == sum_xy_diff)
    {
        if (sum_xy_diff >= 0)
            return 1.0;
        return -1.0;
    }

    /* sqrt(A*B), instead of sqrt(A)*sqrt*(B), avoids round off errors
     * However, it can also lead to overflows...
     *
     * So, we'll sqrt sum_xy_diff, then square it again, to put the errors
     * on an even footing
     */
#if 0
    temp         = sqrt(sum_x2_diff * sum_y2_diff);
#else
    sum_x2_diff  = sqrt(sum_x2_diff);
    sum_y2_diff  = sqrt(sum_y2_diff);
    temp         = sum_x2_diff * sum_y2_diff;

    if (sum_xy_diff >= 0)
    {
        sum_xy_diff = sqrt(sum_xy_diff);
        sum_xy_diff = sum_xy_diff * sum_xy_diff;
    }
    else
    {
        sum_xy_diff = sqrt(-sum_xy_diff);
        sum_xy_diff = -sum_xy_diff * sum_xy_diff;
    }
#endif
  
    if (temp)
    {
      temp = sum_xy_diff / temp;
      
      /* round off errors can lead to slightly greater than 1 */
      if (temp >  1.0)
        return    1.0;
      if (temp < -1.0)
        return   -1.0;
      
      return temp;
    }
  }

  return 0;
}


double calc_pearson_r_weighted_no_missing(double *array1, double *array2,
                                          double *weights, int32_t n)
{
  double x_avg = 0.0, y_avg = 0.0;
  double sum_xy_diff = 0.0, sum_x2_diff = 0.0, sum_y2_diff = 0.0;
  double x_diff, y_diff;
  double temp;
  int32_t i;
  double w_count = 0.0;     /* weighted count, sum of weights */

  for (i = 0; i < n; i++)
  {
    x_avg += weights[i] * array1[i];
    y_avg += weights[i] * array2[i];
    
    w_count += weights[i];
  }
  if (w_count)
  {
    x_avg /= w_count;
    y_avg /= w_count;
  }

  for (i = 0; i < n; i++)
  {
    x_diff = array1[i] - x_avg;
    sum_x2_diff += weights[i] * x_diff * x_diff;

    y_diff = array2[i] - y_avg;
    sum_y2_diff += weights[i] * y_diff * y_diff;

    sum_xy_diff += weights[i] * x_diff * y_diff;
  }

  if (sum_x2_diff && sum_y2_diff)
  {
    /* avoid further round off error entirely, return +/- 1 */
    if (sum_x2_diff == sum_y2_diff && fabs(sum_x2_diff) == sum_xy_diff)
    {
        if (sum_xy_diff >= 0)
            return 1.0;
        return -1.0;
    }

    /* sqrt(A*B), instead of sqrt(A)*sqrt*(B), avoids round off errors
     * However, it can also lead to overflows...
     *
     * So, we'll sqrt sum_xy_diff, then square it again, to put the errors
     * on an even footing
     */
#if 0
    temp         = sqrt(sum_x2_diff * sum_y2_diff);
#else
    sum_x2_diff  = sqrt(sum_x2_diff);
    sum_y2_diff  = sqrt(sum_y2_diff);
    temp         = sum_x2_diff * sum_y2_diff;

    if (sum_xy_diff >= 0)
    {
        sum_xy_diff = sqrt(sum_xy_diff);
        sum_xy_diff = sum_xy_diff * sum_xy_diff;
    }
    else
    {
        sum_xy_diff = sqrt(-sum_xy_diff);
        sum_xy_diff = -sum_xy_diff * sum_xy_diff;
    }
#endif
  
    if (temp)
    {
      temp = sum_xy_diff / temp;
      
      /* round off errors can lead to slightly greater than 1 */
      if (temp >  1.0)
        return    1.0;
      if (temp < -1.0)
        return   -1.0;
      
      return temp;
    }
  }

  return 0;
}


/* cosine between two vectors */
double calc_cosine(double *array1, double *array2, int32_t n)
{
  double sum_xy_diff = 0.0, sum_x2_diff = 0.0, sum_y2_diff = 0.0;
  double temp;
  int32_t i;

  for (i = 0; i < n; i++)
  {
    /* skip missing */
    if (array1[i] != MISSING && array2[i] != MISSING)
    {
      sum_x2_diff += array1[i] * array1[i];
      sum_y2_diff += array2[i] * array2[i];
      sum_xy_diff += array1[i] * array2[i];
    }
  }

  if (sum_x2_diff && sum_y2_diff)
  {
    /* avoid further round off error entirely, return +/- 1 */
    if (sum_x2_diff == sum_y2_diff && fabs(sum_x2_diff) == sum_xy_diff)
    {
        if (sum_xy_diff >= 0)
            return 1.0;
        return -1.0;
    }

    /* sqrt(A*B), instead of sqrt(A)*sqrt*(B), avoids round off errors
     * However, it can also lead to overflows...
     *
     * So, we'll sqrt sum_xy_diff, then square it again, to put the errors
     * on an even footing
     */
#if 0
    temp         = sqrt(sum_x2_diff * sum_y2_diff);
#else
    sum_x2_diff  = sqrt(sum_x2_diff);
    sum_y2_diff  = sqrt(sum_y2_diff);
    temp         = sum_x2_diff * sum_y2_diff;

    if (sum_xy_diff >= 0)
    {
        sum_xy_diff = sqrt(sum_xy_diff);
        sum_xy_diff = sum_xy_diff * sum_xy_diff;
    }
    else
    {
        sum_xy_diff = sqrt(-sum_xy_diff);
        sum_xy_diff = -sum_xy_diff * sum_xy_diff;
    }
#endif
  
    if (temp)
    {
      temp = sum_xy_diff / temp;
      
      /* round off errors can lead to slightly greater than 1 */
      if (temp >  1.0)
        return    1.0;
      if (temp < -1.0)
        return   -1.0;
      
      return temp;
    }
  }

  return 0;
}


/* only use when no missing/weak data */
double calc_cosine_no_missing(double *array1, double *array2, int32_t n)
{
  double sum_xy_diff = 0.0, sum_x2_diff = 0.0, sum_y2_diff = 0.0;
  double temp;
  int32_t i;

  for (i = 0; i < n; i++)
  {
    sum_x2_diff += array1[i] * array1[i];

    sum_y2_diff += array2[i] * array2[i];

    sum_xy_diff += array1[i] * array2[i];
  }

  if (sum_x2_diff && sum_y2_diff)
  {
    /* avoid further round off error entirely, return +/- 1 */
    if (sum_x2_diff == sum_y2_diff && fabs(sum_x2_diff) == sum_xy_diff)
    {
        if (sum_xy_diff >= 0)
            return 1.0;
        return -1.0;
    }

    /* sqrt(A*B), instead of sqrt(A)*sqrt*(B), avoids round off errors
     * However, it can also lead to overflows...
     *
     * So, we'll sqrt sum_xy_diff, then square it again, to put the errors
     * on an even footing
     */
#if 0
    temp         = sqrt(sum_x2_diff * sum_y2_diff);
#else
    sum_x2_diff  = sqrt(sum_x2_diff);
    sum_y2_diff  = sqrt(sum_y2_diff);
    temp         = sum_x2_diff * sum_y2_diff;

    if (sum_xy_diff >= 0)
    {
        sum_xy_diff = sqrt(sum_xy_diff);
        sum_xy_diff = sum_xy_diff * sum_xy_diff;
    }
    else
    {
        sum_xy_diff = sqrt(-sum_xy_diff);
        sum_xy_diff = -sum_xy_diff * sum_xy_diff;
    }
#endif
  
    if (temp)
    {
      temp = sum_xy_diff / temp;
      
      /* round off errors can lead to slightly greater than 1 */
      if (temp >  1.0)
        return    1.0;
      if (temp < -1.0)
        return   -1.0;
      
      return temp;
    }
  }

  return 0;
}


/* cosine between two vectors */
double calc_cosine_weighted(double *array1, double *array2,
                            double *weights, int32_t n)
{
  double sum_xy_diff = 0.0, sum_x2_diff = 0.0, sum_y2_diff = 0.0;
  double temp;
  int32_t i;

  for (i = 0; i < n; i++)
  {
    /* skip missing */
    if (array1[i] != MISSING && array2[i] != MISSING)
    {
      sum_x2_diff += weights[i] * array1[i] * array1[i];
      sum_y2_diff += weights[i] * array2[i] * array2[i];
      sum_xy_diff += weights[i] * array1[i] * array2[i];
    }
  }

  if (sum_x2_diff && sum_y2_diff)
  {
    /* avoid further round off error entirely, return +/- 1 */
    if (sum_x2_diff == sum_y2_diff && fabs(sum_x2_diff) == sum_xy_diff)
    {
        if (sum_xy_diff >= 0)
            return 1.0;
        return -1.0;
    }

    /* sqrt(A*B), instead of sqrt(A)*sqrt*(B), avoids round off errors
     * However, it can also lead to overflows...
     *
     * So, we'll sqrt sum_xy_diff, then square it again, to put the errors
     * on an even footing
     */
#if 0
    temp         = sqrt(sum_x2_diff * sum_y2_diff);
#else
    sum_x2_diff  = sqrt(sum_x2_diff);
    sum_y2_diff  = sqrt(sum_y2_diff);
    temp         = sum_x2_diff * sum_y2_diff;

    if (sum_xy_diff >= 0)
    {
        sum_xy_diff = sqrt(sum_xy_diff);
        sum_xy_diff = sum_xy_diff * sum_xy_diff;
    }
    else
    {
        sum_xy_diff = sqrt(-sum_xy_diff);
        sum_xy_diff = -sum_xy_diff * sum_xy_diff;
    }
#endif
  
    if (temp)
    {
      temp = sum_xy_diff / temp;
      
      /* round off errors can lead to slightly greater than 1 */
      if (temp >  1.0)
        return    1.0;
      if (temp < -1.0)
        return   -1.0;
      
      return temp;
    }
  }

  return 0;
}


/* cosine between two vectors */
/* only use when no missing/weak data */
double calc_cosine_weighted_no_missing(double *array1, double *array2,
                                       double *weights, int32_t n)
{
  double sum_xy_diff = 0.0, sum_x2_diff = 0.0, sum_y2_diff = 0.0;
  double temp;
  int32_t i;

  for (i = 0; i < n; i++)
  {
    sum_x2_diff += weights[i] * array1[i] * array1[i];

    sum_y2_diff += weights[i] * array2[i] * array2[i];

    sum_xy_diff += weights[i] * array1[i] * array2[i];
  }

  if (sum_x2_diff && sum_y2_diff)
  {
    /* avoid further round off error entirely, return +/- 1 */
    if (sum_x2_diff == sum_y2_diff && fabs(sum_x2_diff) == sum_xy_diff)
    {
        if (sum_xy_diff >= 0)
            return 1.0;
        return -1.0;
    }

    /* sqrt(A*B), instead of sqrt(A)*sqrt*(B), avoids round off errors
     * However, it can also lead to overflows...
     *
     * So, we'll sqrt sum_xy_diff, then square it again, to put the errors
     * on an even footing
     */
#if 0
    temp         = sqrt(sum_x2_diff * sum_y2_diff);
#else
    sum_x2_diff  = sqrt(sum_x2_diff);
    sum_y2_diff  = sqrt(sum_y2_diff);
    temp         = sum_x2_diff * sum_y2_diff;

    if (sum_xy_diff >= 0)
    {
        sum_xy_diff = sqrt(sum_xy_diff);
        sum_xy_diff = sum_xy_diff * sum_xy_diff;
    }
    else
    {
        sum_xy_diff = sqrt(-sum_xy_diff);
        sum_xy_diff = -sum_xy_diff * sum_xy_diff;
    }
#endif
  
    if (temp)
    {
      temp = sum_xy_diff / temp;
      
      /* round off errors can lead to slightly greater than 1 */
      if (temp >  1.0)
        return    1.0;
      if (temp < -1.0)
        return   -1.0;
      
      return temp;
    }
  }

  return 0;
}


void flag_missing(double **data_matrix, int32_t num_rows, int32_t num_cols,
                  int ignore_weak_flag)
{
    double *row_ptr;
    double value;
    int32_t row, col;
    
    for (row = 0; row < num_rows; row++)
    {
        row_ptr = data_matrix[row];
        posix_madvise(row_ptr, num_cols * sizeof(double),
                      POSIX_MADV_WILLNEED);
        
        for (col = 0; col < num_cols; col++)
        {
            value = row_ptr[col];
        
            if (value == DBL_MAX ||
                (ignore_weak_flag && fabs(value) < (1.0 - 1E-5)))
            {
                row_ptr[col] = MISSING;
            }
        }
    }
}


void mean_center(double **data_matrix, int32_t num_rows, int32_t num_cols)
{
    double *row_ptr;
    double value;
    double avg;
    int32_t row, col;
    int32_t n;
    
    for (row = 0; row < num_rows; row++)
    {
        row_ptr = data_matrix[row];
        posix_madvise(row_ptr, num_cols * sizeof(double),
                      POSIX_MADV_WILLNEED);
        
        avg = 0.0;
        n = 0;
        for (col = 0; col < num_cols; col++)
        {
            value = row_ptr[col];

            if (value != MISSING)
            {
                avg += value;
                n++;
            }
        }
        
        if (n)
            avg /= n;

        /* subtract the mean */
        if (avg)
        {
            for (col = 0; col < num_cols; col++)
            {
                /* skip missing */
                if (row_ptr[col] != MISSING)
                    row_ptr[col] -= avg;
            }
        }
    }
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
        posix_madvise(row_ptr, num_cols * sizeof(double),
                      POSIX_MADV_WILLNEED);
        
        avg = 0.0;
        n = 0;
        for (col = 0; col < num_cols; col++)
        {
            value = row_ptr[col];

            if (value != MISSING)
            {
                avg += value;
                n++;
            }
        }
        
        if (n)
            avg /= n;

        sd = 0;
        n = 0;
        for (col = 0; col < num_cols; col++)
        {
            value = row_ptr[col];

            if (value != MISSING)
            {
                value -= avg;
                sd += value * value;
                n++;
            }
        }
        
        if (n)
            sd = sqrt(sd / n);

        /* subtract the mean */
        if (avg)
        {
            for (col = 0; col < num_cols; col++)
            {
                if (row_ptr[col] != MISSING)
                    row_ptr[col] -= avg;
            }
        }

        /* scale to unit variance */
        if (sd)
        {
            for (col = 0; col < num_cols; col++)
            {
                if (row_ptr[col] != MISSING)
                    row_ptr[col] /= sd;
            }
        }
    }
}


/* scale so maximum magnitude is unit of 1 */
void scale_unit_max_magnitude(double **data_matrix, int32_t num_rows,
                              int32_t num_cols)
{
    double *row_ptr;
    double value;
    double maxmag;
    int32_t row, col;
    
    for (row = 0; row < num_rows; row++)
    {
        row_ptr = data_matrix[row];
        posix_madvise(row_ptr, num_cols * sizeof(double),
                      POSIX_MADV_WILLNEED);

        maxmag = 0.0;
        for (col = 0; col < num_cols; col++)
        {
            value = row_ptr[col];

            if (value != MISSING)
            {
                value = fabs(value);
                if (value > maxmag)
                    maxmag = value;
            }
        }

        /* scale by maximum magnitude */
        if (maxmag)
        {
            for (col = 0; col < num_cols; col++)
            {
                /* skip missing */
                if (row_ptr[col] != MISSING)
                    row_ptr[col] /= maxmag;
            }
        }
    }
}


/* min --> 0, fraction Nth percentile --> 1 */
void stretch_unit_frac(double **data_matrix, int32_t num_rows,
                       int32_t num_cols, double high_frac)
{
    double  *value_array = NULL;
    double  *row_ptr;
    double   value, value_low, value_high;
    double   min_global  = DBL_MAX;
    double   scale;
    int32_t  row, col;
    int32_t  index_low, index_high;
    
    index_high = high_frac * (num_cols - 1);

    if (num_cols)
    {
        value_array = (double *) malloc(num_cols * sizeof(double));
    }

    /* find global minimum */
    for (row = 0; row < num_rows; row++)
    {
        row_ptr = data_matrix[row];
        posix_madvise(row_ptr, num_cols * sizeof(double),
                      POSIX_MADV_WILLNEED);
        
        for (col = 0; col < num_cols; col++)
        {
            value = row_ptr[col];

            if (value < min_global && value != MISSING)
            {
                min_global = value;
            }
        }
    }
    
    for (row = 0; row < num_rows; row++)
    {
        row_ptr = data_matrix[row];
        posix_madvise(row_ptr, num_cols * sizeof(double),
                      POSIX_MADV_WILLNEED);

        memcpy(value_array, row_ptr, num_cols * sizeof(double));

        /* sort the values from low to high */
        qsort(value_array, num_cols, sizeof(double), cmp_double);
        
        /* find first non-missing value */
        value_low = MISSING;
        index_low = 0;
        for (col = 0; col < num_cols; col++)
        {
            value = value_array[col];

            if (value != MISSING)
            {
                value_low = value;
                index_low = col;
                
                break;
            }
        }
        
        /* entire row is missing, skip to next row */
        if (value_low == MISSING)
            continue;
        
        /* take higher of first non-missing and high frac index */
        col = index_low;
        if (index_high > col)
            col = index_high;
        
        /* find first non-min value starting at fraction high */
        value_high = value_low;
        for (; col < num_cols; col++)
        {
            value = value_array[col];

            if (value != min_global)
            {
                value_high = value;
                
                break;
            }
        }
        
        /* scale, after subtracting out min_global */
        scale = 0.0;
        if (value_high - min_global)
            scale = 1.0 / (value_high - min_global);
        
        posix_madvise(row_ptr, num_cols * sizeof(double),
                      POSIX_MADV_WILLNEED);

        for (col = 0; col < num_cols; col++)
        {
            /* skip missing */
            if (row_ptr[col] != MISSING)
            {
                row_ptr[col] = scale * (row_ptr[col] - min_global);
            }
        }
    }
    
    if (value_array)
        free(value_array);
}


/* scale each sample to the same median */
/* assume unlogged data */
void normalize_cols_median(double **data_matrix, int32_t num_rows,
                           int32_t num_cols)
{
    double *col_value_array   = NULL;    /* store values for a single col */
    double *median_array      = NULL;    /* median values for each col */
    double *median_good_array = NULL;    /* used for median of medians */
    double *scale_array       = NULL;    /* scaling factors to normalize with */
    double *row_ptr;
    double  value;
    double  median, min;
    double  avg_median;
    int32_t row, col, n, n_non_min;
    int32_t n_good_medians;
    int32_t index, index_half, index_next;
    int32_t even_n_flag;
    
    median_array      = (double *) malloc(num_cols * sizeof(double));
    median_good_array = (double *) malloc(num_cols * sizeof(double));
    scale_array       = (double *) malloc(num_cols * sizeof(double));
    col_value_array   = (double *) malloc(num_rows * sizeof(double));
    
    avg_median      = 0.0;
    n_good_medians  = 0;
    
    for (col = 0; col < num_cols; col++)
    {
        min         = DBL_MAX;
        n           = 0;
        even_n_flag = 0;

        for (row = 0; row < num_rows; row++)
        {
            value = data_matrix[row][col];

            if (value != MISSING)
            {
                if (value < min)
                    min = value;

                col_value_array[n++] = value;
            }
        }
        

        median = MISSING;


        if (n)
        {
            /* sort the values from low to high */
            qsort(col_value_array, n, sizeof(double), cmp_double);

            index = n >> 1;
        
            /* odd */
            if (n % 2)
            {
                median      = col_value_array[index];
            }
            /* even */
            /* include check below to not use when exactly half zeroes */
            else
            {
                median      = 0.5 * (col_value_array[index-1] +
                                     col_value_array[index]);
                even_n_flag = 1;
            }
            
            /* HACK -- minimum value is untrustworthy;
             *         could be zero, missing, already imputed, noise, etc.
             *
             *         skip some percentage of non-min points, use that
             */
            if (median == min ||
                (even_n_flag && index && col_value_array[index-1] == min))
            {
                index_half = index;

                /* odd */
                if (even_n_flag == 0)
                {
                    index = index + 1;
                }
                /* even */
                /* the next value after half-way is already index */
            
                /* find first value > min */
                index_next = index;
                for (; index < n; index++)
                {
                    if (col_value_array[index] > min)
                    {
                        median     = col_value_array[index];
                        index_next = index;

                        break;
                    }
                }

                /* number of values > min */
                n_non_min = n - index_next;
                
#if 1
                /* set new index to sliding percentile of non-min;
                 * higher percentile the further past index_half
                 */
                index = index_next + n_non_min *
                        (0.5 * (index_next - index_half) /
                         (double) index_half);
                
                if (index >= n)
                    index = n - 1;

                median = col_value_array[index];
#endif

#if 0
                median = sqrt(median * col_value_array[index]);
#endif

#if 0
                /* If imputed based on min, and the new median is the first
                 * non-zero value, all imputed values will be equal after
                 * median normalization.  Go one higher, so that we have
                 * at least *some* variability in the imputed values.
                 *
                 * Do this even if min=0, since the user may want to impute
                 * the data in such a way later.
                 */
                if (median == col_value_array[index_next])
                {
                    for (index = index + 1; index < n; index++)
                    {
                        if (col_value_array[index] > median)
                        {
                            median = col_value_array[index];

                            break;
                        }
                    }
                }
#endif
            }
            else if (median > 0 && median != MISSING)
            {
                median_good_array[n_good_medians++]  = median;
                avg_median                          += log(median);
            }
        }

        median_array[col] = median;
    }
    
    
    /* no good medians, use all medians instead */
    n = n_good_medians;
    if (1 || n_good_medians == 0)
    {
        avg_median = 0.0;
        n          = 0;
        for (col = 0; col < num_cols; col++)
        {
            value = median_array[col];

            /* skip MISSING regardless of ingore_weak here; bogus value */
            /* take geometric mean, assume unlogged input data */
            if (value > 0.0 && value != MISSING)
            {
                avg_median += log(value);
                n++;
            }
        }
    }
    
    if (n)
    {
        avg_median = exp(avg_median / n);
    }
    

#if 0
    /* use median of good medians */
    if (n_good_medians)
    {
        /* sort the values from low to high */
        qsort(median_good_array, n_good_medians, sizeof(double), cmp_double);

        index      = n_good_medians >> 1;
        avg_median = median_good_array[index];
    }
#endif


    fprintf(stderr, "Representative median:\t%lf\n", avg_median);


    /* calculate scaling factors */
    if (avg_median)
    {
        for (col = 0; col < num_cols; col++)
        {
            /* initialize scaling factor to 1 */
            scale_array[col] = 1.0;


            value = median_array[col];

            /* skip MISSING regardless of ingore_weak here; bogus value */
            if (value && value != MISSING)
                scale_array[col] = avg_median / value;
        }
    }


    /* scale with scaling factors */
    if (avg_median)
    {
        for (row = 0; row < num_rows; row++)
        {
            row_ptr = data_matrix[row];
            posix_madvise(row_ptr, num_cols * sizeof(double),
                          POSIX_MADV_WILLNEED);

            for (col = 0; col < num_cols; col++)
            {
                /* skip missing */
                if (row_ptr[col] != MISSING)
                    row_ptr[col] *= scale_array[col];
            }
        }
    }
    

    if (median_array)
        free(median_array);

    if (median_good_array)
        free(median_good_array);

    if (scale_array)
        free(scale_array);

    if (col_value_array)
        free(col_value_array);
}


/* assume unlogged data, assume zeroes are missing */
void impute_global_min_half(double **data_matrix, int32_t num_rows,
                            int32_t num_cols)
{
    double *row_ptr;
    double  value;
    double  min;
    int32_t row, col;

    min = DBL_MAX;
    for (row = 0; row < num_rows; row++)
    {
        row_ptr = data_matrix[row];
        posix_madvise(row_ptr, num_cols * sizeof(double),
                      POSIX_MADV_WILLNEED);

        for (col = 0; col < num_cols; col++)
        {
            value = row_ptr[col];

            if (value != MISSING)
            {
                if (value > 0.0 && value < min)
                    min = value;
            }
        }
    }
    
    /* halve the global minimum */
    min *= 0.5;
    
    /* impute values */
    for (row = 0; row < num_rows; row++)
    {
        row_ptr = data_matrix[row];
        posix_madvise(row_ptr, num_cols * sizeof(double),
                      POSIX_MADV_WILLNEED);

        for (col = 0; col < num_cols; col++)
        {
            value = row_ptr[col];

            /* assume -0 isn't possible here */
            if (value == 0.0 || value == MISSING)
            {
                row_ptr[col] = min;
            }
        }
    }
}


/* assume unlogged data, assume zeroes are missing */
void impute_col_min_half(double **data_matrix, int32_t num_rows,
                         int32_t num_cols)
{
    double *impute_array = NULL;    /* minimum value for each column */
    double *row_ptr;
    double  value;
    double  min;
    int32_t row, col, n;

    impute_array = (double *) malloc(num_cols * sizeof(double));
    
    for (col = 0; col < num_cols; col++)
    {
        min = DBL_MAX;
        n   = 0;

        for (row = 0; row < num_rows; row++)
        {
            value = data_matrix[row][col];

            if (value != MISSING)
            {
                if (value > 0.0 && value < min)
                    min = value;
            }
        }
        
        impute_array[col] = MISSING;
        if (min < DBL_MAX)
            impute_array[col] = 0.5 * min;
    }

    /* impute values */
    for (row = 0; row < num_rows; row++)
    {
        row_ptr = data_matrix[row];
        posix_madvise(row_ptr, num_cols * sizeof(double),
                      POSIX_MADV_WILLNEED);

        for (col = 0; col < num_cols; col++)
        {
            value = row_ptr[col];

            /* assume -0 isn't possible here */
            if (value == 0.0 || value == MISSING)
                row_ptr[col] = impute_array[col];
        }
    }
    
    if (impute_array)
        free(impute_array);
}


/* assume unlogged data, assume zeroes are missing */
void impute_col_global_min_half(double **data_matrix, int32_t num_rows,
                                int32_t num_cols)
{
    double *impute_array = NULL;    /* minimum value for each column */
    double *row_ptr;
    double  value;
    double  min, min_global;
    int32_t row, col, n;

    impute_array = (double *) malloc(num_cols * sizeof(double));
    
    min_global = DBL_MAX;
    for (col = 0; col < num_cols; col++)
    {
        min = DBL_MAX;
        n   = 0;

        for (row = 0; row < num_rows; row++)
        {
            value = data_matrix[row][col];

            if (value != MISSING)
            {
                if (value > 0.0 && value < min)
                {
                    min = value;

                    if (value < min_global)
                        min_global = value;
                }
            }
        }
        
        impute_array[col] = MISSING;
        if (min < DBL_MAX)
            impute_array[col] = min;
    }
    
    /* half geometric mean of global and col min */
    if (min_global < DBL_MAX)
    {
        for (col = 0; col < num_cols; col++)
        {
            if (impute_array[col] != MISSING)
            {
                impute_array[col] =
                    0.5 * sqrt(min_global * impute_array[col]);
            }
        }
    }

    /* impute values */
    for (row = 0; row < num_rows; row++)
    {
        row_ptr = data_matrix[row];
        posix_madvise(row_ptr, num_cols * sizeof(double),
                      POSIX_MADV_WILLNEED);

        for (col = 0; col < num_cols; col++)
        {
            value = row_ptr[col];

            /* assume -0 isn't possible here */
            if (value == 0.0 || value == MISSING)
                row_ptr[col] = impute_array[col];
        }
    }
    
    if (impute_array)
        free(impute_array);
}


/* assume unlogged data, assume zeroes are missing
 *
 * ENSMUSG00000069049.11 Eif2s3y is good mouse test gene in test dataset,
 * 0% of Female samples have a value
 *
 * Col min should be robust to genes that are high in one set of samples and
 * completely absent in another.
 *
 * However, col min may be affected by scaling up dark samples, leading to
 * overly-high col min for scaled-up dark samples.
 *
 * Row min is not robust to 100%/0% issues like sex-linked genes.
 *
 * However, if row min is less than col min (or less than its 3sd), it should
 * be fine, so we can then average the two together to tamp down any col min
 * that might be overly high.
 */
void impute_row_col(double **data_matrix, int32_t num_rows,
                    int32_t num_cols)
{
    double *col_min_array = NULL;    /* minimum value for each column */
    double *row_min_array = NULL;    /* minimum value for each row */
    double *row_ptr;
    double  value, min_row, min_col;
    double  min_global;
    double  mean_col, sd_col, ub_col, diff;   /* all in log units */
    int32_t row, col, n;

    col_min_array = (double *) malloc(num_cols * sizeof(double));
    row_min_array = (double *) malloc(num_rows * sizeof(double));

    for (col = 0; col < num_cols; col++)
        col_min_array[col] = DBL_MAX;
    
    min_global = DBL_MAX;
    for (row = 0; row < num_rows; row++)
    {
        min_row = DBL_MAX;
        row_ptr = data_matrix[row];
    
        for (col = 0; col < num_cols; col++)
        {
            value = row_ptr[col];

            if (value > 0.0 && value != MISSING)
            {
                if (value < min_row)
                {
                    min_row = value;
                
                    if (value < min_global)
                        min_global = value;
                }
                
                if (value < col_min_array[col])
                    col_min_array[col] = value;
            }
        }
        
        row_min_array[row] = MISSING;
        if (min_row < DBL_MAX)
            row_min_array[row] = min_row;
    }
    
    
    mean_col = 0.0;
    n        = 0.0;
    for (col = 0; col < num_cols; col++)
    {
        value = col_min_array[col];
        
        if (value < DBL_MAX)
        {
            mean_col += log(value);
            n++;
        }
        else
        {
            col_min_array[col] = MISSING;
        }
    }
    if (n)
        mean_col /= n;
    
    sd_col   = 0.0;
    for (col = 0; col < num_cols; col++)
    {
        value = col_min_array[col];
        
        if (value != MISSING)
        {
            diff    = log(value) - mean_col;
            sd_col += diff * diff;
        }
    }
    if (n)
        sd_col = sqrt(sd_col / n);

    ub_col = mean_col + 3.0 * sd_col;    


    if (min_global < DBL_MAX)
        min_global = log(min_global);
    
    
    /* impute values */
    for (row = 0; row < num_rows; row++)
    {
        row_ptr = data_matrix[row];
        posix_madvise(row_ptr, num_cols * sizeof(double),
                      POSIX_MADV_WILLNEED);

        min_row = row_min_array[row];
        if (min_row != MISSING)
            min_row = log(min_row);

        for (col = 0; col < num_cols; col++)
        {
            value = row_ptr[col];

            if (value == 0.0 || value == MISSING)
            {
                min_col = col_min_array[col];
                if (min_col != MISSING)
                    min_col = log(min_col);
            
                if (min_row != MISSING && min_col != MISSING)
                {
                    /* average in both col and row min if within range */
                    if (min_row < min_col || min_row < ub_col)
                    {
                        value = 0.5 * exp(0.5 * (min_global +
                                                 0.5 * (min_row + min_col)));
                    }
                    /* only average in col min */
                    else
                    {
                        value = 0.5 * exp(0.5 * (min_global + min_col));
                    }
                }
                /* only average in col min */
                else if (min_col != MISSING)
                {
                    value = 0.5 * exp(0.5 * (min_global + min_col));
                }
                /* only average in row min */
                else if (min_row != MISSING && min_row < ub_col)
                {
                    value = 0.5 * exp(0.5 * (min_global + min_row));
                }
                /* use only global min */
                else
                {
                    value = 0.5 * exp(min_global);
                }
                
                row_ptr[col] = value;
            }
        }
    }
    
    if (col_min_array)
        free(col_min_array);

    if (row_min_array)
        free(row_min_array);
}


/* estimate limit of detection from col min */
double estimate_lod(double **data_matrix, int32_t num_rows,
                    int32_t num_cols, int32_t floor_lod_flag)
{
    double *col_min_array = NULL;    /* minimum value for each column */
    double *row_ptr;
    double  value;
    double  min_global;
    double  mean_col, sd_col, lb_col, ub_col, diff;   /* all in log units */
    int32_t row, col, n;
    int32_t missing_col_count = 0;
    int32_t missing_col_flag  = 0;

    col_min_array = (double *) malloc(num_cols * sizeof(double));

    for (col = 0; col < num_cols; col++)
        col_min_array[col] = DBL_MAX;
    
    min_global = DBL_MAX;
    for (row = 0; row < num_rows; row++)
    {
        row_ptr = data_matrix[row];
        posix_madvise(row_ptr, num_cols * sizeof(double),
                      POSIX_MADV_WILLNEED);
    
        for (col = 0; col < num_cols; col++)
        {
            value = row_ptr[col];

            if (value > 0.0 && value != MISSING)
            {
                if (value < col_min_array[col])
                {
                    col_min_array[col] = value;

                    if (value < min_global)
                        min_global = value;
                }
            }
        }
    }
    
    
    mean_col = 0.0;
    n        = 0.0;
    for (col = 0; col < num_cols; col++)
    {
        value = col_min_array[col];
        
        if (value < DBL_MAX)
        {
            mean_col += log(value);
            n++;
        }
        else
        {
            col_min_array[col] = MISSING;

            missing_col_flag   = 1;
            missing_col_count++;
        }
    }
    if (n)
        mean_col /= n;
    
    sd_col   = 0.0;
    for (col = 0; col < num_cols; col++)
    {
        value = col_min_array[col];
        
        if (value != MISSING)
        {
            diff    = log(value) - mean_col;
            sd_col += diff * diff;
        }
    }
    if (n)
        sd_col = sqrt(sd_col / n);


    /* upper bound of col mins */
    ub_col = mean_col + 3.0 * sd_col;
    ub_col = exp(ub_col);

    /* lower bound of col mins */
    lb_col = mean_col - 3.0 * sd_col;
    lb_col = exp(lb_col);


    fprintf(stderr, "Noise floor global/lb/mean/ub:  %f  %f  %f  %f\n",
            min_global, lb_col, exp(mean_col), ub_col);
#if 0
    fprintf(stderr, "Noise floor log2/unlog stdev:  %f  %f\n",
            sd_col, exp(sd_col));
#endif


    /* floor the lod lower bound at the global min */
    if (lb_col < min_global)
        lb_col = min_global;


    /* warn on terminal invariant cols */
    if (missing_col_count && missing_col_flag)
    {
        fprintf(stderr, "WARNING -- %d (%f%%) col(s) entirely zero for all rows\n",
            missing_col_count, 
            100.0 * (double) missing_col_count / (double) num_cols);
    }
        
    
    if (col_min_array)
        free(col_min_array);

    /* use the upper bound instead */
    if (floor_lod_flag == 2)
        return ub_col;

    return lb_col;
}


/* assume unlogged data, assume zeroes are missing */
void floor_to_one(double **data_matrix, int32_t num_rows,
                  int32_t num_cols)
{
    double *row_ptr;
    double  value;
    int32_t row, col;

    /* impute values */
    for (row = 0; row < num_rows; row++)
    {
        row_ptr = data_matrix[row];
        posix_madvise(row_ptr, num_cols * sizeof(double),
                      POSIX_MADV_WILLNEED);

        for (col = 0; col < num_cols; col++)
        {
            value = row_ptr[col];

            if (value < 1.0 + 1E-14 || value == MISSING)
            {
                row_ptr[col] = 1.0;
            }
        }
    }
}


/* floor to value, also set missing to value */
void floor_to_value(double **data_matrix, int32_t num_rows,
                    int32_t num_cols, double lod)
{
    double *row_ptr;
    double  value;
    int32_t row, col;
    

    for (row = 0; row < num_rows; row++)
    {
        row_ptr = data_matrix[row];
        posix_madvise(row_ptr, num_cols * sizeof(double),
                      POSIX_MADV_WILLNEED);

        for (col = 0; col < num_cols; col++)
        {
            value = row_ptr[col];

            if (value < lod || value == MISSING)
            {
                row_ptr[col] = lod;
            }
        }
    }
}


void log2_data(double **data_matrix, int32_t num_rows, int32_t num_cols)
{
    double *row_ptr;
    double  value;
    int32_t row, col;
    
    for (row = 0; row < num_rows; row++)
    {
        row_ptr = data_matrix[row];
        posix_madvise(row_ptr, num_cols * sizeof(double),
                      POSIX_MADV_WILLNEED);

        for (col = 0; col < num_cols; col++)
        {
            value = row_ptr[col];

            /* leave MISSING and bad spots as-is */
            if (value == MISSING || value >= DBL_MAX)
                continue;
            
            /* uh oh, negative value or zero, set to missing */
            if (value <= 0.0)
            {
                row_ptr[col] = MISSING;
            }
            /* take the log2 */
            else
            {
                row_ptr[col] = log2(value);
            }
        }
    }
}


void unlog2_data(double **data_matrix, int32_t num_rows, int32_t num_cols)
{
    double *row_ptr;
    double  value;
    int32_t row, col;
    
    for (row = 0; row < num_rows; row++)
    {
        row_ptr = data_matrix[row];
        posix_madvise(row_ptr, num_cols * sizeof(double),
                      POSIX_MADV_WILLNEED);

        for (col = 0; col < num_cols; col++)
        {
            value = row_ptr[col];

            /* leave MISSING and bad spots as-is */
            if (value == MISSING || value >= DBL_MAX)
                continue;
            
            /* uh oh, negative value or zero, set to missing */
            if (value <= 0.0)
            {
                row_ptr[col] = MISSING;
            }
            /* unlog the value */
            else
            {
                row_ptr[col] = pow(2.0, value);
            }
        }
    }
}


double ** filter_rows_by_stats(double  **data_matrix,
                               char   ***return_row_name_array,
                               int32_t  *return_num_rows, int32_t num_cols,
                               struct options *opt)
{
    char    **row_name_array     = *return_row_name_array;
    char    **row_name_array_new = NULL;
    double  **data_matrix_new    = NULL;
    int32_t  *filter_flags_array = NULL;

    double   *row_ptr;
    double    value, diff;
    double    log2_mean, log2_sd;
    double    unlog_mean, unlog_sd, unlog_max;
    double    unlog_mean_plus_2sd, log2_mean_plus_2sd;
    int32_t   num_rows = *return_num_rows;
    int32_t   num_rows_new;
    int32_t   row, col;
    int32_t   n, n_unlog, n_log2;
    int32_t   keep_flag;

    double    unlog_max_signal_cutoff = opt->filter_unlog_max_signal_cutoff;
    double    unlog_mean_cutoff       = opt->filter_unlog_mean_cutoff;
    double    unlog_sd_cutoff         = opt->filter_unlog_sd_cutoff;
    double    log2_mean_cutoff        = opt->filter_log2_mean_cutoff;
    double    log2_sd_cutoff          = opt->filter_log2_sd_cutoff;

    double    unlog_mean_plus_2sd_cutoff =
                  unlog_mean_cutoff + 2.0 * unlog_sd_cutoff;
    double    log2_mean_plus_2sd_cutoff =
                  log2_mean_cutoff + 2.0 * log2_sd_cutoff;

    filter_flags_array = (int32_t *) calloc(num_rows, sizeof(int32_t));
    
    /* calculate average and sd for each row */
    num_rows_new = 0;
    for (row = 0; row < num_rows; row++)
    {
        row_ptr = data_matrix[row];
        posix_madvise(row_ptr, num_cols * sizeof(double),
                      POSIX_MADV_WILLNEED);
        
        unlog_max  = -DBL_MAX;
        unlog_mean = log2_mean = 0.0;
        unlog_sd   = log2_sd   = 0.0;
        n_unlog    = n_log2    = 0;

        for (col = 0; col < num_cols; col++)
        {
            value = row_ptr[col];

            if (value == MISSING)
                continue;

            if (value > unlog_max)
                unlog_max = value;

            unlog_mean += value;
            n_unlog++;
            
            if (value >= 0)
            {
                log2_mean += log2(value);
                n_log2++;
            }
        }
        
        if (n_unlog)
            unlog_mean /= n_unlog;
        if (n_log2)
            log2_mean  /= n_log2;

        for (col = 0; col < num_cols; col++)
        {
            value = row_ptr[col];

            if (value == MISSING)
                continue;

            diff      = value - unlog_mean;
            unlog_sd += diff * diff;

            if (value >= 0)
            {
                diff     = log2(value) - log2_mean;
                log2_sd += diff * diff;
            }
        }
        
        if (n_unlog)
            unlog_sd = sqrt(unlog_sd / n_unlog);

        if (n_log2)
            log2_sd = sqrt(log2_sd / n_log2);


        unlog_mean_plus_2sd = unlog_mean + 2.0 * unlog_sd;
        log2_mean_plus_2sd  = log2_mean  + 2.0 * log2_sd;


        keep_flag = 1;

        if (unlog_max_signal_cutoff &&
            unlog_max < unlog_max_signal_cutoff - 1E-14)
        {
            keep_flag = 0;
        }
        if (unlog_mean_cutoff &&
            unlog_mean < unlog_mean_cutoff - 1E-14)
        {
            keep_flag = 0;
        }
        if (unlog_sd_cutoff &&
            unlog_sd < unlog_sd_cutoff - 1E-14)
        {
            keep_flag = 0;
        }
        if (log2_mean_cutoff &&
            log2_mean < log2_mean_cutoff - 1E-14)
        {
            keep_flag = 0;
        }
        if (log2_sd_cutoff &&
            log2_sd < log2_sd_cutoff - 1E-14)
        {
            keep_flag = 0;
        }

        if (unlog_mean_cutoff && unlog_sd_cutoff &&
            unlog_mean_plus_2sd < unlog_mean_plus_2sd_cutoff - 1E-14)
        {
            keep_flag = 0;
        }
        if (log2_mean_cutoff && log2_sd_cutoff &&
            log2_mean_plus_2sd < log2_mean_plus_2sd_cutoff - 1E-14)
        {
            keep_flag = 0;
        }


        if (keep_flag)
        {
            filter_flags_array[row] = 1;
            num_rows_new++;
        }
    }
    

    /* nothing was filtered, return without doing anything */
    if (num_rows == num_rows_new)
        return data_matrix;
    
    
    /* allocate new filtered arrays */
    row_name_array_new = (char **)   calloc(num_rows_new, sizeof(char   *));
    data_matrix_new    = (double **) calloc(num_rows_new, sizeof(double *));


    /* filter into new arrays */
    n = 0;
    for (row = 0; row < num_rows; row++)
    {
        if (filter_flags_array[row])
        {
            row_name_array_new[n] = row_name_array[row];
            data_matrix_new[n]    = data_matrix[row];
            n++;
        }
        else if (data_matrix[row])
        {
            free(data_matrix[row]);
        }
    }
    

    /* free old arrays */
    if (row_name_array)
        free(row_name_array);

    if (data_matrix)
        free(data_matrix);

    /* clean up */
    if (filter_flags_array)
        free(filter_flags_array);


    /* set new number of rows, etc. */
    *return_num_rows       = num_rows_new;
    *return_row_name_array = row_name_array_new;

    return data_matrix_new;
}


double ** filter_rows_by_present(double  **data_matrix,
                                 char   ***return_row_name_array,
                                 int32_t  *return_num_rows, int32_t num_cols,
                                 struct options *opt)
{
    char    **row_name_array     = *return_row_name_array;
    char    **row_name_array_new = NULL;
    double  **data_matrix_new    = NULL;
    int32_t  *filter_flags_array = NULL;

    double   *row_ptr;
    double    mag_cutoff = opt->filter_present_mag;
    double    global_min = DBL_MAX;
    double    value;
    int32_t   num_rows = *return_num_rows;
    int32_t   num_rows_new;
    int32_t   row, col;
    int32_t   n, n_present;
    int32_t   keep_flag;
    int32_t   count_passed = 0;

    double    fraction_present_cutoff = opt->filter_present_cutoff;
    double    fraction_present;


    /* nothing was filtered, return without doing anything */
    if (fraction_present_cutoff < 1E-14)
        return data_matrix;


    /* fudge magnitude cutoff for floating point roundoff error;
     * also floor negative cutoffs to zero
     */
    mag_cutoff += 1E-14;
    if (mag_cutoff < 0.0)
        mag_cutoff = 0.0;


    filter_flags_array = (int32_t *) calloc(num_rows, sizeof(int32_t));
    
    num_rows_new = 0;
    for (row = 0; row < num_rows; row++)
    {
        row_ptr = data_matrix[row];
        posix_madvise(row_ptr, num_cols * sizeof(double),
                      POSIX_MADV_WILLNEED);
        
        n_present = 0;
        for (col = 0; col < num_cols; col++)
        {
            value = row_ptr[col];
            
            if (value != MISSING)
            {
                if (value < global_min)
                    global_min = value;

                if (value > mag_cutoff)
                    n_present++;
            }
        }

        fraction_present = (double) n_present / (double) num_cols;
        keep_flag        = 1;

        if (fraction_present < fraction_present_cutoff - 1E-14)
        {
            keep_flag = 0;
        }

        if (keep_flag)
        {
            filter_flags_array[row] = 1;
            num_rows_new++;
        }
        
        count_passed += n_present;
    }
    
    
    /* HACK -- data may have been imputed, try global min for cutoff */
    if (count_passed == num_rows * num_cols &&
        opt->filter_present_mag == 0.0 &&
        global_min != 0.0 && global_min != DBL_MAX)
    {
        memset(filter_flags_array, 0, num_rows * sizeof(int32_t));
    
        mag_cutoff   = global_min;
        count_passed = 0;
        num_rows_new = 0;

        for (row = 0; row < num_rows; row++)
        {
            row_ptr = data_matrix[row];
            posix_madvise(row_ptr, num_cols * sizeof(double),
                          POSIX_MADV_WILLNEED);
            
            n_present = 0;
            for (col = 0; col < num_cols; col++)
            {
                value = row_ptr[col];
                
                if (value > mag_cutoff && value != MISSING)
                    n_present++;
            }

            fraction_present = (double) n_present / (double) num_cols;
            keep_flag        = 1;

            if (fraction_present < fraction_present_cutoff - 1E-14)
            {
                keep_flag = 0;
            }

            if (keep_flag)
            {
                filter_flags_array[row] = 1;
                num_rows_new++;
            }
            
            count_passed += n_present;
        }
    }

    /* nothing was filtered, return without doing anything */
    if (num_rows == num_rows_new)
        return data_matrix;
    
    
    /* allocate new filtered arrays */
    row_name_array_new = (char **)   calloc(num_rows_new, sizeof(char   *));
    data_matrix_new    = (double **) calloc(num_rows_new, sizeof(double *));


    /* filter into new arrays */
    n = 0;
    for (row = 0; row < num_rows; row++)
    {
        if (filter_flags_array[row])
        {
            row_name_array_new[n] = row_name_array[row];
            data_matrix_new[n]    = data_matrix[row];

            n++;
        }
        else if (data_matrix[row])
        {
            free(data_matrix[row]);
        }
    }
    

    /* free old arrays */
    if (row_name_array)
        free(row_name_array);

    if (data_matrix)
        free(data_matrix);

    /* clean up */
    if (filter_flags_array)
        free(filter_flags_array);


    /* set new number of rows, etc. */
    *return_num_rows       = num_rows_new;
    *return_row_name_array = row_name_array_new;

    return data_matrix_new;
}


void scale_ten_thousand(double **data_matrix, int32_t num_rows, int32_t num_cols)
{
    double *row_ptr;
    double value;
    double max;
    int32_t row, col;
    
    /* find max */
    for (row = 0; row < num_rows; row++)
    {
        row_ptr = data_matrix[row];
        posix_madvise(row_ptr, num_cols * sizeof(double),
                      POSIX_MADV_WILLNEED);
        
        max = -DBL_MAX;

        for (col = 0; col < num_cols; col++)
        {
            value = row_ptr[col];

            if (value > max)
                max = value;
        }
        
        /* scale to 10000 */
        if (max > -DBL_MAX)
        {
            for (col = 0; col < num_cols; col++)
            {
                /* skip bad spots */
                if (row_ptr[col] < DBL_MAX)
                    row_ptr[col] = 10000.0 * row_ptr[col] / max;
            }
        }
    }
}


void *fill_dist_matrix_thread(void *passed_ptr)
{
    struct dist_thread_data *tdata_ptr =
        (struct dist_thread_data *) passed_ptr;

    int *t_flags_array = tdata_ptr->t_flags_array;
    int  t_index       = tdata_ptr->t_index;
    int  t_flag        = 0;
    pthread_mutex_t *mutex_flags_ptr = &tdata_ptr->mutex_flags;
    pthread_mutex_t *mutex_vars_ptr  = &tdata_ptr->mutex_vars;
    pthread_mutex_t *mutex_wait_ptr  = &tdata_ptr->mutex_wait;
    pthread_cond_t  *cond_wait_ptr   = &tdata_ptr->cond_wait;
    struct timespec  ns_to_wait;

    double **data_matrix     = tdata_ptr->data_matrix;
    double **dist_matrix     = tdata_ptr->dist_matrix;
    double  *weights         = tdata_ptr->weights;
    double  *flip_row1       = tdata_ptr->flip_row1;
    int32_t *row_good_counts = tdata_ptr->row_good_counts;
    int32_t *x_coord_array   = tdata_ptr->x_coord_array;
    int32_t *y_coord_array   = tdata_ptr->y_coord_array;
    int32_t  num_cols        = tdata_ptr->num_cols;
    struct   options *opt    = tdata_ptr->opt;
    int32_t  row2_inc        = tdata_ptr->row2_inc;

    double  *row_ptr1;
    int32_t  n1;
    int32_t  row1;
    int32_t  row2_start;

    double  *row_ptr2;
    double   max_dist = DBL_MIN;
    uint64_t n_scanned;

    double   value1, value2;
    double   euc2;
    double   dist, diff, pearson;
    double   similarity;
    double   w_count;
    int32_t  row2, col;
    int32_t  n, n2, n_max, n_min;
    int32_t  a, b, c, d, p;

    char rmsd_vs_euclidian_flag = opt->rmsd_vs_euclidian_flag;
    char pearson_flag           = opt->pearson_flag;
    char cosine_flag            = opt->cosine_flag;
    char heuristic_flag         = opt->heuristic_flag;
    char geomean_flag           = opt->geomean_flag;
    char similarity_flag        = opt->similarity_flag;
    char absolute_flag          = opt->absolute_flag;
    char mostly_trig_flag       = opt->mostly_trig_flag;

    double minkowski            = opt->minkowski;


    /* https://stackoverflow.com/questions/72145863/
     * 84 ns on Windows, 18 ns on Linux
     * 25 should be fine for Linux, but we'll use 100 just to be safe
     */
    ns_to_wait.tv_sec  = 0;
    ns_to_wait.tv_nsec = 100;


    /* flags:
     *   0  nothing to do
     *   1  just finished
     *   2  something to do
     *   3  signal to terminate
     */

    while(1)
    {
      /* call is to a real thread, not an unthreaded call */
      if (t_index)
      {
#if 0
          /* paranoia, in case it didn't get signalled properly earlier */
          pthread_mutex_lock(&mutex_wait_main);
          if (n_threads_left <= 0)
              pthread_cond_signal(&cond_wait_main);
          pthread_mutex_unlock(&mutex_wait_main);
#endif


          /* sleep until there is something to do */
          pthread_mutex_lock(mutex_wait_ptr);
          pthread_mutex_lock(mutex_flags_ptr);
          t_flag = t_flags_array[t_index];
          pthread_mutex_unlock(mutex_flags_ptr);
          if (t_flag < 2)
              pthread_cond_wait(cond_wait_ptr, mutex_wait_ptr);
          pthread_mutex_unlock(mutex_wait_ptr);

          pthread_mutex_lock(mutex_flags_ptr);
          t_flag = t_flags_array[t_index];
          pthread_mutex_unlock(mutex_flags_ptr);


          /* fallback in case of spurious wakeup */
          while (t_flag < 2)
          {
              /* Small numbers of threads can deadlock, so wait a bit.
               * We are highly unlikely to ever get here anyways, so we
               * aren't worried about performance at this point.
               */
              nanosleep(&ns_to_wait, NULL);

              pthread_mutex_lock(mutex_flags_ptr);
              t_flag = t_flags_array[t_index];
              pthread_mutex_unlock(mutex_flags_ptr);
          }


          /* terminate the thread */
          if (t_flag == 3)
          {
              pthread_mutex_destroy(mutex_flags_ptr);
              pthread_mutex_destroy(mutex_vars_ptr);
              pthread_mutex_destroy(mutex_wait_ptr);
              pthread_cond_destroy(cond_wait_ptr);
              pthread_exit(NULL);
          }
      }


      /* do new work */
    
      pthread_mutex_lock(mutex_vars_ptr);
      row_ptr1   = tdata_ptr->row_ptr1;
      n1         = tdata_ptr->n1;
      row1       = tdata_ptr->row1;
      row2_start = tdata_ptr->row2_start;
      pthread_mutex_unlock(mutex_vars_ptr);

      n_scanned  = 0;
      for (row2 = row2_start; row2 < row1; row2 += row2_inc)
      {
          row_ptr2 = data_matrix[row2];
          posix_madvise(row_ptr2, num_cols * sizeof(double),
                        POSIX_MADV_WILLNEED);
          posix_madvise(row_ptr1, num_cols * sizeof(double),
                        POSIX_MADV_WILLNEED);
          
          n2 = row_good_counts[row2];
          
          dist    = 0;
          w_count = 0.0;

          /* no missing data */
          n = 0;
          if (n1 == num_cols && n2 == num_cols)
          {
              n = n1;
              
              /* not necessarily Euclidian, could be any l-norm */
              if (opt->calc_euclidian_flag)
              {
                  if (opt->cityblock_flag)
                  {
                      for (col = 0; col < num_cols; col++)
                      {
                          diff = fabs(row_ptr1[col] - row_ptr2[col]);

                          if (!weights)
                          {
                              dist    += diff;
                          }
                          else
                          {
                              dist    += weights[col] * diff;
                              w_count += weights[col];
                          }
                      }
                  }
                  /* l-#-norm */
                  else if (minkowski)
                  {
                      if (minkowski == 1.5)
                      {
                          for (col = 0; col < num_cols; col++)
                          {
                              diff = fabs(row_ptr1[col] - row_ptr2[col]);

                              if (!weights)
                              {
                                  dist    += diff * sqrt(diff);
                              }
                              else
                              {
                                  dist    += weights[col] * diff * sqrt(diff);
                                  w_count += weights[col];
                              }
                          }
                      }
                      else if (minkowski == 1.75)
                      {
                          double temp;
                      
                          for (col = 0; col < num_cols; col++)
                          {
                              diff = fabs(row_ptr1[col] - row_ptr2[col]);
                              temp = sqrt(diff);

                              if (!weights)
                              {
                                  dist    += diff * temp * sqrt(temp);
                              }
                              else
                              {
                                  dist    += weights[col] * diff *
                                             temp * sqrt(temp);
                                  w_count += weights[col];
                              }
                          }
                      }
                      else if (minkowski == 1.25)
                      {
                          for (col = 0; col < num_cols; col++)
                          {
                              diff = fabs(row_ptr1[col] - row_ptr2[col]);

                              if (!weights)
                              {
                                  dist    += diff * sqrt(sqrt(diff));
                              }
                              else
                              {
                                  dist    += weights[col] * diff *
                                             sqrt(sqrt(diff));
                                  w_count += weights[col];
                              }
                          }
                      }
                      else if (minkowski == 2.0)
                      {
                          for (col = 0; col < num_cols; col++)
                          {
                              diff = row_ptr1[col] - row_ptr2[col];

                              if (!weights)
                              {
                                  dist    += diff * diff;
                              }
                              else
                              {
                                  dist    += weights[col] * diff * diff;
                                  w_count += weights[col];
                              }
                          }
                      }
                      else if (minkowski == 1.0)
                      {
                          for (col = 0; col < num_cols; col++)
                          {
                              diff = fabs(row_ptr1[col] - row_ptr2[col]);

                              if (!weights)
                              {
                                  dist    += diff;
                              }
                              else
                              {
                                  dist    += weights[col] * diff;
                                  w_count += weights[col];
                              }
                          }
                      }
                      else
                      {
                          for (col = 0; col < num_cols; col++)
                          {
                              diff = fabs(row_ptr1[col] - row_ptr2[col]);

                              if (!weights)
                              {
                                  dist    += pow(diff, minkowski);
                              }
                              else
                              {
                                  dist    += weights[col] * pow(diff,minkowski);
                                  w_count += weights[col];
                              }
                          }
                      }
                  }
                  /* Euclidian */
                  else
                  {
                      for (col = 0; col < num_cols; col++)
                      {
                          diff = row_ptr1[col] - row_ptr2[col];

                          if (!weights)
                          {
                              dist    += diff * diff;
                          }
                          else
                          {
                              dist    += weights[col] * diff * diff;
                              w_count += weights[col];
                          }
                      }
                  }
              }
          }
          /* missing data, but both have at least *some* data */
          else if (n1 && n2)
          {
              for (col = 0; col < num_cols; col++)
              {
                  value1 = row_ptr1[col];
                  value2 = row_ptr2[col];
                  
                  if (value1 == MISSING || value2 == MISSING)
                  {
                      continue;
                  }

                  n++;
                  
                  if (opt->calc_euclidian_flag)
                  {
                      if (opt->cityblock_flag)
                      {
                          diff = fabs(value1 - value2);

                          if (!weights)
                          {
                              dist    += diff;
                          }
                          else
                          {
                              dist    += weights[col] * diff;
                              w_count += weights[col];
                          }
                      }
                      else if (minkowski)
                      {
                          diff = fabs(value1 - value2);

                          if (minkowski == 1.5)
                          {
                              if (!weights)
                              {
                                  dist    += diff * sqrt(diff);
                              }
                              else
                              {
                                  dist    += weights[col] * diff * sqrt(diff);
                                  w_count += weights[col];
                              }
                          }
                          else if (minkowski == 1.75)
                          {
                              double temp;
                              
                              temp = sqrt(diff);
                          
                              if (!weights)
                              {
                                  dist    += diff * temp * sqrt(temp);
                              }
                              else
                              {
                                  dist    += weights[col] * diff *
                                             temp * sqrt(temp);
                                  w_count += weights[col];
                              }
                          }
                          else if (minkowski == 1.25)
                          {
                              if (!weights)
                              {
                                  dist    += diff * sqrt(sqrt(diff));
                              }
                              else
                              {
                                  dist    += weights[col] * diff *
                                             sqrt(sqrt(diff));
                                  w_count += weights[col];
                              }
                          }
                          else if (minkowski == 2.0)
                          {
                              if (!weights)
                              {
                                  dist    += diff * diff;
                              }
                              else
                              {
                                  dist    += weights[col] * diff * diff;
                                  w_count += weights[col];
                              }
                          }
                          else if (minkowski == 1.0)
                          {
                              if (!weights)
                              {
                                  dist    += diff;
                              }
                              else
                              {
                                  dist    += weights[col] * diff;
                                  w_count += weights[col];
                              }
                          }
                          else
                          {
                              if (!weights)
                              {
                                  dist    += pow(diff, minkowski);
                              }
                              else
                              {
                                  dist    += weights[col] * pow(diff,minkowski);
                                  w_count += weights[col];
                              }
                          }
                      }
                      /* Euclidian */
                      else
                      {
                          diff = value1 - value2;

                          if (!weights)
                          {
                              dist    += diff * diff;
                          }
                          else
                          {
                              dist    += weights[col] * diff * diff;
                              w_count += weights[col];
                          }
                      }
                  }
              }
          }

          /* no points in common */
          if (n == 0)
          {
              dist_matrix[row1][row2] = DBL_MAX;
              n_scanned++;

              continue;
          }
          
          /* store sum of squared distances for later optimization */
          euc2 = dist;
          
          if (!weights)
              w_count = n;

          /* RMAD */
          if (opt->cityblock_flag)
          {
              dist = dist / w_count;
          }
          else if (minkowski)
          {
              /* x^(2/3) = (x^1/3) * (x^1/3) */
              if (minkowski == 1.5)
              {
                  /* NOTE -- cbrt() is slower than pow(1.0 / 3.0),
                   *         but should, in theory, be more accurate
                   *
                   * since this isn't a rate limiting step, use cbrt()
                   */
                  dist  = cbrt(dist / w_count);
                  dist *= dist;
              }
              /* RMSD */
              else if (minkowski == 2.0)
              {
                  dist = sqrt(dist / w_count);
              }
              /* RMAD */
              else if (minkowski == 1.0)
              {
                  dist = dist / w_count;
              }
              else
              {
                  dist = pow(dist / w_count, 1.0 / minkowski);
              }
          }
          /* RMSD */
          else if (rmsd_vs_euclidian_flag)
          {
              dist = sqrt(dist / w_count);
          }
          
          /* Scale up the sum like R does, rather than RMSD.
           * This doesn't really matter to us, since it just globally
           *  scales RMSD by sqrt(num_cols).
           */
          else
          {
              /* City Block */
              if (opt->cityblock_flag)
              {
                  dist = dist * (double) num_cols / w_count;
              }
              else if (minkowski)
              {
                  /* x^(2/3) = (x^1/3) * (x^1/3) */
                  if (minkowski == 1.5)
                  {
                      dist  = cbrt(dist * (double) num_cols / w_count);
                      dist *= dist;
                  }
                  /* Euclidian */
                  else if (minkowski == 2.0)
                  {
                      dist = sqrt(dist * (double) num_cols / w_count);
                  }
                  /* City Block */
                  else if (minkowski == 1.0)
                  {
                      dist = dist * (double) num_cols / w_count;
                  }
                  else
                  {
                      dist = pow(dist * (double) num_cols / w_count,
                             1.0 / minkowski);
                  }
              }
              /* Euclidian */
              else
              {
                  dist = sqrt(dist * (double) num_cols / w_count);
              }
          }


          /* no Euclidian/RMSD calculation */
          if (opt->calc_euclidian_flag == 0 &&
              opt->cityblock_flag      == 0 &&
              opt->minkowski           == 0.0)
          {
              dist_matrix[row1][row2] = DBL_MAX;
          }

          
          dist_matrix[row1][row2] = dist;

          /* shut -Wall up, gcc isn't smart enough
           *
           * If I truly have messed up, initializing to DBL_MAX should yield
           * garbage results and it should be immediately noticable.
           */
          pearson = DBL_MAX;

          if (pearson_flag || cosine_flag)
          {
              /* unweighted */
              if (!weights)
              {
                  /* optimize if no missing data */
                  if (n == num_cols && opt->mean_center_flag &&
                      opt->transpose_last_flag == 0)
                  {
                      /* Transform Euclidian distance if mean centered uv.
                       *
                       * When vectors are L2-normalized, |vector| = 1.
                       *
                       * If the data is scaled to unit variance, we know
                       * that the resulting vector length is whichever
                       * n was used in the unit vector standard deviation
                       * calculations (in this case, n, since we're using
                       * population standard deviation in this program).
                       *
                       * Rather than scale the input data to unit length,
                       * we can just down-scale the sum of squared
                       * differences by the n we used in the unit variance
                       * calculations.  Or, in other words, assuming we
                       * used population standard deviations, d = RMSD.
                       *
                       *     d = sqrt(sum_squared_diffs / n)
                       *     d^2 = 2(1 - r)
                       *     r = 1 - d^2 / 2
                       *
                       * If sample standard deviations (n-1) are used
                       * instead of population standard deviations (n),
                       * then d needs to be calculated with n-1 instead
                       * of n. In the code below, that would be (n+n-2).
                       */
                      if (opt->unit_variance_flag)
                      {
                          pearson = 1.0 - euc2 / (n+n);
                      }
                      /* centered cosine is same as pearson */
                      else
                      {
                          pearson = calc_cosine_no_missing(row_ptr1,
                                                           row_ptr2,
                                                           num_cols);
                      }
                  }
                  else if (cosine_flag)
                  {
                      if (n == num_cols)
                      {
                          pearson = calc_cosine_no_missing(row_ptr1,
                                                           row_ptr2,
                                                           num_cols);
                      }
                      else
                      {
                          pearson = calc_cosine(row_ptr1, row_ptr2, num_cols);
                      }
                  }
                  /* Pearson correlation */
                  else
                  {
                      if (n == num_cols)
                      {
                          pearson = calc_pearson_r_no_missing(row_ptr1,
                                                              row_ptr2,
                                                              num_cols);
                      }
                      else
                      {
                          pearson = calc_pearson_r(row_ptr1, row_ptr2,
                                                   num_cols);
                      }
                  }
              }
              /* weighted */
              else
              {
                  /* optimize if no missing data */
                  if (n == num_cols && opt->mean_center_flag &&
                      opt->transpose_last_flag == 0)
                  {
                      /* centered cosine is same as pearson */
                      pearson = calc_cosine_weighted_no_missing(row_ptr1,
                                    row_ptr2, weights, num_cols);
                  }
                  else if (cosine_flag)
                  {
                      if (n == num_cols)
                      {
                          pearson =
                              calc_cosine_weighted_no_missing(row_ptr1,
                                  row_ptr2, weights, num_cols);
                      }
                      else
                      {
                          pearson =
                              calc_cosine_weighted(row_ptr1,
                                  row_ptr2, weights, num_cols);
                      }
                  }
                  /* Pearson correlation */
                  else
                  {
                      if (n == num_cols)
                      {
                          pearson =
                              calc_pearson_r_weighted_no_missing(row_ptr1,
                                  row_ptr2, weights, num_cols);
                      }
                      else
                      {
                          pearson =
                              calc_pearson_r_weighted(row_ptr1,
                                  row_ptr2, weights, num_cols);
                      }
                  }
              }
          }
          

          if ((pearson_flag || cosine_flag) &&
              heuristic_flag == 0 && geomean_flag == 0)
          {
              /* transform r into a metric distance [sqrt(0.5 * (1-r))]
               */
              /* see:
               *   https://doi.org/10.1101/582106
               *   Chen et al. 2019
               *   "On triangular Inequalities of correlation-based
               *    distances for gene expression profiles"
               *
               *  and
               *
               *   Stijn van Dongen and Anton J. Enright 2012
               *   "Metric distance derived from cosine similarity and
               *    Pearson and Spearman correlations"
               *
               * sqrt(1 - |r|) is metric, and works better for biological
               * clustering than sqrt(1 - r^2).  My own experience agrees
               * that it works better as well.
               */
              
              /* random vectors length 429 cause r --> 0.0
               * due to double precision limits
               */

              if (absolute_flag)
              {
                  pearson = 1.0 - fabs(pearson);

                  /* deal with potential NaN-inducing roundoff error */
                  if (pearson < 0.0)
                      pearson = 0.0;

                  pearson = sqrt(pearson);
              }
              else
              {
                  pearson = 1.0 - pearson;

                  /* deal with potential NaN-inducing roundoff error */
                  if (pearson < 0.0)
                      pearson = 0.0;

                  pearson = sqrt(0.5 * pearson);
              }

              if (mostly_trig_flag == 0)
                  dist_matrix[row1][row2] = pearson;
              else
                  dist_matrix[row1][row2] = 0.001*dist_matrix[row1][row2] +
                                            0.999*pearson;
          }
          
          /* flip one of the vectors, average the two distances if better */
          /* this should never happen, given current command line options */
          else if (heuristic_flag && pearson_flag == 0 && cosine_flag == 0)
          {
              dist = 0;
              n = 0;
              for (col = 0; col < num_cols; col++)
              {
                  value1 = flip_row1[col];
                  value2 = row_ptr2[col];
                  
                  if (value1 != MISSING && value2 != MISSING)
                  {
                      diff = value1 - value2;
                      dist += diff * diff;
                      n++;
                  }
              }

              if (n)
              {
                  /* RMSD */
                  if (rmsd_vs_euclidian_flag)
                  {
                      dist = sqrt(dist / n);
                  }
                  /* Scale up the sum like R does, rather than RMSD.
                   * This doesn't really matter, since it just globally
                   *  scales RMSD by sqrt(num_cols).
                   */
                  else
                  {
                      dist = sqrt(dist * (double) num_cols / (double) n);
                  }
              }
              
              /* average the two distances if flipped is better */
              if (dist < dist_matrix[row1][row2])
              {
                  dist_matrix[row1][row2] =
                      0.5 * (dist_matrix[row1][row2] + dist);
              }
          }
          /* combine Euclidian with trigonometric distance */
          else if ((heuristic_flag || geomean_flag) &&
                   (pearson_flag || cosine_flag))
          {
              /* always use absolute value for --heuristic */
              /* respect --signed --absolute flags with --geomean */
              if (heuristic_flag || absolute_flag)
              {
                  pearson = 1.0 - fabs(pearson);

                  /* deal with potential NaN-inducing roundoff error */
                  if (pearson < 0.0)
                      pearson = 0.0;

                  pearson = sqrt(pearson);
              }
              else
              {
                  pearson = 1.0 - pearson;

                  /* deal with potential NaN-inducing roundoff error */
                  if (pearson < 0.0)
                      pearson = 0.0;

                  pearson = sqrt(0.5 * pearson);
              }
              
#if 1
              dist_matrix[row1][row2] =
                  sqrt(dist_matrix[row1][row2] * (1.0 + pearson));
#endif

#if 0
              dist_matrix[row1][row2] =
                  sqrt(dist_matrix[row1][row2] * pearson);
#endif

#if 0
              dist_matrix[row1][row2] =
                  dist_matrix[row1][row2] * (1.0 + pearson);
#endif
          }

          
          /* Todeschini, Consonni, et al. J Chem Inf Model 2012 */
          similarity = 1.0;
          if (similarity_flag)
          {
              p = num_cols;
              a = n;

              /* abc = row_good_counts[row1] + row_good_counts[row2] - n;
               * d = p - abc;
               * b = abc - row_good_counts[row2];
               * c = abc - row_good_counts[row1];
               */
              b = row_good_counts[row1] - a;    /* 2x2 row sum minus a */
              c = row_good_counts[row2] - a;    /* 2x2 col sum minus a */
              d = p - a - b - c;                /* 2x2 total - (a+b+c) */
              
              n_max = row_good_counts[row1];
              n_min = row_good_counts[row2];
              if (row_good_counts[row2] > n_max)
              {
                  n_max = row_good_counts[row2];
                  n_min = row_good_counts[row1];
              }
              
#if 0
              /* Jaccard distance = 1 - J */
              similarity = 1.0L - (double) a / (a + b + c);
              
              /* Jaccard distance -> zero for 100% similar, so add +1 */
              similarity += 1.0L;
#endif

#if 0
              /* Jaccard-Tanimoto similarity */
              similarity = (double) a / (a + b + c);
#endif

#if 0
              /* #42 CT4: Jaccard-like similarity, with logs */
              similarity = log(a + 1) / log(a + b + c + 1);
#endif

#if 0
              /* #43 CT5: ranges from -1,1 */
              /* useless, since if one vectors is fully present,
               * the similarity = 0
               */
              similarity = (log(1 + a*d) - log(1 + b*c)) /
                           log(1 + 0.25 * p*p);
#endif

#if 0
              /* #38 HL: Harris-Lahey (1978) */
              /* this behaves like a distance, not a similarity */
              /* exhibits very poor disimilarity behavior */
              similarity = (double) a * (d+d + b + c) / (double) (2 * (a + b + c));
              if (b + c + d)
                  similarity += (double) d * (a+a + b + c) / (double) (2 * (b + c + d));
#endif

#if 0
              /* #39 CT1 */
              similarity = log(1 + a + d) / log(1 + p);
#endif

#if 0
              /* Kulczynski-2 */
              /* arithmetic mean conditional probability */
              similarity = 0.5 * ((double) a / (a + b) +
                                  (double) a / (a + c));
#endif

#if 0
              /* #22 Fossum in Holiday et al. (2002) */
              /* this behaves like a distance, not a similarity */
              /* exhibits very poor disimilarity behavior */
              similarity = p * (a - 0.5) * (a - 0.5) /
                           (double) ((a+b) * (a+c));
#endif

#if 0
              /* #33 Sorgenfrei */
              /* use log form to avoid overflows */
              similarity = (double) (a*a) / ((a + b) * (a + c));
#endif

#if 0
              /* #6 Forbes */
              /* results in poor clustering when dividing */
              similarity = (double) (p*a) / ((a + b) * (a + c));
#endif

#if 0
              /* Sorensen-Dice */
              /* similarity = (double) (n) / (n_max + n_min); */
              similarity = (double) (a+a) / (2*a + b + c);
#endif

#if 0
              /* Sokal-Sneath */
              similarity = (double) (3*a) / (3*a + b + c);
#endif

#if 0
              /* simple fraction of num_cols */
              /* does not work as well as Ochiai */
              similarity = (double) a / (double) p;
#endif

#if 0
              /* I've found that log(n) works well for estimating the
               * "importance" of a correlation, given the vector length.
               *
               * Scale by log(num_cols) to normalize it.
               *
               * Deal with potential log(1) = 0 by flooring
               */
              /* this works pretty well, but is too lax compared to Ochiai
               */
              if (n >= 2)
                  similarity = log(n) / log(p);
              else
                  similarity = log(2.0) / log(p);
#endif

#if 1
              /* Ochiai coefficient (binary cosine) */
              /* geometric mean conditional probability */
              similarity = (double) a / sqrt((a + b) * (a + c));
#endif

#if 0
fprintf(stderr, "%f   %03d %03d   %03d %03d %03d %03d\n",
      similarity,
      row_good_counts[row1], row_good_counts[row2],
      a, b, c, d);
#endif
          
/*                dist_matrix[row1][row2] *= sqrt(n_max / n); */
/*                dist_matrix[row1][row2] *= log(n_max + 1) /
                                           log(n + 1);
*/

#if 1
              /* Dividing by Ochiai similarity produces ill-behvaed
               * trees for auto-clustering.  Taking the sqrt() first
               * results in better clusters.
               *
               * Dividing by log10(n_max * similarity + 1.0) doesn't
               * work so well.
               */
              dist_matrix[row1][row2] = dist_matrix[row1][row2] /
                                        sqrt(similarity);
#endif

#if 0
              /* Transform Ochiai the same as regular cosine/pearson
               * similarity, so that it should wind up a metric distance.
               * Take geometric mean the same as we do the geomean
               * heuristic when calculating the euclidian+correlation
               * part of the distance.
               */
              dist_matrix[row1][row2] =
                  sqrt(dist_matrix[row1][row2] *
                       (1.0 + sqrt(1.0 - similarity)));
#endif


#if 0
              /* Transform Ochiai similarly to metric cosine,
               * then take the arithmetic mean.
               *
               * This produces ill-behaved clusters for auto-clustering,
               * even though the heatmaps themselves look reasonable.
               * Removing the sqrt() doesn't help.
               */
              dist_matrix[row1][row2] = 0.5 *
                  (dist_matrix[row1][row2] +
                   sqrt(1.0 - similarity));
#endif

#if 0
              /* geometric mean, add 1 to similarity to avoid zero */
              dist_matrix[row1][row2] = sqrt(dist_matrix[row1][row2] *
                                        (2.0 - similarity));
#endif

#if 0
              dist_matrix[row1][row2] = dist_matrix[row1][row2] +
                                        (1.0 - similarity);
#endif
#if 0
              dist_matrix[row1][row2] = dist_matrix[row1][row2] +
                                        0.5 * (1.0 - similarity);
#endif

#if 0
              /* transform Ochiai (binary cosine) into a metric distance */
              similarity = sqrt(1.0 - similarity);

/*
              dist_matrix[row1][row2] = sqrt(dist_matrix[row1][row2] *
                                        (1.0 + similarity));
*/

              dist_matrix[row1][row2] = dist_matrix[row1][row2] +
                                        similarity;
#endif
          }
          
          /* fudge distance based on spatial proximity */
          if (opt->spatial_flag)
          {
            int32_t x1, x2, y1, y2;

            x1    = x_coord_array[row1];
            x2    = x_coord_array[row2];
            y1    = y_coord_array[row1];
            y2    = y_coord_array[row2];
            
            /* only calculate on valid coordinates */
            if (x1 >= 0 && x2 >= 0 && y1 >= 0 && y2 >= 0)
            {
              /* Euclidian distance between 2 points */
              diff  = x2 - x1;
              dist  = diff * diff;
              diff  = y2 - y1;
              dist += diff * diff;
              dist  = sqrt(dist);

/* not strong enough */
#if 0
              /* exponentiate based on proximity */
              dist_matrix[row1][row2] =
                  pow(dist_matrix[row1][row2], 1.0 + 0.5 / pow(2, dist));
#endif

/* works really well */
#if 1
              /* exponentiate based on proximity */
              dist_matrix[row1][row2] =
                  pow(dist_matrix[row1][row2], 1.0 + 1.0 / pow(2, dist));
#endif

/* too strong */
#if 0
              /* exponentiate based on proximity */
              dist_matrix[row1][row2] =
                  pow(dist_matrix[row1][row2], 1.0 + 2.0 / pow(2, dist));
#endif

#if 0
              /* exponentiate based on proximity */
              dist_matrix[row1][row2] =
                  pow(dist_matrix[row1][row2], 1.0 + 1.0 / (dist*dist));
#endif
#if 0
              /* exponentiate based on proximity */
              dist_matrix[row1][row2] =
                  pow(dist_matrix[row1][row2], 1.0 + 1.0 / dist);
#endif
#if 0
              /* exponentiate based on proximity */
              dist_matrix[row1][row2] =
                  pow(dist_matrix[row1][row2], 1.0 + 2.0 / dist);
#endif
#if 0
              /* exponentiate based on proximity */
              dist_matrix[row1][row2] =
                  pow(dist_matrix[row1][row2], 1.0 + sqrt(1.0 / dist));
#endif
            }
          }
          
          /* store max RMSD for later */
          if (dist_matrix[row1][row2] > max_dist)
              max_dist = dist_matrix[row1][row2];

          n_scanned++;
      }


      /* store values we need to return */
      pthread_mutex_lock(mutex_vars_ptr);
      tdata_ptr->n_scanned = n_scanned;
      tdata_ptr->max_dist  = max_dist;
      pthread_mutex_unlock(mutex_vars_ptr);
      
      /* finished current work */
      pthread_mutex_lock(mutex_flags_ptr);
      t_flags_array[t_index] = 1;
      pthread_mutex_unlock(mutex_flags_ptr);
      t_flag = 1;

      /* return as normal function, since it wasn't launched as thread */
      if (t_index == 0)
          return(NULL);

      /* wake up main thread if all threads are now finished */
      pthread_mutex_lock(&mutex_wait_main);
      if (--n_threads_left <= 0)
          pthread_cond_signal(&cond_wait_main);
      pthread_mutex_unlock(&mutex_wait_main);
    }

    /* we should never get here, but put these here just in case */
    pthread_mutex_destroy(mutex_flags_ptr);
    pthread_mutex_destroy(mutex_vars_ptr);
    pthread_mutex_destroy(mutex_wait_ptr);
    pthread_cond_destroy(cond_wait_ptr);
    pthread_exit(NULL);
}


void fill_dist_matrix(double **data_matrix,
                      double **dist_matrix,
                      double  *weights,
                      int32_t num_rows, int32_t num_cols,
                      int32_t *x_coord_array,
                      int32_t *y_coord_array,
                      struct options *opt)
{
    int       t_flags_array[MAX_THREADS] = { 0 };
    int       pt_ret[MAX_THREADS];
    pthread_t threads[MAX_THREADS];
    int32_t   max_threads = opt->threads;
    int32_t   cur_threads = max_threads;
    int32_t   t, count_finished;
    int32_t   n_left;

    struct    dist_thread_data *tdata_array = NULL;
    struct    timespec ns_to_wait;

    int32_t  *row_good_counts = NULL;
    double   *row_ptr1;
    double   *flip_row1 = NULL;
    double    value1;
    double    max_dist = DBL_MIN, worst_dist;
    int32_t   row1, row2, col;
    int32_t   n, n1;
    int32_t   never_missing_flag = 1;
    
    char      pearson_flag    = opt->pearson_flag;
    char      cosine_flag     = opt->cosine_flag;
    char      heuristic_flag  = opt->heuristic_flag;
    char      similarity_flag = opt->similarity_flag;
    
    /* start and end times for system clock
     * we're using 1 second resolution due to ANSI-compatibility
     * we could use clock_gettime() for high-precision, but that is POSIX
     */
    time_t start_time, end_time;
    double delta_time, time_left;
    uint64_t n_scanned, n_left_to_go;
    char time_str1[20], time_str2[20];    /* printf can't re-eval one str */

    /* sanity check max threads */
    if (max_threads > MAX_THREADS)
        max_threads = MAX_THREADS;

    /* allocate thread data */
    tdata_array = (struct dist_thread_data *)
        calloc(max_threads, sizeof(struct dist_thread_data));


    if (heuristic_flag && pearson_flag == 0 && cosine_flag == 0)
        flip_row1 = (double *) calloc(num_cols, sizeof(double));

    /* always count how many good values we have */
    if (1 || similarity_flag)
    {
        row_good_counts = (int32_t *) calloc(num_rows, sizeof(int32_t));
        
        for (row1 = 0; row1 < num_rows; row1++)
        {
            row_ptr1 = data_matrix[row1];
            posix_madvise(row_ptr1, num_cols * sizeof(double),
                          POSIX_MADV_WILLNEED);
            
            n = 0;
            for (col = 0; col < num_cols; col++)
            {
                value1 = row_ptr1[col];
                
                if (value1 != MISSING)
                    n++;
            }
            
            row_good_counts[row1] = n;
        }
    }


    pthread_mutex_init(&mutex_wait_main, NULL);
    pthread_cond_init(&cond_wait_main,   NULL);

    ns_to_wait.tv_sec  = 0;
    ns_to_wait.tv_nsec = 100;

    /* initialize constant part of thread data */
    for (t = 0; t < max_threads; t++)
    {
        pthread_mutex_init(&tdata_array[t].mutex_flags, NULL);
        pthread_mutex_init(&tdata_array[t].mutex_vars,  NULL);
        pthread_mutex_init(&tdata_array[t].mutex_wait,  NULL);
        pthread_cond_init(&tdata_array[t].cond_wait,  NULL);

        tdata_array[t].t_flags_array   = t_flags_array;
        tdata_array[t].t_index         = t;

        tdata_array[t].data_matrix     = data_matrix;
        tdata_array[t].dist_matrix     = dist_matrix;
        tdata_array[t].weights         = weights;
        tdata_array[t].flip_row1       = flip_row1;
        tdata_array[t].row_good_counts = row_good_counts;
        tdata_array[t].x_coord_array   = x_coord_array;
        tdata_array[t].y_coord_array   = y_coord_array;
        tdata_array[t].num_cols        = num_cols;
        tdata_array[t].row2_inc        = max_threads;
        tdata_array[t].opt             = opt;
    }

    /* launch each worker thread, which will wait for data */
    for (t = 1; t < max_threads; t++)
    {
        pt_ret[t] = pthread_create(&threads[t],
                        NULL,
                        fill_dist_matrix_thread,
                        (void *) &tdata_array[t]);
    }

    
    fprintf(stderr, "NumRowsLeft: %d\n", num_rows);
    n_scanned = 0;

    start_time = time(NULL);
    for (row1 = num_rows - 1; row1 >= 0; --row1)
    {
        if (n_scanned && (row1 + 1) % (max_threads * 100) == 0)
        {
            end_time = time(NULL);
            delta_time = difftime(end_time, start_time);
            
            n_left_to_go = 0.5 * (row1 + 1) * row1;

            /* accumulate at least 15 seconds before estimating */
            if (delta_time >= 15.0)
            {
                time_left = (double) n_left_to_go * delta_time /
                            (double) n_scanned;

                fprintf(stderr, "NumRowsLeft: %d   TimeElapsed: %s   TimeLeft: %s\n",
                    row1 + 1,
                    seconds_to_str(time_str1, 20, delta_time),
                    seconds_to_str(time_str2, 20, time_left));
            }
            else
            {
                fprintf(stderr, "NumRowsLeft: %d   TimeElapsed: %s   TimeLeft: %s\n",
                    row1 + 1,
                    seconds_to_str(time_str1, 20, delta_time),
                    "estimating...");
            }
        }

        /* allocate new distance matrix row */
        dist_matrix[row1] = (double *) calloc(row1+1, sizeof(double));
        row_ptr1 = data_matrix[row1];
        
        n1 = row_good_counts[row1];

        /* dataset has some missing data, note this for later */
        if (n1 != num_cols)
            never_missing_flag = 0;
        
        /* flip vector for anticorrelation heuristic */
        if (heuristic_flag && pearson_flag == 0 && cosine_flag == 0)
        {
            double min=MISSING, max=MISSING, mid=MISSING;
            
            /* find min/max */
            for (col = 0; col < num_cols; col++)
            {
                if (row_ptr1[col] == MISSING)
                {
                    continue;
                }
                
                /* if min hasn't been stored, max hasn't been either */
                if (min == MISSING)
                {
                    min = row_ptr1[col];
                    max = row_ptr1[col];
                    
                    continue;
                }
                
                if (row_ptr1[col] < min)
                    min = row_ptr1[col];
                if (row_ptr1[col] > max)
                    max = row_ptr1[col];
            }
            
            mid = MISSING;
            if (min != MISSING)
            {
                mid = 0.5 * (max - min) + min;
            }
            
            /* flip vector1 */
            memcpy(flip_row1, row_ptr1, num_cols * sizeof(double));
            if (mid != MISSING)
            {
                for (col = 0; col < num_cols; col++)
                {
                    if (flip_row1[col] == MISSING)
                    {
                        continue;
                    }

                    flip_row1[col] = -(flip_row1[col] - mid) + mid;
                }
            }
        }

        
        /* initialize the rows to scan
         *
         * This loop, in particular, highlights the dangers of not using
         * mutexes, since incorrect results can sometimes occur if everything
         * isn't properly mutexed.
         */
        pthread_mutex_lock(&mutex_wait_main);
        n_threads_left = max_threads - 1;
        if (row1 < max_threads)
            n_threads_left = row1 - 1;
        if (n_threads_left < 0)
            n_threads_left = 0;
        n_left = n_threads_left;
        pthread_mutex_unlock(&mutex_wait_main);
        for (t = 1; t < max_threads && t < row1; t++)
        {
            pthread_mutex_lock(&tdata_array[t].mutex_vars);
            tdata_array[t].row_ptr1   = row_ptr1;
            tdata_array[t].row1       = row1;
            tdata_array[t].n1         = n1;
            tdata_array[t].row2_start = t;
            tdata_array[t].n_scanned  = 0;
            pthread_mutex_unlock(&tdata_array[t].mutex_vars);
            
            /* signal worker thread to process data */
            pthread_mutex_lock(&tdata_array[t].mutex_flags);
            t_flags_array[t] = 2;
            pthread_mutex_unlock(&tdata_array[t].mutex_flags);

            /* wake up worker thread */
            pthread_mutex_lock(&tdata_array[t].mutex_wait);
            pthread_cond_signal(&tdata_array[t].cond_wait);
            pthread_mutex_unlock(&tdata_array[t].mutex_wait);
        }
        cur_threads = t;


        /* call the worker function in this main thread while the others
         * are finishing
         */
        t = 0;
        tdata_array[t].row_ptr1   = row_ptr1;
        tdata_array[t].row1       = row1;
        tdata_array[t].n1         = n1;
        tdata_array[t].row2_start = t;
        tdata_array[t].n_scanned  = 0;
        t_flags_array[t]          = 2;
        fill_dist_matrix_thread((void *) &tdata_array[t]);


        /* sleep until other threads finish */
        pthread_mutex_lock(&mutex_wait_main);
        if (n_threads_left > 0)
            pthread_cond_wait(&cond_wait_main, &mutex_wait_main);
        pthread_mutex_unlock(&mutex_wait_main);

        /* nastiness in case there is a spurious wakeup */
        pthread_mutex_lock(&mutex_wait_main);
        n_left = n_threads_left;
        pthread_mutex_unlock(&mutex_wait_main);
        if (n_left > 0)
        {
            fprintf(stderr, "Recovering from spurious main thread wakeup:  %d\n",
                n_left);

            /* wait for other threads to finish */
            do
            {
                count_finished = 1;
                for (t = 1; t < cur_threads; t++)
                {
                    pthread_mutex_lock(&tdata_array[t].mutex_flags);
                    if (t_flags_array[t] == 1)
                    {
                        count_finished++;
                    }
                    pthread_mutex_unlock(&tdata_array[t].mutex_flags);
                }
                
                /* Small numbers of threads can deadlock, so wait a bit.
                 * We are highly unlikely to ever get here anyways, so we
                 * aren't worried about performance at this point.
                 */
                nanosleep(&ns_to_wait, NULL);
            } while (count_finished < cur_threads);
        }

        /* iterate through return values */
        for (t = 0; t < cur_threads; t++)
        {
            /* flag thread to await new data */
            pthread_mutex_lock(&tdata_array[t].mutex_flags);
            t_flags_array[t] = 0;
            pthread_mutex_unlock(&tdata_array[t].mutex_flags);
        
            pthread_mutex_lock(&tdata_array[t].mutex_vars);
            n_scanned += tdata_array[t].n_scanned;
            if (tdata_array[t].max_dist > max_dist)
            {
                max_dist = tdata_array[t].max_dist;
            }
            pthread_mutex_unlock(&tdata_array[t].mutex_vars);
        }
        
        /* free data matrix row to minimize memory footprint */
        if (data_matrix[row1])
        {
#if 0
            posix_madvise(data_matrix[row1], num_cols * sizeof(double),
                          POSIX_MADV_DONTNEED);
#endif

            free(data_matrix[row1]);

#if 0
            /* release data off top of heap */
            malloc_trim(0);
#endif
            
            /* we're never going to intentionally use it again, but set it
             * to NULL to make sure that we crash if we accidentally do,
             * rather than risk bogus results we don't know are bogus
             */
            data_matrix[row1] = NULL;
        }
    }
    if (data_matrix)
        free(data_matrix);
    

    /* set bad dists to 2x max dist */
    if (never_missing_flag == 0)
    {
        worst_dist = max_dist + max_dist;
        for (row1 = 1; row1 < num_rows; row1++)
        {
            row_ptr1 = dist_matrix[row1];
            posix_madvise(row_ptr1, row1 * sizeof(double), POSIX_MADV_WILLNEED);

            for (row2 = 0; row2 < row1; row2++)
            {
                if (row_ptr1[row2] == DBL_MAX)
                {
                    row_ptr1[row2] = worst_dist;
                }
            }
        }
    }


    /* signal all the threads to terminate */
    for (t = 1; t < max_threads; t++)
    {
        pthread_mutex_lock(&tdata_array[t].mutex_flags);
        t_flags_array[t] = 3;
        pthread_mutex_unlock(&tdata_array[t].mutex_flags);
        
        pthread_mutex_lock(&tdata_array[t].mutex_wait);
        pthread_cond_signal(&tdata_array[t].cond_wait);
        pthread_mutex_unlock(&tdata_array[t].mutex_wait);
    }
    /* wait for all threads to finish */
    for (t = 1; t < max_threads; t++)
    {
        pthread_join(threads[t], NULL);
    }
    
    
    pthread_mutex_destroy(&mutex_wait_main);
    pthread_cond_destroy(&cond_wait_main);
    

    end_time = time(NULL);
    delta_time = difftime(end_time, start_time);
    fprintf(stderr, "Time spent calculating distance matrix:\t%s\n",
            seconds_to_str(time_str1, 20, delta_time));


    if (tdata_array)
        free(tdata_array);
    
    if (flip_row1)
        free(flip_row1);
    if (row_good_counts)
        free(row_good_counts);
    
    return;
}


void raise_dist_matrix_to_power(double **dist_matrix,
                                int32_t num_rows, int32_t num_cols,
                                double power)
{
    double  *dptr;
    int32_t  row, col;
    
    /* do nothing */
    if (power == 0.0 || power == 1.0)
        return;

    fprintf(stderr, "Exponentiating distances:  X^%g\n", power);
    
    /* square it */
    if (power == 2.0)
    {
        for (row = 0; row < num_rows; row++)
        {
            dptr = dist_matrix[row];
        
            for (col = 0; col < row; col++)
            {
                if (dptr[col] != DBL_MAX)
                    dptr[col] *= dptr[col];
            }
        }
    }
    /* multiply by sqrt */
    else if (power == 1.5)
    {
        for (row = 0; row < num_rows; row++)
        {
            dptr = dist_matrix[row];
        
            for (col = 0; col < row; col++)
            {
                if (dptr[col] != DBL_MAX)
                    dptr[col] *= sqrt(dptr[col]);
            }
        }
    }
    else if (power == 1.75)
    {
        double temp;
    
        for (row = 0; row < num_rows; row++)
        {
            dptr = dist_matrix[row];
        
            for (col = 0; col < row; col++)
            {
                if (dptr[col] != DBL_MAX)
                {
                    temp = sqrt(dptr[col]);
                
                    dptr[col] *= temp * sqrt(temp);
                }
            }
        }
    }
    else if (power == 1.25)
    {
        for (row = 0; row < num_rows; row++)
        {
            dptr = dist_matrix[row];
        
            for (col = 0; col < row; col++)
            {
                if (dptr[col] != DBL_MAX)
                    dptr[col] *= sqrt(sqrt(dptr[col]));
            }
        }
    }
    else
    {
        for (row = 0; row < num_rows; row++)
        {
            dptr = dist_matrix[row];
        
            for (col = 0; col < row; col++)
            {
                if (dptr[col] != DBL_MAX)
                    dptr[col] = pow(dptr[col], power);
            }
        }
    }
}


void print_data_matrix(double **data_matrix,
                       char **row_name_array, char **col_name_array,
                       int32_t num_rows, int32_t num_cols)
{
    double  *row_ptr;
    char    str_14g[50], str_06f[50];
    int32_t len_14g, len_06f;
    int32_t i, j;

    /* header */
    printf("%s", "DataMatrix");
    for (i = 0; i < num_cols; i++)
    {
        printf("\t%s", col_name_array[i]);
    }
    printf("\n");
    
    for (i = 0; i < num_rows; i++)
    {
        printf("%s", row_name_array[i]);
        
        row_ptr = data_matrix[i];
        
        for (j = 0; j < num_cols; j++)
        {
            /* leave missing or bad spots as blank */
            if (row_ptr[j] == MISSING)
            {
                printf("\t");
            }
            else
            {
                len_14g = snprintf(str_14g, 100 * sizeof(char),
                                   "%.14g", row_ptr[j]);
                len_06f = snprintf(str_06f, 100 * sizeof(char),
                                   "%.6f",  row_ptr[j]);

                if (len_06f < len_14g)
                    printf("\t%s", str_06f);
                else
                    printf("\t%s", str_14g);
            }
        }
        printf("\n");
        
        free(row_ptr);
        data_matrix[i] = NULL;
    }

    free(data_matrix);
}


void print_dist_matrix(double **dist_matrix, char **row_name_array,
                       int32_t num_rows)
{
    double *row_ptr;
    int32_t i, j;

    if (dist_matrix)
    {
        for (i = 0; i < num_rows; i++)
        {
            printf("%s", row_name_array[i]);

            row_ptr = dist_matrix[i];

            for (j = 0; j <= i; j++)
            {
                printf("\t");

                if (row_ptr[j] != MISSING)
                {
                    printf("%.13g", row_ptr[j]);
                }
            }
            printf("\n");

            free(row_ptr);
            dist_matrix[i] = NULL;
        }
    
        free(dist_matrix);
    }
}


int32_t main(int32_t argc, char *argv[])
{
    double          **data_matrix     = NULL;
    double          **dist_matrix     = NULL;
    double           *weights         = NULL;
    double           *dptr;
    double            weight;
    double            lod;
    char            **row_name_array     = NULL;
    char            **col_name_array     = NULL;
    char            **filename_array     = NULL;
    char             *filename           = NULL;
    char             *newick_str         = NULL;
    char             *tree_string        = NULL;    /* original from file */
    char             *str_ptr;
    int32_t          *file_col_counts    = NULL;
    int32_t          *x_coord_array      = NULL;    /* X for spatial data */
    int32_t          *y_coord_array      = NULL;    /* Y for spatial data */
    int32_t           num_rows  = 0;
    int32_t           num_cols  = 0;
    int32_t           num_cols_chunk = 0;    /* #cols in file last read in */
    int32_t           num_files = 0;
    int32_t           i, col, max_num_cols_chunk;
    int32_t           tree_return_value = 0;
    int32_t           n_clusters;
    int32_t           x, y;

    /* tree weight stuff */
    struct tree_node **node_ptr_array = NULL;
    struct tree_node **leaf_ptr_array = NULL;
    struct tree_node  *tree_root;
    double             eff_count;
    int32_t            num_nodes;
    int32_t            num_leaves;

    struct options opt;
    int32_t            took_log2_flag = 0;


    /* set various non-default malloc options */
#ifdef __GLIBC__

#if (__GLIBC__ >= 2 && __GLIBC_MINOR__ >= 15)
    /* fprintf(stderr, "glibc malloc(): setting M_ARENA_MAX=1\n"); */
    mallopt(M_ARENA_MAX, 1);   /* we're not worried about thread contention */
#endif

#if (__GLIBC__ >= 2 && __GLIBC_MINOR__ >=  3)
    /* fprintf(stderr, "glibc malloc(): setting M_MXFAST=0\n"); */
    mallopt(M_MXFAST, 0);      /* reduce memory footprint and fragmentation */
#endif

#if 1
    /* size of pointer or double = 8 bytes */
    /* fprintf(stderr, "glibc malloc(): setting M_TRIM_THRESHOLD=8\n"); */
    mallopt(M_TRIM_THRESHOLD, 8);
#endif

#if 0
    /* setting M_TRIM_THRESHOLD disables dynamic M_MMAP_THRESHOLD
     *  adjustment as well, so if M_TRIM_THRESHOLD is set,
     *  M_MMAP_THRESHOLD defaults to 128 KiB anyways
     */
    fprintf(stderr, "glibc malloc(): setting M_MMAP_THRESHOLD=131072\n");
    mallopt(M_MMAP_THRESHOLD, 131072);     /* faster, less footprint */
#endif

#if 1
    /* default of 65536 may be too small
     *
     * 2147483648 acts like it may be disabling MMAP allocations,
     * potentially represented by signed int32 that goes negative at 2^31
     */
    /* fprintf(stderr, "glibc malloc(): setting M_MMAP_MAX=2147483647  (2^31 - 1)\n"); */
    mallopt(M_MMAP_MAX, 2147483647);
#endif

#if 0
    /* doesn't appear to help in terms of max resident memory */
    fprintf(stderr, "glibc malloc(): setting M_TOP_PAD=0\n");
    mallopt(M_TOP_PAD, 0);
#endif

#endif


    /* initialize options */
    memset(&opt, 0, sizeof(struct options));
    opt.calc_euclidian_flag    = 1;    /* default to RMSD */
    opt.rmsd_vs_euclidian_flag = 1;    /* default to RMSD */
    opt.ties_fewest_flag       = 1;
    opt.threads                = 1;


    if (argc > 1)
    {
        filename_array = (char **) calloc(argc, sizeof(char *));

        for (i = 1; i < argc; i++)
        {
            if (strncmp(argv[i], "--", 2) == 0)
            {
                if (strcmp(argv[i], "--ignore-weak") == 0)
                {
                    opt.ignore_weak_flag = 1;
                }
                else if (strcmp(argv[i], "--absolute") == 0)
                {
                    opt.absolute_flag = 1;
                }
                else if (strcmp(argv[i], "--signed") == 0)
                {
                    opt.absolute_flag = 0;
                }
                else if (strcmp(argv[i], "--pearson") == 0)
                {
                    opt.pearson_flag = 1;
                    opt.cosine_flag  = 0;
                    opt.set_any_distance_flag = 1;
                }
                else if (strcmp(argv[i], "--cosine") == 0)
                {
                    opt.cosine_flag  = 1;
                    opt.pearson_flag = 0;
                    opt.set_any_distance_flag = 1;
                }
                else if (strcmp(argv[i], "--euclidian") == 0)
                {
                    opt.rmsd_vs_euclidian_flag = 0;
                    opt.cityblock_flag = 0;
                    opt.pearson_flag = 0;
                    opt.cosine_flag  = 0;
                    opt.set_any_distance_flag = 1;
                    opt.calc_euclidian_flag = 1;
                }
                else if (strcmp(argv[i], "--rmsd") == 0)
                {
                    opt.rmsd_vs_euclidian_flag = 1;
                    opt.cityblock_flag = 0;
                    opt.pearson_flag = 0;
                    opt.cosine_flag  = 0;
                    opt.set_any_distance_flag = 1;
                    opt.calc_euclidian_flag = 1;
                }
                else if (strcmp(argv[i], "--cityblock") == 0)
                {
                    opt.rmsd_vs_euclidian_flag = 0;
                    opt.cityblock_flag = 1;
                    opt.pearson_flag = 0;
                    opt.cosine_flag  = 0;
                    opt.set_any_distance_flag = 1;
                    opt.calc_euclidian_flag = 1;
                }
                /* l-N-norm distances */
                else if (strncmp(argv[i], "--minkowski=",
                         strlen("--minkowski=")) == 0 &&
                         strlen(argv[i]) > strlen("--minkowski="))
                {
                    opt.minkowski =
                        atof(argv[i] + strlen("--minkowski="));

                    opt.rmsd_vs_euclidian_flag = 0;
                    opt.cityblock_flag = 0;
                    opt.pearson_flag = 0;
                    opt.cosine_flag  = 0;
                    opt.set_any_distance_flag = 1;
                    opt.calc_euclidian_flag = 1;
                }
                /* old minowski typo, retain for backwards compatability */
                else if (strncmp(argv[i], "--minowski=",
                         strlen("--minowski=")) == 0 &&
                         strlen(argv[i]) > strlen("--minowski="))
                {
                    opt.minkowski =
                        atof(argv[i] + strlen("--minowski="));

                    opt.rmsd_vs_euclidian_flag = 0;
                    opt.cityblock_flag = 0;
                    opt.pearson_flag = 0;
                    opt.cosine_flag  = 0;
                    opt.set_any_distance_flag = 1;
                    opt.calc_euclidian_flag = 1;
                }
                else if (strncmp(argv[i], "--distpow=",
                         strlen("--distpow=")) == 0 &&
                         strlen(argv[i]) > strlen("--distpow="))
                {
                    opt.dist_pow =
                        atof(argv[i] + strlen("--distpow="));
                }
                else if (strcmp(argv[i], "--depower") == 0)
                {
                    opt.depower_flag = 1;
                }
                else if (strcmp(argv[i], "--include-weak") == 0)
                {
                    opt.ignore_weak_flag = 0;
                }
                else if (strcmp(argv[i], "--mean-center") == 0)
                {
                    opt.mean_center_flag = 1;
                }
                else if (strcmp(argv[i], "--unit-max-mag") == 0)
                {
                    opt.maxmag_flag = 1;
                }
                else if (strcmp(argv[i], "--unit-variance") == 0)
                {
                    opt.unit_variance_flag = 1;
                    opt.mean_center_flag   = 1;
                }
                else if (strcmp(argv[i], "--scale-10000") == 0)
                {
                    opt.scale_ten_thousand_flag = 1;
                    opt.mean_center_flag = 0;
                    opt.unit_variance_flag = 0;
                }
                else if (strcmp(argv[i], "--unaltered") == 0)
                {
                    opt.mean_center_flag   = 0;
                    opt.unit_variance_flag = 0;
                }
                else if (strcmp(argv[i], "--heuristic") == 0)
                {
                    opt.geomean_flag     = 0;
                    opt.heuristic_flag   = 1;
                    opt.mostly_trig_flag = 0;
                    opt.calc_euclidian_flag = 1;
                }
                else if (strcmp(argv[i], "--mostly-trig") == 0)
                {
                    opt.geomean_flag     = 0;
                    opt.heuristic_flag   = 0;
                    opt.mostly_trig_flag = 1;
                    opt.calc_euclidian_flag = 1;
                }
                else if (strcmp(argv[i], "--geomean") == 0)
                {
                    opt.geomean_flag     = 1;
                    opt.heuristic_flag   = 0;
                    opt.mostly_trig_flag = 0;
                    opt.calc_euclidian_flag = 1;
                }
                else if (strcmp(argv[i], "--spatial") == 0)
                {
                    opt.spatial_flag     = 1;
                }
                else if (strcmp(argv[i], "--similarity") == 0)
                {
                    opt.similarity_flag = 1;
                }
                else if (strcmp(argv[i], "--read-dists") == 0)
                {
                    opt.dists_from_file_flag  = 1;
                    opt.build_tree_flag       = 1;
                    opt.tree_from_file_flag   = 0;
                    opt.dump_data_matrix_flag = 0;
                }
                else if (strcmp(argv[i], "--read-tree") == 0)
                {
                    opt.tree_from_file_flag   = 1;
                    opt.dists_from_file_flag  = 0;
                    opt.build_tree_flag       = 0;
                    opt.dump_data_matrix_flag = 0;
                }
                else if (strcmp(argv[i], "--dump-matrix") == 0 ||
                         strcmp(argv[i], "--output-data") == 0)
                {
                    opt.dump_data_matrix_flag = 1;
                    opt.dists_from_file_flag  = 0;
                    opt.build_tree_flag       = 0;
                    opt.tree_from_file_flag   = 0;
                    opt.clusters_flag         = 0;
                    opt.tree_weights_flag     = 0;
                }
                else if (strcmp(argv[i], "--leaf-weights") == 0)
                {
                    opt.tree_weights_flag     = 1;
                    opt.tree_defaults_flag    = 1;
                    opt.dump_data_matrix_flag = 0;
                    opt.clusters_flag         = 0;
                }
                else if (strcmp(argv[i], "--leaf-custom") == 0)
                {
                    opt.tree_weights_flag     = 1;
                    opt.tree_defaults_flag    = 0;
                    opt.dump_data_matrix_flag = 0;
                    opt.clusters_flag         = 0;
                }
                else if (strcmp(argv[i], "--clusters") == 0)
                {
                    opt.clusters_flag         = 1;
                    opt.tree_weights_flag     = 0;
                    opt.dump_data_matrix_flag = 0;
                }
                /* don't penalize small clusters */
                else if (strcmp(argv[i], "--cluster-no-merge") == 0)
                {
                    opt.clusters_no_merge_flag = 1;
                }
                /* merge small clusters and penalize them */
                else if (strcmp(argv[i], "--cluster-merge") == 0)
                {
                    opt.clusters_no_merge_flag = 0;
                }

                else if (strcmp(argv[i], "--log2") == 0)
                {
                    opt.log2_prior_flag       = 1;
                }
                else if (strcmp(argv[i], "--unlog2") == 0)
                {
                    opt.unlog2_last_flag      = 1;
                }
                else if (strcmp(argv[i], "--norm-median") == 0)
                {
                    opt.norm_median_flag      = 1;
                }
                else if (strcmp(argv[i], "--impute-col-min") == 0)
                {
                    opt.impute_col_min_half_flag            = 1;
                    opt.impute_global_min_half_flag         = 0;
                    opt.impute_col_global_min_half_flag     = 0;
                    opt.impute_row_col_global_min_half_flag = 0;
                }
                else if (strcmp(argv[i], "--impute-global") == 0)
                {
                    opt.impute_global_min_half_flag         = 1;
                    opt.impute_col_min_half_flag            = 0;
                    opt.impute_col_global_min_half_flag     = 0;
                    opt.impute_row_col_global_min_half_flag = 0;
                }
                else if (strcmp(argv[i], "--impute-col-geo") == 0)
                {
                    opt.impute_global_min_half_flag         = 0;
                    opt.impute_col_min_half_flag            = 0;
                    opt.impute_col_global_min_half_flag     = 1;
                    opt.impute_row_col_global_min_half_flag = 0;
                }
                else if (strcmp(argv[i], "--impute-row-col") == 0)
                {
                    opt.impute_global_min_half_flag         = 0;
                    opt.impute_col_min_half_flag            = 0;
                    opt.impute_col_global_min_half_flag     = 0;
                    opt.impute_row_col_global_min_half_flag = 1;
                }
                else if (strcmp(argv[i], "--floor-to-one") == 0)
                {
                    opt.floor_to_one_flag = 1;
                }
                else if (strcmp(argv[i], "--floor-to-lod") == 0)
                {
                    opt.floor_to_lod_flag = 1;
                }
                else if (strcmp(argv[i], "--floor-lod-ub") == 0)
                {
                    opt.floor_to_lod_flag = 2;   /* use upper bound instead */
                }
                else if (strncmp(argv[i], "--floor-value=",
                         strlen("--floor-value=")) == 0 &&
                         strlen(argv[i]) > strlen("--floor-value="))
                {
                    opt.floor_to_value    =
                        atof(argv[i] + strlen("--floor-value="));
                }
                else if (strcmp(argv[i], "--impute-later") == 0)
                {
                    opt.impute_later_flag = 1;
                }

                else if (strcmp(argv[i], "--transpose-first") == 0)
                {
                    opt.transpose_first_flag  = 1;
                }
                else if (strcmp(argv[i], "--transpose-last") == 0)
                {
                    opt.transpose_last_flag   = 1;
                }

                else if (strncmp(argv[i], "--nclusters=",
                         strlen("--nclusters=")) == 0 &&
                         strlen(argv[i]) > strlen("--nclusters="))
                {
                    opt.target_n_clusters     =
                        atol(argv[i] + strlen("--nclusters="));
                }
                else if (strncmp(argv[i], "--n-clusters=",
                         strlen("--n-clusters=")) == 0 &&
                         strlen(argv[i]) > strlen("--n-clusters="))
                {
                    opt.target_n_clusters     =
                        atol(argv[i] + strlen("--n-clusters="));
                }
                /* also accept --target-clusters= */
                else if (strncmp(argv[i], "--target-clusters=",
                         strlen("--target-clusters=")) == 0 &&
                         strlen(argv[i]) > strlen("--target-clusters="))
                {
                    opt.target_n_clusters     =
                        atol(argv[i] + strlen("--target-clusters="));
                }
                else if (strncmp(argv[i], "--threads=",
                         strlen("--threads=")) == 0 &&
                         strlen(argv[i]) > strlen("--threads="))
                {
                    opt.threads               =
                        atol(argv[i] + strlen("--threads="));
                }
                else if (strncmp(argv[i], "--filter-log2-sd=",
                         strlen("--filter-log2-sd=")) == 0 &&
                         strlen(argv[i]) > strlen("--filter-log2-sd="))
                {
                    opt.filter_log2_sd_cutoff      =
                        atof(argv[i] + strlen("--filter-log2-sd="));
                    opt.filter_log2_flag = 1;
                }
                else if (strncmp(argv[i], "--filter-log2-mean=",
                         strlen("--filter-log2-mean=")) == 0 &&
                         strlen(argv[i]) > strlen("--filter-log2-mean="))
                {
                    opt.filter_log2_mean_cutoff      =
                        atof(argv[i] + strlen("--filter-log2-mean="));
                    opt.filter_log2_flag = 1;
                }
                else if (strncmp(argv[i], "--filter-unlog-sd=",
                         strlen("--filter-unlog-sd=")) == 0 &&
                         strlen(argv[i]) > strlen("--filter-unlog-sd="))
                {
                    opt.filter_unlog_sd_cutoff      =
                        atof(argv[i] + strlen("--filter-unlog-sd="));
                    opt.filter_unlog_flag = 1;
                }
                else if (strncmp(argv[i], "--filter-unlog-max=",
                         strlen("--filter-unlog-max=")) == 0 &&
                         strlen(argv[i]) > strlen("--filter-unlog-max="))
                {
                    opt.filter_unlog_max_signal_cutoff      =
                        atof(argv[i] + strlen("--filter-unlog-max="));
                    opt.filter_unlog_flag = 1;
                }
                else if (strncmp(argv[i], "--filter-unlog-mean=",
                         strlen("--filter-unlog-mean=")) == 0 &&
                         strlen(argv[i]) > strlen("--filter-unlog-mean="))
                {
                    opt.filter_unlog_mean_cutoff      =
                        atof(argv[i] + strlen("--filter-unlog-mean="));
                    opt.filter_unlog_flag = 1;
                }
                else if (strncmp(argv[i], "--filter-present=N",
                         strlen("--filter-present=")) == 0 &&
                         strlen(argv[i]) > strlen("--filter-present="))
                {
                    opt.filter_present_cutoff      =
                        atof(argv[i] + strlen("--filter-present="));

                    /* up to the user to format options correctly,
                     * no error checking for proper format
                     */
                    if (strstr(argv[i], ",") &&
                        strlen(strstr(argv[i], ",")) > 1)
                    {
                        opt.filter_present_mag =
                            atof(strstr(argv[i], ",") + 1);
                    }
                }
                else if (strncmp(argv[i], "--unit-stretch=N",
                         strlen("--unit-stretch=")) == 0 &&
                         strlen(argv[i]) > strlen("--unit-stretch="))
                {
                    opt.unit_stretch =
                        atof(argv[i] + strlen("--unit-stretch="));
                }

                else if (strncmp(argv[i], "--weight-file=",
                         strlen("--weight-file=")) == 0 &&
                         strlen(argv[i]) > strlen("--weight-file="))
                {
                    opt.weight_filename =
                        argv[i] + strlen("--weight-file=");
                }
                
                else if (strcmp(argv[i], "--upgma") == 0)
                    opt.linkage_method = UPGMA;
                else if (strcmp(argv[i], "--wpgma") == 0)
                    opt.linkage_method = WPGMA;
                else if (strcmp(argv[i], "--slink") == 0)
                    opt.linkage_method = SLINK;
                else if (strcmp(argv[i], "--clink") == 0)
                    opt.linkage_method = CLINK;
                else if (strcmp(argv[i], "--wardu") == 0)
                    opt.linkage_method = WARDU;
                else if (strcmp(argv[i], "--ward2") == 0)
                    opt.linkage_method = WARD2;
                else if (strcmp(argv[i], "--upgmcu") == 0)
                    opt.linkage_method = UPGMCU;
                else if (strcmp(argv[i], "--upgmc2") == 0)
                    opt.linkage_method = UPGMC2;
                else if (strcmp(argv[i], "--wpgmcu") == 0)
                    opt.linkage_method = WPGMCU;
                else if (strcmp(argv[i], "--wpgmc2") == 0)
                    opt.linkage_method = WPGMC2;

                else if (strcmp(argv[i], "--tree-flip-size") == 0)
                {
                    opt.tree_flip_size_flag = 1;
                    opt.tree_flip_edge_flag = 0;
                }
                else if (strcmp(argv[i], "--tree-flip-edge") == 0)
                {
                    opt.tree_flip_edge_flag = 1;
                    opt.tree_flip_size_flag = 0;
                }

                else if (strcmp(argv[i], "--ties-fewest") == 0)
                {
                    opt.ties_fewest_flag = 1;
                    opt.ties_random_flag = 0;
                }
                else if (strcmp(argv[i], "--ties-order") == 0)
                {
                    opt.ties_fewest_flag = 0;
                    opt.ties_random_flag = 0;
                }
                else if (strcmp(argv[i], "--ties-random") == 0)
                {
                    opt.ties_fewest_flag = 0;
                    opt.ties_random_flag = 1;
                }

                else if (strcmp(argv[i], "--multi-line") == 0)
                    opt.one_line_flag = 0;
                else if (strcmp(argv[i], "--single-line") == 0)
                    opt.one_line_flag = 1;

                else
                {
                    printf("Usage: hcdist [OPTIONS] input_text_file [#jumbles [integer seed]]\n");
                    printf("\n");
                    printf("If multiple file names are given, they *cannot* be all digits.\n");
                    printf("Multiple filenames are currently only used for calculation of\n");
                    printf("weighted \"co-clustering\" distances.\n");
                    printf("\n");
                    printf("If a 2nd non-option is an integer, it is taken as the number of \"jumbles\",\n");
                    printf("the number of trees to generate from randomly reordered copies of the\n");
                    printf("original distance matrix.  If --ties-random is also specified, the\n");
                    printf("distance matrix is copied as-is, instead of jumbled, since jumbling would\n");
                    printf("not increase the randomness of the resulting trees, so it is not worth the\n");
                    printf("additional CPU time spent jumbling the order.\n");
                    printf("\n");
                    printf("If an additional non-option integer is given, it is used as the PRNG seed.\n");
                    printf("If no seed is given, or is 0, seed is set to current system time.\n");
                    printf("\n");
                    printf("  Metric options:\n");
                    printf("    --absolute         use |r| when calculating Pearson and cosine metrics\n");
                    printf("    --cityblock        Manhattan distance, divided by n (MAD instead of RM#D)\n");
                    printf("    --cosine           sqrt((1 - r)/2), satisfies triangle inequality\n");
                    printf("    --euclidian        Euclidian distance\n");
                    printf("    --geomean          geometric mean of RM#D and correlation\n");
                    printf("                       --cosine (default) or --pearson select correlation metric\n");
                    printf("                       --rmsd --euclidian --cityblock --minkowski options honored\n");
                    printf("    --heuristic        geometric mean of RM#D and (1 + absolute correlation)\n");
                    printf("                       --cosine (default) or --pearson select correlation metric\n");
                    printf("                       --rmsd --euclidian --cityblock --minkowski options honored\n");
                    printf("                       recommended: --heuristic --cosine --unit-variance\n");
                    printf("    --minkowski=p      Minkowski distances (1: City Block, 2: Euclidian, etc.)\n");
                    printf("    --mostly-trig      if only pearson/cosine, use 999:1 Trig:Euclidian\n");
                    printf("    --pearson          sqrt((1 - r)/2), satisfies triangle inequality\n");
                    printf("    --rmsd             root mean squared deviation (default)\n");
                    printf("\n");
                    printf("    --signed           keep sign of Pearson and cosine metrics (default)\n");
                    printf("    --similarity       divide distance by sqrt(Ochiai) present/absent similarity\n");
                    printf("                       most useful with --ignore-weak\n");
                    printf("    --spatial          factor local x#y# distances into distance calculations\n");
                    printf("    --weight-file=file weights for rows during distance calculations\n");
                    printf("\n");
                    printf("  Pre-processing options:\n");
                    printf("    --depower          de-exponentiate tree dists generated with --distpow\n");
                    printf("    --distpow=#        raise all distances to # power\n");
                    printf("    --floor-lod-ub     floor with lod upper bound instead of lower bound\n");
                    printf("    --floor-to-lod     floor values/missing to limit of detection after norm\n");
                    printf("    --floor-to-one     floor missing and values < 1 to 1\n");
                    printf("    --floor-value=X    floor missing and values < X to X\n");
                    printf("    --ignore-weak      ignore  values < 1E-5\n");
                    printf("    --impute-col-geo   impute 0/missing w/ half geomean of col and global min\n");
                    printf("    --impute-col-min   impute 0/missing w/ half column minimum\n");
                    printf("    --impute-global    impute 0/missing w/ half global minimum\n");
                    printf("    --impute-later     impute/floor after norm/filtering, instead of before\n");
                    printf("    --include-weak     include values < 1E-5 (default)\n");
                    printf("    --log2             log2 transform input prior to \"most\" calculations\n");
                    printf("    --mean-center      center vectors on their mean\n");
                    printf("    --norm-median      scale each column to have the same median\n");
                    printf("    --scale-10000      scale each row to max=10000\n");
                    printf("    --transpose-first  transpose input data before anything else\n");
                    printf("    --transpose-last   transpose input data after all other pre-processing\n");
                    printf("    --unaltered        do not transform vectors prior to calculations (default)\n");
                    printf("    --unit-max-mag     max magnitude --> 1, does not force mean-centered\n");
                    printf("    --unit-stretch=N   global min --> 0, row fraction Nth value --> 1\n");
                    printf("    --unit-variance    mean-center and scale vectors to unit variance\n");
                    printf("    --unlog2           un-log2 transform input after all pre-processing\n");
                    printf("\n");
                    printf("  Filtering options:  (before mean-center/uv, median norm, etc.)\n");
                    printf("    --filter-log2-mean=N     remove rows w/ mean  log2 values < N\n");
                    printf("    --filter-log2-sd=N       remove rows w/ stdev log2 values < N\n");
                    printf("    --filter-present=N,M     remove rows w/ fraction N above signal M\n");
                    printf("                             (,M optional; M is zero if not specified)\n");
                    printf("    --filter-unlog-mean=N    remove rows w/ mean  values      < N\n");
                    printf("    --filter-unlog-sd=N      remove rows w/ stdev values      < N\n");
                    printf("    --filter-unlog-max=N     remove rows w/ max   value       < N\n");
                    printf("\n");
                    printf("  Tree building options:\n");
                    printf("    --slink            single linkage\n");
                    printf("    --clink            complete linkage\n");
                    printf("    --upgma            average linkage (default)\n");
                    printf("    --wpgma            McQuitty's weighted average linkage\n");
                    printf("\n");
                    printf("    --upgmc2           centroid linkage,               for unsquared distances\n");
                    printf("    --ward2            Ward's minimum variance method, for unsquared distances\n");
                    printf("    --wpgmc2           Gower's median linkage,         for unsquared distances\n");
                    printf("\n");
                    printf("    --upgmcu           centroid linkage, user must pre-exponentiate distances\n");
                    printf("    --wardu            Ward's linkage,   user must pre-exponentiate distances\n");
                    printf("    --wpgmcu           median linkage,   user must pre-exponentiate distances\n");
                    printf("\n");
                    printf("    --ties-order       break ties using only nearest input/jumbled order\n");
                    printf("    --ties-fewest      break ties first using fewest descendants (default)\n");
                    printf("    --ties-random      break ties randomly\n");
                    printf("\n");
                    printf("  Tree clustering and STDOUT output options:\n");
                    printf("    --cluster-merge    merge and penalize overly small clusters (default)\n");
                    printf("    --cluster-no-merge do not merge or penalize overly small clusters\n");
                    printf("    --clusters         print clusters to STDOUT instead of tree(s)\n");
                    printf("    --leaf-custom      print leaf weights calculated from tree, no defaults\n");
                    printf("    --leaf-weights     print leaf weights calculated from tree, force defaults:\n");
                    printf("                         (--absolute --cosine --mean-center --upgma)\n");
                    printf("    --multi-line       output one line for each leaf in the tree (default)\n");
                    printf("    --nclusters=N      target N clusters for clustering, instead of auto\n");
                    printf("    --single-line      output a single long line per tree\n");
                    printf("    --tree-flip-size   reorder larger child branches first\n");
                    printf("    --tree-flip-edge   reorder shorter parent->child edges first\n");
                    printf("\n");
                    printf("  General input/output and threading options:\n");
                    printf("    --output-data      print (transformed) input data matrix to STDOUT\n");
                    printf("    --read-dists       read distance matrix from file/stdin\n");
                    printf("    --read-tree        read tree from file/stdin\n");
                    printf("    --threads=N        use N threads for all vs. all calculations (default: 1)\n");

                    free_filenames(filename_array, num_files);
                    exit(1);
                }
            }
            else
            {
                /* store first filename, for tree dist matrix */
                /* initialize filename_array */
                if (filename == NULL)
                {
                    filename = argv[i];
                    filename_array[num_files++] = strdup(argv[i]);
                }
                else
                {
                    /* TODO -- fix lazy tree jumble/seed parsing */
                    if (is_all_digits(argv[i]) == 0)
                    {
                        filename_array[num_files++] = strdup(argv[i]);
                    }
                }
            }
        }
    }
    

    /* sanity check on number of threads */
    if (opt.threads <= 0)
    {
        opt.threads = 1;
    }

    
    /* fill in missing mixed Euclidian / trig options */
    /* default to cosine, if no trig option specified */
    if (opt.heuristic_flag || opt.mostly_trig_flag ||
        opt.geomean_flag)
    {
        opt.calc_euclidian_flag = 1;
    
        if (opt.cosine_flag == 0 && opt.pearson_flag == 0)
        {
            opt.cosine_flag = 1;
        }
    }
    /* we might use the Euclidian --> cosine optimization */
    else if (opt.mean_center_flag && opt.unit_variance_flag &&
             (opt.cosine_flag || opt.pearson_flag) &&
             !(num_files > 1 || opt.weight_filename ||
               opt.transpose_last_flag))
    {
        opt.calc_euclidian_flag = 1;
    }
    /* purely pearson/cosine, disable Euclidian/RMSD calculations */
    else if (opt.cosine_flag || opt.pearson_flag)
    {
        opt.calc_euclidian_flag = 0;
    }
    
    /* specifying a linkage method triggers tree building */
    if (opt.linkage_method)
        opt.build_tree_flag = 1;

    /* build tree for tree weights if we haven't read one in */
    if (opt.tree_weights_flag)
        opt.build_tree_flag = 1;

    
    /* sanity check, disable tree building if read from file,
     * or if no tree need be built
     */
    if (opt.tree_from_file_flag || opt.dump_data_matrix_flag)
        opt.build_tree_flag = 0;


    
    /* fill empty filename with stdin character - */
    if (num_files == 0)
    {
        filename_array              = (char **) realloc(filename_array,
                                                        1 * sizeof(char *));
        filename_array[num_files++] = strdup("-");
    }

    if (opt.unit_variance_flag)
        opt.mean_center_flag = 1;

    /* set default heuristic similarity metric, if none already specified */
    if (opt.heuristic_flag && !opt.set_any_distance_flag)
    {
        opt.cosine_flag  = 1;
        opt.pearson_flag = 0;
    }
    /* set default geomean similarity metric, if none already specified */
    if (opt.geomean_flag && opt.pearson_flag == 0 && opt.cosine_flag == 0)
    {
        opt.cosine_flag  = 1;
        opt.pearson_flag = 0;
    }
    
    /* set defaults for tree weights absolute value calculations */
    if (opt.tree_defaults_flag)
    {
        opt.cosine_flag      = 1;
        opt.mean_center_flag = 1;
        opt.absolute_flag    = 1;
        opt.linkage_method   = UPGMA;
    }

    /* row weights can't currently be used with multiple input files
     * figure out what to do, and warn, about this later
     */
    if (opt.weight_filename && num_files > 1)
        opt.weight_filename = NULL;

    /* HACK -- sort of size first, before sorting on edge length */
    if (opt.tree_flip_edge_flag)
    {
        /* flip_nodes_multi() is coded to sort size before sorting edges;
         * enabling size effectively breaks edge length ties using size
         */
        opt.tree_flip_size_flag = 1;
    }

    
    if (opt.dists_from_file_flag == 0 && opt.tree_from_file_flag == 0)
    {
      /* allocate file row count array */
      file_col_counts = (int32_t *) calloc(num_files, sizeof(int32_t));

      for (i = 0; i < num_files; i++)
      {
        data_matrix = read_data_matrix(filename_array[i], data_matrix,
                                       &row_name_array, &col_name_array,
                                       &num_rows, &num_cols, &num_cols_chunk);

        if (data_matrix == NULL)
           exit(1);

        /* transpose the matrix before doing anything else */
        if (data_matrix && opt.transpose_first_flag)
        {
            data_matrix = transpose_matrix(data_matrix,
                                       &row_name_array, &col_name_array,
                                       &num_rows, &num_cols, &num_cols_chunk);
        }

        file_col_counts[i] = num_cols_chunk;
      }
      

      /* weight columns */
      if (num_files > 1)
      {
        weights = (double *) calloc(num_cols, sizeof(double));;
        dptr    = weights;

        /* find largest dataset */
        max_num_cols_chunk = file_col_counts[0];
        for (i = 1; i < num_files; i++)
            if (file_col_counts[i] > max_num_cols_chunk)
                max_num_cols_chunk = file_col_counts[i];

        for (i = 0; i < num_files; i++)
        {
            /* up-weight under-represented datasets */
            weight = (double) max_num_cols_chunk / file_col_counts[i];

            for (col = 0; col < num_cols_chunk; col++)
                *dptr++ = weight;
        }
      }
      

      /* flag missing data, convert DBL_MAX to MISSING */
      flag_missing(data_matrix, num_rows, num_cols, opt.ignore_weak_flag);


      /* estimate lod prior to filtering, imputing, etc. */
      lod = 0.0;
      if (opt.floor_to_lod_flag)
          lod = estimate_lod(data_matrix, num_rows, num_cols,
                             opt.floor_to_lod_flag);


      /* filter by missing data prior to imputing */
      if (opt.filter_present_cutoff)
      {
        data_matrix = filter_rows_by_present(data_matrix,
                                            &row_name_array,
                                            &num_rows, num_cols, &opt);

        fprintf(stderr, "NumFilteredRows:\t%d\n", (int) num_rows);
      }


      /* impute before filtering */
      if (opt.impute_later_flag == 0)
      {
          if (opt.impute_col_min_half_flag)
          {
              impute_col_min_half(data_matrix, num_rows, num_cols);
          }
          if (opt.impute_global_min_half_flag)
          {
              impute_global_min_half(data_matrix, num_rows, num_cols);
          }
          if (opt.impute_col_global_min_half_flag)
          {
              impute_col_global_min_half(data_matrix, num_rows, num_cols);
          }
          if (opt.impute_row_col_global_min_half_flag)
          {
              impute_row_col(data_matrix, num_rows, num_cols);
          }
          if (opt.floor_to_one_flag)
          {
              floor_to_one(data_matrix, num_rows, num_cols);
          }
          if (opt.floor_to_value != 0.0)
          {
              floor_to_value(data_matrix, num_rows, num_cols,
                             opt.floor_to_value);
          }
      }


      /* filter before normalization */
      if (opt.filter_log2_flag || opt.filter_unlog_flag)
      {
        data_matrix = filter_rows_by_stats(data_matrix,
                                       &row_name_array, &num_rows, num_cols,
                                       &opt);

        fprintf(stderr, "NumFilteredRows:\t%d\n", (int) num_rows);
      }

      if (opt.norm_median_flag)
        normalize_cols_median(data_matrix, num_rows, num_cols);

      if (opt.scale_ten_thousand_flag)
        scale_ten_thousand(data_matrix, num_rows, num_cols);
      
      if (opt.maxmag_flag)
        scale_unit_max_magnitude(data_matrix, num_rows, num_cols);

      if (opt.unit_stretch)
        stretch_unit_frac(data_matrix, num_rows, num_cols, opt.unit_stretch);


      /* impute/floor after filtering, normalization, etc. */
      if (opt.impute_later_flag)
      {
          if (opt.impute_col_min_half_flag)
          {
              impute_col_min_half(data_matrix, num_rows, num_cols);
          }
          if (opt.impute_global_min_half_flag)
          {
              impute_global_min_half(data_matrix, num_rows, num_cols);
          }
          if (opt.impute_col_global_min_half_flag)
          {
              impute_col_global_min_half(data_matrix, num_rows, num_cols);
          }
          if (opt.impute_row_col_global_min_half_flag)
          {
              impute_row_col(data_matrix, num_rows, num_cols);
          }
          if (opt.floor_to_one_flag)
          {
              floor_to_one(data_matrix, num_rows, num_cols);
          }
          if (opt.floor_to_value != 0.0)
          {
              floor_to_value(data_matrix, num_rows, num_cols,
                             opt.floor_to_value);
          }
      }
      
      
      if (opt.floor_to_lod_flag)
          floor_to_value(data_matrix, num_rows, num_cols, lod);


      /* log2 after normalization, since normalization expects unlogged */
      
      /* after normalization */
      if (opt.unit_variance_flag)
      {
        if (opt.log2_prior_flag)
        {
          log2_data(data_matrix, num_rows, num_cols);
          took_log2_flag = 1;
        }

        mean_center_uv(data_matrix, num_rows, num_cols);
      }
      /* after normalization */
      else if (opt.mean_center_flag)
      {
        if (opt.log2_prior_flag)
        {
          log2_data(data_matrix, num_rows, num_cols);
          took_log2_flag = 1;
        }

        mean_center(data_matrix, num_rows, num_cols);
      }
      
      /* log2 the data */
      else if (opt.log2_prior_flag && opt.unlog2_last_flag == 0)
      {
        log2_data(data_matrix, num_rows, num_cols);
        took_log2_flag = 1;
      }


      /* un-log2 after everything else */
      if (opt.unlog2_last_flag &&
          (opt.log2_prior_flag == 0 || took_log2_flag))
      {
        unlog2_data(data_matrix, num_rows, num_cols);
      }


      /* transpose data after pre-processing */
      if (opt.transpose_last_flag)
      {
        data_matrix = transpose_matrix(data_matrix,
                                       &row_name_array, &col_name_array,
                                       &num_rows, &num_cols, &num_cols_chunk);
      }


      /* weight cols from file
       * must come after we're done transposing
       */
      if (opt.weight_filename)
      {
          weights = (double *) read_weights_for_data(opt.weight_filename,
                                                     weights,
                                                     col_name_array,
                                                     num_cols);
          if (weights)
          {
              fprintf(stderr,
                      "finished reading column weights from file: %s\n",
                      opt.weight_filename);
          }
          else
          {
              fprintf(stderr, "no valid column weights found in file: %s\n",
                  opt.weight_filename);
          }
      }


      /* names must contain x#y# format coordinates; y#x# should also work */
      if (opt.spatial_flag && row_name_array)
      {
        /* allocate x,y coordinate arrays */
        x_coord_array = (int32_t *) malloc(num_rows * sizeof(int32_t));
        y_coord_array = (int32_t *) malloc(num_rows * sizeof(int32_t));
      
        for (i = 0; i < num_rows; i++)
        {
          x = y = -1;

          str_ptr = strstr(row_name_array[i], "x");
          if (str_ptr && strlen(str_ptr) > 1)
          {
            x = atol(str_ptr + 1);
          }

          str_ptr = strstr(row_name_array[i], "y");
          if (str_ptr && strlen(str_ptr) > 1)
          {
            y = atol(str_ptr + 1);
          }
          
          x_coord_array[i] = x;
          y_coord_array[i] = y;
        }
      }


      if (opt.dump_data_matrix_flag)
      {
        /* NOTE -- data matrix is freed within the funciton */
        print_data_matrix(data_matrix, row_name_array, col_name_array,
                          num_rows, num_cols);
        
        data_matrix = NULL;
      }

      /* free col names */
      if (col_name_array)
      {
        for (i = 0; i < num_cols; i++)
            if (col_name_array[i])
                free(col_name_array[i]);

        free(col_name_array);
      }
      col_name_array = NULL;

      if (opt.dump_data_matrix_flag == 0)
      {
          /* allocate distance matrix pointers */
          if (num_rows && num_cols)
            dist_matrix = (double **) calloc(num_rows, sizeof(double *));

          /* NOTE -- data matrix is freed within the funciton */
          fill_dist_matrix(data_matrix, dist_matrix, weights,
                           num_rows, num_cols,
                           x_coord_array, y_coord_array,
                           &opt);
      }

      /* data matrix was already freed earlier */
      data_matrix = NULL;
    }
    
    
    /* exponentiate distances */
    if (dist_matrix && opt.dist_pow)
        raise_dist_matrix_to_power(dist_matrix, num_rows, num_cols,
                                   opt.dist_pow);
    

    if (opt.tree_from_file_flag)
    {
        tree_string = read_in_tree_string(filename_array[0]);

        if (tree_string == NULL)
            exit(1);

        tree_root   = create_tree_from_string(tree_string,
                                              &node_ptr_array, &num_nodes,
                                              &leaf_ptr_array, &num_leaves);
        
        if (tree_string)
        {
            free(tree_string);
            tree_string = NULL;
        }

#if 0
        /* DISABLE_ME -- just for testing */
        opt.tree_flip_size_flag = 0;
        opt.tree_flip_edge_flag = 0;
        opt.tree_flip_avg_flag  = 1;
        for (i = 0; i < num_nodes; i++)
        {
            node_ptr_array[i]->edge_length =
                pow(node_ptr_array[i]->edge_length, 1.0 / 1.5);
        }
#endif
        
        bless_tree(tree_root);
        calc_node_clevels(leaf_ptr_array, num_leaves);

#if 0
        /* DISABLE_ME -- just for testing */
        calc_node_avg_lengths(node_ptr_array, num_nodes,
                              leaf_ptr_array, num_leaves);
#endif

        /* de-exponentiate if requested */
        if (opt.dist_pow && opt.depower_flag)
        {
            depower_tree(node_ptr_array, num_nodes, &opt);
        }

        /* flip nodes to reorder on size or edge length */
        if (opt.tree_flip_size_flag || opt.tree_flip_edge_flag ||
            opt.tree_flip_avg_flag)
        {
            flip_nodes_multi(node_ptr_array, num_nodes, &opt);
        }

        if (opt.tree_weights_flag)
        {
#if WEIGHT_TYPE == 1
            eff_count = calc_tree_weights_gsc(node_ptr_array, num_nodes,
                                      leaf_ptr_array, num_leaves, 1);
#elif WEIGHT_TYPE == 2
            eff_count = calc_tree_weights_experimental(node_ptr_array,
                                      num_nodes,
                                      leaf_ptr_array, num_leaves, 1);
#else
            eff_count = calc_tree_weights_clustalw(node_ptr_array, num_nodes,
                                      leaf_ptr_array, num_leaves, 1);
#endif

            fprintf(stderr, "Effective count:\t%0.8g\n", eff_count);
        }
        else
        {
            /* WARNING -- make sure we're passing same cutoffs as we
             *            do in tree.c
             */
            scan_clusters(tree_root, node_ptr_array, num_nodes,
                          0.025, 5, 1.0, 1.0, opt.target_n_clusters,
                          &n_clusters, &opt, 1);

            /* output clusters */
            if (opt.clusters_flag)
            {
                output_clusters_multi(tree_root, 1, NULL);
            }
            /* output tree */
            else
            {
                /* re-create output string, with cluster nodes flagged */
                if (newick_str) free(newick_str);

                newick_str = create_newick_string_multi(tree_root,
                                                        opt.one_line_flag, 0);

                if (newick_str)
                    printf("%s", newick_str);
            }

            /* print effective row count */
#if WEIGHT_TYPE == 1
            fprintf(stderr, "Effective count:\t%0.8g\n",
                    calc_tree_weights_gsc(node_ptr_array, num_nodes,
                                          leaf_ptr_array, num_leaves, 0));
#elif WEIGHT_TYPE == 2
            fprintf(stderr, "Effective count:\t%0.8g\n",
                    calc_tree_weights_experimental(node_ptr_array, num_nodes,
                                          leaf_ptr_array, num_leaves, 0));
#else
            fprintf(stderr, "Effective count:\t%0.8g\n",
                    calc_tree_weights_clustalw(node_ptr_array, num_nodes,
                                          leaf_ptr_array, num_leaves, 0));
#endif

            /* don't free newick string before this without setting it to NULL,
             * otherwise, free_tree_stuff() frees it again
             */
            free_tree_stuff(newick_str, node_ptr_array, num_nodes,
                            leaf_ptr_array, num_leaves, 1);
        }
    }


    /* free some arrays prior to tree building to save memory */
    free_filenames(filename_array, num_files);
    filename_array = NULL;
    
    if (file_col_counts)
        free(file_col_counts);
    
    if (weights)
        free(weights);


    /* build tree */
    if (opt.build_tree_flag)
    {
        tree_return_value = launch_tree(dist_matrix, row_name_array, num_rows,
                                        argc, argv, &opt);

        /* some arrays were already freed during tree building */
        dist_matrix    = NULL;
        row_name_array = NULL;
    }
    /* print distance matrix */
    else if (opt.tree_from_file_flag == 0)
    {
        /* NOTE -- distance matrix is freed within the function */
        if (opt.dump_data_matrix_flag == 0)
        {
            print_dist_matrix(dist_matrix, row_name_array, num_rows);

            dist_matrix = NULL;
        }
    }


    /* free dist matrix */
    if (dist_matrix)
    {
        for (i = 0; i < num_rows; i++)
            if (dist_matrix[i])
                free(dist_matrix[i]);

        free(dist_matrix);
    }

    /* free name arrays */
    if (row_name_array)
    {
        for (i = 0; i < num_rows; i++)
            if (row_name_array[i])
                free(row_name_array[i]);

        free(row_name_array);
    }

    /* free data matrix */
    if (data_matrix)
    {
        for (i = 0; i < num_rows; i++)
            if (data_matrix[i])
                free(data_matrix[i]);

        free(data_matrix);
    }

    /* free col names */
    if (col_name_array)
    {
        for (i = 0; i < num_cols; i++)
            if (col_name_array[i])
                free(col_name_array[i]);

        free(col_name_array);
    }
    
    /* free x,y coordinates */
    if (x_coord_array)
        free(x_coord_array);
    if (y_coord_array)
        free(y_coord_array);


    if (opt.build_tree_flag && tree_return_value < 0)
        return tree_return_value;
    
    pthread_exit(NULL);
}
