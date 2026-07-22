#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

extern char * fgets_strip_realloc(char **return_string,
                                  int32_t *return_max_length,
                                  FILE *infile);

extern int32_t split_tabs(char *string, char ***fields,
                          int32_t *return_max_field);

extern int is_all_digits(char *string);

extern double ** read_data_matrix(char     *filename,
                                  double  **data_matrix,
                                  char   ***return_row_name_array,
                                  char   ***return_col_name_array,
                                  int32_t  *return_num_rows,
                                  int32_t  *return_num_cols,
                                  int32_t  *return_num_cols_chunk);
extern double ** transpose_matrix(double  **data_matrix,
                                  char   ***return_row_name_array,
                                  char   ***return_col_name_array,
                                  int32_t  *return_num_rows,
                                  int32_t  *return_num_cols,
                                  int32_t  *return_num_cols_chunk);

extern int32_t cmp_name_index_by_name(const void *keyval, const void *datum);
char * seconds_to_str(char *string, int32_t buf_len, double seconds);
