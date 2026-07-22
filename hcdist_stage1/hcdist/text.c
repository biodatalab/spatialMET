#include <string.h>
#include "text.h"
#include "hcdist.h"

#define MEM_OVERHEAD 1.01    /* speed hack -- overallocate to avoid reallocs */


int32_t cmp_name_index_by_name(const void *keyval, const void *datum)
{
    struct name_index_pair *ptr1, *ptr2;
    int    compare = 0;
    
    ptr1 = (struct name_index_pair *) keyval;
    ptr2 = (struct name_index_pair *) datum;

    if (ptr1->name && ptr2->name) compare = strcmp(ptr1->name, ptr2->name);
    if (compare) return(compare);

#if 0
    /* must leave this out for use with bsearch */
    if (ptr1->index < ptr2->index) return -1;
    if (ptr1->index > ptr2->index) return  1;
#endif
    
    return 0;
}


/* realloc input string, return new max array length (including terminal null)
 * handles \r\n \n \r, including mixes of EOL characters within same file
 * strips EOL from end of string
 */
char * fgets_strip_realloc(char **return_string, int32_t *return_max_length,
                           FILE *infile)
{
    char c;
    char *string          = *return_string;
    int32_t length        = 0;
    int32_t total_length;
    int32_t max_length    = *return_max_length;
    char old_c            = '\0';
    int32_t anything_flag = 0;

    while((c = fgetc(infile)) != EOF)
    {
        anything_flag = 1;
    
        /* EOL: \n or \r\n */
        if (c == '\n')
        {
            /* MSDOS, get rid of the previously stored \r */
            if (old_c == '\r')
            {
                string[length - 1] = '\0';
            }

            old_c = c;
            
            break;
        }
        /* EOL: \r */
        /* may be a Mac text line, back up a character */
        else if (old_c == '\r')
        {
            fseek(infile, -1 * sizeof(char), SEEK_CUR);

            break;
        }
        
        old_c = c;
    
        total_length = length + 2;    /* increment, plus terminal null */

        if (total_length > max_length)
        {
            max_length = MEM_OVERHEAD * total_length;
            string     = (char *) realloc(string, max_length * sizeof(char));
        }

        string[length++] = c;
    }
    
    /* check for dangling \r from reading in Mac lines */
    if (length && string[length-1] == '\r')
    {
        string[length-1] = '\0';
    }
    
    if (length == 0)
    {
        if (!anything_flag)
        {
            return NULL;
        }
            
        if (1 > max_length)
        {
            max_length = 1;
            string     = (char *) realloc(string, sizeof(char));
        }
    }

    string[length]     = '\0';
    
    *return_string     = string;
    *return_max_length = max_length;

    return string;
}


/* replace tabs with nulls, fill array of field pointers,
 * return number of fields
 * WARNING -- clobbers tabs in original input string
 */
int32_t split_tabs(char *string, char ***fields, int32_t *return_max_field)
{
    char *cptr, *sptr;
    int32_t count     = 0;
    int32_t max_field = *return_max_field;
    
    sptr = string;
    for (cptr = string; *cptr; cptr++)
    {
        if (*cptr == '\t')
        {
            count++;

            if (count > max_field)
            {
                max_field = MEM_OVERHEAD * count;
                *fields   = realloc(*fields, max_field * sizeof(char *));
            }

            (*fields)[count-1] = sptr;
            sptr = cptr + 1;
            *cptr = '\0';
        }
    }
    
    /* final field */
    count++;
    if (count > max_field)
    {
        max_field = MEM_OVERHEAD * count;
        *fields = realloc(*fields, max_field * sizeof(char *));
    }
    (*fields)[count-1] = sptr;
    
    *return_max_field = max_field;
    
    return count;
}


int is_all_digits(char *string)
{
    int length = 0;
    int num_digits = 0;
    char *sptr;
    
    if (string)
       length = strlen(string);
    else
       return 0;
    
    sptr = string;
    while (*sptr)
    {
        if (isdigit(*sptr++))
            num_digits++;
        else
            break;
    }
    
    if (num_digits == length)
        return 1;
    
    return 0;
}


/*
 * Input matrix format is tab-delimited full distance matrix, but with
 *  labels added as the first column.
 * There can be no blank/extra rows or columns
 *
 * return pointer to distance matrix
 *
 * 2020-10-23  support concatenating multiple files by calling sequentially
 * 2020-11-10  deal with blank lines; deal with trailing blank columns
 */
double ** read_data_matrix(char     *filename,
                           double  **data_matrix,
                           char   ***return_row_name_array,
                           char   ***return_col_name_array,
                           int32_t  *return_num_rows,
                           int32_t  *return_num_cols,
                           int32_t  *return_num_cols_chunk)
{
    FILE     *infile;
    char     *string         = NULL;
    char    **fields         = NULL;
    int32_t   max_string_len = 0;
    int32_t   num_fields     = 0;
    int32_t   max_num_fields = 0;
    int32_t   max_num_rows   = *return_num_rows;
    char     *buffer         = NULL;
    char     *sptr;
    int32_t   line_num       = 0;

    double   *dptr           = NULL;
    char    **row_name_array = *return_row_name_array;
    char    **col_name_array = *return_col_name_array;
    int32_t   orig_num_rows  = *return_num_rows;
    int32_t   orig_num_cols  = *return_num_cols;
    int32_t   num_rows       = 0;
    int32_t   num_cols       = 0;
    int32_t   max_filled_col = -1;
    int32_t   temp_int32_t;
    int32_t   row;
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
        /* skip entirely blank lines */
        temp_int32_t = 0;
        for (sptr = string; *sptr; sptr++)
        {
            if (!isspace(*sptr))
            {
                temp_int32_t = 1;
                break;
            }
        }
        if (temp_int32_t == 0)
            continue;


        num_fields = split_tabs(string, &fields, &max_num_fields);

        /* header line */
        if (line_num == 0)
        {
            num_cols = num_fields - 1;
            
            if (num_cols)
                col_name_array = (char **) realloc(col_name_array,
                                     (orig_num_cols + num_cols) *
                                     sizeof(char *));
            
            /* null out newly added columns */
            for (i = 0; i < num_cols; i++)
                col_name_array[orig_num_cols + i] = NULL;
            
            /* store column names */
            for (i = 1; i < num_fields; i++)
            {
                col_name_array[orig_num_cols + i-1] = strdup(fields[i]);
            }
            
            /* scan for blank trailing columns */
            for (i = 1; i < num_fields; i++)
            {
                temp_int32_t = 0;
                for (sptr = fields[i]; *sptr; sptr++)
                {
                    if (!isspace(*sptr))
                    {
                        temp_int32_t = 1;
                        break;
                    }
                }
                if (temp_int32_t && i-1 > max_filled_col)
                    max_filled_col = i-1;
            }

            line_num++;
            
            continue;
        }
        
        row = num_rows++;

        /* allocate new row pointers if needed */
        if (num_rows > max_num_rows)
        {
            temp_int32_t = max_num_rows;
            max_num_rows = MEM_OVERHEAD * num_rows;

            /* new row pointers in data matrix */
            data_matrix    = (double **) realloc(data_matrix,
                                             max_num_rows * sizeof(double *));
            /* new row pointers in row name matrix */
            row_name_array = (char   **) realloc(row_name_array,
                                             max_num_rows * sizeof(char *));

            /* zero out new row pointers */
            memset(data_matrix    + temp_int32_t, 0,
                   (max_num_rows - temp_int32_t) * sizeof(double *));
            memset(row_name_array + temp_int32_t, 0,
                   (max_num_rows - temp_int32_t) * sizeof(char *));
        }
        
        /* reallocate row */
        data_matrix[row] = (double *) realloc(data_matrix[row],
                                              (orig_num_cols + num_cols) *
                                              sizeof(double));
        /* zero out new data */
        memset(data_matrix[row] + orig_num_cols, 0,
               num_cols * sizeof(double));
        
        /* store row name */
        if (row_name_array[row])
            free(row_name_array[row]);
        row_name_array[row] = strdup(fields[0]);
        
        /* store data in data matrix */
        dptr = data_matrix[row] + orig_num_cols;
        for (i = 1; i < num_fields; i++)
        {
            *dptr = atof(fields[i]);

    	    /* Inf and NaN */
    	    if (!isfinite(*dptr))
    	    {
                *dptr = DBL_MAX;
    	    }
    	    /* anything else that isn't [0-9.+-]
    	     * " 0" would be a false-positive, but I am assuming that
    	     *  leading/trailing spaces have been stripped from input file
    	     * 0xDEADBEEF would be a false-negative?
    	     */
    	    if (*dptr == 0)
    	    {
    	        if (!isdigit(fields[i][0]) &&
    	            fields[i][0] != '-' &&
    	            fields[i][0] != '+' &&
    	            fields[i][0] != '.')
    	        {
                    *dptr = DBL_MAX;
    	        }
    	    }
    	    
    	    if (*dptr != DBL_MAX && i-1 > max_filled_col)
    	        max_filled_col = i-1;
            
            dptr++;
        }
        
        line_num++;
    }
    fclose(infile);
    
    /* set num_cols to last /data/ column containing data */
    num_cols = 0;
    if (max_filled_col >= 0)
        num_cols = max_filled_col + 1;

    /* fprintf(stderr, "Rows: %d\tCols: %d\n", num_rows, num_cols); */
    
    /* HACK -- sanity check if we're reading in a second file */
    if (orig_num_rows && num_rows != orig_num_rows)
    {
        fprintf(stderr,
                "ABORT -- number of rows different between files: %s\n",
                filename);

        exit(2);
    }

    if (buffer) free(buffer);
    if (string) free(string);
    if (fields) free(fields);

    *return_row_name_array = row_name_array;
    *return_col_name_array = col_name_array;
    *return_num_rows       = num_rows;
    *return_num_cols       = orig_num_cols + num_cols;
    *return_num_cols_chunk = num_cols;
    
    return data_matrix;
}


double ** transpose_matrix(double  **data_matrix_orig,
                           char   ***return_row_name_array,
                           char   ***return_col_name_array,
                           int32_t  *return_num_rows,
                           int32_t  *return_num_cols,
                           int32_t  *return_num_cols_chunk)
{
    double  **data_matrix_new = NULL;
    char    **row_name_array  = *return_row_name_array;
    char    **col_name_array  = *return_col_name_array;
    int32_t   num_rows        = *return_num_rows;
    int32_t   num_cols        = *return_num_cols;
    int32_t   row, col;
    
    
    /* allocate new data matrix */
    data_matrix_new = (double **) malloc(num_cols * sizeof(double *));
    for (col = 0; col < num_cols; col++)
        data_matrix_new[col] = (double *) malloc(num_rows * sizeof(double));


    /* transpose original matrix into new matrix */
    for (row = 0; row < num_rows; row++)
        for (col = 0; col < num_cols; col++)
            data_matrix_new[col][row] = data_matrix_orig[row][col];


    /* free original data matrix */
    if (data_matrix_orig)
    {
        for (row = 0; row < num_rows; row++)
            if (data_matrix_orig[row])
                free(data_matrix_orig[row]);

        free(data_matrix_orig);
    }


    /* swap everything else */
    *return_row_name_array = col_name_array;
    *return_col_name_array = row_name_array;
    *return_num_rows       = num_cols;
    *return_num_cols       = num_rows;
    *return_num_cols_chunk = num_rows;
    
    return data_matrix_new;
}


/* input file is 2-column tab delimited, contains a header line */
double * read_weights_for_data(char     *filename,
                               double   *weight_array,
                               char    **data_names,
                               int32_t   num_points)
{
    FILE     *infile;
    char     *string         = NULL;
    char    **fields         = NULL;
    int32_t   max_string_len = 0;
    int32_t   num_fields     = 0;
    int32_t   max_num_fields = 0;
    char     *buffer         = NULL;
    char     *sptr;
    int32_t   line_num       = 0;

    int32_t   temp_int32_t;
    int32_t   i;

    struct name_index_pair *name_index_pairs = NULL;
    struct name_index_pair *pair_ptr;
    struct name_index_pair  query_pair;
    double value;


    /* allocate my own i/o buffer, since we can't trust the system and/or
     * compiler to allocate a decently large one...
     */
    
    buffer = (char *) malloc(1000000 * sizeof(char));

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


    /* initialize the weights */
    if (weight_array)
    {
        for (i = 0; i < num_points; i++)
            weight_array[i] = 1.0;
    }


    while(fgets_strip_realloc(&string, &max_string_len, infile))
    {
        /* skip entirely blank lines */
        temp_int32_t = 0;
        for (sptr = string; *sptr; sptr++)
        {
            if (!isspace(*sptr))
            {
                temp_int32_t = 1;
                break;
            }
        }
        if (temp_int32_t == 0)
            continue;


        num_fields = split_tabs(string, &fields, &max_num_fields);

        /* header line */
        if (line_num == 0)
        {
            line_num++;

            continue;
        }
        
        if (num_fields >= 2)
        {
            value = atof(fields[1]);

            /* Inf and NaN */
            if (!isfinite(value))
            {
                continue;
            }

            /* anything else that isn't [0-9.+-]
             * " 0" would be a false-positive, but I am assuming that
             *  leading/trailing spaces have been stripped from input file
             * 0xDEADBEEF would be a false-negative?
             */
            if (value == 0.0)
            {
                if (!isdigit(fields[i][0]) &&
                    fields[i][0] != '-' &&
                    fields[i][0] != '+' &&
                    fields[i][0] != '.')
                {
                    continue;
                }
            }


            /* initialize the query */
            query_pair.name = fields[0];

            /* search by name */
            pair_ptr = bsearch(&query_pair, name_index_pairs, num_points,
                               sizeof(struct name_index_pair),
                               cmp_name_index_by_name);

            /* found it */
            if (pair_ptr)
            {
                /* allocate array if needed */
                if (weight_array == NULL)
                {
                    /* new row pointers in data matrix */
                    weight_array =
                        (double *) malloc(num_points * sizeof(double));

                    /* initialize the weights */
                    for (i = 0; i < num_points; i++)
                        weight_array[i] = 1.0;
                }

                weight_array[pair_ptr->index] = value;
            }
        }
        
        line_num++;
    }
    fclose(infile);

    
    if (buffer) free(buffer);
    if (string) free(string);
    if (fields) free(fields);
    
    if (name_index_pairs)
        free(name_index_pairs);


    return weight_array;
}


/* write ddd:hh:mm:ss string
 *
 * buf_len = maximum length of resulting string, including terminal null byte
 *
 * return address of original string
 */
char * seconds_to_str(char *string, int32_t buf_len, double seconds)
{
    double days = 0.0, hours = 0.0, minutes = 0.0;
    
    /* round to nearest non-negative second */
    if (seconds < 0.0) seconds = 0.0;
    seconds = floor(seconds + 0.5);

    if (seconds >= 60)
    {
        minutes = floor(seconds / 60.0);
        seconds = seconds - 60.0 * minutes;
    }

    if (minutes >= 60)
    {
        hours   = floor(minutes / 60.0);
        minutes = minutes - 60.0 * hours;
    }

    if (hours >= 24)
    {
        days   = floor(hours / 24.0);
        hours  = hours    - 24.0 * days;
    }
    
    if (days)
    {
        snprintf(string, buf_len, "%03dd:%02dh:%02dm:%02ds",
                 (int) days, (int) hours, (int) minutes, (int) seconds);
    }
    else if (hours)
    {
        snprintf(string, buf_len, "%02dh:%02dm:%02ds",
                 (int) hours, (int) minutes, (int) seconds);
    }
    else
    {
        snprintf(string, buf_len, "%02dm:%02ds",
                 (int) minutes, (int) seconds);
    }
    
    return string;
}
