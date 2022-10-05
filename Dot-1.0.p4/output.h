#ifndef OUTPUT_H_
#define OUTPUT_H_


#include <stdio.h>
#include <string.h>
#include "algorithm.h"

typedef enum {dot, bed} output_format_t;

struct outfile {
  output_format_t filetype;
  FILE *pf;
};

struct outfile *output_create (char *filename);
void output_destroy (struct outfile *out);
void print_header (struct outfile *out);
void print_TRs_list_toFile (struct outfile *out, char* label, char* sequence, TRs_Result_Bundle* trs_bundle);
void dump_TR_elem (struct outfile *out, char* sequence, TRs_result_t* tr_result, unsigned short int* motif_lengths, char *outStr, char *tmpbuff, size_t n);

/*
 * Param: (int* array, int length)
 *
 * Print 'int_Array' made of 'length' values
 */
void print_usintArray(unsigned short int* array, int length);

#endif
