#ifndef DOT_MATRIX_H_
#define DOT_MATRIX_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"
#include "fileManager.h"

/*** MASK values **/
#define UNCHECKED '-'
#define CHECKED 'c'

struct dot_matrix {
  long int sequence_len;
  char* sequence;
  char IDseq[MAX_LABEL_LENGTH];
  MATCH_ARRAY_TYPE *match_arrays[MAXIND_ARRAY]; /*  Not really necessary  */
  MATCH_ARRAY_TYPE **pointers_sequence;
  char* mask;
};

void dot_free (struct dot_matrix* dm);
struct dot_matrix *__dot_create_matrix (struct sequence_t *seq);
MATCH_ARRAY_TYPE *__dot_create_line (struct dot_matrix *dm, MATCH_ARRAY_TYPE **wm, int ch);
struct dot_matrix *dot_init (struct sequence_t *seq, MATCH_ARRAY_TYPE **wm);

#endif /* STRINGMANAGING_STRING_MANAGER_H_ */

