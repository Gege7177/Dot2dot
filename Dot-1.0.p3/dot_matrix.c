#include "dot_matrix.h"

void dot_free (struct dot_matrix* dm) {
  int i;
  free (dm->sequence);
  free (dm->pointers_sequence);
  free (dm->mask);
  for (i = 0; i < MAXIND_ARRAY; i ++) {
    if (dm->match_arrays[i] != NULL) 
      free (dm->match_arrays[i]);
  }
  free (dm);
}

struct dot_matrix *__dot_create_matrix (struct sequence_t *seq) {
  struct dot_matrix *dm = NULL;
  int i;
  dm = (struct dot_matrix *) malloc (sizeof(struct dot_matrix));
  if (dm == NULL) {
    perror("Malloc Error in __dot_create_matrix ();\n");
    return NULL;
  }
  /*  XXX - Copying seems an unreasonable behaviour, but makes modules independent  */
  dm->sequence = (char *) malloc ((seq->sequence_size + 1) * sizeof (char));  /*  Including \0  */
  if (dm->sequence == NULL) {
    free (dm);
    perror("Malloc Error in __dot_create_matrix ();\n");
    return NULL;
  }
  strcpy (dm->sequence, seq->sequence);
  dm->sequence_len = seq->sequence_size;
  
  strcpy (dm->IDseq, seq->label);  /*  Sequence name  */
  for (i = 0; i < MAXIND_ARRAY; i ++) dm->match_arrays[i] = NULL;
  dm->pointers_sequence = (MATCH_ARRAY_TYPE**) malloc(dm->sequence_len * sizeof(MATCH_ARRAY_TYPE*));
  if (dm->pointers_sequence == NULL) { 
    free (dm);
    free (dm->sequence);
    perror("Malloc Error in __dot_create_matrix ();\n");
    return NULL;
  }
  dm->mask = (char *) malloc ((dm->sequence_len + 1) * sizeof (char));
  if (dm->mask == NULL) { 
    free (dm);
    free (dm->sequence);
    free (dm->pointers_sequence);
    perror("Malloc Error in __dot_create_matrix ();\n");
    return NULL;
  }
  for(i = 0; i < dm->sequence_len; i ++) {
    dm->mask[i] = UNCHECKED;
  }
  dm->mask[i]='\0';
  return dm;
}

/*  accept ch as int and not char just to avoid a compiler warning  */
MATCH_ARRAY_TYPE *__dot_create_line (struct dot_matrix *dm, MATCH_ARRAY_TYPE **wm, int ch) {
  MATCH_ARRAY_TYPE *line;
  int i;
  line = (MATCH_ARRAY_TYPE*) malloc ((dm->sequence_len + 1) * sizeof(MATCH_ARRAY_TYPE));
  if (line == NULL) return NULL;
  for (i = 0; i < dm->sequence_len; i ++)
    line[i] = wm[ch][(int) dm->sequence[i]];
  line[i] = 0;
  return line;
}

struct dot_matrix *dot_init (struct sequence_t *seq, MATCH_ARRAY_TYPE **wm) {
  int i, sequence_i;
  struct dot_matrix *dm = __dot_create_matrix (seq);
  if (dm == NULL) return NULL;
  for (i = 0; i < dm->sequence_len; i ++) {
    sequence_i = (int) dm->sequence[i];
    if (dm->match_arrays[sequence_i] == NULL) {
      dm->match_arrays[sequence_i] = __dot_create_line (dm, wm, (int) sequence_i);
      if (dm->match_arrays[sequence_i] == NULL) {
	perror("Malloc Error in __dot_create_matrix ();\n");
	dot_free (dm);
	return NULL;
      }
    }
    dm->pointers_sequence[i] = dm->match_arrays[sequence_i];
  }
  /*
  int j,k;
  for (j = 0; j < dm->sequence_len; j ++) {
    for (k = 0; k < dm->sequence_len; k ++) {
      printf ("%f ",dm->pointers_sequence[j][k]);
    }
  }
  return NULL;
  */
  return dm;
}
