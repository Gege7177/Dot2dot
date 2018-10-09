#ifndef __COMMON__
#define __COMMON__

#include <stdbool.h>  /*  Makes _Bool and true/false available  */

#define MATCH_ARRAY_TYPE float
#define MAXIND_ARRAY (256 * sizeof (char))
#define MAX_PARAM_LEN 1024
#define MAX_LABEL_LENGTH  65


/*  Native type for filtering  */
typedef enum {NONE, THRESHOLD, LIGHT, FAIR, HEAVY} strategy_t;
/*  Native typer for file manager  */
typedef enum {UNKNOWN, FASTA, FASTQ} seqfile_t;

#endif /* __COMMON__ */
