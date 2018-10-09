#ifndef _CONFIG_READER_H_
#define _CONFIG_READER_H_

#include "common.h"
#include "utils.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>

#define MAX_FILELINE 256
/*  #define MAX_PARAM_LEN 1024  defined in common.h  */

#define FILT_FLAG 32768
#define TOLL_FLAG 16384
#define MTRL_FLAG 8192
#define MPUR_FLAG 4096
#define OVER_FLAG 2048
#define J_FLAG 1024
#define S_FLAG 512
#define C_FLAG 256
#define N_FLAG 128
#define X_FLAG 64
#define F_FLAG 32
#define O_FLAG 16
#define U_FLAG 8
#define I_FLAG 4
#define G_FLAG 2
#define T_FLAG 1


struct config {
  unsigned short int flags;
  int nvalue; /* min window length  */
  int xvalue; /* max window length  */
  int jvalue; /* jumps allowed  */
  int gvalue; /* gaps allowed  */
  int thread_max;
  char svalue[MAX_PARAM_LEN]; /**/
  char cvalue[MAX_PARAM_LEN];  /**/
  char output_filename[MAX_PARAM_LEN];
  MATCH_ARRAY_TYPE fvalue;
  strategy_t filter_type;
  float tollerance;
  _Bool allow_overlap;
  unsigned int min_TR_len; 
  float min_purity; 
  _Bool verbose;
};

struct paramList {
  char param[MAX_FILELINE];
  char value[MAX_FILELINE];
  struct paramList *next;
};

typedef struct paramList* paramList_pointer;


struct config *command_line_parser (int argc, char* argv[]);
void param_list_parser (struct paramList *par_list, struct config *cfg);

struct paramList *loadConfigFromFile (char *filename);
void free_paramList (struct paramList* pl);
MATCH_ARRAY_TYPE **read_weights_matrix (struct paramList *param);
void free_weights_matrix (MATCH_ARRAY_TYPE **m);

#endif /* _CONFIG_READER_H_ */
