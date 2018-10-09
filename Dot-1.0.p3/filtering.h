#ifndef FILTERING_H_
#define FILTERING_H_

#include "common.h"
#include "algorithm.h"
#include "configReader.h"

typedef struct _tandemR_list_elem TRs_list;

/*  typedef enum {NONE, THRESHOLD, LIGHT, FAIR, HEAVY} strategy_t;  Declared in common.h */

struct filter_object {
  unsigned int clunk_max_TR_len;    /* Update for each clunk with the longer */
  float clunk_purity; /* Update for each clunk with the purest */      
  unsigned int clunk_min_mot_len;   /* Update for each clunk with the lower motif length */
  unsigned int clunk_max_mot_len;   /* Update for each clunk with the higher motif length */
  unsigned long int clunk_init;  /* Array position of the first TR in the clunk */
  unsigned long int clunk_end;  /*  Array position of the last TR the clunk */
  unsigned int rightmost;
  _Bool finish;
  TRs_result_t *trlist;  /*  Pointer to the vector with results  */
  unsigned long int tot_trs;
};

struct clunk_stats {
  unsigned int clunk_max_TR_len;
  float clunk_purity;
  unsigned int clunk_min_mot_len;
  unsigned int clunk_max_mot_len;
  unsigned long int clunk_init;
  unsigned long int clunk_end;
};

struct filter_object *init_filter_obj (TRs_result_t *trlist, unsigned long int tot_trs);
struct clunk_stats *getNextContiguousSetOfTRs (struct filter_object *fObj, struct clunk_stats *);
void applyFilter (struct filter_object *fObj, struct clunk_stats *stats, struct config *params);
void filter_lightFilter (struct filter_object *fObj, struct clunk_stats *stats);
void filter_fairFilter (struct filter_object *fObj, struct clunk_stats *stats);
void filter_heavyFilter (struct filter_object *fObj, struct clunk_stats *stats, struct config *params);
void filter_overlapReducerFilter (struct filter_object *fObj, struct clunk_stats *stats, struct config *params);
void filter_init_trlist (TRs_result_t *trlist, unsigned long int tot_trs);
void filter_thresholdFilter (TRs_result_t *trlist,  unsigned long int tot_trs, struct config *params);
void filter (TRs_Result_Bundle *results, struct config *params);


#endif /* FILTERING_H_ */
