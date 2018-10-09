#include "filtering.h"

struct filter_object *init_filter_obj (TRs_result_t *trlist, unsigned long int tot_trs) {
  struct filter_object *fObj;
  unsigned long int elem_pos;
  fObj = (struct filter_object *) malloc (sizeof (struct filter_object));
  if (fObj == NULL) {
    perror ("Error allocating memory while filtering - not applyed\n");
    return NULL;
  }
  /*   Initialize fObj  */
  fObj->clunk_init = 0;
  fObj->clunk_end = 0;
  fObj->rightmost = 0;
  fObj->clunk_max_TR_len = 0;
  fObj->clunk_purity = 0;
  fObj->clunk_min_mot_len = 0;
  fObj->clunk_max_mot_len = 0;
  fObj->finish = false;
  fObj->trlist = trlist;  /*  Pointer to the vector with results  */
  fObj->tot_trs = tot_trs;   /*  Number of elements in the vector  */
  for (elem_pos = 0; elem_pos < tot_trs; elem_pos ++) {
    if (trlist[elem_pos].valid_TR == true) {
      fObj->clunk_init = elem_pos;  /*  Posizione nella lista del primo TR del clunk */
      fObj->clunk_end = elem_pos;  /*  Posizione nella lista del primo TR del clunk */
      fObj->rightmost = trlist[elem_pos].origin_position + trlist[elem_pos].full_length;  /*  The rightmost coordinate found till now */
      /*  Update statistics  */
      fObj->clunk_max_TR_len = trlist[elem_pos].full_length;
      fObj->clunk_purity = trlist[elem_pos].purity_percentage;
      fObj->clunk_min_mot_len = trlist[elem_pos].period;
      fObj->clunk_max_mot_len = trlist[elem_pos].period;
      break;
    }
  }
  if (fObj->clunk_max_TR_len == 0) /*  No valid TRs in list  */
    fObj->finish = true;  

  return fObj;
}

struct clunk_stats *getNextContiguousSetOfTRs (struct filter_object *fObj, struct clunk_stats *stats) {
  TRs_result_t *trlist = fObj->trlist; /*  Pointer to the results vector */
  unsigned long int i;
  if (fObj->finish == true) {
    if (stats != NULL) free (stats); /*  Release the memory of the stats  */
    return NULL;
  }
  if (stats == NULL) {   /*  if null allocate memory, otherwise just reuse it  */
    stats = (struct clunk_stats *) malloc (sizeof (struct clunk_stats));
    if (stats == NULL) {
      perror ("Memory allocation failure during filtering - results partially filtered\n");
      return NULL;
    }
  }
  for (i = fObj->clunk_end + 1; i < fObj->tot_trs; i ++) {
    if (trlist[i].valid_TR == false) continue;  /*  Ignore it  */

    if (trlist[i].origin_position < fObj->rightmost) { /* It is overlapping  */
      fObj->clunk_end = i; /* store position in the vector */
      if (fObj->rightmost < trlist[i].origin_position + trlist[i].full_length)
	  fObj->rightmost = trlist[i].origin_position + trlist[i].full_length;
	/*  Update statistics  */
	if (trlist[i].full_length > fObj->clunk_max_TR_len)
	  fObj->clunk_max_TR_len = trlist[i].full_length;
	if (trlist[i].purity_percentage > fObj->clunk_purity) 
	  fObj->clunk_purity = trlist[i].purity_percentage;
	if (trlist[i].period < fObj->clunk_min_mot_len)
	  fObj->clunk_min_mot_len = trlist[i].period;
	if (trlist[i].period > fObj->clunk_max_mot_len)
	  fObj->clunk_max_mot_len = trlist[i].period;
    }
    else {   /*  This is a new clunk  */
      stats->clunk_init = fObj->clunk_init;
      stats->clunk_end = fObj->clunk_end;
      stats->clunk_max_TR_len = fObj->clunk_max_TR_len;
      stats->clunk_purity = fObj->clunk_purity;
      stats->clunk_min_mot_len = fObj->clunk_min_mot_len;
      stats->clunk_max_mot_len = fObj->clunk_max_mot_len;

      fObj->clunk_init = i;  /*  Posizione nella lista del primo TR del clunk */
      fObj->clunk_end = i;  /*  Posizione nella lista del primo TR del clunk */
      fObj->rightmost = trlist[i].origin_position + trlist[i].full_length;  /*  The rightmost coordinate found till now */
      /*  Update statistics  */
      fObj->clunk_max_TR_len = trlist[i].full_length;
      fObj->clunk_purity = trlist[i].purity_percentage;
      fObj->clunk_min_mot_len = trlist[i].period;
      fObj->clunk_max_mot_len = trlist[i].period;
      return stats;
    }
  }
  fObj->finish = true;
  stats->clunk_init = fObj->clunk_init;
  stats->clunk_end = fObj->clunk_end;
  stats->clunk_max_TR_len = fObj->clunk_max_TR_len;
  stats->clunk_purity = fObj->clunk_purity;
  stats->clunk_min_mot_len = fObj->clunk_min_mot_len;
  stats->clunk_max_mot_len = fObj->clunk_max_mot_len;
  return stats;
}

void applyFilter (struct filter_object *fObj, struct clunk_stats *stats, struct config *params) {
#if DEBUG
  unsigned long int i;
  TRs_result_t *trlist = fObj->trlist; /*  Pointer to the results vector */
#endif
  if (params->filter_type == LIGHT) filter_lightFilter (fObj, stats);
  if (params->filter_type == FAIR) filter_fairFilter (fObj, stats);
  if (params->filter_type == HEAVY) filter_heavyFilter (fObj, stats, params);
  if (params->allow_overlap == false) filter_overlapReducerFilter (fObj, stats, params);
#if DEBUG
  printf ("maxTR %u pur %.3f minMot %u maxMot %u\n",
	  stats->clunk_max_TR_len, stats->clunk_purity,
	  stats->clunk_min_mot_len, stats->clunk_max_mot_len);
  for (i = stats->clunk_init; i <= stats->clunk_end; i ++) {
    if (trlist[i].valid_TR == true) 
      printf ("%ld\t%ld\t%d\t%d\t%.3f\t%.3f\tok\n", trlist[i].origin_position + 1,
	      trlist[i].origin_position + trlist[i].full_length,
	      trlist[i].copy_number,trlist[i].period, trlist[i].purity_percentage,
	      trlist[i].stats);
    else  
      printf ("%ld\t%ld\t%d\t%d\t%.3f\t%.3f\tno\n", trlist[i].origin_position + 1,
	      trlist[i].origin_position + trlist[i].full_length,
	      trlist[i].copy_number,trlist[i].period, trlist[i].purity_percentage,
	      trlist[i].stats);
  }
  fflush (stdout);
#endif
}

void filter_lightFilter (struct filter_object *fObj, struct clunk_stats *stats) {
  TRs_result_t *trlist = fObj->trlist; /*  Pointer to the results vector */
  unsigned long int i;
  for (i = stats->clunk_init; i <= stats->clunk_end; i ++) {
    if ((trlist[i].full_length != stats->clunk_max_TR_len) &&
	(trlist[i].purity_percentage != stats->clunk_purity) &&
	(trlist[i].period != stats->clunk_min_mot_len))                               
      trlist[i].valid_TR = false;
  }
}

void filter_fairFilter (struct filter_object *fObj, struct clunk_stats *stats) {
  TRs_result_t *trlist = fObj->trlist; /*  Pointer to the results vector */
  unsigned long int i;
  unsigned char maxval = 0;  /*  char perche' sono tirchio  */
  /*   Collect statistics  */
  for (i = stats->clunk_init; i <= stats->clunk_end; i ++) {
    trlist[i].stats = 0;
    if (trlist[i].valid_TR == false) continue;

    if (trlist[i].full_length == stats->clunk_max_TR_len) trlist[i].stats ++;
    if (trlist[i].purity_percentage == stats->clunk_purity) trlist[i].stats ++;
    if (trlist[i].period == stats->clunk_min_mot_len) trlist[i].stats ++;
    if (trlist[i].stats > maxval) maxval = trlist[i].stats;
  }
  /*   filter based on statistics  */
  for (i = stats->clunk_init; i <= stats->clunk_end; i ++) {
    if (trlist[i].stats < maxval) trlist[i].valid_TR = false;
  }
}

void filter_heavyFilter (struct filter_object *fObj, struct clunk_stats *stats, struct config *params) {
  TRs_result_t *trlist = fObj->trlist; /*  Pointer to the results vector */
  unsigned long int i;
  float maxval = 0;  
  /*   Collect statistics  */
  for (i = stats->clunk_init; i <= stats->clunk_end; i ++) {
    trlist[i].stats = 0;
    if (trlist[i].valid_TR == false) continue;

    trlist[i].stats = (float) trlist[i].full_length / (float) stats->clunk_max_TR_len;
    trlist[i].stats += trlist[i].purity_percentage;
    /*trlist[i].stats += (1 - ((float) trlist[i].period / (float) stats->clunk_max_mot_len));*/
    if (trlist[i].stats > maxval) maxval = trlist[i].stats;
  }

  /*   filter based on statistics  */
  for (i = stats->clunk_init; i <= stats->clunk_end; i ++) {
    if (trlist[i].stats < maxval - params->tollerance) trlist[i].valid_TR = false;
  }
}

void filter_overlapReducerFilter (struct filter_object *fObj, struct clunk_stats *stats, struct config *params) {
  TRs_result_t *trlist = fObj->trlist; /*  Pointer to the results vector */
  unsigned long int i, next;
  i = stats->clunk_init;
  next = stats->clunk_init;

  /*   XXX  - keep this as a while and not for  */
  while (i < stats->clunk_end) {   /*  Last element is igneored because no one is ahead  */
    if (trlist[i].valid_TR == false) {
      i ++;
      next = i;
      continue;
    }
    if (next == i) next ++; /*  Next in comparison, may be far ahead  */
    while (trlist[next].valid_TR == false && next != stats->clunk_end) {
      next ++;
    }
    if (next == stats->clunk_end && trlist[next].valid_TR == false) {
      i ++;
      next = i;
      continue;
    }
    /*  Non overlapping, thus stop  */
    if (trlist[i].origin_position + trlist[i].full_length < trlist[next].origin_position) {
      i = next;
      continue;
    }
    if (trlist[i].purity_percentage > trlist[next].purity_percentage) {
      trlist[next].valid_TR = false;
      continue;
    }
    if (trlist[i].purity_percentage < trlist[next].purity_percentage) {
      trlist[i].valid_TR = false;
      i = next;
      continue;
    }

    if (trlist[i].stats > trlist[next].stats) {
      trlist[next].valid_TR = false;
      continue;
    }
    if (trlist[i].stats < trlist[next].stats) {
      trlist[i].valid_TR = false;
      i = next;
      continue;
    }

    if (trlist[i].full_length >= trlist[next].full_length) {  /*  Here breaks the tie  */
      trlist[next].valid_TR = false;
    }
    else {
      trlist[i].valid_TR = false;
      i = next;
    }
  }
}

void filter_init_trlist (TRs_result_t *trlist, unsigned long int tot_trs) {
  unsigned long int i;
  for (i = 0; i < tot_trs; i ++) {
    trlist[i].valid_TR = true;
    trlist[i].stats = 0;
  }
}

void filter_thresholdFilter (TRs_result_t *trlist,  unsigned long int tot_trs, struct config *params) {
  unsigned long int i;
  for (i = 0; i < tot_trs; i ++) {
    if ((trlist[i].purity_percentage < params->min_purity) || (trlist[i].partial_length < params->min_TR_len))
      trlist[i].valid_TR = false;
    else
      trlist[i].valid_TR = true;
    trlist[i].stats = 0;
  }
}

void filter (TRs_Result_Bundle *results, struct config *params) {
  struct filter_object *fObj = NULL;
  struct clunk_stats *stats = NULL;
  TRs_result_t *trlist = results->TRs_found;
  unsigned long int tot_trs = results->trs_found_offset;
    
  if (params->filter_type == NONE) {  /*  No filtering - set values for valid_TR to true  */
    filter_init_trlist (trlist, tot_trs);
    return;
  }
  filter_thresholdFilter (trlist, tot_trs, params);

  if (params->filter_type == THRESHOLD) return;  /*  filter out only too small TRs  */
  fObj = init_filter_obj (trlist, tot_trs);
  if (fObj == NULL) return;

  
  while ((stats = getNextContiguousSetOfTRs (fObj, stats)) != NULL) {
    applyFilter (fObj, stats, params);
  }
  
  free (fObj);  /*  Release the memory of fObj  */
}
