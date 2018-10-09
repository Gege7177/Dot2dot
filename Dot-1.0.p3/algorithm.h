#ifndef ALGORITHM_H_
#define ALGORITHM_H_

#include "common.h"
#include "dot_matrix.h"
#include "threads_manager.h"
#include "configReader.h" 


#define MAX(x,y)	(((x)<(y)) ? (y) : (x))
#define MIN(x,y)	(((x)>(y)) ? (y) : (x))
/*** findTandemRepeats function RESULTS ***/
#define NO_TR_FOUND 0
#define TR_FOUND 1

/* Sizes for thread's TRs array... */
#define RESIZE_TRS_AMOUNT 100
#define RESIZE_MOTIFS_AMOUNT 1000
/* and for local TRs array  */
#define RESIZE_TR_BUNDLE_AMOUNT 10
#define RESIZE_TR_MOTIFS_AMOUNT 100

/*****************************************************************************\
******************************** TRS_RESULT STRUCT ****************************
\*****************************************************************************/

/*single TR data*/
struct _TRs_result_t {
  long int origin_position;
  int insertions_count;
  int full_length;
  int partial_length;
  short int copy_number;
  short int period;
  _Bool valid_TR;  /*  Used for filtering  */
  float stats;  /*  Used for filtering  */
  MATCH_ARRAY_TYPE purity_percentage;
  unsigned long int motif_start_index; /* start index of motifs in motif_legths array of trs_result_bundle  */
  unsigned short int motifs_number; /* copy_number + insertion groups + final group */
};

typedef struct _TRs_result_t TRs_result_t;

/*
 * Reset TRs_result_t struct
 *
 * return: void
 */
void reset_TRs_result(TRs_result_t* tr_result);

/*
 * Copy TRs_result_t struct main in dest
 *
 * return: void
 */
void copy_TR_struct(TRs_result_t* main, TRs_result_t* dest);

/*****************************************************************************\
*******************************************************************************
\*****************************************************************************/

/*****************************************************************************\
**************************** TRS_RESULT_BUNDLE STRUCT *************************
\*****************************************************************************/

struct _trs_result_bundle {
  TRs_result_t* TRs_found;  /* final TRs array  */
  unsigned long int TR_array_size; /* dynamic size of TRs array (bytes)  */
  unsigned long int max_trs_number;
  /* next index available TRs position in array  */
  unsigned long int trs_found_offset;  /* XXXX  */
  /*  Motif structs  */
  unsigned short int *motif_lengths;
  unsigned long int motif_lengths_size;
  unsigned long int motif_lengths_offset;
  unsigned long int max_motifs_number;
};
typedef struct _trs_result_bundle TRs_Result_Bundle;

/*
 * Initialize TRs_Result_Bundle struct
 *
 * return: Pointer to inizialized structure
 */
TRs_Result_Bundle* init_TRs_Bundle(int init_n_elem_trs, int init_n_elem_motifs);
/*
 * Copy TRs_result_Bundle struct main in dest
 *
 * return: 0 for success, otherwise 1
 */
short int copy_TRs_Bundle(TRs_Result_Bundle* main, TRs_Result_Bundle* dest);

/*
 * Check and Resize TRs array (TRs_found) of TRs_Result_Bundle* tr_result
 *
 * return: 0 for success, otherwise 1
 */
short int motifs_amount_check_and_resize(TRs_Result_Bundle* tr_result, int init_n_elem_motifs);

/*
 * Check and Resize motif lengths array (motif_lengths) of TRs_Result_Bundle* obj
 *
 * return: 0 for success, otherwise 1
 */
short int trs_amount_check_and_resize(TRs_Result_Bundle* obj, int init_n_elem_trs);

/*
 * Insert 1 TRs_result_t (saved in trs_bundle_single) in TRs_Result_Bundle* trs_bundle
 * 
 * return: 0 for success, otherwise 1
 */
short int insert_TRresult_inBundle(TRs_Result_Bundle* trs_bundle, TRs_Result_Bundle* trs_bundle_single, int resize_n_elem_trs, int resize_n_elem_motifs);

/*
 * Insert 1 motif length in TRs_Result_Bundle* trs_bundle
 * 
 * return: 0 for success, otherwise 1
 */
short int insert_TRmotif_inTRresult(TRs_Result_Bundle* tr_result, short int length, int init_n_elem_motifs);

/*
 * Deallocate TRs_Result_Bundle* trs_bundle
 * 
 * return: void
 */
void destroy_TRs_Bundle(TRs_Result_Bundle** trs_bundle);

/*
 * Reset TRs_Result_Bundle* trs_bundle (No free used)
 * 
 * return: void
 */
void reset_TRs_Bundle(TRs_Result_Bundle* trs_bundle);

/*****************************************************************************\
*******************************************************************************
\*****************************************************************************/

/*****************************************************************************\
****************************** RESULT_FINDTR STRUCT ***************************
\*****************************************************************************/

/*
 * Function findTandemRepeats return this structure to make the extern cycle handle window's movement
 * 'result_code' (TR_FOUND or NO_TR_FOUND)
 */
struct _result_findTR {
	int result_code;
	TRs_Result_Bundle* resulted_TR;
};
typedef struct _result_findTR result_findTR;

result_findTR* init_result_struct();
void reset_result_struct(result_findTR* rs);
void destroy_result_struct(result_findTR** rs);

/*****************************************************************************\
*******************************************************************************
\*****************************************************************************/

/*****************************************************************************\
**************************** DOT_THREAD_INPUT STRUCT **************************
\*****************************************************************************/

/*Thread input data structure*/
struct _Dot_Thread_input {
  char* sequence;  /* input sequence  */
  char* IDSeq;  /*    sequence name  */
  TRs_Result_Bundle *thread_TRs_bundle;
  struct dot_matrix *matrix;
  struct config *config_params;
  struct filemanager *file_manager;
  MATCH_ARRAY_TYPE **weight_matrix;
  struct outfile *output;
  short int t_id;
};
typedef struct _Dot_Thread_input Dot_Thread_input;

Dot_Thread_input* dot_Thread_obj_init(struct config *cp, MATCH_ARRAY_TYPE **wm, struct filemanager *fm, struct outfile* output, short int t_id );
void destroy_dot_Thread_obj(Dot_Thread_input** obj);
void reset_dot_Thread_obj(Dot_Thread_input* obj);


/*****************************************************************************\
***************************** ALGORITHM FUNCTIONS *****************************
\*****************************************************************************/

char* copy_string(char* seq);

/*
 * Param: (int start_index, int length, char* seq)
 *
 * Copy a part of 'seq' starting from 'start_index' and made of 'length' chars
 * seq[start_index, start_index+length]
 *
 * return: copied sequence (already allocated in memory)
 */
char* copy_seqPart(int start_index, int length, char* seq) ;

/*
 * Param: (int window_length, int window_index, int current_match_index, pointers_array pointers)
 *
 * Get the sum of 'window_length' values in the current string (current_match_index), using the arrays' pointers
 * referring to the nucleobase in window which start from 'window_index'
 */
MATCH_ARRAY_TYPE getSumFromMatchValues(int window_length, int window_index, int current_match_index, MATCH_ARRAY_TYPE** pointers);


/*
 * Initialize result_findTR* structure and return it.
 * result_findTR* used by findTandemRepeats function
 */
result_findTR *init_results_struct ();


/*
 * Param: (int window_length, int window_index, matrix_struct m, MATCH_ARRAY_TYPE minThreshold, int max_jumps)
 *
 * Find Tandem Repeats of 'window_length' chars referring to window_index. It finds a values sum from m->pointers_array which reach at least 'minThreshold'.
 * If max_jumps was higher then 0 this function could jump forward on the sequence (m->sequence) to find a Tandem Repeat stopped by some Insertions.
 *
 */
int findTandemRepeats(int window_length, int window_index, struct dot_matrix * m, MATCH_ARRAY_TYPE minThreshold, int max_jumps, result_findTR* rs);

/*
 * Param: (TRs_result_t* prec, TRs_result_t* last)
 *
 * Check if 'last' is included in 'prec' tandem
 * 
 * return: 1 if included, else 0.
 */
int isLastIntersected(TRs_result_t* prec, TRs_result_t* last);

/*
 * Param: (TRs_result_t* prec, TRs_result_t* last)
 *
 * Check if 'last' is intersected and not Included in 'prec' tandem
 * 
 * return: 1 if intersected, else 0.
 */
int isLastIncluded(TRs_result_t* prec, TRs_result_t* last);

/*
 * Param: (TRs_result_t** TRs_bundle_list, TRs_result_t** last_tandem_found, int biggest_full_length, float tollerance )
 *
 * Filter the TRs_result_bundle TRs_bundle by choosing only one among the TRs through score, partial_length and purity percentage filter.
 * 
 * return: void
*/
void expansion_filter(TRs_Result_Bundle* TRs_bundle, TRs_Result_Bundle * last_tandem_found, int biggest_full_length, float tollerance);

/*
 * Param: (Dot_Thread_input* param);
 *
 * Uses the structure above to start calling 'findTandemRepeats' function to find repeated sequences
 *
 */
int start_TRs_search(Dot_Thread_input* param);

/*****************************************************************************\
*******************************************************************************
\*****************************************************************************/

#endif /* ALGORITHM_H_ */
