#include "algorithm.h"

/*****************************************************************************\
******************************** TRS_RESULT STRUCT ****************************
\*****************************************************************************/

void reset_TRs_result(TRs_result_t* tr_result) {   /*  XXX  replace me with memset  */
	if (tr_result != NULL) {
		tr_result->origin_position = 0;
		tr_result->insertions_count = 0;
		tr_result->full_length = 0;
		tr_result->partial_length = 0;
		tr_result->copy_number = 0;
		tr_result->period = 0;
		tr_result->valid_TR = true;  /*  Used for filtering  */
		tr_result->stats = 0;  /*  Used for filtering  */
		tr_result->purity_percentage = 0;
		tr_result->motif_start_index = 0;
		tr_result->motifs_number = 0;
	}
}


void copy_TR_struct(TRs_result_t* main, TRs_result_t* dest) { /*  XXX  replace me with memcpy  */	

		dest->origin_position = main->origin_position;
		dest->insertions_count = main->insertions_count;
		dest->full_length = main->full_length;
		dest->partial_length = main->partial_length;
		dest->copy_number = main->copy_number;
		dest->period = main->period;
		dest->valid_TR = main->valid_TR;  /*  Used for filtering  */
		dest->stats = main->stats;  /*  Used for filtering  */
		dest->purity_percentage = main->purity_percentage;
		dest->motif_start_index = main->motif_start_index;
		dest->motifs_number = main->motifs_number;
}

/*****************************************************************************\
*******************************************************************************
\*****************************************************************************/

/*****************************************************************************\
**************************** TRS_RESULT_BUNDLE STRUCT *************************
\*****************************************************************************/

TRs_Result_Bundle* init_TRs_Bundle(int init_n_elem_trs, int init_n_elem_motifs) {

	TRs_Result_Bundle* trs_b;
	unsigned long int size;

	if ( ( trs_b = (TRs_Result_Bundle*) malloc( sizeof(TRs_Result_Bundle) ) ) == NULL) {
			perror("Malloc error for trs_b in init_TRs_Bundle\n");
			return NULL;
	};

	/*  Allocate TRS bundle  */
	size = init_n_elem_trs * sizeof(TRs_result_t);
	if ((trs_b->TRs_found = (TRs_result_t*) malloc( size )) == NULL) {
			perror("Malloc error for trs_b->TRs_found in init_TRs_Bundle\n");
			return NULL;
	};
	trs_b->TR_array_size = size;
	trs_b->max_trs_number = init_n_elem_trs;
	trs_b->trs_found_offset = 0;
	/*
	for ( i = 0; i < trs_b->max_trs_number; ++i ) {
		reset_TRs_result(&(trs_b->TRs_found[i]));		
	}
	*/
	/* Allocate motifs array  */
	size = init_n_elem_motifs * sizeof(unsigned short int);
	if ((trs_b->motif_lengths = (unsigned short int*) malloc( size )) == NULL) {
			perror("Malloc error for trs_b->motif_lengths in init_TRs_Bundle\n");
			return NULL;
	};
	trs_b->motif_lengths_size = size;
	trs_b->max_motifs_number = init_n_elem_motifs;
	trs_b->motif_lengths_offset = 0;
	/*
	for ( i = 0; i < trs_b->max_motifs_number; ++i ) {
		trs_b->motif_lengths[i] = 0;		
	}
	*/

	return trs_b;

}

short int copy_TRs_Bundle(TRs_Result_Bundle* main, TRs_Result_Bundle* dest) {
	
	unsigned long int i;

	if ( dest->max_trs_number < main->max_trs_number ) {
	  /*printf("TRS COPY NUMBER (size lower):%u\n", dest->copy_number); */
		
		if ((dest->TRs_found = (TRs_result_t*) realloc(dest->TRs_found, main->TR_array_size)) == NULL) {
			perror("Realloc error for dest->TRs_found in copy_TRs_Bundle\n");
			return 1;
		};

		dest->TR_array_size = main->TR_array_size;
		dest->max_trs_number = main->max_trs_number;
		dest->trs_found_offset = main->trs_found_offset;

		for ( i = 0; i < main->trs_found_offset ; ++i ) {
			copy_TR_struct(&(main->TRs_found[i]), &(dest->TRs_found[i]));
		}

	} 
	else {/*Copy only the content as the size is ok (Equal or bigger)  */
	  /*printf("TRS COPY NUMBER (size equal or bigger):%u\n", dest->copy_number);*/
		
		dest->trs_found_offset = main->trs_found_offset;
		for ( i = 0; i < main->trs_found_offset ; ++i ) {
			copy_TR_struct(&(main->TRs_found[i]), &(dest->TRs_found[i]));
		}
		/*
		for ( i = dest->trs_found_offset; i < dest->max_trs_number ; ++i ) {
			reset_TRs_result(&(dest->TRs_found[i]));
		}
		*/
	}

	if ( dest->max_motifs_number < main->max_motifs_number ) {
	  /*printf("TRS COPY NUMBER (size lower):%u\n", dest->copy_number);*/
		
		if ((dest->motif_lengths = (unsigned short int *) realloc(dest->motif_lengths, main->motif_lengths_size)) == NULL) {
			perror("Realloc error for dest->motif_lengths in copy_TR_struct\n");
			return 1;
		};

		dest->motif_lengths_size = main->motif_lengths_size;
		dest->max_motifs_number = main->max_motifs_number;
		dest->motif_lengths_offset = main->motif_lengths_offset;

		for ( i = 0; i < main->motif_lengths_offset ; ++i ) {
			dest->motif_lengths[i] = main->motif_lengths[i];
		}

	} 
	else {/*Copy only the content as the size is ok (Equal or bigger)*/
		/*printf("TRS COPY NUMBER (size equal or bigger):%u\n", dest->copy_number);*/
		
		dest->motif_lengths_offset = main->motif_lengths_offset;
		for ( i = 0; i < main->motif_lengths_offset ; ++i ) {
			dest->motif_lengths[i] = main->motif_lengths[i];
		}
		/*
		for ( i = dest->motif_lengths_offset; i < dest->max_motifs_number ; ++i ) {
			dest->motif_lengths[i] = 0;
		}
		*/
	}

	return 0;
}

short int motifs_amount_check_and_resize(TRs_Result_Bundle* trs_bundle, int init_n_elem_motifs) {
  if ( trs_bundle->motif_lengths_offset == trs_bundle->max_motifs_number ) {
    if ((trs_bundle->motif_lengths = (unsigned short int *) realloc(trs_bundle->motif_lengths, trs_bundle->motif_lengths_size + init_n_elem_motifs * sizeof(unsigned short int))) == NULL) {
      perror("Realloc error for trs_bundle->motif_lengths in motifs_amount_check_resize\n");
      return 1;
    };
    trs_bundle->motif_lengths_size = trs_bundle->motif_lengths_size + init_n_elem_motifs * sizeof(unsigned short int);
    trs_bundle->max_motifs_number += init_n_elem_motifs;
  }
  return 0;
}

short int trs_amount_check_and_resize(TRs_Result_Bundle* obj, int init_n_elem_trs) {
  if ( obj->trs_found_offset == obj->max_trs_number ) {
    if ((obj->TRs_found = (TRs_result_t*) realloc(obj->TRs_found, obj->TR_array_size + init_n_elem_trs * sizeof(TRs_result_t))) == NULL) {
      perror("Realloc error for obj->TRs_found in trs_amount_check_and_resize\n");
      return 1;
    };
    obj->TR_array_size = obj->TR_array_size + init_n_elem_trs * sizeof(TRs_result_t);
    obj->max_trs_number += init_n_elem_trs;
  }
  return 0;
}

short int insert_TRresult_inBundle(TRs_Result_Bundle* trs_bundle, TRs_Result_Bundle* trs_bundle_single, int resize_n_elem_trs, int resize_n_elem_motifs) {

	unsigned long int j;
	/*printf("TRS_OFFSET %u\n", trs_bundle_single->trs_found_offset);*/
	if (trs_bundle_single->trs_found_offset == 1) {
	
		if (trs_amount_check_and_resize(trs_bundle, resize_n_elem_trs)) {
			perror("Error in insert_TRresult_inBundle\n");
			return 1;
		}
		/*  Copy the TR_found  */
		copy_TR_struct(&(trs_bundle_single->TRs_found[0]), &(trs_bundle->TRs_found[trs_bundle->trs_found_offset]));
		/*  Copy the start index of the motif lengths  */
		trs_bundle->TRs_found[trs_bundle->trs_found_offset].motif_start_index =	trs_bundle->motif_lengths_offset;				
		trs_bundle->trs_found_offset++;

		/*  Copy the motifs lengths  */
		for ( j = 0; j < trs_bundle_single->motif_lengths_offset ; ++j) {
			insert_TRmotif_inTRresult(trs_bundle, trs_bundle_single->motif_lengths[j], resize_n_elem_motifs);
		}
		
	} else {
		printf("Error in insert_TRresult_inBundle: Must have 1 result\n");
		return 1;
	}
	return 0;
}

short int insert_TRmotif_inTRresult(TRs_Result_Bundle* tr_result, short int length, int init_n_elem_motifs) {

	if ( motifs_amount_check_and_resize(tr_result, init_n_elem_motifs) ) {
		perror("Error in insert_TRmotif_inBundle\n");
		return 1;
	}
	tr_result->motif_lengths[tr_result->motif_lengths_offset] = length;
	tr_result->motif_lengths_offset++;


	return 0;
}

void destroy_TRs_Bundle(TRs_Result_Bundle** trs_bundle) {

	if (*trs_bundle != NULL) {

		free((*trs_bundle)->motif_lengths);		
		free((*trs_bundle)->TRs_found);
		free(*trs_bundle);
		*trs_bundle = NULL;
	}
}

void reset_TRs_Bundle(TRs_Result_Bundle* trs_bundle) {
	if (trs_bundle != NULL) {
		/*
		for ( i = 0; i < trs_bundle->trs_found_offset; ++i ) {
			reset_TRs_result(&(trs_bundle->TRs_found[i]));
		}
		*/
		trs_bundle->trs_found_offset = 0;
		/*
		for ( i = 0; i < trs_bundle->motif_lengths_offset; ++i ) {
			trs_bundle->motif_lengths[i] = 0;
		}
		*/
		trs_bundle->motif_lengths_offset = 0;

	}
}


/*****************************************************************************\
*******************************************************************************
\*****************************************************************************/

/*****************************************************************************\
****************************** RESULT_FINDTR STRUCT ***************************
\*****************************************************************************/

result_findTR* init_result_struct() {

	result_findTR *rs;
	
	if ((rs = (result_findTR*) malloc(sizeof(result_findTR))) == NULL) {
			perror("Malloc Error in init_results_struct()\n");
			return NULL;
	}
	rs->result_code = NO_TR_FOUND;

	if ((rs->resulted_TR = init_TRs_Bundle( 1, RESIZE_TR_MOTIFS_AMOUNT)) == NULL) {
			perror("Malloc Error in init_results_struct()\n");
			return NULL;
	}	
	
	return rs;
}

void reset_result_struct(result_findTR* rs) {

	if (rs != NULL) {
		rs->result_code = NO_TR_FOUND;
		reset_TRs_Bundle(rs->resulted_TR);
	}
}

void destroy_result_struct(result_findTR** rs) {

	destroy_TRs_Bundle(&((*rs)->resulted_TR));
	free(*rs);
	*rs = NULL;
}

/*****************************************************************************\
*******************************************************************************
\*****************************************************************************/

/*****************************************************************************\
**************************** DOT_THREAD_INPUT STRUCT **************************
\*****************************************************************************/

Dot_Thread_input* dot_Thread_obj_init(struct config *cp, MATCH_ARRAY_TYPE **wm, struct filemanager *fm, struct outfile* output, short int t_id ) {
	
	Dot_Thread_input* dot_Thread_input;

	if ((dot_Thread_input = (Dot_Thread_input*) malloc(sizeof(Dot_Thread_input))) == NULL) {
		perror("Error in allocating memory for dot_Thread_input in dot_Thread_obj_init\n");
		return NULL;
	}
	
	/*  input sequence  */
	dot_Thread_input->sequence = NULL;
	/*  sequence name  */
	dot_Thread_input->IDSeq = NULL;
	dot_Thread_input->thread_TRs_bundle = init_TRs_Bundle(RESIZE_TRS_AMOUNT, RESIZE_MOTIFS_AMOUNT);	
	dot_Thread_input->matrix = NULL;
	dot_Thread_input->config_params = cp;
	dot_Thread_input->file_manager = fm;
	dot_Thread_input->weight_matrix = wm;
	dot_Thread_input->output = output;
	dot_Thread_input->t_id = t_id;

	return dot_Thread_input;
}

void destroy_dot_Thread_obj(Dot_Thread_input** obj) {

	
	destroy_TRs_Bundle(&((*obj)->thread_TRs_bundle));
	/* freed at the end  */
	(*obj)->matrix = NULL;
	(*obj)->config_params = NULL;
	(*obj)->weight_matrix = NULL;
	
	free(*obj);
	*obj = NULL;
}

void reset_dot_Thread_obj(Dot_Thread_input* obj) {

	free(obj->sequence);
	free(obj->IDSeq);
	reset_TRs_Bundle(obj->thread_TRs_bundle);
	obj->matrix = NULL;

}

/*****************************************************************************\
*******************************************************************************
\*****************************************************************************/

/*****************************************************************************\
*******************************************************************************
\*****************************************************************************/

char* copy_string(char* seq) {
	int length, i=0;
	char* s; 
	
	if (seq == NULL) return NULL;
	for (length = 0; seq[length] != '\0'; length ++);
	/*  XXX - perche' + 1 ?? */
	if ((s = malloc(length * sizeof(char) + 1))==NULL) {
		perror("Malloc Error in copy_string()\n");
		return NULL;
	}
	/*This cycle copies seq into s */
	for (i = 0; i < length; i ++) {
	  s[i] = seq[i];
	}
	s[i] = '\0';
	return s;
}

char* copy_seqPart(int start_index, int length, char* seq) {
	int seq_ind=start_index, i=0;
	char* s; 
	
	if ((s = malloc(length * sizeof(char)+1))==NULL) {
		perror("Malloc Error in copy_seqPart()\n");
		return NULL;
	}
	/*This cycle copies and return seq from start_index to start_index+length */
	while (i<length) {
		s[i] = seq[seq_ind];
		seq_ind++;
		i++;
	}
	s[i] = '\0';
	return s;
}

/*****************************************************************************\
*******************************************************************************
\*****************************************************************************/

/*****************************************************************************\
***************************** ALGORITHM FUNCTIONS *****************************
\*****************************************************************************/


int getFinalTandemLength (int window_index, int current_match_index, int window_length, int max_length, MATCH_ARRAY_TYPE** pointers) {
	
	int w_i=window_index, value=0, c_i=current_match_index, final_length=0;
	while(final_length <= max_length) {
		value = pointers[w_i][c_i];
		if (value == 1) {
			final_length ++;
			w_i ++;
			c_i ++;
		} else {
			break;
		}		
	}
	return final_length;
}
MATCH_ARRAY_TYPE getValueFromMatchArrays(int char_index, int match_char_index, MATCH_ARRAY_TYPE** pointers) {
	/*return the match value in the array pointed by 'char_index' in position match_char_index*/
	return pointers[char_index][match_char_index];
}

MATCH_ARRAY_TYPE getSumFromMatchValues(int window_length, int window_index, int current_match_index, MATCH_ARRAY_TYPE** pointers) {
	int i=window_index, j=current_match_index;
	MATCH_ARRAY_TYPE sum=0;
	/*  XXX -  Secondo me qui c'e' un bug  */
	while(j<current_match_index+window_length) {
		sum += getValueFromMatchArrays(i, j, pointers);
		i++;
		j++;
	}
	/*Sum of match values in current window*/
	return sum;
}


int findTandemRepeats(int window_length, int window_index, struct dot_matrix *m, MATCH_ARRAY_TYPE minThreshold, int max_jumps, result_findTR* rs) {
	int current_match_index=window_index+window_length, i=0, 
		final_tandem_length=0, max_length=0, insertions_count=0,
		current_insertion_len=0, match_stop=0, jump_limit_modifier=0, 
		first_useless=0, jumps_limit=max_jumps, jumps_bundle_allowed=window_length/2;
	unsigned short motifs_number=0;
	short copy_number=0;
	MATCH_ARRAY_TYPE purity_sum = window_length;
	MATCH_ARRAY_TYPE sum;

	if ( insert_TRmotif_inTRresult(rs->resulted_TR, window_length, RESIZE_TR_MOTIFS_AMOUNT) ) {
		perror("Error in findTandemRepeats: first motif (window_index)\n");
		return 1;
	}
	motifs_number++;
	copy_number++;

	while (current_match_index <= m->sequence_len-window_length) {
		sum = getSumFromMatchValues(window_length, window_index, current_match_index, m->pointers_sequence);

    	if ((sum >= minThreshold) && (sum <= window_length)) {

			
			if (match_stop) { /* save skipped sequence as insertion */

				current_insertion_len = current_match_index-first_useless;				
				
				if ( insert_TRmotif_inTRresult(rs->resulted_TR, current_insertion_len, RESIZE_TR_MOTIFS_AMOUNT) ) {
					perror("Error in findTandemRepeats: match_stop block\n");
					return 1;
				}
				motifs_number++;

				if (!jump_limit_modifier) { 
					jump_limit_modifier=1; 
					jumps_limit=current_insertion_len; 
				}
				match_stop=0;
				i = 0;
				insertions_count += current_insertion_len;
				first_useless = 0;
				jumps_bundle_allowed--;
			}

			if (rs->result_code != TR_FOUND) { rs->result_code = TR_FOUND; }
			/*I use always the first because in this function I find 1 TR  */
			copy_number++;
			
			m->mask[current_match_index]=CHECKED;
			purity_sum += sum;


			if ( insert_TRmotif_inTRresult(rs->resulted_TR, window_length, RESIZE_TR_MOTIFS_AMOUNT) ) {
				perror("Error in findTandemRepeats:  TR motif\n");
				return 1;
			}
			motifs_number++;
			
			current_match_index += window_length;
		} else { /* sum doesn't match */
		
			if (!match_stop) { first_useless = current_match_index; match_stop=1; } /* if true it would be the first stop*/
			if ((rs->result_code != TR_FOUND) && (i>=max_jumps)) { /* if there were no Tandem and i CAN NOT jump, stop finding */
				rs->result_code = NO_TR_FOUND;
			} else {/* if there was at least one Tandem found and I can make some skips more */
				if ((i < jumps_limit) && (jumps_bundle_allowed > 0)) { /*if CAN jump, make 1 step forward and continue the cycle*/
					i++;
					current_match_index++;					
					continue; 
				}
			}						
			break;
		}
	}		
	
	/* create the element TandemRepeats to insert in Tandem List */
	if (copy_number > 1) {
				
		/* Look for part of cut tandem (include decimal exponent) HAMMING DISTANCE */
		current_match_index -= i;	
		max_length = (current_match_index+window_length > m->sequence_len) ? (m->sequence_len-current_match_index) : (window_length-1);
		final_tandem_length=getFinalTandemLength (window_index, current_match_index, window_length, max_length, m->pointers_sequence);	
		if (final_tandem_length > 0) {

			purity_sum += final_tandem_length;

			if ( insert_TRmotif_inTRresult(rs->resulted_TR, final_tandem_length, RESIZE_TR_MOTIFS_AMOUNT) ) {
				perror("Error in findTandemRepeats:  final tandem block\n");
				return 1;
			}
			motifs_number++;			
		}
		
		/* +1 because I add the count of the first window_copy as well */
		rs->resulted_TR->TRs_found[0].copy_number = copy_number;
		rs->resulted_TR->TRs_found[0].full_length = (window_length * copy_number) + insertions_count + final_tandem_length;
		rs->resulted_TR->TRs_found[0].partial_length = (window_length * copy_number) + final_tandem_length;
		rs->resulted_TR->TRs_found[0].purity_percentage = (MATCH_ARRAY_TYPE) purity_sum/(rs->resulted_TR->TRs_found[0].full_length);
		rs->resulted_TR->TRs_found[0].origin_position = window_index;
		rs->resulted_TR->TRs_found[0].period = window_length;	
		rs->resulted_TR->TRs_found[0].stats = 0;
		rs->resulted_TR->TRs_found[0].valid_TR = true;			
		rs->resulted_TR->TRs_found[0].insertions_count = insertions_count;
		rs->resulted_TR->TRs_found[0].motif_start_index = 0;
		rs->resulted_TR->TRs_found[0].motifs_number = motifs_number;
		/*The bundle has one result */
		rs->resulted_TR->trs_found_offset = 1;
	} else {
		/* NO TR FOUND */
		/* Leave the result_struct unmodified  */
	        /* It is resetted out of this function */
	}
	return 0;
}

int isLastIntersected(TRs_result_t* prec, TRs_result_t* last) {
	int last_index_t1=prec->origin_position+prec->full_length-1, last_index_t2=last->origin_position+last->full_length-1;
	
	return ((last_index_t1 > last->origin_position) && (last_index_t1 < last_index_t2));
}
int isLastIncluded(TRs_result_t* prec, TRs_result_t* last) {
	int last_index_t1=prec->origin_position+prec->full_length-1, last_index_t2=last->origin_position+last->full_length-1;
	
	return ((last_index_t1 > last->origin_position) && (last_index_t1 >= last_index_t2));
}

void expansion_filter(TRs_Result_Bundle* TRs_bundle, TRs_Result_Bundle* last_tandem_found, int biggest_full_length, float tollerance) {
	
	unsigned long int i, tot_trs_left, j;
	unsigned short int motifs_number;
	float current_score=0, max_score=0;
	TRs_result_t* results, *last;

	if (TRs_bundle != NULL) {

		results = TRs_bundle->TRs_found;	
		tot_trs_left = TRs_bundle->trs_found_offset;
		
		if (TRs_bundle->trs_found_offset > 0) {
			
#ifdef DEBUG_ALG
			printf("\t|-- Exp Filter --> There is one or more TRs: %d\n", tot_trs_left);		
#endif	
			
			for ( i = 0; i < TRs_bundle->trs_found_offset; ++i ) {
				current_score = (float) results[i].full_length/(float)biggest_full_length + results[i].purity_percentage;
				if (current_score > max_score) max_score = current_score;
			}
			
			/*First filter by score and partial length*/
			last = &results[0];
			for ( i = 0; i < TRs_bundle->trs_found_offset; ++i ) {
				current_score = (float) results[i].full_length/(float)biggest_full_length + results[i].purity_percentage;
				/*printf("CURR SCORE AFTER: %.3f\n", current_score);*/
				if (current_score >= max_score - tollerance) {
					last = &results[i];
				} else if ((results[i].partial_length == (last)->partial_length) && (results[i].purity_percentage > (last)->purity_percentage)) {	
					last = &results[i];
				} else {
					/*Leave the same result of before
					and reset the current one*/
					reset_TRs_result(&results[i]);
					tot_trs_left--;
				}
			}
			/* Cycle all the array because results could be between NULL pointers */
			/*Second, filter by purity */
			/*Start from the first valid*/
#ifdef DEBUG_ALG
			printf("\t|-- Exp Filter --> There are %d TRs left after filtering by score and partial length\n", tot_trs_left);			
#endif
			i = 0;
			while ( (results[i].copy_number == 0) && ( i < TRs_bundle->trs_found_offset ) ) ++i;
			last = &results[i];
			++i;

			/*Still there are more than 1 results*/
			if (tot_trs_left > 1) {
				while( i < TRs_bundle->trs_found_offset ) {
					if (results[i].copy_number == 0) { ++i; continue; }
					if (results[i].purity_percentage > (last)->purity_percentage) {
						reset_TRs_result(last);
						last = &results[i];
					} else {
						reset_TRs_result(&results[i]);
					}
					++i;				
				}
			}
			/*Copy last into last_tandem_found getting the motifs too*/
			copy_TR_struct(last, &(last_tandem_found->TRs_found[0]));
			last_tandem_found->TRs_found[0].motif_start_index = 0; /*always zero at this point  */
			last_tandem_found->trs_found_offset = 1;
			last_tandem_found->motif_lengths_offset = 0; /*always zero at this point*/
			motifs_number = last->motifs_number;
			for ( j = last->motif_start_index; motifs_number > 0; ++j, --motifs_number ) {
				insert_TRmotif_inTRresult(last_tandem_found, TRs_bundle->motif_lengths[j], RESIZE_TR_MOTIFS_AMOUNT);
			}			
		} else {
#ifdef DEBUG_ALG
			printf("\t|-- Exp Filter --> NO TRs in bundle\n");		
#endif
		}
	} 
}

void print_usintArray_new(unsigned short * array,int init, int length) {
	int i;
	for(i=init; i<init+length; i++) {
		printf("%hu ", array[i]);
	}
	printf("\n");
}

int start_TRs_search (Dot_Thread_input* param) {
	int window_length, window_length_max, window_index, window_end, max_length_flag, 
		jumps_limit=0, gaps_limit=0, biggest_full_length=0;
	MATCH_ARRAY_TYPE minThresholdByPercentage, minThresholdByGaps, minThreshold;
	struct dot_matrix *m = param->matrix;
	struct config *cfg = param->config_params;
	result_findTR *current_result;	
	TRs_Result_Bundle *trs_current_bundle, *last_tandem_found, *previous_window_tandem;;

	if ((cfg->fvalue < 0) || (cfg->fvalue > 1)) {
		printf ("MinMatch ranges in [0,1]\n");
		return 1;
	}
	if (cfg->gvalue < 0) { 
		printf ("max_gaps cannot be negative\n"); 
		return 1; 
	} 
	else {
		gaps_limit=cfg->gvalue;
	}

	if (cfg->xvalue < 1) { /* max_length not inserted in command line */
		max_length_flag = 0;
	} 
	else { /* max_length inserted */
		if (cfg->xvalue < cfg->nvalue) { /* if lower than min_length */
			printf ("max_length < min_length not allowed\n");
			return 1;
		}
		max_length_flag = 1;
	}
	if (cfg->nvalue < 1) { 
		printf ("motif length cannot be lower than 1\n");
		return 1;
	} 
	else { 
		window_length = cfg->nvalue; 
	}

	window_index=0;   
	window_end=m->sequence_len; 


	if (cfg->jvalue < 0) {
		printf ("MaxInsert cannot be negative\n");
		return 1;
	}

	if ((current_result = init_result_struct()) == NULL) {
		perror ("Error in start_TRs_search() for current_result initialization\n");
		return 1;
	}

	if (( trs_current_bundle = init_TRs_Bundle(RESIZE_TR_BUNDLE_AMOUNT, RESIZE_TR_MOTIFS_AMOUNT) ) == NULL) {
			perror("Error in initializing trs_current_bundle\n");
			return 1;
	}

	if (( previous_window_tandem = init_TRs_Bundle( 1, RESIZE_TR_MOTIFS_AMOUNT ) ) == NULL) {
			perror("Error in initializing previous_window_tandem\n");
			return 1;
	}

	if (( last_tandem_found = init_TRs_Bundle( 1, RESIZE_TR_MOTIFS_AMOUNT ) ) == NULL) {
			perror("Error in initializing trs_current_bundle\n");
			return 1;
	}

	/*  Da qui  */
	while ( window_index < window_end ) {
#ifdef DEBUG_ALG
		printf("TRS WINDOW INDEX STARTED:%d\n", window_index);
#endif
		if (max_length_flag && ( cfg->xvalue <= ((m->sequence_len - window_index)/2) )) {
			window_length_max = cfg->xvalue;			
		} else {
			window_length_max = (m->sequence_len - window_index)/2;			
		}		
		
		
		biggest_full_length = 0;
		while( window_length <= window_length_max ) {			
#ifdef DEBUG_ALG
			printf("\tTRS WINDOW LENGTH STARTED:%d\n", window_length);
#endif					
			minThresholdByPercentage = (MATCH_ARRAY_TYPE) cfg->fvalue * window_length;
			minThresholdByGaps = (gaps_limit >= window_length) ? ((MATCH_ARRAY_TYPE) window_length-1) : ((MATCH_ARRAY_TYPE) window_length - gaps_limit);					
			if ((cfg->fvalue == 1) && (gaps_limit > 0)) { /* only gaps parameter */
				minThreshold = minThresholdByGaps;
			} else if ((gaps_limit == 0) && (cfg->fvalue < 1)) { /* only percentage parameter */
				minThreshold = minThresholdByPercentage;
			} else { /* both or no parameter: percentage and gaps */
				minThreshold = MAX(minThresholdByPercentage, minThresholdByGaps); 
			}

#ifdef DEBUG_ALG
			printf("\tMINTHRESHOLD %f\n", minThreshold);
#endif	
			/* Jump limit is at most window_length or cfg->jvalue */
			if (cfg->jvalue > window_length) { jumps_limit = window_length; } else { jumps_limit = cfg->jvalue; }

			if (findTandemRepeats(window_length, window_index, m, minThreshold, jumps_limit, current_result)) {
				printf ("Error in start_TRs_search() for current_result\n");
				return 1;
			}
			switch (current_result->result_code) {
				case (NO_TR_FOUND) : {
					window_length++;
					break;
				}
				case (TR_FOUND) : { 
					window_length ++;
					/* Put all the TRs with the same window_index in the TRs_Bundle_list and check which is the highest length */
					if (current_result->resulted_TR->TRs_found[0].full_length > biggest_full_length) biggest_full_length = current_result->resulted_TR->TRs_found[0].full_length;
					if (insert_TRresult_inBundle(trs_current_bundle, current_result->resulted_TR, RESIZE_TR_BUNDLE_AMOUNT, RESIZE_TR_MOTIFS_AMOUNT)) {
						perror("Error in inserting the current result in the TRs bundle\n");
						return 1;
					}
#ifdef DEBUG_ALG
					printf("\tTRFOUND at %d INSERTED\n", window_index);
#endif
					break;
				}
				default : {
					printf ("Unknown case!\n");
					return 1;
				}
			}
			/* Reset current_result  */
			reset_result_struct(current_result);
			
		}
		
		/* Apply Expansion Filter to the TRs list of the same zone (same Window_index)*/
		/*last_tandem_found is NULL or has a valid TR to compare to the previous one */
		expansion_filter(trs_current_bundle, last_tandem_found, biggest_full_length, cfg->tollerance);
#ifdef DEBUG_ALG
		printf("\tTRS_BUNDLE of %d FILTERED\n", window_index);
#endif
		
		if ( last_tandem_found->motif_lengths_offset == 0 ) {/* Tandem not found */
			window_index++; 
		} else {
			
		  /*Insert the tandem as it is the first found	 */
			if (previous_window_tandem->motif_lengths_offset == 0) {
#ifdef DEBUG_ALG
				printf("\tTRS WINDOW INDEX PREVIUOS WITH 0 COPY NUMBER:%d\n", window_index);
#endif
				if (insert_TRresult_inBundle(param->thread_TRs_bundle , last_tandem_found, RESIZE_TRS_AMOUNT, RESIZE_MOTIFS_AMOUNT)) {
					perror("Error in inserting the current result in the Thread TRs bundle\n");
					return 1;
				};
				if ( copy_TRs_Bundle(last_tandem_found, previous_window_tandem) ) {
					perror("Error in copying last_tandem_found to previous_window_tandem in start_TRs_search\n");
					return 1;
				};
				window_index++;
			} else {
				
				if (isLastIncluded(&(previous_window_tandem->TRs_found[0]), &(last_tandem_found->TRs_found[0]))) {
					/* included Tandems */
#ifdef DEBUG_ALG
					printf("\tTRS WINDOW INDEX INCLUDED:%d\n", window_index);
#endif
					if (last_tandem_found->TRs_found[0].purity_percentage > previous_window_tandem->TRs_found[0].purity_percentage/*(last_tandem_found->period > previous_window_tandem->period)*/) { /* if false it has found an included TR which is the same but fragmented */
						switch (m->mask[window_index]) {
							case (UNCHECKED) : {
								/* tandem has not been checked before */ 				
								if (insert_TRresult_inBundle(param->thread_TRs_bundle , last_tandem_found, RESIZE_TRS_AMOUNT, RESIZE_MOTIFS_AMOUNT)) {
									perror("Error in inserting the current result in the Thread TRs bundle\n");
									return 1;
								};
								window_index++;
								break;											
							}
							case (CHECKED) : {
								/* it is a part of an other tandem found before */
								window_index++;
								break;
							}
							default : { break; }
						}
					} else {
						window_index++;
					}
				} else { /* intersected Tandems */		
#ifdef DEBUG_ALG
					printf("\tTRS WINDOW INDEX INTERSECTED:%d\n", window_index);				
#endif
					if (insert_TRresult_inBundle(param->thread_TRs_bundle , last_tandem_found, RESIZE_TRS_AMOUNT, RESIZE_MOTIFS_AMOUNT)) {
						perror("Error in inserting the current result in the Thread TRs bundle\n");
						return 1;
					};
					/*reset_TRs_Bundle(previous_window_tandem);*/
					if ( copy_TRs_Bundle(last_tandem_found, previous_window_tandem) ) {
						perror("Error in copying last_tandem_found to previous_window_tandem in start_TRs_search\n");
						return 1;
					}
					window_index++;
				}
			}
		}
		reset_TRs_Bundle(last_tandem_found);
		reset_TRs_Bundle(trs_current_bundle);
		window_length=cfg->nvalue;
	
	}
	destroy_TRs_Bundle(&last_tandem_found);
	destroy_TRs_Bundle(&previous_window_tandem);
	destroy_result_struct(&current_result);
	destroy_TRs_Bundle(&trs_current_bundle);
	return 0;
}

/*****************************************************************************\
*******************************************************************************
\*****************************************************************************/

