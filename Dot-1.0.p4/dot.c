#include <stdio.h>
#include "threads_manager.h" /* includes mark_time.h */
#include "algorithm.h"
#include "mark_time.h"
#include "utils.h"
#include "fileManager.h"
#include "dot_matrix.h"
#include "output.h"
#include "filtering.h"

_Bool verbose_output;

void* dot_thread_fn (void* args);

void* dot_thread_fn (void* args) {
    Dot_Thread_input* param = (Dot_Thread_input*) args;
    int res; 
    struct dot_matrix *dm;
    struct sequence_t *seq; /*   contains information on the last loaded seq  */

    /*tid = pthread_self();*/
    /*printf("thread %u will perform algorithm run attempt\n", tid);*/
    /* algorithm run attempt */	
    seq = NULL;
    do {
        
        pt_m_lock (&mtx_getseq);
        while ( param->t_id != read_seq_global ) {
            pthread_cond_wait(&cond_read_seq_global, &mtx_getseq);
        }
        seq = filemanager_next_seq (param->file_manager, seq);
        read_seq_global = ( read_seq_global + 1 ) % param->config_params->thread_max;
        if (verbose_output) {
            if (seq != NULL) {
                fprintf (stderr,"Processing %s length: %ld bp\n",seq->label, seq->sequence_size);
                fflush (stderr);
            }
        }
        if ( seq == NULL ) {
            /*printf("No more sequence to get from file, Thread %u\n", tid);*/
            pthread_cond_broadcast(&cond_read_seq_global);            
            pt_m_unlock(&mtx_getseq);
            continue;
        } 
        else {

          if ((param->sequence = copy_string(seq->sequence)) == NULL) {
              perror("Error in copying sequence in dot_thread_fn\n");
              pthread_cond_broadcast(&cond_read_seq_global);              
              pt_m_unlock(&mtx_getseq);
              return (void*) 1;
          };
          if ((param->IDSeq = copy_string(seq->label)) == NULL) {
              perror("Error in copying label in dot_thread_fn\n");
              pthread_cond_broadcast(&cond_read_seq_global);              
              pt_m_unlock(&mtx_getseq);
              return (void*) 1;
          };
          pthread_cond_broadcast(&cond_read_seq_global);          
          pt_m_unlock(&mtx_getseq);

          dm = dot_init (seq, param->weight_matrix);
          if (dm == NULL) {
              perror("Error in creating dot_matrix\n");
              return (void*) 1; 
          }
          param->matrix = dm;
          res = start_TRs_search(param);
          if (res != 0) {
              perror("Something was wrong searching TRs\n");
              return (void*) 1; 
          }

            

	  /* Print TRs on file  */
          pt_m_lock(&mtx_wrtfile);
        
          while ( param->t_id != write_file_global ) {
            pthread_cond_wait(&cond_write_file_global, &mtx_wrtfile);
          }
#ifdef DEBUG_THREAD
          printf("THREAD ID %u -- write_file_global %u -- LABEL %s\n", param->t_id, write_file_global, param->IDSeq);        
          printf("THREAD %u IS WRITING...\n", param->t_id);        
#endif          
          if (param->thread_TRs_bundle->trs_found_offset > 0) {
            /* FINAL LIST FILTERING */
            filter (param->thread_TRs_bundle, param->config_params);
            print_TRs_list_toFile (param->output , param->IDSeq, param->sequence, param->thread_TRs_bundle);
          } else {
            printf("THREAD ID %u: No Tandem Repeats found in %s\n", param->t_id, param->IDSeq);        
          }
          write_file_global = ( write_file_global + 1 ) % param->config_params->thread_max;
          pthread_cond_broadcast(&cond_write_file_global);
                    
          pt_m_unlock(&mtx_wrtfile);
          
          reset_dot_Thread_obj(param);
          dot_free(dm);
        }    

    } while (seq != NULL);

    if ((param->file_manager)->finish != true) {
        perror ("Error reading input file\n");
        return (void*) 1;
    }
    return (void*) 0;
}

int main (int argc, char* argv[]) {
  int i;
  /*long int timer;*/
  Thread_List_Elem* last_inThreadList=NULL, *p;
  int err;
  pthread_t tid;
  Dot_Thread_input* search_param;
  struct paramList *par_list=NULL;
  struct config *cfg;
  MATCH_ARRAY_TYPE **wm;
  struct filemanager *fm;  /*  fasta/fastq file maneger  */
  struct outfile *output;

  
  cfg = command_line_parser (argc, argv);

  verbose_output = cfg->verbose;  /*  rely on external variable  */

  if ((cfg->flags & U_FLAG) != 0) {
    print_usage();
    exit(EXIT_SUCCESS);
  }
  if ((cfg->flags & (S_FLAG + C_FLAG)) != S_FLAG + C_FLAG) {
    printf("Missing REQUIRED files: SequenceFile or ConfigurationFile\n");
    print_usage();
    exit (EXIT_FAILURE);
  }
  
  if ((par_list = loadConfigFromFile (cfg->cvalue)) == NULL) {
    perror ("opening configuration file\n");
    exit (EXIT_FAILURE);
  }

  param_list_parser (par_list, cfg); /*  Update cfg  */
  wm = read_weights_matrix (par_list);
  free_paramList (par_list);   /*  From here it is no longer used  */
  if (wm == NULL) exit (EXIT_FAILURE);
  /*
    printf("NVALUE: %s\n", cfg->output_filename);
    printf("XVALUE: %d\n", cfg->xvalue);
    printf("FVALUE: %f\n", cfg->fvalue);
    printf("JVALUE: %d\n", cfg->jvalue);
    printf("GVALUE: %d\n", cfg->gvalue);
  */

  
  /*  Integrity check - must be moved in the appropr. place */
  if ((cfg->xvalue < cfg->nvalue) && (cfg->xvalue > 0)) {
    printf("max_length must be higher than min_length\n");
    print_usage();
    free_weights_matrix (wm);
    free (cfg);  /* free config params struct  */
    exit(EXIT_FAILURE);
  }

  if ((cfg->fvalue <= 0) || (cfg->fvalue > 1)) {
    printf("MinMatch parameter must be included in (0,1]\n");
    print_usage();
    free_weights_matrix (wm);
    free (cfg);  /* free config params struct  */
    exit(EXIT_FAILURE);
  }	
  /*  Initialize the file manager  */
  fm = filemanager_init (cfg->svalue);
  if (fm == NULL) {
    free_weights_matrix (wm);
    free (cfg);  /* free config params struct  */
    exit (EXIT_FAILURE);
  }

  pthread_mutex_init (&mtx_getseq, NULL);
  pthread_cond_init(&cond_read_seq_global, NULL);
  
  pthread_mutex_init (&mtx_wrtfile, NULL);
  pthread_cond_init(&cond_write_file_global, NULL);

    
  /*  Output file not specifyed  */
  if ((cfg->flags & O_FLAG) == 0)
    output = output_create (NULL);
  else 
    output = output_create (cfg->output_filename);

  print_header (output);
  
  /* Decide who is the next thread that can read a seq or write. Starts from the first thread.  */
  write_file_global = 0;
  read_seq_global = 0;

  /* threads run attempt */

  for (i = 0; i < cfg->thread_max; i++) {
    Dot_Thread_input* t_param;
    Thread_List_Elem* new_thread;
    /* I need i to handle write requests on file   */
    if ((t_param = dot_Thread_obj_init(cfg, wm, fm, output, i)) ==NULL) { 
      printf("Error in main() for t_param pointer\n"); 
      free_weights_matrix (wm);
      free (cfg);  /* free config params struct  */
      filemanager_destroy (fm);  /*  Releasefile managemant resources  */
      exit (EXIT_FAILURE); 
    }


    if ((err=pthread_create(&tid, NULL, &dot_thread_fn, (void*) t_param)) == 0 ) { 	
      if ((new_thread = (Thread_List_Elem*) malloc(sizeof (Thread_List_Elem))) == NULL) {
        perror("Error in insert_threadInList() for new_thread element\n");
        return 1;
      }
      new_thread->thread_id = tid;
      new_thread->thread_info = (void*) t_param;
      new_thread->next=NULL;
      insert_threadInList(&threads_list, &last_inThreadList, new_thread);
      /*printf("Thread %u inserted!\n",new_thread->thread_id);*/
    } 
    else { 
      check_Thread_Error(err); 
      destroy_dot_Thread_obj(&t_param);
    }
  }

  p=threads_list;
  while (p != NULL) {
    pthread_join (p->thread_id, (void*) &(p->status));
    if ( ((int)p->status) == 1 ) 
      printf ("Thread creation failure...the computation should still be safe\n");
    /* free the struct except for dot_matrix, file_manager, weight_matrix */
    search_param = (Dot_Thread_input*) p->thread_info;
    
    destroy_dot_Thread_obj(&search_param);
    /*printf("THREAD_STATUS: %d\n", p->status);*/
    p=p->next;
  }
  remove_thread (&threads_list, 0);
  
  output_destroy (output);
  filemanager_destroy  (fm);
  free_weights_matrix (wm);
  free (cfg);  /* free config params struct  */
  exit (EXIT_SUCCESS);
}
