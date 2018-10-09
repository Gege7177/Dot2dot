#ifndef THREADSMANAGING_THREADS_MANAGER_H_
#define THREADSMANAGING_THREADS_MANAGER_H_

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>
#include <sys/types.h>
#include <string.h>
#include <errno.h>

/**
*	Status = 0 -> Thread da eliminare; = 1 -> Thread attivo;
*	
*/
struct _thread_list_elem {

	pthread_t thread_id;
	int status;
	void* thread_info;
	struct _thread_list_elem* next;

};
typedef struct _thread_list_elem Thread_List_Elem;

Thread_List_Elem* threads_list;

pthread_mutex_t mtx_getseq;
pthread_cond_t cond_read_seq_global;
short int read_seq_global;

pthread_mutex_t mtx_wrtfile;
pthread_cond_t cond_write_file_global;
short int write_file_global;

int insert_threadInList(Thread_List_Elem** list_head, Thread_List_Elem** last_inList, Thread_List_Elem* new_thread);

Thread_List_Elem* remove_thread(Thread_List_Elem** list_head, pthread_t t_id);

int pt_m_lock (pthread_mutex_t *mtx);
int pt_m_unlock (pthread_mutex_t *mtx);
void check_Thread_Error(int error_code);

#endif /* THREADSMANAGING_THREADS_MANAGER_H_ */
