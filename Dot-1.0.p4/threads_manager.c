#include "threads_manager.h"

int insert_threadInList(Thread_List_Elem** list_head, Thread_List_Elem** last_inList, Thread_List_Elem* new_thread) {	
	Thread_List_Elem* curr;
	
	if (*list_head == NULL) {
		*list_head=new_thread;		
		if (last_inList != NULL) *last_inList=new_thread;
	} else {
		if (last_inList != NULL) curr = *last_inList; else curr = *list_head;
		while (curr->next != NULL) {
			curr = curr->next;
		}
		curr->next = new_thread;
		if (last_inList != NULL) *last_inList = new_thread;
	}
	
	return 0;
}

Thread_List_Elem* remove_thread(Thread_List_Elem** list_head, pthread_t t_id) {

	if (t_id != 0) { /* remove t_id thread */
		if (*list_head == NULL) return NULL;
		else {
			if ((*list_head)->thread_id == t_id) {
				Thread_List_Elem* succ = (*list_head)->next;
				free(*list_head);
				/* printf("thread %u removed\n", t_id); */
				*list_head = succ;
				return *list_head;
			}
			else {
				(*list_head)->next = remove_thread(&((*list_head)->next), t_id);
				return *list_head;
			}
		}
	} else { /* Remove all threads */
		if (*list_head == NULL) return NULL;
		else {
			(*list_head)->next = remove_thread(&((*list_head)->next), t_id);
			/* printf("thread %u removed\n", (*list_head)->thread_id); */
			free(*list_head);
			return NULL;
		}
	}
}

int pt_m_lock (pthread_mutex_t *mtx) {
  int err;
  if ((err = pthread_mutex_lock(mtx)) != 0) {
    perror ("lock\n");
    pthread_exit (NULL);
  } 
  return err;
}

int pt_m_unlock (pthread_mutex_t *mtx) {
  int err;
  if ((err=pthread_mutex_unlock(mtx)) != 0) {
    perror("unlock\n");
    pthread_exit (NULL);
  }
  return err;
}

void check_Thread_Error(int error_code) {
  switch (error_code) {
  case (EAGAIN) : {
    printf ("System lacked the necessary memory resource!\n");
    break;
  }
  case (EINVAL) : {
    printf ("Attribute Error!\n");
    break;
  }
  default : {
    printf ("Unknown Error!\n");
    break;
  }
  }
}

