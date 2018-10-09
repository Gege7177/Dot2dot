
#ifndef MARK_TIME_C_
#define MARK_TIME_C_


#include <stdlib.h>
#include <sys/time.h>
#include <stdbool.h>  /*  Makes _Bool and true/false available  */

#define PRINT_TIME(str,x) printf("%s in: %.3f seconds\n", str,(float) (x)/1000000);

/*
 * This time function must be called 2 times.
 * First: to start the timer.
 * Second: to stop the timer and to get the value returned which means the time elapsed in milliseconds
 *
 */
long int Mark_Time (void);

#endif /* MARK_TIME_C_ */
