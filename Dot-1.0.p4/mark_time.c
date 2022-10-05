
#include "mark_time.h"

long int Mark_Time (void) {
   static struct timeval TimeStart;
   struct timeval TimeEnd;
   static char set = 0;
   if (set == 0) {
       set ++;
       gettimeofday (&TimeStart, NULL);
       return 0;
   }
   else {
       set = 0;
       gettimeofday (&TimeEnd, NULL);
       return ((TimeEnd.tv_sec - TimeStart.tv_sec) * 1000000) +
               (TimeEnd.tv_usec - TimeStart.tv_usec);
   }
}

