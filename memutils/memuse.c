#ifdef __sgi
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>

#define RUSAGE_SELF      0         /* calling process */
#define RUSAGE_CHILDREN  -1        /* terminated child processes */

int memuse_(void)
{
 struct rusage my_rusage;
 int iret;

 my_rusage.ru_maxrss = 0;
 iret = getrusage(RUSAGE_SELF,&my_rusage);
 iret = (int) my_rusage.ru_maxrss;
 return iret;
 /* printf("Max memuse in KB is %d \n",iret); */
}
#endif
