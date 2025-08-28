/***********************************************************************
 *                             Apache License 2.0
 *
 * This file is part of the GFDL Flexible Modeling System (FMS).
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * FMS is distributed in the hope that it will be useful, but WITHOUT
 * WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the License for the specific language
 * governing permissions and limitations under the License.
 ***********************************************************************/

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <unistd.h>
#include <sched.h>
#include <errno.h>
#include <sys/resource.h>
#include <sys/syscall.h>
#ifdef __APPLE__
#include <pthread.h>
#endif
/** \endcond
 */
// skips doc parsing for includes and license

/**
 * \addtogroup affinity
 * \@{
 */

/**
 * gettid function for systems that do not have this function (e.g. on Mac OS.)
 */
#ifndef HAVE_GETTID
static pid_t gettid(void)
{
#if defined(__APPLE__)
  uint64_t tid64;
  pthread_threadid_np(NULL, &tid64);
  pid_t tid = (pid_t)tid64;
#else
  pid_t tid = syscall(__NR_gettid);
#endif
  return tid;
}
#endif

/**
 * Returns this thread's CPU affinity, if bound to a single core,
 * or else -1.
 */
int get_cpu_affinity(void)
{
#ifdef HAVE_SCHED_GETAFFINITY
  cpu_set_t coremask;           /* core affinity mask */

  CPU_ZERO(&coremask);
  if (sched_getaffinity(gettid(),sizeof(cpu_set_t),&coremask) != 0) {
    fprintf(stderr,"Unable to get thread %d affinity. %s\n",gettid(),strerror(errno));
  }

  int cpu;
  for (cpu=0;cpu < CPU_SETSIZE;cpu++) {
    if (CPU_ISSET(cpu,&coremask)) {
      return cpu;
    }
  }
#endif
  return -1;
}

/**
 * Returns this groups CPUSET
 * and also the CPUSET size or -1 (in case of a storage error)
 */
int get_cpuset(int fsz, int *output, int pe, _Bool debug)
{
#ifdef HAVE_SCHED_GETAFFINITY
  cpu_set_t coremask; /* core affinity mask */

  CPU_ZERO(&coremask);
  if (sched_getaffinity(gettid(),sizeof(cpu_set_t),&coremask) != 0) {
    fprintf(stderr,"Unable to get thread %d affinity. %s\n",gettid(),strerror(errno));
  }

  int  cpu;
  int  count;

  if (debug) {
    for (cpu=0;cpu < CPU_SETSIZE;cpu++) {
      if (CPU_ISSET(cpu,&coremask)) {
        printf("=> get_cpuset - pe %d: %d\n",pe, cpu);
      }
    }
  }

  count = 0;
  for (cpu=0;cpu < CPU_SETSIZE;cpu++) {
    if (CPU_ISSET(cpu,&coremask)) {
      if (count > fsz) {
        return -1;
      }
      output[count] = cpu;
      count ++;
    }
  }
  return count;
#else
  return fsz;
#endif
}

/**
 * Set CPU affinity to one core.
 */
int set_cpu_affinity(int cpu)
{
#ifdef HAVE_SCHED_GETAFFINITY
  cpu_set_t coremask; /* core affinity mask */

  CPU_ZERO(&coremask);
  CPU_SET(cpu,&coremask);
  if (sched_setaffinity(gettid(),sizeof(cpu_set_t),&coremask) != 0) {
    return -1;
  }
#endif
  return 0;
}
///@}
