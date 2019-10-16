/***********************************************************************
 *                   GNU Lesser General Public License
 *
 * This file is part of the GFDL Flexible Modeling System (FMS).
 *
 * FMS is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * FMS is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
 **********************************************************************/
#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sched.h>
#include <errno.h>
#include <sys/resource.h>
#include <sys/syscall.h>

#ifndef __APPLE__
static pid_t gettid(void)
{
  return syscall(__NR_gettid);
}
#endif

/*
 * Returns this thread's CPU affinity, if bound to a single core,
 * or else -1.
 */
int get_cpu_affinity(void)
{
#ifndef __APPLE__
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

int get_cpu_affinity_(void) { return get_cpu_affinity(); }      /* Fortran interface */


/*
 * Returns this groups CPUSET
 * and also the CPUSET size or -1 (in case of a storage error)
 */
int get_cpuset(int fsz, int *output, int pe, _Bool debug)
{
#ifndef __APPLE__
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

int get_cpuset_(int *fsz, int *output, int *pe, _Bool *debug) { return get_cpuset(*fsz, output, *pe, *debug); } /* Fortran interface */


/*
 * Set CPU affinity to one core.
 */
int set_cpu_affinity(int cpu)
{
#ifndef __APPLE__
  cpu_set_t coremask; /* core affinity mask */

  CPU_ZERO(&coremask);
  CPU_SET(cpu,&coremask);
  if (sched_setaffinity(gettid(),sizeof(cpu_set_t),&coremask) != 0) {
    return -1;
  }
#endif
  return 0;
}

int set_cpu_affinity_(int *cpu) { return set_cpu_affinity(*cpu); } /* Fortran interface */
