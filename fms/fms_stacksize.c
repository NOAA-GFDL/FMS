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

#include <sys/resource.h>

/*
 * Set the stack size limit to its maximum permissible value
 */

void maximize_system_stacksize_limit()
{
  struct rlimit stacksize;

  getrlimit(RLIMIT_STACK, &stacksize);
  stacksize.rlim_cur = stacksize.rlim_max;
  setrlimit(RLIMIT_STACK, &stacksize);
}
