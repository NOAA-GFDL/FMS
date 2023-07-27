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
