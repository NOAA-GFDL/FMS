/************************************************************************
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
 ************************************************************************/

#include <stdio.h>

#define GETPEAKRSS_FC FC_FUNC (getpeakrss, GETPEAKRSS)
#ifdef __cplusplus
extern "C" /* prevent C++ name mangling */
#endif

size_t GETPEAKRSS_FC();

int main()
{
  size_t maxrss = GETPEAKRSS_FC();
  unsigned short return_val = 0;

  /* This test is mostly dubious, as maxrss can never be less than zero */
  if (maxrss <= 0) {
    return_val = 1;
  }

  return return_val;
}
