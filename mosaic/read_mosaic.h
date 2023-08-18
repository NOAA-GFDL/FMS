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
#ifndef GET_GLOBAL_AREA_
#define GET_GLOBAL_AREA_

/* get global area */


#ifdef OVERLOAD_R4

float get_global_area(void);

#else

double get_global_area(void);

#endif

#ifdef OVERLOAD_R4

float get_global_area_(void);

#else

double get_global_area_(void);

#endif  /* OVERLOAD_R4 */

#endif
