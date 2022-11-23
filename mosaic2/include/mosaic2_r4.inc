!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
!> @file
!> @brief include file for mosaic22_mod to generate subroutines/functions for r4_kind arguments

!> @addtogroup mosaic2_mod
!> @{

#undef  FMS_MOS_KIND_
#define FMS_MOS_KIND_ r4_kind

#undef  GET_MOSAIC_XGRID_
#define GET_MOSAIC_XGRID_ get_mosaic_xgrid_r4

#undef  CALC_MOSAIC_GRID_AREA_
#define CALC_MOSAIC_GRID_AREA_ calc_mosaic_grid_area_r4

#undef  CALC_MOSAIC_GRID_GREAT_CIRCLE_AREA_
#define CALC_MOSAIC_GRID_GREAT_CIRCLE_AREA_ calc_mosaic_grid_great_circle_area_r4

#undef  IS_INSIDE_POLYGON_
#define IS_INSIDE_POLYGON_ is_inside_polygon_r4

#include <mosaic2.fh>
!> @}
