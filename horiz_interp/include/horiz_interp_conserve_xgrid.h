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
#ifndef HORIZ_INTERP_CREATE_XGRID_H_
#define HORIZ_INTERP_CREATE_XGRID_H_

#ifndef MAXXGRID
#define MAXXGRID 1e6
#endif

#define AREA_RATIO_THRESH (1.e-6)
#define MASK_THRESH (0.5)
#define MAX_V 8

int create_xgrid_1dx2d_order1(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out, const double *lon_in,
                              const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *i_out,
                              int *j_out, double *xgrid_area);

int create_xgrid_1dx2d_order1_(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
                               const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                               const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out, double *xgrid_area);

int create_xgrid_1dx2d_order2(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
                              const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
                              double *xgrid_area, double *xgrid_clon, double *xgrid_clat);

int create_xgrid_1dx2d_order2_(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
                               const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                               const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
                               double *xgrid_area, double *xgrid_clon, double *xgrid_clat);

int create_xgrid_2dx1d_order1(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out, const double *lon_in,
                              const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *i_out,
                              int *j_out, double *xgrid_area);

int create_xgrid_2dx1d_order1_(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
                               const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                               const double *mask_in, int *i_in, int *j_in, int *i_out,
                               int *j_out, double *xgrid_area);

int create_xgrid_2dx1d_order2(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
                              const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
                              double *xgrid_area, double *xgrid_clon, double *xgrid_clat);

int create_xgrid_2dx1d_order2_(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
                               const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                               const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
                               double *xgrid_area, double *xgrid_clon, double *xgrid_clat);

int create_xgrid_2dx2d_order1(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
                              const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *i_out,
                              int *j_out, double *xgrid_area);

int create_xgrid_2dx2d_order2(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
                              const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
                              double *xgrid_area, double *xgrid_clon, double *xgrid_clat);

int create_xgrid_great_circle(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
                              const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
                              double *xgrid_area, double *xgrid_clon, double *xgrid_clat);

int create_xgrid_1dx2d_order1_ug(const int *nlon_in, const int *nlat_in, const int *npts_out, const double *lon_in,
                              const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *l_out, double *xgrid_area);

int create_xgrid_great_circle_ug(const int *nlon_in, const int *nlat_in, const int *npts_out,
                               const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                               const double *mask_in, int *i_in, int *j_in, int *l_out,
                               double *xgrid_area, double *xgrid_clon, double *xgrid_clat);

int create_xgrid_2dx2d_order1_(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
                               const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                               const double *mask_in, int *i_in, int *j_in, int *i_out,
                               int *j_out, double *xgrid_area);

int create_xgrid_2dx2d_order2_(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
                               const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                               const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
                               double *xgrid_area, double *xgrid_clon, double *xgrid_clat);

int create_xgrid_1dx2d_order1_ug_(const int *nlon_in, const int *nlat_in, const int *npts_out, const double *lon_in,
                              const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *l_out, double *xgrid_area);

int create_xgrid_great_circle_ug_(const int *nlon_in, const int *nlat_in, const int *npts_out,
                               const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                               const double *mask_in, int *i_in, int *j_in, int *l_out,
                               double *xgrid_area, double *xgrid_clon, double *xgrid_clat);

int get_maxxgrid(void);
int get_maxxgrid_(void);

#endif
