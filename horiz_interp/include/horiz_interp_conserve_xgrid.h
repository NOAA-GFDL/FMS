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
