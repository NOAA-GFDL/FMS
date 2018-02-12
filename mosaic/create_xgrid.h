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
#ifndef CREATE_XGRID_H_
#define CREATE_XGRID_H_

#ifndef MAXXGRID
#define MAXXGRID 1e6
#endif

#define MV 50
/* this value is small compare to earth area */

double grid_box_radius(const double *x, const double *y, const double *z, int n);

double dist_between_boxes(const double *x1, const double *y1, const double *z1, int n1,
                          const double *x2, const double *y2, const double *z2, int n2);

int inside_edge(double x0, double y0, double x1, double y1, double x, double y);

int line_intersect_2D_3D(double *a1, double *a2, double *q1, double *q2, double *q3,
                         double *intersect, double *u_a, double *u_q, int *inbound);

double poly_ctrlon(const double lon[], const double lat[], int n, double clon);

double poly_ctrlat(const double lon[], const double lat[], int n);

double box_ctrlon(double ll_lon, double ll_lat, double ur_lon, double ur_lat, double clon);

double box_ctrlat(double ll_lon, double ll_lat, double ur_lon, double ur_lat);

int get_maxxgrid(void);

int get_maxxgrid_(void);

void get_grid_area(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area);

void get_grid_great_circle_area(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area);

void get_grid_area_dimensionless(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area);

void get_grid_area_no_adjust(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area);

int clip(const double lon_in[], const double lat_in[], int n_in, double ll_lon, double ll_lat,
         double ur_lon, double ur_lat, double lon_out[], double lat_out[]);

int clip_2dx2d(const double lon1_in[], const double lat1_in[], int n1_in,
               const double lon2_in[], const double lat2_in[], int n2_in,
               double lon_out[], double lat_out[]);

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

int clip_2dx2d_great_circle(const double x1_in[], const double y1_in[], const double z1_in[], int n1_in,
                            const double x2_in[], const double y2_in[], const double z2_in [], int n2_in,
                            double x_out[], double y_out[], double z_out[]);

int create_xgrid_great_circle(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
                              const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
                              double *xgrid_area, double *xgrid_clon, double *xgrid_clat);

void get_grid_area_ug(const int *npts, const double *lon, const double *lat, double *area);
int create_xgrid_1dx2d_order1_ug(const int *nlon_in, const int *nlat_in, const int *npts_out, const double *lon_in,
                              const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *l_out, double *xgrid_area);
void get_grid_great_circle_area_ug(const int *npts, const double *lon, const double *lat, double *area);
int create_xgrid_great_circle_ug(const int *nlon_in, const int *nlat_in, const int *npts_out,
                               const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                               const double *mask_in, int *i_in, int *j_in, int *l_out,
                               double *xgrid_area, double *xgrid_clon, double *xgrid_clat);

#ifndef __AIX
void get_grid_area_(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area);

void get_grid_great_circle_area_(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area);

int create_xgrid_2dx2d_order1_(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
                               const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                               const double *mask_in, int *i_in, int *j_in, int *i_out,
                               int *j_out, double *xgrid_area);

int create_xgrid_2dx2d_order2_(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
                               const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                               const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
                               double *xgrid_area, double *xgrid_clon, double *xgrid_clat);
void get_grid_area_ug_(const int *npts, const double *lon, const double *lat, double *area);
int create_xgrid_1dx2d_order1_ug_(const int *nlon_in, const int *nlat_in, const int *npts_out, const double *lon_in,
                              const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *l_out, double *xgrid_area);
void get_grid_great_circle_area_ug_(const int *npts, const double *lon, const double *lat, double *area);
int create_xgrid_great_circle_ug_(const int *nlon_in, const int *nlat_in, const int *npts_out,
                               const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                               const double *mask_in, int *i_in, int *j_in, int *l_out,
                               double *xgrid_area, double *xgrid_clon, double *xgrid_clat);

#endif

#endif
