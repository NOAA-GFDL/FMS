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
/***********************************************************************
                      mosaic_util.h
    This header file provide some utilities routine that will be used in many tools.

    contact: Zhi.Liang@noaa.gov
***********************************************************************/
#ifndef GRID_UTILS_H_
#define GRID_UTILS_H_

#define TOLERANCE (1.e-6)
#ifndef RANGE_CHECK_CRITERIA
#define RANGE_CHECK_CRITERIA 0.05
#endif

#define MV 50

#define min(a,b) (a<b ? a:b)
#define max(a,b) (a>b ? a:b)
#define SMALL_VALUE ( 1.e-10 )

void error_handler(const char *msg);

int lon_fix(double *x, double *y, int n_in, double tlon);

double minval_double(int size, const double *data);

double maxval_double(int size, const double *data);

double avgval_double(int size, const double *data);

void latlon2xyz(int size, const double *lon, const double *lat, double *x, double *y, double *z);

void xyz2latlon(int size, const double *x, const double *y, const double *z, double *lon, double *lat);

int delete_vtx(double x[], double y[], int n, int n_del);

int insert_vtx(double x[], double y[], int n, int n_ins, double lon_in, double lat_in);

int fix_lon(double lon[], double lat[], int n, double tlon);

double great_circle_distance(double *p1, double *p2);

void vect_cross(const double *p1, const double *p2, double *e );

double spherical_angle(const double *v1, const double *v2, const double *v3);

double great_circle_area(int n, const double *x, const double *y, const double *z);

double * cross(const double *p1, const double *p2);

double dot(const double *p1, const double *p2);

void normalize_vect(double *e);

void unit_vect_latlon(int size, const double *lon, const double *lat, double *vlon, double *vlat);

int intersect_tri_with_line(const double *plane, const double *l1, const double *l2, double *p,
                            double *t);

int invert_matrix_3x3(long double m[], long double m_inv[]);

void mult(long double m[], long double v[], long double out_v[]);

double metric(const double *p);

int inside_a_polygon( double *lon1, double *lat1, int *npts, double *lon2, double *lat2);

int samePoint(double x1, double y1, double z1, double x2, double y2, double z2);

int inside_a_polygon_(double *lon1, double *lat1, int *npts, double *lon2, double *lat2);

int inside_edge(double x0, double y0, double x1, double y1, double x, double y);

int line_intersect_2D_3D(double *a1, double *a2, double *q1, double *q2, double *q3,
                         double *intersect, double *u_a, double *u_q, int *inbound);

double poly_ctrlon(const double lon[], const double lat[], int n, double clon);

double poly_ctrlat(const double lon[], const double lat[], int n);

int get_maxxgrid(void);

int get_maxxgrid_(void);

double get_global_area(void);

double get_global_area_(void);

double poly_area(const double lon[], const double lat[], int n);

double poly_area_dimensionless(const double x[], const double y[], int n);

double spherical_excess_area(const double* p_ll, const double* p_ul,
                             const double* p_lr, const double* p_ur, double radius);

void get_grid_area(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area);

void get_grid_great_circle_area(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area);

void get_grid_area_no_adjust(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area);

int clip(const double lon_in[], const double lat_in[], int n_in, double ll_lon, double ll_lat,
         double ur_lon, double ur_lat, double lon_out[], double lat_out[]);

int clip_2dx2d(const double lon1_in[], const double lat1_in[], int n1_in,
               const double lon2_in[], const double lat2_in[], int n2_in,
               double lon_out[], double lat_out[]);

int clip_2dx2d_great_circle(const double x1_in[], const double y1_in[], const double z1_in[], int n1_in,
                            const double x2_in[], const double y2_in[], const double z2_in [], int n2_in,
                            double x_out[], double y_out[], double z_out[]);

void get_grid_area_ug(const int *npts, const double *lon, const double *lat, double *area);

void get_grid_great_circle_area_ug(const int *npts, const double *lon, const double *lat, double *area);

void get_grid_area_(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area);

void get_grid_great_circle_area_(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area);

void get_grid_area_ug_(const int *npts, const double *lon, const double *lat, double *area);

void get_grid_great_circle_area_ug_(const int *npts, const double *lon, const double *lat, double *area);

#endif
