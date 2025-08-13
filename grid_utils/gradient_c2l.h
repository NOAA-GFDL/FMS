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
#ifndef GRADIENT_H_
#define GRADIENT_H_

void a2b_ord2(int nx, int ny, const double *qin, const double *edge_w, const double *edge_e,
              const double *edge_s, const double *edge_n, double *qout,
              int on_west_edge, int on_east_edge, int on_south_edge, int on_north_edge);

void grad_c2l(const int *nlon, const int *nlat, const double *pin, const double *dx, const double *dy, const double *area,
              const double *edge_w, const double *edge_e, const double *edge_s, const double *edge_n,
              const double *en_n, const double *en_e, const double *vlon, const double *vlat,
              double *grad_x, double *grad_y, const int *on_west_edge, const int *on_east_edge,
              const int *on_south_edge, const int *on_north_edge);

void grad_c2l_(const int *nlon, const int *nlat, const double *pin, const double *dx, const double *dy, const double *area,
               const double *edge_w, const double *edge_e, const double *edge_s, const double *edge_n,
               const double *en_n, const double *en_e, const double *vlon, const double *vlat,
               double *grad_x, double *grad_y, const int *on_west_edge, const int *on_east_edge,
               const int *on_south_edge, const int *on_north_edge);

void calc_c2l_grid_info(int *nx_pt, int *ny_pt, const double *xt, const double *yt, const double *xc, const double *yc,
                        double *dx, double *dy, double *area, double *edge_w, double *edge_e, double *edge_s,
                        double *edge_n, double *en_n, double *en_e, double *vlon, double *vlat,
                        int *on_west_edge, int *on_east_edge, int *on_south_edge, int *on_north_edge);

void calc_c2l_grid_info_(int *nx_pt, int *ny_pt, const double *xt, const double *yt, const double *xc, const double *yc,
                         double *dx, double *dy, double *area, double *edge_w, double *edge_e, double *edge_s,
                         double *edge_n, double *en_n, double *en_e, double *vlon, double *vlat,
                         int *on_west_edge, int *on_east_edge, int *on_south_edge, int *on_north_edge);

void get_edge(int nx, int ny, const double *lont, const double *latt,
              const double *lonc, const double *latc, double *edge_w, double *edge_e, double *edge_s, double *edge_n,
              int on_west_edge, int on_east_edge, int on_south_edge, int on_north_edge );

void mid_pt_sphere(const double *p1, const double *p2, double *pm);

void mid_pt3_cart(const double *p1, const double *p2, double *e);

#endif
