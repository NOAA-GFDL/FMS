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

// redefine functions for fortran usage (ie. append underscore)
// TODO proper fortran bind(C) routines
#define get_maxxgrid get_maxxgrid_
#define create_xgrid_1dx2d_order1 create_xgrid_1dx2d_order1_
#define create_xgrid_2dx1d_order1 create_xgrid_2dx1d_order1_ 
#define create_xgrid_2dx2d_order1 create_xgrid_2dx2d_order1_
#define create_xgrid_2dx2d_order2 create_xgrid_2dx2d_order2_
#define create_xgrid_great_circle create_xgrid_great_circle_

int create_xgrid_1dx2d_order1(const int nlon_in, const int nlat_in, const int nlon_out, const int nlat_out, const double *lon_in,
                              const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *i_out,
                              int *j_out, double *xgrid_area);

int create_xgrid_1dx2d_order2(const int nlon_in, const int nlat_in, const int nlon_out, const int nlat_out,
                              const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
                              double *xgrid_area, double *xgrid_clon, double *xgrid_clat);

int create_xgrid_2dx1d_order1(const int nlon_in, const int nlat_in, const int nlon_out, const int nlat_out, const double *lon_in,
                              const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *i_out,
                              int *j_out, double *xgrid_area);

int create_xgrid_2dx1d_order2(const int nlon_in, const int nlat_in, const int nlon_out, const int nlat_out,
                              const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
                              double *xgrid_area, double *xgrid_clon, double *xgrid_clat);

int create_xgrid_2dx2d_order1(const int nlon_input_cells, const int nlat_input_cells,
                              const int nlon_output_cells, const int nlat_output_cells,
                              const double *input_grid_lon, const double *input_grid_lat,
                              const double *output_grid_lon, const double *output_grid_lat,
                              const double *skip_input_cells, const double *skip_output_cells,
                              int *i_in, int *j_in, int *i_out,
                              int *j_out, double *xgrid_area);

int create_xgrid_2dx2d_order2(const int nlon_in, const int nlat_in, const int nlon_out, const int nlat_out,
                              const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
                              double *xgrid_area, double *xgrid_clon, double *xgrid_clat);

int create_xgrid_great_circle(const int nlon_in, const int nlat_in, const int nlon_out, const int nlat_out,
                              const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
                              double *xgrid_area, double *xgrid_clon, double *xgrid_clat);

int create_xgrid_1dx2d_order1_ug(const int nlon_in, const int nlat_in, const int npts_out, const double *lon_in,
                              const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *l_out, double *xgrid_area);

int create_xgrid_great_circle_ug(const int nlon_in, const int nlat_in, const int npts_out,
                               const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                               const double *mask_in, int *i_in, int *j_in, int *l_out,
                               double *xgrid_area, double *xgrid_clon, double *xgrid_clat);

int get_maxxgrid(void);

int block_setup(int ** i_in, int ** j_in, int ** i_out, int ** j_out,
                 double** xgrid_area,
                 int **pi_in, int **pj_in, int **pi_out, int **pj_out,
                 double **pxgrid_area,
                 int **istart2, int **iend2, int **pstart, int **pnxgrid,
                 const int nx2, const int ny2);

void compute_output_cell_bounds(const int nx_output_cells, const int ny_output_cells,
                                const int nx_output_points,
                                const double *output_grid_lon, const double *output_grid_lat,
                                double *lon_out_min_list, double *lon_out_max_list,
                                double *lat_out_min_list, double *lat_out_max_list,
                                double *lon_out_avg, int *n2_list,
                                double *lon_out_list, double *lat_out_list);


void clip_and_calc_2dx2d_order1(int curr_thread,
                                int nx_input_cells, int ny_input_cells, int nx_input_points,
                                const double *input_grid_lon, const double *input_grid_lat,
                                const double *skip_input_cells,
                                int nx_output_cells,
                                const double *lat_out_min_list, const double *lat_out_max_list,
                                const int *n2_list,
                                const double *lon_out_list, const double *lat_out_list,
                                const double *lon_out_min_list, const double *lon_out_max_list,
                                const double *lon_out_avg,
                                const double *area_in, const double *area_out,
                                double *pxgrid_area, int *pnxgrid, int *pi_in, int *pj_in, int *pi_out, int *pj_out, int *pstart,
                                const int *istart2, const int *iend2,
                                const double *skip_output_cells, int nthreads);

#endif
