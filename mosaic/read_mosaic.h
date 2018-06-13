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
#ifndef READ_MOSAIC_H_
#define READ_MOSAIC_H_

/* netcdf helpers */
/* perhaps should consider making static, or breaking out into seperate file,
   some of these names (field_exist) could pollute namespace... */

void handle_netcdf_error(const char *msg, int status );

void get_file_dir(const char *file, char *dir);

int field_exist(const char* file, const char *name);

int get_dimlen(const char* file, const char *name);

void get_string_data_level(const char *file, const char *name, char *data, const unsigned int* level);

void get_var_data(const char *file, const char *name, void *data);

void get_var_data_region(const char *file, const char *name, const size_t *start, const size_t *nread, void *data);

void get_string_data(const char *file, const char *name, char *data);

void get_var_text_att(const char *file, const char *name, const char *attname, char *att);
/* end netcdf helpers */

int read_mosaic_xgrid_size( const char *xgrid_file );

#ifdef OVERLOAD_R4

void read_mosaic_xgrid_order1(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2, float *area );

void read_mosaic_xgrid_order1_region(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2, float *area, int *isc, int *iec );

void read_mosaic_xgrid_order2(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2,
                              float *area, float *di, float *dj );

float get_global_area(void);

#else

void read_mosaic_xgrid_order1(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2, double *area );

void read_mosaic_xgrid_order1_region(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2, double *area, int *isc, int *iec );

void read_mosaic_xgrid_order2(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2,
                              double *area, double *di, double *dj );

double get_global_area(void);

#endif

int read_mosaic_ntiles(const char *mosaic_file);

int read_mosaic_ncontacts(const char *mosaic_file);

void read_mosaic_grid_sizes(const char *mosaic_file, int *nx, int *ny);

void read_mosaic_contact(const char *mosaic_file, int *tile1, int *tile2, int *istart1, int *iend1,
                         int *jstart1, int *jend1, int *istart2, int *iend2, int *jstart2, int *jend2);

int transfer_to_model_index(int istart_in, int iend_in, int *istart_out, int *iend_out, int refine_ratio);

void read_mosaic_grid_data(const char *mosaic_file, const char *name, int nx, int ny,
                           double *data, unsigned int level, int ioff, int joff);


#ifndef __AIX

void read_mosaic_contact_(const char *mosaic_file, int *tile1, int *tile2, int *istart1, int *iend1,
                          int *jstart1, int *jend1, int *istart2, int *iend2, int *jstart2, int *jend2);

int read_mosaic_xgrid_size_( const char *xgrid_file );

int read_mosaic_ntiles_(const char *mosaic_file);

int read_mosaic_ncontacts_(const char *mosaic_file);

void read_mosaic_grid_sizes_(const char *mosaic_file, int *nx, int *ny);

#ifdef OVERLOAD_R4

float get_global_area_(void);

void read_mosaic_xgrid_order1_(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2, float *area );

void read_mosaic_xgrid_order1_region_(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2, float *area, int *isc, int *iec );

void read_mosaic_xgrid_order2_(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2, float *area, float *di, float *dj );

#else

double get_global_area_(void);

void read_mosaic_xgrid_order1_(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2, double *area );

void read_mosaic_xgrid_order1_region_(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2, double *area, int *isc, int *iec );

void read_mosaic_xgrid_order2_(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2, double *area, double *di, double *dj );

#endif  /* OVERLOAD_R4 */

#endif  /* AIX */


#endif
