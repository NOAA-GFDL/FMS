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
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "grid_utils.h"
#include "tree_utils.h"
#include "horiz_interp_conserve_xgrid.h"
#include "constant.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

/** \file
 *  \ingroup mosaic
 *  \brief Grid creation and calculation functions for use in @ref mosaic_mod
 * /

/*******************************************************************************
  void create_xgrid_1dx2d_order1
  This routine generate exchange grids between two grids for the first order
  conservative interpolation. nlon_in,nlat_in,nlon_out,nlat_out are the size of the grid cell
  and lon_in,lat_in are 1-D grid bounds, lon_out,lat_out are geographic grid location of grid cell bounds.
*******************************************************************************/
int create_xgrid_1dx2d_order1_(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
                               const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                               const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out, double *xgrid_area)
{
  int nxgrid;

  nxgrid = create_xgrid_1dx2d_order1(nlon_in, nlat_in, nlon_out, nlat_out, lon_in, lat_in, lon_out, lat_out, mask_in,
                                     i_in, j_in, i_out, j_out, xgrid_area);
  return nxgrid;

}

int create_xgrid_1dx2d_order1(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out, const double *lon_in,
                              const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *i_out,
                              int *j_out, double *xgrid_area)
{

  int nx1, ny1, nx2, ny2, nx1p, nx2p;
  int i1, j1, i2, j2, nxgrid;
  double ll_lon, ll_lat, ur_lon, ur_lat, x_in[MV], y_in[MV], x_out[MV], y_out[MV];
  double *area_in, *area_out, min_area;
  double *tmpx, *tmpy;

  nx1 = *nlon_in;
  ny1 = *nlat_in;
  nx2 = *nlon_out;
  ny2 = *nlat_out;

  nxgrid = 0;
  nx1p = nx1 + 1;
  nx2p = nx2 + 1;

  area_in = (double *)malloc(nx1*ny1*sizeof(double));
  area_out = (double *)malloc(nx2*ny2*sizeof(double));
  tmpx = (double *)malloc((nx1+1)*(ny1+1)*sizeof(double));
  tmpy = (double *)malloc((nx1+1)*(ny1+1)*sizeof(double));
  for(j1=0; j1<=ny1; j1++) for(i1=0; i1<=nx1; i1++) {
      tmpx[j1*nx1p+i1] = lon_in[i1];
      tmpy[j1*nx1p+i1] = lat_in[j1];
    }
  /* This is just a temporary fix to solve the issue that there is one point in zonal direction */
  if(nx1 > 1)
    get_grid_area(nlon_in, nlat_in, tmpx, tmpy, area_in);
  else
    get_grid_area_no_adjust(nlon_in, nlat_in, tmpx, tmpy, area_in);

  get_grid_area(nlon_out, nlat_out, lon_out, lat_out, area_out);
  free(tmpx);
  free(tmpy);

  for(j1=0; j1<ny1; j1++) for(i1=0; i1<nx1; i1++) if( mask_in[j1*nx1+i1] > MASK_THRESH ) {

        ll_lon = lon_in[i1];   ll_lat = lat_in[j1];
        ur_lon = lon_in[i1+1]; ur_lat = lat_in[j1+1];
        for(j2=0; j2<ny2; j2++) for(i2=0; i2<nx2; i2++) {
            int n_in, n_out;
            double Xarea;

            y_in[0] = lat_out[j2*nx2p+i2];
            y_in[1] = lat_out[j2*nx2p+i2+1];
            y_in[2] = lat_out[(j2+1)*nx2p+i2+1];
            y_in[3] = lat_out[(j2+1)*nx2p+i2];
            if (  (y_in[0]<=ll_lat) && (y_in[1]<=ll_lat)
                  && (y_in[2]<=ll_lat) && (y_in[3]<=ll_lat) ) continue;
            if (  (y_in[0]>=ur_lat) && (y_in[1]>=ur_lat)
                  && (y_in[2]>=ur_lat) && (y_in[3]>=ur_lat) ) continue;

            x_in[0] = lon_out[j2*nx2p+i2];
            x_in[1] = lon_out[j2*nx2p+i2+1];
            x_in[2] = lon_out[(j2+1)*nx2p+i2+1];
            x_in[3] = lon_out[(j2+1)*nx2p+i2];
            n_in = fix_lon(x_in, y_in, 4, (ll_lon+ur_lon)/2);

            if ( (n_out = clip ( x_in, y_in, n_in, ll_lon, ll_lat, ur_lon, ur_lat, x_out, y_out )) > 0 ) {
              Xarea = poly_area (x_out, y_out, n_out ) * mask_in[j1*nx1+i1];
              min_area = min(area_in[j1*nx1+i1], area_out[j2*nx2+i2]);
              if( Xarea/min_area > AREA_RATIO_THRESH ) {
                xgrid_area[nxgrid] = Xarea;
                i_in[nxgrid]    = i1;
                j_in[nxgrid]    = j1;
                i_out[nxgrid]   = i2;
                j_out[nxgrid]   = j2;
                ++nxgrid;
                if(nxgrid > MAXXGRID) error_handler("nxgrid is greater than MAXXGRID, increase MAXXGRID");
              }
            }
          }
      }

  free(area_in);
  free(area_out);

  return nxgrid;

} /* create_xgrid_1dx2d_order1 */


/*******************************************************************************
  void create_xgrid_1dx2d_order1_ug
  This routine generate exchange grids between two grids for the first order
  conservative interpolation. nlon_in,nlat_in,nlon_out,nlat_out are the size of the grid cell
  and lon_in,lat_in are 1-D grid bounds, lon_out,lat_out are geographic grid location of grid cell bounds.
*******************************************************************************/
int create_xgrid_1dx2d_order1_ug_(const int *nlon_in, const int *nlat_in, const int *npts_out,
                               const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                               const double *mask_in, int *i_in, int *j_in, int *l_out, double *xgrid_area)
{
  int nxgrid;

  nxgrid = create_xgrid_1dx2d_order1_ug(nlon_in, nlat_in, npts_out, lon_in, lat_in, lon_out, lat_out, mask_in,
                                     i_in, j_in, l_out, xgrid_area);
  return nxgrid;

}

int create_xgrid_1dx2d_order1_ug(const int *nlon_in, const int *nlat_in, const int *npts_out, const double *lon_in,
                              const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *l_out, double *xgrid_area)
{

  int nx1, ny1, nx1p, nv, npts2;
  int i1, j1, l2, nxgrid;
  double ll_lon, ll_lat, ur_lon, ur_lat, x_in[MV], y_in[MV], x_out[MV], y_out[MV];
  double *area_in, *area_out, min_area;
  double *tmpx, *tmpy;

  nx1 = *nlon_in;
  ny1 = *nlat_in;
  nv  = 4;
  npts2 = *npts_out;

  nxgrid = 0;
  nx1p = nx1 + 1;

  area_in = (double *)malloc(nx1*ny1*sizeof(double));
  area_out = (double *)malloc(npts2*sizeof(double));
  tmpx = (double *)malloc((nx1+1)*(ny1+1)*sizeof(double));
  tmpy = (double *)malloc((nx1+1)*(ny1+1)*sizeof(double));
  for(j1=0; j1<=ny1; j1++) for(i1=0; i1<=nx1; i1++) {
      tmpx[j1*nx1p+i1] = lon_in[i1];
      tmpy[j1*nx1p+i1] = lat_in[j1];
    }
  /* This is just a temporary fix to solve the issue that there is one point in zonal direction */
  if(nx1 > 1)
    get_grid_area(nlon_in, nlat_in, tmpx, tmpy, area_in);
  else
    get_grid_area_no_adjust(nlon_in, nlat_in, tmpx, tmpy, area_in);

  get_grid_area_ug(npts_out, lon_out, lat_out, area_out);
  free(tmpx);
  free(tmpy);

  for(j1=0; j1<ny1; j1++) for(i1=0; i1<nx1; i1++) if( mask_in[j1*nx1+i1] > MASK_THRESH ) {

        ll_lon = lon_in[i1];   ll_lat = lat_in[j1];
        ur_lon = lon_in[i1+1]; ur_lat = lat_in[j1+1];
        for(l2=0; l2<npts2; l2++) {
            int n_in, n_out;
            double Xarea;

            y_in[0] = lat_out[l2*nv];
            y_in[1] = lat_out[l2*nv+1];
            y_in[2] = lat_out[l2*nv+2];
            y_in[3] = lat_out[l2*nv+3];
            if (  (y_in[0]<=ll_lat) && (y_in[1]<=ll_lat)
                  && (y_in[2]<=ll_lat) && (y_in[3]<=ll_lat) ) continue;
            if (  (y_in[0]>=ur_lat) && (y_in[1]>=ur_lat)
                  && (y_in[2]>=ur_lat) && (y_in[3]>=ur_lat) ) continue;

            x_in[0] = lon_out[l2*nv];
            x_in[1] = lon_out[l2*nv+1];
            x_in[2] = lon_out[l2*nv+2];
            x_in[3] = lon_out[l2*nv+3];
            n_in = fix_lon(x_in, y_in, 4, (ll_lon+ur_lon)/2);

            if ( (n_out = clip ( x_in, y_in, n_in, ll_lon, ll_lat, ur_lon, ur_lat, x_out, y_out )) > 0 ) {
              Xarea = poly_area (x_out, y_out, n_out ) * mask_in[j1*nx1+i1];
              min_area = min(area_in[j1*nx1+i1], area_out[l2]);
              if( Xarea/min_area > AREA_RATIO_THRESH ) {
                xgrid_area[nxgrid] = Xarea;
                i_in[nxgrid]    = i1;
                j_in[nxgrid]    = j1;
                l_out[nxgrid]   = l2;
                ++nxgrid;
                if(nxgrid > MAXXGRID) error_handler("nxgrid is greater than MAXXGRID, increase MAXXGRID");
              }
            }
          }
      }

  free(area_in);
  free(area_out);

  return nxgrid;

} /* create_xgrid_1dx2d_order1_ug */

/********************************************************************************
  void create_xgrid_1dx2d_order2
  This routine generate exchange grids between two grids for the second order
  conservative interpolation. nlon_in,nlat_in,nlon_out,nlat_out are the size of the grid cell
  and lon_in,lat_in are 1-D grid bounds, lon_out,lat_out are geographic grid location of grid cell bounds.
********************************************************************************/
int create_xgrid_1dx2d_order2_(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
                               const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                               const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
                               double *xgrid_area, double *xgrid_clon, double *xgrid_clat)
{
  int nxgrid;
  nxgrid = create_xgrid_1dx2d_order2(nlon_in, nlat_in, nlon_out, nlat_out, lon_in, lat_in, lon_out, lat_out, mask_in, i_in,
                                     j_in, i_out, j_out, xgrid_area, xgrid_clon, xgrid_clat);
  return nxgrid;

}
int create_xgrid_1dx2d_order2(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
                              const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
                              double *xgrid_area, double *xgrid_clon, double *xgrid_clat)
{

  int nx1, ny1, nx2, ny2, nx1p, nx2p;
  int i1, j1, i2, j2, nxgrid;
  double ll_lon, ll_lat, ur_lon, ur_lat, x_in[MV], y_in[MV], x_out[MV], y_out[MV];
  double *area_in, *area_out, min_area;
  double *tmpx, *tmpy;

  nx1 = *nlon_in;
  ny1 = *nlat_in;
  nx2 = *nlon_out;
  ny2 = *nlat_out;

  nxgrid = 0;
  nx1p = nx1 + 1;
  nx2p = nx2 + 1;

  area_in      = (double *)malloc(nx1*ny1*sizeof(double));
  area_out     = (double *)malloc(nx2*ny2*sizeof(double));
  tmpx = (double *)malloc((nx1+1)*(ny1+1)*sizeof(double));
  tmpy = (double *)malloc((nx1+1)*(ny1+1)*sizeof(double));
  for(j1=0; j1<=ny1; j1++) for(i1=0; i1<=nx1; i1++) {
      tmpx[j1*nx1p+i1] = lon_in[i1];
      tmpy[j1*nx1p+i1] = lat_in[j1];
    }
  get_grid_area(nlon_in, nlat_in, tmpx, tmpy, area_in);
  get_grid_area(nlon_out, nlat_out, lon_out, lat_out, area_out);
  free(tmpx);
  free(tmpy);

  for(j1=0; j1<ny1; j1++) for(i1=0; i1<nx1; i1++) if( mask_in[j1*nx1+i1] > MASK_THRESH ) {

        ll_lon = lon_in[i1];   ll_lat = lat_in[j1];
        ur_lon = lon_in[i1+1]; ur_lat = lat_in[j1+1];
        for(j2=0; j2<ny2; j2++) for(i2=0; i2<nx2; i2++) {
            int n_in, n_out;
            double xarea, lon_in_avg;

            y_in[0] = lat_out[j2*nx2p+i2];
            y_in[1] = lat_out[j2*nx2p+i2+1];
            y_in[2] = lat_out[(j2+1)*nx2p+i2+1];
            y_in[3] = lat_out[(j2+1)*nx2p+i2];
            if (  (y_in[0]<=ll_lat) && (y_in[1]<=ll_lat)
                  && (y_in[2]<=ll_lat) && (y_in[3]<=ll_lat) ) continue;
            if (  (y_in[0]>=ur_lat) && (y_in[1]>=ur_lat)
                  && (y_in[2]>=ur_lat) && (y_in[3]>=ur_lat) ) continue;

            x_in[0] = lon_out[j2*nx2p+i2];
            x_in[1] = lon_out[j2*nx2p+i2+1];
            x_in[2] = lon_out[(j2+1)*nx2p+i2+1];
            x_in[3] = lon_out[(j2+1)*nx2p+i2];
            n_in = fix_lon(x_in, y_in, 4, (ll_lon+ur_lon)/2);
            lon_in_avg = avgval_double(n_in, x_in);

            if (  (n_out = clip ( x_in, y_in, n_in, ll_lon, ll_lat, ur_lon, ur_lat, x_out, y_out )) > 0 ) {
              xarea = poly_area (x_out, y_out, n_out ) * mask_in[j1*nx1+i1];
              min_area = min(area_in[j1*nx1+i1], area_out[j2*nx2+i2]);
              if(xarea/min_area > AREA_RATIO_THRESH ) {
                xgrid_area[nxgrid] = xarea;
                xgrid_clon[nxgrid] = poly_ctrlon(x_out, y_out, n_out, lon_in_avg);
                xgrid_clat[nxgrid] = poly_ctrlat (x_out, y_out, n_out );
                i_in[nxgrid]    = i1;
                j_in[nxgrid]    = j1;
                i_out[nxgrid]   = i2;
                j_out[nxgrid]   = j2;
                ++nxgrid;
                if(nxgrid > MAXXGRID) error_handler("nxgrid is greater than MAXXGRID, increase MAXXGRID");
              }
            }
          }
      }
  free(area_in);
  free(area_out);

  return nxgrid;

} /* create_xgrid_1dx2d_order2 */

/*******************************************************************************
  void create_xgrid_2dx1d_order1
  This routine generate exchange grids between two grids for the first order
  conservative interpolation. nlon_in,nlat_in,nlon_out,nlat_out are the size of the grid cell
  and lon_out,lat_out are 1-D grid bounds, lon_in,lat_in are geographic grid location of grid cell bounds.
  mask is on grid lon_in/lat_in.
*******************************************************************************/
int create_xgrid_2dx1d_order1_(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
                               const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                               const double *mask_in, int *i_in, int *j_in, int *i_out,
                               int *j_out, double *xgrid_area)
{
  int nxgrid;

  nxgrid = create_xgrid_2dx1d_order1(nlon_in, nlat_in, nlon_out, nlat_out, lon_in, lat_in, lon_out, lat_out, mask_in,
                                     i_in, j_in, i_out, j_out, xgrid_area);
  return nxgrid;

}
int create_xgrid_2dx1d_order1(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out, const double *lon_in,
                              const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *i_out,
                              int *j_out, double *xgrid_area)
{

  int nx1, ny1, nx2, ny2, nx1p, nx2p;
  int i1, j1, i2, j2, nxgrid;
  double ll_lon, ll_lat, ur_lon, ur_lat, x_in[MV], y_in[MV], x_out[MV], y_out[MV];
  double *area_in, *area_out, min_area;
  double *tmpx, *tmpy;
  int n_in, n_out;
  double Xarea;


  nx1 = *nlon_in;
  ny1 = *nlat_in;
  nx2 = *nlon_out;
  ny2 = *nlat_out;

  nxgrid = 0;
  nx1p = nx1 + 1;
  nx2p = nx2 + 1;
  area_in = (double *)malloc(nx1*ny1*sizeof(double));
  area_out = (double *)malloc(nx2*ny2*sizeof(double));
  tmpx = (double *)malloc((nx2+1)*(ny2+1)*sizeof(double));
  tmpy = (double *)malloc((nx2+1)*(ny2+1)*sizeof(double));
  for(j2=0; j2<=ny2; j2++) for(i2=0; i2<=nx2; i2++) {
      tmpx[j2*nx2p+i2] = lon_out[i2];
      tmpy[j2*nx2p+i2] = lat_out[j2];
    }
  get_grid_area(nlon_in, nlat_in, lon_in, lat_in, area_in);
  get_grid_area(nlon_out, nlat_out, tmpx, tmpy, area_out);

  free(tmpx);
  free(tmpy);

  for(j2=0; j2<ny2; j2++) for(i2=0; i2<nx2; i2++) {

      ll_lon = lon_out[i2];   ll_lat = lat_out[j2];
      ur_lon = lon_out[i2+1]; ur_lat = lat_out[j2+1];
      for(j1=0; j1<ny1; j1++) for(i1=0; i1<nx1; i1++) if( mask_in[j1*nx1+i1] > MASK_THRESH ) {

            y_in[0] = lat_in[j1*nx1p+i1];
            y_in[1] = lat_in[j1*nx1p+i1+1];
            y_in[2] = lat_in[(j1+1)*nx1p+i1+1];
            y_in[3] = lat_in[(j1+1)*nx1p+i1];
            if (  (y_in[0]<=ll_lat) && (y_in[1]<=ll_lat)
                  && (y_in[2]<=ll_lat) && (y_in[3]<=ll_lat) ) continue;
            if (  (y_in[0]>=ur_lat) && (y_in[1]>=ur_lat)
                  && (y_in[2]>=ur_lat) && (y_in[3]>=ur_lat) ) continue;

            x_in[0] = lon_in[j1*nx1p+i1];
            x_in[1] = lon_in[j1*nx1p+i1+1];
            x_in[2] = lon_in[(j1+1)*nx1p+i1+1];
            x_in[3] = lon_in[(j1+1)*nx1p+i1];

            n_in = fix_lon(x_in, y_in, 4, (ll_lon+ur_lon)/2);

            if ( (n_out = clip ( x_in, y_in, n_in, ll_lon, ll_lat, ur_lon, ur_lat, x_out, y_out )) > 0 ) {
              Xarea = poly_area ( x_out, y_out, n_out ) * mask_in[j1*nx1+i1];
              min_area = min(area_in[j1*nx1+i1], area_out[j2*nx2+i2]);
              if( Xarea/min_area > AREA_RATIO_THRESH ) {
                xgrid_area[nxgrid] = Xarea;
                i_in[nxgrid]    = i1;
                j_in[nxgrid]    = j1;
                i_out[nxgrid]   = i2;
                j_out[nxgrid]   = j2;
                ++nxgrid;
                if(nxgrid > MAXXGRID) error_handler("nxgrid is greater than MAXXGRID, increase MAXXGRID");
              }
            }
          }
    }

  free(area_in);
  free(area_out);

  return nxgrid;

} /* create_xgrid_2dx1d_order1 */


/********************************************************************************
  void create_xgrid_2dx1d_order2
  This routine generate exchange grids between two grids for the second order
  conservative interpolation. nlon_in,nlat_in,nlon_out,nlat_out are the size of the grid cell
  and lon_out,lat_out are 1-D grid bounds, lon_in,lat_in are geographic grid location of grid cell bounds.
  mask is on grid lon_in/lat_in.
********************************************************************************/
int create_xgrid_2dx1d_order2_(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
                               const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                               const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
                               double *xgrid_area, double *xgrid_clon, double *xgrid_clat)
{
  int nxgrid;
  nxgrid = create_xgrid_2dx1d_order2(nlon_in, nlat_in, nlon_out, nlat_out, lon_in, lat_in, lon_out, lat_out, mask_in, i_in,
                                     j_in, i_out, j_out, xgrid_area, xgrid_clon, xgrid_clat);
  return nxgrid;

}

int create_xgrid_2dx1d_order2(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
                              const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
                              double *xgrid_area, double *xgrid_clon, double *xgrid_clat)
{

  int nx1, ny1, nx2, ny2, nx1p, nx2p;
  int i1, j1, i2, j2, nxgrid;
  double ll_lon, ll_lat, ur_lon, ur_lat, x_in[MV], y_in[MV], x_out[MV], y_out[MV];
  double *tmpx, *tmpy;
  double *area_in, *area_out, min_area;
  double  lon_in_avg;
  int n_in, n_out;
  double xarea;


  nx1 = *nlon_in;
  ny1 = *nlat_in;
  nx2 = *nlon_out;
  ny2 = *nlat_out;

  nxgrid = 0;
  nx1p = nx1 + 1;
  nx2p = nx2 + 1;

  area_in      = (double *)malloc(nx1*ny1*sizeof(double));
  area_out     = (double *)malloc(nx2*ny2*sizeof(double));
  tmpx = (double *)malloc((nx2+1)*(ny2+1)*sizeof(double));
  tmpy = (double *)malloc((nx2+1)*(ny2+1)*sizeof(double));
  for(j2=0; j2<=ny2; j2++) for(i2=0; i2<=nx2; i2++) {
      tmpx[j2*nx2p+i2] = lon_out[i2];
      tmpy[j2*nx2p+i2] = lat_out[j2];
    }
  get_grid_area(nlon_in, nlat_in, lon_in, lat_in, area_in);
  get_grid_area(nlon_out, nlat_out, tmpx, tmpy, area_out);

  free(tmpx);
  free(tmpy);

  for(j2=0; j2<ny2; j2++) for(i2=0; i2<nx2; i2++) {

      ll_lon = lon_out[i2];   ll_lat = lat_out[j2];
      ur_lon = lon_out[i2+1]; ur_lat = lat_out[j2+1];
      for(j1=0; j1<ny1; j1++) for(i1=0; i1<nx1; i1++) if( mask_in[j1*nx1+i1] > MASK_THRESH ) {

            y_in[0] = lat_in[j1*nx1p+i1];
            y_in[1] = lat_in[j1*nx1p+i1+1];
            y_in[2] = lat_in[(j1+1)*nx1p+i1+1];
            y_in[3] = lat_in[(j1+1)*nx1p+i1];
            if (  (y_in[0]<=ll_lat) && (y_in[1]<=ll_lat)
                  && (y_in[2]<=ll_lat) && (y_in[3]<=ll_lat) ) continue;
            if (  (y_in[0]>=ur_lat) && (y_in[1]>=ur_lat)
                  && (y_in[2]>=ur_lat) && (y_in[3]>=ur_lat) ) continue;

            x_in[0] = lon_in[j1*nx1p+i1];
            x_in[1] = lon_in[j1*nx1p+i1+1];
            x_in[2] = lon_in[(j1+1)*nx1p+i1+1];
            x_in[3] = lon_in[(j1+1)*nx1p+i1];

            n_in = fix_lon(x_in, y_in, 4, (ll_lon+ur_lon)/2);
            lon_in_avg = avgval_double(n_in, x_in);

            if (  (n_out = clip ( x_in, y_in, n_in, ll_lon, ll_lat, ur_lon, ur_lat, x_out, y_out )) > 0 ) {
              xarea = poly_area (x_out, y_out, n_out ) * mask_in[j1*nx1+i1];
              min_area = min(area_in[j1*nx1+i1], area_out[j2*nx2+i2]);
              if(xarea/min_area > AREA_RATIO_THRESH ) {
                xgrid_area[nxgrid] = xarea;
                xgrid_clon[nxgrid] = poly_ctrlon(x_out, y_out, n_out, lon_in_avg);
                xgrid_clat[nxgrid] = poly_ctrlat (x_out, y_out, n_out );
                i_in[nxgrid]  = i1;
                j_in[nxgrid]  = j1;
                i_out[nxgrid] = i2;
                j_out[nxgrid] = j2;
                ++nxgrid;
                if(nxgrid > MAXXGRID) error_handler("nxgrid is greater than MAXXGRID, increase MAXXGRID");
              }
            }
          }
    }

  free(area_in);
  free(area_out);

  return nxgrid;

} /* create_xgrid_2dx1d_order2 */

/*******************************************************************************
  void create_xgrid_2DX2D_order1
  This routine generate exchange grids between two grids for the first order
  conservative interpolation. nlon_in,nlat_in,nlon_out,nlat_out are the size of the grid cell
  and lon_in,lat_in, lon_out,lat_out are geographic grid location of grid cell bounds.
  mask is on grid lon_in/lat_in.
*******************************************************************************/
int create_xgrid_2dx2d_order1_(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
                               const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                               const double *mask_in, int *i_in, int *j_in, int *i_out,
                               int *j_out, double *xgrid_area)
{
  int nxgrid;

  nxgrid = create_xgrid_2dx2d_order1(nlon_in, nlat_in, nlon_out, nlat_out, lon_in, lat_in, lon_out, lat_out, mask_in,
                                     i_in, j_in, i_out, j_out, xgrid_area);
  return nxgrid;

}
int create_xgrid_2dx2d_order1(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
                              const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *i_out,
                              int *j_out, double *xgrid_area)
{

  int nx1, nx2, ny1, ny2, nx1p, nx2p, nxgrid;
  double *area_in, *area_out;
  int nblocks =1;
  int *istart2=NULL, *iend2=NULL;
  int npts_left, nblks_left, pos, m, npts_my, ij;
  double *lon_out_min_list,*lon_out_max_list,*lon_out_avg,*lat_out_min_list,*lat_out_max_list;
  double *lon_out_list, *lat_out_list;
  int *pnxgrid=NULL, *pstart;
  int *pi_in=NULL, *pj_in=NULL, *pi_out=NULL, *pj_out=NULL;
  double *pxgrid_area=NULL;
  int    *n2_list;
  int nthreads, nxgrid_block_max;

  nx1 = *nlon_in;
  ny1 = *nlat_in;
  nx2 = *nlon_out;
  ny2 = *nlat_out;
  nx1p = nx1 + 1;
  nx2p = nx2 + 1;

  area_in  = (double *)malloc(nx1*ny1*sizeof(double));
  area_out = (double *)malloc(nx2*ny2*sizeof(double));
  get_grid_area(nlon_in, nlat_in, lon_in, lat_in, area_in);
  get_grid_area(nlon_out, nlat_out, lon_out, lat_out, area_out);

  nthreads = 1;
#if defined(_OPENMP)
#pragma omp parallel
  nthreads = omp_get_num_threads();
#endif

  nblocks = nthreads;

  istart2 = (int *)malloc(nblocks*sizeof(int));
  iend2 = (int *)malloc(nblocks*sizeof(int));

  pstart = (int *)malloc(nblocks*sizeof(int));
  pnxgrid = (int *)malloc(nblocks*sizeof(int));

  nxgrid_block_max = MAXXGRID/nblocks;

  for(m=0; m<nblocks; m++) {
    pnxgrid[m] = 0;
    pstart[m] = m*nxgrid_block_max;
  }

  if(nblocks == 1) {
    pi_in = i_in;
    pj_in = j_in;
    pi_out = i_out;
    pj_out = j_out;
    pxgrid_area = xgrid_area;
  }
  else {
    pi_in = (int *)malloc(MAXXGRID*sizeof(int));
    pj_in = (int *)malloc(MAXXGRID*sizeof(int));
    pi_out = (int *)malloc(MAXXGRID*sizeof(int));
    pj_out = (int *)malloc(MAXXGRID*sizeof(int));
    pxgrid_area = (double *)malloc(MAXXGRID*sizeof(double));
  }

  npts_left = nx2*ny2;
  nblks_left = nblocks;
  pos = 0;
  for(m=0; m<nblocks; m++) {
    istart2[m] = pos;
    npts_my = npts_left/nblks_left;
    iend2[m] = istart2[m] + npts_my - 1;
    pos = iend2[m] + 1;
    npts_left -= npts_my;
    nblks_left--;
  }

  lon_out_min_list = (double *)malloc(nx2*ny2*sizeof(double));
  lon_out_max_list = (double *)malloc(nx2*ny2*sizeof(double));
  lat_out_min_list = (double *)malloc(nx2*ny2*sizeof(double));
  lat_out_max_list = (double *)malloc(nx2*ny2*sizeof(double));
  lon_out_avg = (double *)malloc(nx2*ny2*sizeof(double));
  n2_list     = (int *)malloc(nx2*ny2*sizeof(int));
  lon_out_list = (double *)malloc(MAX_V*nx2*ny2*sizeof(double));
  lat_out_list = (double *)malloc(MAX_V*nx2*ny2*sizeof(double));
#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(nx2,ny2,nx2p,lon_out,lat_out,lat_out_min_list, \
                                              lat_out_max_list,lon_out_min_list,lon_out_max_list, \
                                              lon_out_avg,n2_list,lon_out_list,lat_out_list)
#endif
  for(ij=0; ij<nx2*ny2; ij++){
    int i2, j2, n, n0, n1, n2, n3, n2_in, l;
    double x2_in[MV], y2_in[MV];
    i2 = ij%nx2;
    j2 = ij/nx2;
    n = j2*nx2+i2;
    n0 = j2*nx2p+i2; n1 = j2*nx2p+i2+1;
    n2 = (j2+1)*nx2p+i2+1; n3 = (j2+1)*nx2p+i2;
    x2_in[0] = lon_out[n0]; y2_in[0] = lat_out[n0];
    x2_in[1] = lon_out[n1]; y2_in[1] = lat_out[n1];
    x2_in[2] = lon_out[n2]; y2_in[2] = lat_out[n2];
    x2_in[3] = lon_out[n3]; y2_in[3] = lat_out[n3];

    lat_out_min_list[n] = minval_double(4, y2_in);
    lat_out_max_list[n] = maxval_double(4, y2_in);
    n2_in = fix_lon(x2_in, y2_in, 4, M_PI);
    if(n2_in > MAX_V) error_handler("create_xgrid.c: n2_in is greater than MAX_V");
    lon_out_min_list[n] = minval_double(n2_in, x2_in);
    lon_out_max_list[n] = maxval_double(n2_in, x2_in);
    lon_out_avg[n] = avgval_double(n2_in, x2_in);
    n2_list[n] = n2_in;
    for(l=0; l<n2_in; l++) {
      lon_out_list[n*MAX_V+l] = x2_in[l];
      lat_out_list[n*MAX_V+l] = y2_in[l];
    }
  }

  nxgrid = 0;

#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(nblocks,nx1,ny1,nx1p,mask_in,lon_in,lat_in, \
                                              istart2,iend2,nx2,lat_out_min_list,lat_out_max_list, \
                                              n2_list,lon_out_list,lat_out_list,lon_out_min_list, \
                                              lon_out_max_list,lon_out_avg,area_in,area_out, \
                                              pxgrid_area,pnxgrid,pi_in,pj_in,pi_out,pj_out,pstart,nthreads)
#endif
  for(m=0; m<nblocks; m++) {
    int i1, j1, ij;
    for(j1=0; j1<ny1; j1++) for(i1=0; i1<nx1; i1++) if( mask_in[j1*nx1+i1] > MASK_THRESH ) {
          int n0, n1, n2, n3, l,n1_in;
          double lat_in_min,lat_in_max,lon_in_min,lon_in_max,lon_in_avg;
          double x1_in[MV], y1_in[MV], x_out[MV], y_out[MV];

          n0 = j1*nx1p+i1;       n1 = j1*nx1p+i1+1;
          n2 = (j1+1)*nx1p+i1+1; n3 = (j1+1)*nx1p+i1;
          x1_in[0] = lon_in[n0]; y1_in[0] = lat_in[n0];
          x1_in[1] = lon_in[n1]; y1_in[1] = lat_in[n1];
          x1_in[2] = lon_in[n2]; y1_in[2] = lat_in[n2];
          x1_in[3] = lon_in[n3]; y1_in[3] = lat_in[n3];
          lat_in_min = minval_double(4, y1_in);
          lat_in_max = maxval_double(4, y1_in);
          n1_in = fix_lon(x1_in, y1_in, 4, M_PI);
          lon_in_min = minval_double(n1_in, x1_in);
          lon_in_max = maxval_double(n1_in, x1_in);
          lon_in_avg = avgval_double(n1_in, x1_in);
          for(ij=istart2[m]; ij<=iend2[m]; ij++) {
            int n_out, i2, j2, n2_in;
            double xarea, dx, lon_out_min, lon_out_max;
            double x2_in[MAX_V], y2_in[MAX_V];

            i2 = ij%nx2;
            j2 = ij/nx2;

            if(lat_out_min_list[ij] >= lat_in_max || lat_out_max_list[ij] <= lat_in_min ) continue;
            /* adjust x2_in according to lon_in_avg*/
            n2_in = n2_list[ij];
            for(l=0; l<n2_in; l++) {
              x2_in[l] = lon_out_list[ij*MAX_V+l];
              y2_in[l] = lat_out_list[ij*MAX_V+l];
            }
            lon_out_min = lon_out_min_list[ij];
            lon_out_max = lon_out_max_list[ij];
            dx = lon_out_avg[ij] - lon_in_avg;
            if(dx < -M_PI ) {
              lon_out_min += TPI;
              lon_out_max += TPI;
              for (l=0; l<n2_in; l++) x2_in[l] += TPI;
            }
            else if (dx >  M_PI) {
              lon_out_min -= TPI;
              lon_out_max -= TPI;
              for (l=0; l<n2_in; l++) x2_in[l] -= TPI;
            }

            /* x2_in should in the same range as x1_in after lon_fix, so no need to
               consider cyclic condition
            */
            if(lon_out_min >= lon_in_max || lon_out_max <= lon_in_min ) continue;
            if (  (n_out = clip_2dx2d( x1_in, y1_in, n1_in, x2_in, y2_in, n2_in, x_out, y_out )) > 0) {
              double min_area;
              int    nn;
              xarea = poly_area (x_out, y_out, n_out ) * mask_in[j1*nx1+i1];
              min_area = min(area_in[j1*nx1+i1], area_out[j2*nx2+i2]);
              if( xarea/min_area > AREA_RATIO_THRESH ) {
                pnxgrid[m]++;
                if(pnxgrid[m]>= MAXXGRID/nthreads)
                  error_handler("nxgrid is greater than MAXXGRID/nthreads, increase MAXXGRID, decrease nthreads, or increase number of MPI ranks");
                nn = pstart[m] + pnxgrid[m]-1;

                pxgrid_area[nn] = xarea;
                pi_in[nn]       = i1;
                pj_in[nn]       = j1;
                pi_out[nn]      = i2;
                pj_out[nn]      = j2;
              }

            }

          }
        }
  }

  /*copy data if nblocks > 1 */
  if(nblocks == 1) {
    nxgrid = pnxgrid[0];
    pi_in = NULL;
    pj_in = NULL;
    pi_out = NULL;
    pj_out = NULL;
    pxgrid_area = NULL;
  }
  else {
    int nn, i;
    nxgrid = 0;
    for(m=0; m<nblocks; m++) {
      for(i=0; i<pnxgrid[m]; i++) {
        nn = pstart[m] + i;
        i_in[nxgrid] = pi_in[nn];
        j_in[nxgrid] = pj_in[nn];
        i_out[nxgrid] = pi_out[nn];
        j_out[nxgrid] = pj_out[nn];
        xgrid_area[nxgrid] = pxgrid_area[nn];
        nxgrid++;
      }
    }
    free(pi_in);
    free(pj_in);
    free(pi_out);
    free(pj_out);
    free(pxgrid_area);
  }

  free(area_in);
  free(area_out);
  free(lon_out_min_list);
  free(lon_out_max_list);
  free(lat_out_min_list);
  free(lat_out_max_list);
  free(lon_out_avg);
  free(n2_list);
  free(lon_out_list);
  free(lat_out_list);

  return nxgrid;

}/* get_xgrid_2Dx2D_order1 */

/********************************************************************************
  void create_xgrid_2dx1d_order2
  This routine generate exchange grids between two grids for the second order
  conservative interpolation. nlon_in,nlat_in,nlon_out,nlat_out are the size of the grid cell
  and lon_in,lat_in, lon_out,lat_out are geographic grid location of grid cell bounds.
  mask is on grid lon_in/lat_in.
********************************************************************************/
int create_xgrid_2dx2d_order2_(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
                               const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                               const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
                               double *xgrid_area, double *xgrid_clon, double *xgrid_clat)
{
  int nxgrid;
  nxgrid = create_xgrid_2dx2d_order2(nlon_in, nlat_in, nlon_out, nlat_out, lon_in, lat_in, lon_out, lat_out, mask_in, i_in,
                                     j_in, i_out, j_out, xgrid_area, xgrid_clon, xgrid_clat);
  return nxgrid;

}
int create_xgrid_2dx2d_order2(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
                              const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
                              double *xgrid_area, double *xgrid_clon, double *xgrid_clat)
{

  int nx1, nx2, ny1, ny2, nx1p, nx2p, nxgrid;
  double *area_in, *area_out;
  int nblocks =1;
  int *istart2=NULL, *iend2=NULL;
  int npts_left, nblks_left, pos, m, npts_my, ij;
  double *lon_out_min_list,*lon_out_max_list,*lon_out_avg,*lat_out_min_list,*lat_out_max_list;
  double *lon_out_list, *lat_out_list;
  int *pnxgrid=NULL, *pstart;
  int *pi_in=NULL, *pj_in=NULL, *pi_out=NULL, *pj_out=NULL;
  double *pxgrid_area=NULL, *pxgrid_clon=NULL, *pxgrid_clat=NULL;
  int    *n2_list;
  int nthreads, nxgrid_block_max;

  nx1 = *nlon_in;
  ny1 = *nlat_in;
  nx2 = *nlon_out;
  ny2 = *nlat_out;
  nx1p = nx1 + 1;
  nx2p = nx2 + 1;

  area_in  = (double *)malloc(nx1*ny1*sizeof(double));
  area_out = (double *)malloc(nx2*ny2*sizeof(double));
  get_grid_area(nlon_in, nlat_in, lon_in, lat_in, area_in);
  get_grid_area(nlon_out, nlat_out, lon_out, lat_out, area_out);

  nthreads = 1;
#if defined(_OPENMP)
#pragma omp parallel
  nthreads = omp_get_num_threads();
#endif

  nblocks = nthreads;

  istart2 = (int *)malloc(nblocks*sizeof(int));
  iend2 = (int *)malloc(nblocks*sizeof(int));

  pstart = (int *)malloc(nblocks*sizeof(int));
  pnxgrid = (int *)malloc(nblocks*sizeof(int));

  nxgrid_block_max = MAXXGRID/nblocks;

  for(m=0; m<nblocks; m++) {
    pnxgrid[m] = 0;
    pstart[m] = m*nxgrid_block_max;
  }

  if(nblocks == 1) {
    pi_in = i_in;
    pj_in = j_in;
    pi_out = i_out;
    pj_out = j_out;
    pxgrid_area = xgrid_area;
    pxgrid_clon = xgrid_clon;
    pxgrid_clat = xgrid_clat;
  }
  else {
    pi_in = (int *)malloc(MAXXGRID*sizeof(int));
    pj_in = (int *)malloc(MAXXGRID*sizeof(int));
    pi_out = (int *)malloc(MAXXGRID*sizeof(int));
    pj_out = (int *)malloc(MAXXGRID*sizeof(int));
    pxgrid_area = (double *)malloc(MAXXGRID*sizeof(double));
    pxgrid_clon = (double *)malloc(MAXXGRID*sizeof(double));
    pxgrid_clat = (double *)malloc(MAXXGRID*sizeof(double));
  }

  npts_left = nx2*ny2;
  nblks_left = nblocks;
  pos = 0;
  for(m=0; m<nblocks; m++) {
    istart2[m] = pos;
    npts_my = npts_left/nblks_left;
    iend2[m] = istart2[m] + npts_my - 1;
    pos = iend2[m] + 1;
    npts_left -= npts_my;
    nblks_left--;
  }

  lon_out_min_list = (double *)malloc(nx2*ny2*sizeof(double));
  lon_out_max_list = (double *)malloc(nx2*ny2*sizeof(double));
  lat_out_min_list = (double *)malloc(nx2*ny2*sizeof(double));
  lat_out_max_list = (double *)malloc(nx2*ny2*sizeof(double));
  lon_out_avg = (double *)malloc(nx2*ny2*sizeof(double));
  n2_list     = (int *)malloc(nx2*ny2*sizeof(int));
  lon_out_list = (double *)malloc(MAX_V*nx2*ny2*sizeof(double));
  lat_out_list = (double *)malloc(MAX_V*nx2*ny2*sizeof(double));
#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(nx2,ny2,nx2p,lon_out,lat_out,lat_out_min_list, \
                                              lat_out_max_list,lon_out_min_list,lon_out_max_list, \
                                              lon_out_avg,n2_list,lon_out_list,lat_out_list)
#endif
  for(ij=0; ij<nx2*ny2; ij++){
    int i2, j2, n, n0, n1, n2, n3, n2_in, l;
    double x2_in[MV], y2_in[MV];
    i2 = ij%nx2;
    j2 = ij/nx2;
    n = j2*nx2+i2;
    n0 = j2*nx2p+i2; n1 = j2*nx2p+i2+1;
    n2 = (j2+1)*nx2p+i2+1; n3 = (j2+1)*nx2p+i2;
    x2_in[0] = lon_out[n0]; y2_in[0] = lat_out[n0];
    x2_in[1] = lon_out[n1]; y2_in[1] = lat_out[n1];
    x2_in[2] = lon_out[n2]; y2_in[2] = lat_out[n2];
    x2_in[3] = lon_out[n3]; y2_in[3] = lat_out[n3];

    lat_out_min_list[n] = minval_double(4, y2_in);
    lat_out_max_list[n] = maxval_double(4, y2_in);
    n2_in = fix_lon(x2_in, y2_in, 4, M_PI);
    if(n2_in > MAX_V) error_handler("create_xgrid.c: n2_in is greater than MAX_V");
    lon_out_min_list[n] = minval_double(n2_in, x2_in);
    lon_out_max_list[n] = maxval_double(n2_in, x2_in);
    lon_out_avg[n] = avgval_double(n2_in, x2_in);
    n2_list[n] = n2_in;
    for(l=0; l<n2_in; l++) {
      lon_out_list[n*MAX_V+l] = x2_in[l];
      lat_out_list[n*MAX_V+l] = y2_in[l];
    }
  }

  nxgrid = 0;

#if defined(_OPENMP)
#pragma omp parallel for default(none) shared(nblocks,nx1,ny1,nx1p,mask_in,lon_in,lat_in, \
                                              istart2,iend2,nx2,lat_out_min_list,lat_out_max_list, \
                                              n2_list,lon_out_list,lat_out_list,lon_out_min_list, \
                                              lon_out_max_list,lon_out_avg,area_in,area_out, \
                                              pxgrid_area,pnxgrid,pxgrid_clon,pxgrid_clat,pi_in, \
                                              pj_in,pi_out,pj_out,pstart,nthreads)
#endif
  for(m=0; m<nblocks; m++) {
    int i1, j1, ij;
    for(j1=0; j1<ny1; j1++) for(i1=0; i1<nx1; i1++) if( mask_in[j1*nx1+i1] > MASK_THRESH ) {
          int n0, n1, n2, n3, l,n1_in;
          double lat_in_min,lat_in_max,lon_in_min,lon_in_max,lon_in_avg;
          double x1_in[MV], y1_in[MV], x_out[MV], y_out[MV];

          n0 = j1*nx1p+i1;       n1 = j1*nx1p+i1+1;
          n2 = (j1+1)*nx1p+i1+1; n3 = (j1+1)*nx1p+i1;
          x1_in[0] = lon_in[n0]; y1_in[0] = lat_in[n0];
          x1_in[1] = lon_in[n1]; y1_in[1] = lat_in[n1];
          x1_in[2] = lon_in[n2]; y1_in[2] = lat_in[n2];
          x1_in[3] = lon_in[n3]; y1_in[3] = lat_in[n3];
          lat_in_min = minval_double(4, y1_in);
          lat_in_max = maxval_double(4, y1_in);
          n1_in = fix_lon(x1_in, y1_in, 4, M_PI);
          lon_in_min = minval_double(n1_in, x1_in);
          lon_in_max = maxval_double(n1_in, x1_in);
          lon_in_avg = avgval_double(n1_in, x1_in);
          for(ij=istart2[m]; ij<=iend2[m]; ij++) {
            int n_out, i2, j2, n2_in;
            double xarea, dx, lon_out_min, lon_out_max;
            double x2_in[MAX_V], y2_in[MAX_V];

            i2 = ij%nx2;
            j2 = ij/nx2;

            if(lat_out_min_list[ij] >= lat_in_max || lat_out_max_list[ij] <= lat_in_min ) continue;
            /* adjust x2_in according to lon_in_avg*/
            n2_in = n2_list[ij];
            for(l=0; l<n2_in; l++) {
              x2_in[l] = lon_out_list[ij*MAX_V+l];
              y2_in[l] = lat_out_list[ij*MAX_V+l];
            }
            lon_out_min = lon_out_min_list[ij];
            lon_out_max = lon_out_max_list[ij];
            dx = lon_out_avg[ij] - lon_in_avg;
            if(dx < -M_PI ) {
              lon_out_min += TPI;
              lon_out_max += TPI;
              for (l=0; l<n2_in; l++) x2_in[l] += TPI;
            }
            else if (dx >  M_PI) {
              lon_out_min -= TPI;
              lon_out_max -= TPI;
              for (l=0; l<n2_in; l++) x2_in[l] -= TPI;
            }

            /* x2_in should in the same range as x1_in after lon_fix, so no need to
               consider cyclic condition
            */
            if(lon_out_min >= lon_in_max || lon_out_max <= lon_in_min ) continue;
            if (  (n_out = clip_2dx2d( x1_in, y1_in, n1_in, x2_in, y2_in, n2_in, x_out, y_out )) > 0) {
              double min_area;
              int nn;
              xarea = poly_area (x_out, y_out, n_out ) * mask_in[j1*nx1+i1];
              min_area = min(area_in[j1*nx1+i1], area_out[j2*nx2+i2]);
              if( xarea/min_area > AREA_RATIO_THRESH ) {
                pnxgrid[m]++;
                if(pnxgrid[m]>= MAXXGRID/nthreads)
                  error_handler("nxgrid is greater than MAXXGRID/nthreads, increase MAXXGRID, decrease nthreads, or increase number of MPI ranks");
                nn = pstart[m] + pnxgrid[m]-1;
                pxgrid_area[nn] = xarea;
                pxgrid_clon[nn] = poly_ctrlon(x_out, y_out, n_out, lon_in_avg);
                pxgrid_clat[nn] = poly_ctrlat (x_out, y_out, n_out );
                pi_in[nn]       = i1;
                pj_in[nn]       = j1;
                pi_out[nn]      = i2;
                pj_out[nn]      = j2;
              }
            }
          }
        }
  }

  /*copy data if nblocks > 1 */
  if(nblocks == 1) {
    nxgrid = pnxgrid[0];
    pi_in = NULL;
    pj_in = NULL;
    pi_out = NULL;
    pj_out = NULL;
    pxgrid_area = NULL;
    pxgrid_clon = NULL;
    pxgrid_clat = NULL;
  }
  else {
    int nn, i;
    nxgrid = 0;
    for(m=0; m<nblocks; m++) {
      for(i=0; i<pnxgrid[m]; i++) {
        nn = pstart[m] + i;
        i_in[nxgrid] = pi_in[nn];
        j_in[nxgrid] = pj_in[nn];
        i_out[nxgrid] = pi_out[nn];
        j_out[nxgrid] = pj_out[nn];
        xgrid_area[nxgrid] = pxgrid_area[nn];
        xgrid_clon[nxgrid] = pxgrid_clon[nn];
        xgrid_clat[nxgrid] = pxgrid_clat[nn];
        nxgrid++;
      }
    }
    free(pi_in);
    free(pj_in);
    free(pi_out);
    free(pj_out);
    free(pxgrid_area);
    free(pxgrid_clon);
    free(pxgrid_clat);
  }

  free(area_in);
  free(area_out);
  free(lon_out_min_list);
  free(lon_out_max_list);
  free(lat_out_min_list);
  free(lat_out_max_list);
  free(lon_out_avg);
  free(n2_list);
  free(lon_out_list);
  free(lat_out_list);

  return nxgrid;

}/* get_xgrid_2Dx2D_order2 */

int create_xgrid_great_circle_(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
                               const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                               const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
                               double *xgrid_area, double *xgrid_clon, double *xgrid_clat)
{
  int nxgrid;
  nxgrid = create_xgrid_great_circle(nlon_in, nlat_in, nlon_out, nlat_out, lon_in, lat_in, lon_out, lat_out,
                                     mask_in, i_in, j_in, i_out, j_out, xgrid_area, xgrid_clon, xgrid_clat);

  return nxgrid;
}

int create_xgrid_great_circle(const int *nlon_in, const int *nlat_in, const int *nlon_out, const int *nlat_out,
                              const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *i_out, int *j_out,
                              double *xgrid_area, double *xgrid_clon, double *xgrid_clat)
{

  int nx1, nx2, ny1, ny2, nx1p, nx2p, ny1p, ny2p, nxgrid, n1_in, n2_in;
  int n0, n1, n2, n3, i1, j1, i2, j2;
  double x1_in[MV], y1_in[MV], z1_in[MV];
  double x2_in[MV], y2_in[MV], z2_in[MV];
  double x_out[MV], y_out[MV], z_out[MV];
  double *x1=NULL, *y1=NULL, *z1=NULL;
  double *x2=NULL, *y2=NULL, *z2=NULL;

  double *area1, *area2, min_area;

  nx1 = *nlon_in;
  ny1 = *nlat_in;
  nx2 = *nlon_out;
  ny2 = *nlat_out;
  nxgrid = 0;
  nx1p = nx1 + 1;
  nx2p = nx2 + 1;
  ny1p = ny1 + 1;
  ny2p = ny2 + 1;

  /* first convert lon-lat to cartesian coordinates */
  x1 = (double *)malloc(nx1p*ny1p*sizeof(double));
  y1 = (double *)malloc(nx1p*ny1p*sizeof(double));
  z1 = (double *)malloc(nx1p*ny1p*sizeof(double));
  x2 = (double *)malloc(nx2p*ny2p*sizeof(double));
  y2 = (double *)malloc(nx2p*ny2p*sizeof(double));
  z2 = (double *)malloc(nx2p*ny2p*sizeof(double));

  latlon2xyz(nx1p*ny1p, lon_in, lat_in, x1, y1, z1);
  latlon2xyz(nx2p*ny2p, lon_out, lat_out, x2, y2, z2);

  area1  = (double *)malloc(nx1*ny1*sizeof(double));
  area2 = (double *)malloc(nx2*ny2*sizeof(double));
  get_grid_great_circle_area(nlon_in, nlat_in, lon_in, lat_in, area1);
  get_grid_great_circle_area(nlon_out, nlat_out, lon_out, lat_out, area2);
  n1_in = 4;
  n2_in = 4;

  for(j1=0; j1<ny1; j1++) for(i1=0; i1<nx1; i1++) if( mask_in[j1*nx1+i1] > MASK_THRESH ) {
        /* clockwise */
        n0 = j1*nx1p+i1;       n1 = (j1+1)*nx1p+i1;
        n2 = (j1+1)*nx1p+i1+1; n3 = j1*nx1p+i1+1;
        x1_in[0] = x1[n0]; y1_in[0] = y1[n0]; z1_in[0] = z1[n0];
        x1_in[1] = x1[n1]; y1_in[1] = y1[n1]; z1_in[1] = z1[n1];
        x1_in[2] = x1[n2]; y1_in[2] = y1[n2]; z1_in[2] = z1[n2];
        x1_in[3] = x1[n3]; y1_in[3] = y1[n3]; z1_in[3] = z1[n3];

        for(j2=0; j2<ny2; j2++) for(i2=0; i2<nx2; i2++) {
            int n_out;
            double xarea;

            n0 = j2*nx2p+i2;       n1 = (j2+1)*nx2p+i2;
            n2 = (j2+1)*nx2p+i2+1; n3 = j2*nx2p+i2+1;
            x2_in[0] = x2[n0]; y2_in[0] = y2[n0]; z2_in[0] = z2[n0];
            x2_in[1] = x2[n1]; y2_in[1] = y2[n1]; z2_in[1] = z2[n1];
            x2_in[2] = x2[n2]; y2_in[2] = y2[n2]; z2_in[2] = z2[n2];
            x2_in[3] = x2[n3]; y2_in[3] = y2[n3]; z2_in[3] = z2[n3];

            if (  (n_out = clip_2dx2d_great_circle( x1_in, y1_in, z1_in, n1_in, x2_in, y2_in, z2_in, n2_in,
                                                    x_out, y_out, z_out)) > 0) {
              xarea = great_circle_area ( n_out, x_out, y_out, z_out ) * mask_in[j1*nx1+i1];
              min_area = min(area1[j1*nx1+i1], area2[j2*nx2+i2]);
              if( xarea/min_area > AREA_RATIO_THRESH ) {
                xgrid_area[nxgrid] = xarea;
                xgrid_clon[nxgrid] = 0; /*z1l: will be developed very soon */
                xgrid_clat[nxgrid] = 0;
                i_in[nxgrid]       = i1;
                j_in[nxgrid]       = j1;
                i_out[nxgrid]      = i2;
                j_out[nxgrid]      = j2;
                ++nxgrid;
                if(nxgrid > MAXXGRID) error_handler("nxgrid is greater than MAXXGRID, increase MAXXGRID");
              }
            }
          }
      }


  free(area1);
  free(area2);

  free(x1);
  free(y1);
  free(z1);
  free(x2);
  free(y2);
  free(z2);

  return nxgrid;

}/* create_xgrid_great_circle */

int create_xgrid_great_circle_ug_(const int *nlon_in, const int *nlat_in, const int *npts_out,
                               const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                               const double *mask_in, int *i_in, int *j_in, int *l_out,
                               double *xgrid_area, double *xgrid_clon, double *xgrid_clat)
{
  int nxgrid;
  nxgrid = create_xgrid_great_circle_ug(nlon_in, nlat_in, npts_out, lon_in, lat_in, lon_out, lat_out,
                                     mask_in, i_in, j_in, l_out, xgrid_area, xgrid_clon, xgrid_clat);

  return nxgrid;
}

int create_xgrid_great_circle_ug(const int *nlon_in, const int *nlat_in, const int *npts_out,
                              const double *lon_in, const double *lat_in, const double *lon_out, const double *lat_out,
                              const double *mask_in, int *i_in, int *j_in, int *l_out,
                              double *xgrid_area, double *xgrid_clon, double *xgrid_clat)
{

  int nx1, ny1, npts2, nx1p, ny1p, nxgrid, n1_in, n2_in, nv;
  int n0, n1, n2, n3, i1, j1, l2;
  double x1_in[MV], y1_in[MV], z1_in[MV];
  double x2_in[MV], y2_in[MV], z2_in[MV];
  double x_out[MV], y_out[MV], z_out[MV];
  double *x1=NULL, *y1=NULL, *z1=NULL;
  double *x2=NULL, *y2=NULL, *z2=NULL;

  double *area1, *area2, min_area;

  nx1 = *nlon_in;
  ny1 = *nlat_in;
  nv = 4;
  npts2 = *npts_out;
  nxgrid = 0;
  nx1p = nx1 + 1;
  ny1p = ny1 + 1;

  /* first convert lon-lat to cartesian coordinates */
  x1 = (double *)malloc(nx1p*ny1p*sizeof(double));
  y1 = (double *)malloc(nx1p*ny1p*sizeof(double));
  z1 = (double *)malloc(nx1p*ny1p*sizeof(double));
  x2 = (double *)malloc(npts2*nv*sizeof(double));
  y2 = (double *)malloc(npts2*nv*sizeof(double));
  z2 = (double *)malloc(npts2*nv*sizeof(double));

  latlon2xyz(nx1p*ny1p, lon_in, lat_in, x1, y1, z1);
  latlon2xyz(npts2*nv, lon_out, lat_out, x2, y2, z2);

  area1  = (double *)malloc(nx1*ny1*sizeof(double));
  area2 = (double *)malloc(npts2*sizeof(double));
  get_grid_great_circle_area(nlon_in, nlat_in, lon_in, lat_in, area1);
  get_grid_great_circle_area_ug(npts_out, lon_out, lat_out, area2);
  n1_in = 4;
  n2_in = 4;

  for(j1=0; j1<ny1; j1++) for(i1=0; i1<nx1; i1++) if( mask_in[j1*nx1+i1] > MASK_THRESH ) {
        /* clockwise */
        n0 = j1*nx1p+i1;       n1 = (j1+1)*nx1p+i1;
        n2 = (j1+1)*nx1p+i1+1; n3 = j1*nx1p+i1+1;
        x1_in[0] = x1[n0]; y1_in[0] = y1[n0]; z1_in[0] = z1[n0];
        x1_in[1] = x1[n1]; y1_in[1] = y1[n1]; z1_in[1] = z1[n1];
        x1_in[2] = x1[n2]; y1_in[2] = y1[n2]; z1_in[2] = z1[n2];
        x1_in[3] = x1[n3]; y1_in[3] = y1[n3]; z1_in[3] = z1[n3];

        for(l2=0; l2<npts2; l2++) {
            int n_out;
            double xarea;

            n0 = l2*nv;   n1 = l2*nv+1;
            n2 = l2*nv+2; n3 = l2*nv+3;
            x2_in[0] = x2[n0]; y2_in[0] = y2[n0]; z2_in[0] = z2[n0];
            x2_in[1] = x2[n1]; y2_in[1] = y2[n1]; z2_in[1] = z2[n1];
            x2_in[2] = x2[n2]; y2_in[2] = y2[n2]; z2_in[2] = z2[n2];
            x2_in[3] = x2[n3]; y2_in[3] = y2[n3]; z2_in[3] = z2[n3];

            if (  (n_out = clip_2dx2d_great_circle( x1_in, y1_in, z1_in, n1_in, x2_in, y2_in, z2_in, n2_in,
                                                    x_out, y_out, z_out)) > 0) {
              xarea = great_circle_area ( n_out, x_out, y_out, z_out ) * mask_in[j1*nx1+i1];
              min_area = min(area1[j1*nx1+i1], area2[l2]);
              if( xarea/min_area > AREA_RATIO_THRESH ) {
                xgrid_area[nxgrid] = xarea;
                xgrid_clon[nxgrid] = 0; /*z1l: will be developed very soon */
                xgrid_clat[nxgrid] = 0;
                i_in[nxgrid]       = i1;
                j_in[nxgrid]       = j1;
                l_out[nxgrid]      = l2;
                ++nxgrid;
                if(nxgrid > MAXXGRID) error_handler("nxgrid is greater than MAXXGRID, increase MAXXGRID");
              }
            }
          }
      }


  free(area1);
  free(area2);

  free(x1);
  free(y1);
  free(z1);
  free(x2);
  free(y2);
  free(z2);

  return nxgrid;

}/* create_xgrid_great_circle_ug */

/*******************************************************************************
  int get_maxxgrid
  return constants MAXXGRID.
*******************************************************************************/
int get_maxxgrid(void)
{
  return MAXXGRID;
}

int get_maxxgrid_(void)
{
  return get_maxxgrid();
}
