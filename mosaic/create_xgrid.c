#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mosaic_util.h"
#include "create_xgrid.h"


const double MASK_THRESH = 1.0e-8;
const double EPSLN = 1.0e-10;
double grid_box_radius(const double *x, const double *y, const double *z, int n);
double dist_between_boxes(const double *x1, const double *y1, const double *z1, int n1,
			  const double *x2, const double *y2, const double *z2, int n2);

/*******************************************************************************
  int get_maxxgrid
  return constants MAXXGRID.
*******************************************************************************/
int get_maxxgrid(void)  
{
  return MAXXGRID;
}

#ifndef __AIX
int get_maxxgrid_(void)
{
  return get_maxxgrid();
}
#endif

/*******************************************************************************
void get_grid_area(const int *nlon, const int *nlat, const double *lon, const double *lat, const double *area)
  return the grid area.
*******************************************************************************/
#ifndef __AIX
void get_grid_area_(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area)
{
  get_grid_area(nlon, nlat, lon, lat, area);
}
#endif

void get_grid_area(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area)
{
  int nx, ny, nxp, i, j, n_in;
  double x_in[20], y_in[20];
  
  nx = *nlon;
  ny = *nlat;
  nxp = nx + 1;
  
  for(j=0; j<ny; j++) for(i=0; i < nx; i++) {
    x_in[0] = lon[j*nxp+i];
    x_in[1] = lon[j*nxp+i+1];
    x_in[2] = lon[(j+1)*nxp+i+1];
    x_in[3] = lon[(j+1)*nxp+i];
    y_in[0] = lat[j*nxp+i];
    y_in[1] = lat[j*nxp+i+1];
    y_in[2] = lat[(j+1)*nxp+i+1];
    y_in[3] = lat[(j+1)*nxp+i];
    n_in = fix_lon(x_in, y_in, 4, M_PI);    
    area[j*nx+i] = poly_area(x_in, y_in, n_in);
  }

};  /* get_grid_area */


/*******************************************************************************
  void create_xgrid_1dx2d_order1
  This routine generate exchange grids between two grids for the first order
  conservative interpolation. nlon1,nlat1,nlon2,nlat2 are the size of the grid cell
  and lon1,lat1 are 1-D grid bounds, lon2,lat2 are geographic grid location of grid cell bounds.
  mask is on grid lon1/lat1. 
*******************************************************************************/
#ifndef __AIX
int create_xgrid_1dx2d_order1_(const int *nlon1, const int *nlat1, const int *nlon2, const int *nlat2,
			       const double *lon1, const double *lat1, const double *lon2, const double *lat2,
			       const double *mask1, int *grid1_i, int *grid1_j, int *grid2_i,
			       int *grid2_j, double *xgrid_area)
{
  int nxgrid;
  
  nxgrid = create_xgrid_1dx2d_order1(nlon1, nlat1, nlon2, nlat2, lon1, lat1, lon2, lat2, mask1,
			       grid1_i, grid1_j, grid2_i, grid2_j, xgrid_area);
  return nxgrid;
    
};  
#endif
int create_xgrid_1dx2d_order1(const int *nlon1, const int *nlat1, const int *nlon2, const int *nlat2, const double *lon1,
			      const double *lat1, const double *lon2, const double *lat2,
			      const double *mask1, int *grid1_i, int *grid1_j, int *grid2_i,
			      int *grid2_j, double *xgrid_area)
{

  int nx1, ny1, nx2, ny2, nx1p, nx2p;
  int i1, j1, i2, j2, nxgrid;
  double ll_lon, ll_lat, ur_lon, ur_lat, lon_in[MV], lat_in[MV], lon_out[MV], lat_out[MV];
  
  nx1 = *nlon1;
  ny1 = *nlat1;
  nx2 = *nlon2;
  ny2 = *nlat2;

  nxgrid = 0;
  nx1p = nx1 + 1;
  nx2p = nx2 + 1;

  for(j1=0; j1<ny1; j1++) for(i1=0; i1<nx1; i1++) if( mask1[j1*nx1+i1] > MASK_THRESH ) {

    ll_lon = lon1[i1];   ll_lat = lat1[j1];
    ur_lon = lon1[i1+1]; ur_lat = lat1[j1+1];
    for(j2=0; j2<ny2; j2++) for(i2=0; i2<nx2; i2++) {
      int n_in, n_out;
      double Xarea;
      
      lat_in[0] = lat2[j2*nx2p+i2];
      lat_in[1] = lat2[j2*nx2p+i2+1];
      lat_in[2] = lat2[(j2+1)*nx2p+i2+1];
      lat_in[3] = lat2[(j2+1)*nx2p+i2];
      if (  (lat_in[0]<=ll_lat) && (lat_in[1]<=ll_lat)
	    && (lat_in[2]<=ll_lat) && (lat_in[3]<=ll_lat) ) continue;
      if (  (lat_in[0]>=ur_lat) && (lat_in[1]>=ur_lat)
	    && (lat_in[2]>=ur_lat) && (lat_in[3]>=ur_lat) ) continue;

      lon_in[0] = lon2[j2*nx2p+i2];
      lon_in[1] = lon2[j2*nx2p+i2+1];
      lon_in[2] = lon2[(j2+1)*nx2p+i2+1];
      lon_in[3] = lon2[(j2+1)*nx2p+i2];

      n_in = fix_lon(lon_in, lat_in, 4, (ll_lon+ur_lon)/2);
      
      if ( (n_out = clip ( lon_in, lat_in, n_in, ll_lon, ll_lat, ur_lon, ur_lat, lon_out, lat_out )) > 0 ) {
	Xarea = poly_area ( lon_out, lat_out, n_out ) * mask1[j1*nx1+i1];
	if( Xarea > AREA_THRESH ) {
      	  xgrid_area[nxgrid] = Xarea;
	  grid1_i[nxgrid]    = i1;
	  grid1_j[nxgrid]    = j1;
	  grid2_i[nxgrid]    = i2;
	  grid2_j[nxgrid]    = j2;
	  ++nxgrid;
	  if(nxgrid > MAXXGRID) error_handler("nxgrid is greater than MAXXGRID, increase MAXXGRID");
	}
      }
    }
  }

  return nxgrid;
  
}; /* create_xgrid_1dx2d_order1 */


/********************************************************************************
  void create_xgrid_1dx2d_order2
  This routine generate exchange grids between two grids for the second order
  conservative interpolation. nlon1,nlat1,nlon2,nlat2 are the size of the grid cell
  and lon1,lat1 are 1-D grid bounds, lon2,lat2 are geographic grid location of grid cell bounds.
  mask is on grid lon1/lat1. 
********************************************************************************/
#ifndef __AIX
int create_xgrid_1dx2d_order2_(const int *nlon1, const int *nlat1, const int *nlon2, const int *nlat2,
			 const double *lon1, const double *lat1, const double *lon2, const double *lat2,
			 const double *mask1, int *grid1_i, int *grid1_j, int *grid2_i, int *grid2_j,
			 double *xgrid_area, double *grid1_di, double *grid1_dj, double *grid2_di, double *grid2_dj)
{
  int nxgrid;
  nxgrid = create_xgrid_1dx2d_order2(nlon1, nlat1, nlon2, nlat2, lon1, lat1, lon2, lat2, mask1, grid1_i,
                                     grid1_j, grid2_i, grid2_j, xgrid_area, grid1_di, grid1_dj, grid2_di, grid2_dj);
  return nxgrid;

};
#endif
int create_xgrid_1dx2d_order2(const int *nlon1, const int *nlat1, const int *nlon2, const int *nlat2,
			const double *lon1, const double *lat1, const double *lon2, const double *lat2,
			const double *mask1, int *grid1_i, int *grid1_j, int *grid2_i, int *grid2_j,
			double *xgrid_area, double *grid1_di, double *grid1_dj, double *grid2_di, double *grid2_dj)
{

  int nx1, ny1, nx2, ny2, nx1p, nx2p;
  int i1, j1, i2, j2, nxgrid;
  double ll_lon, ll_lat, ur_lon, ur_lat, lon_in[MV], lat_in[MV], lon_out[MV], lat_out[MV];
  double *ctrlon1, *ctrlat1, ctrlon2, *ctrlat2;
  
  nx1 = *nlon1;
  ny1 = *nlat1;
  nx2 = *nlon2;
  ny2 = *nlat2;

  nxgrid = 0;
  nx1p = nx1 + 1;
  nx2p = nx2 + 1;

  ctrlon1    = (double *)malloc(nx1*ny1*sizeof(double));
  ctrlat1    = (double *)malloc(nx1*ny1*sizeof(double));
  ctrlat2    = (double *)malloc(nx2*ny2*sizeof(double));
  
  for(j1=0; j1<ny1; j1++) {
    for(i1=0; i1<nx1; i1++) {
      ll_lon = lon1[i1]; ll_lat = lat1[j1];
      ur_lon = lon1[i1+1]; ur_lat = lat1[j1+1];
      ctrlon1[j1*nx1+i1] = box_ctrlon(ll_lon, ll_lat, ur_lon, ur_lat, 0.5*(ll_lon+ur_lon) );
      ctrlat1[j1*nx1+i1] = box_ctrlat(ll_lon, ll_lat, ur_lon, ur_lat );
    }
  }

  for(j2=0; j2<ny2; j2++) {
    for(i2=0; i2<nx2; i2++) {
      int n_in;
      lon_in[0] = lon2[j2*nx2p+i2];
      lon_in[1] = lon2[j2*nx2p+i2+1];
      lon_in[2] = lon2[(j2+1)*nx2p+i2+1];
      lon_in[3] = lon2[(j2+1)*nx2p+i2];
      lat_in[0] = lat2[j2*nx2p+i2];
      lat_in[1] = lat2[j2*nx2p+i2+1];
      lat_in[2] = lat2[(j2+1)*nx2p+i2+1];
      lat_in[3] = lat2[(j2+1)*nx2p+i2];
      n_in = fix_lon(lon_in, lat_in, 4, M_PI);
      ctrlat2[j2*nx2+i2]    = poly_ctrlat(lon_in, lat_in, n_in);
    }
  }

  for(j1=0; j1<ny1; j1++) for(i1=0; i1<nx1; i1++) if( mask1[j1*nx1+i1] > MASK_THRESH ) {

    ll_lon = lon1[i1];   ll_lat = lat1[j1];
    ur_lon = lon1[i1+1]; ur_lat = lat1[j1+1];
    for(j2=0; j2<ny2; j2++) for(i2=0; i2<nx2; i2++) {
      int n_in, n_out;
      double xarea, xctrlon, xctrlat;
      
      lat_in[0] = lat2[j2*nx2p+i2];
      lat_in[1] = lat2[j2*nx2p+i2+1];
      lat_in[2] = lat2[(j2+1)*nx2p+i2+1];
      lat_in[3] = lat2[(j2+1)*nx2p+i2];
      if (  (lat_in[0]<=ll_lat) && (lat_in[1]<=ll_lat)
	    && (lat_in[2]<=ll_lat) && (lat_in[3]<=ll_lat) ) continue;
      if (  (lat_in[0]>=ur_lat) && (lat_in[1]>=ur_lat)
	    && (lat_in[2]>=ur_lat) && (lat_in[3]>=ur_lat) ) continue;

      lon_in[0] = lon2[j2*nx2p+i2];
      lon_in[1] = lon2[j2*nx2p+i2+1];
      lon_in[2] = lon2[(j2+1)*nx2p+i2+1];
      lon_in[3] = lon2[(j2+1)*nx2p+i2];

      n_in = fix_lon(lon_in, lat_in, 4, (ll_lon+ur_lon)/2);
      
      if (  (n_out = clip ( lon_in, lat_in, n_in, ll_lon, ll_lat, ur_lon, ur_lat, lon_out, lat_out )) > 0 ) {
	xarea = poly_area (lon_out, lat_out, n_out ) * mask1[j1*nx1+i1];
	if(xarea > AREA_THRESH ) {	  
	  xgrid_area[nxgrid] = xarea;
	  grid1_i[nxgrid]    = i1;
	  grid1_j[nxgrid]    = j1;
	  grid2_i[nxgrid]    = i2;
	  grid2_j[nxgrid]    = j2;
	  ctrlon2 = poly_ctrlon(lon_in, lat_in, n_in, 0.5*(ll_lon+ur_lon));
	  xctrlon = poly_ctrlon ( lon_out, lat_out, n_out, 0.5*(ll_lon+ur_lon));
	  xctrlat = poly_ctrlat ( lon_out, lat_out, n_out );
	  grid1_di[nxgrid] = xctrlon - ctrlon1[j1*nx1+i1];
	  grid1_dj[nxgrid] = xctrlat - ctrlat1[j1*nx1+i1];
	  grid2_di[nxgrid] = xctrlon - ctrlon2;
	  grid2_dj[nxgrid] = xctrlat - ctrlat2[j2*nx2+i2];
	  ++nxgrid;
	  if(nxgrid > MAXXGRID) error_handler("nxgrid is greater than MAXXGRID, increase MAXXGRID");
	}
      }
    }
  }

  free(ctrlon1);
  free(ctrlat1);
  free(ctrlat2);
  
  return nxgrid;
  
}; /* create_xgrid_1dx2d_order2 */

/*******************************************************************************
  void create_xgrid_2dx1d_order1
  This routine generate exchange grids between two grids for the first order
  conservative interpolation. nlon1,nlat1,nlon2,nlat2 are the size of the grid cell
  and lon2,lat2 are 1-D grid bounds, lon1,lat1 are geographic grid location of grid cell bounds.
  mask is on grid lon1/lat1. 
*******************************************************************************/
#ifndef __AIX
int create_xgrid_2dx1d_order1_(const int *nlon1, const int *nlat1, const int *nlon2, const int *nlat2,
			       const double *lon1, const double *lat1, const double *lon2, const double *lat2,
			       const double *mask1, int *grid1_i, int *grid1_j, int *grid2_i,
			       int *grid2_j, double *xgrid_area)
{
  int nxgrid;
  
  nxgrid = create_xgrid_2dx1d_order1(nlon1, nlat1, nlon2, nlat2, lon1, lat1, lon2, lat2, mask1,
			       grid1_i, grid1_j, grid2_i, grid2_j, xgrid_area);
  return nxgrid;
    
};  
#endif
int create_xgrid_2dx1d_order1(const int *nlon1, const int *nlat1, const int *nlon2, const int *nlat2, const double *lon1,
			      const double *lat1, const double *lon2, const double *lat2,
			      const double *mask1, int *grid1_i, int *grid1_j, int *grid2_i,
			      int *grid2_j, double *xgrid_area)
{

  int nx1, ny1, nx2, ny2, nx1p, nx2p;
  int i1, j1, i2, j2, nxgrid;
  double ll_lon, ll_lat, ur_lon, ur_lat, lon_in[MV], lat_in[MV], lon_out[MV], lat_out[MV];
  
  nx1 = *nlon1;
  ny1 = *nlat1;
  nx2 = *nlon2;
  ny2 = *nlat2;

  nxgrid = 0;
  nx1p = nx1 + 1;
  nx2p = nx2 + 1;

  for(j2=0; j2<ny2; j2++) for(i2=0; i2<nx2; i2++) {

    ll_lon = lon2[i2];   ll_lat = lat2[j2];
    ur_lon = lon2[i2+1]; ur_lat = lat2[j2+1];
    for(j1=0; j1<ny1; j1++) for(i1=0; i1<nx1; i1++) if( mask1[j1*nx1+i1] > MASK_THRESH ) {
      int n_in, n_out;
      double Xarea;
      
      lat_in[0] = lat1[j1*nx1p+i1];
      lat_in[1] = lat1[j1*nx1p+i1+1];
      lat_in[2] = lat1[(j1+1)*nx1p+i1+1];
      lat_in[3] = lat1[(j1+1)*nx1p+i1];
      if (  (lat_in[0]<=ll_lat) && (lat_in[1]<=ll_lat)
	    && (lat_in[2]<=ll_lat) && (lat_in[3]<=ll_lat) ) continue;
      if (  (lat_in[0]>=ur_lat) && (lat_in[1]>=ur_lat)
	    && (lat_in[2]>=ur_lat) && (lat_in[3]>=ur_lat) ) continue;

      lon_in[0] = lon1[j1*nx1p+i1];
      lon_in[1] = lon1[j1*nx1p+i1+1];
      lon_in[2] = lon1[(j1+1)*nx1p+i1+1];
      lon_in[3] = lon1[(j1+1)*nx1p+i1];

      n_in = fix_lon(lon_in, lat_in, 4, (ll_lon+ur_lon)/2);
      
      if ( (n_out = clip ( lon_in, lat_in, n_in, ll_lon, ll_lat, ur_lon, ur_lat, lon_out, lat_out )) > 0 ) {
	Xarea = poly_area ( lon_out, lat_out, n_out ) * mask1[j1*nx1+i1];
	if( Xarea > AREA_THRESH ) {
      	  xgrid_area[nxgrid] = Xarea;
	  grid1_i[nxgrid]    = i1;
	  grid1_j[nxgrid]    = j1;
	  grid2_i[nxgrid]    = i2;
	  grid2_j[nxgrid]    = j2;
	  ++nxgrid;
	  if(nxgrid > MAXXGRID) error_handler("nxgrid is greater than MAXXGRID, increase MAXXGRID");
	}
      }
    }
  }

  return nxgrid;
  
}; /* create_xgrid_2dx1d_order1 */


/********************************************************************************
  void create_xgrid_2dx1d_order2
  This routine generate exchange grids between two grids for the second order
  conservative interpolation. nlon1,nlat1,nlon2,nlat2 are the size of the grid cell
  and lon2,lat2 are 1-D grid bounds, lon1,lat1 are geographic grid location of grid cell bounds.
  mask is on grid lon1/lat1. 
********************************************************************************/
#ifndef __AIX
int create_xgrid_2dx1d_order2_(const int *nlon1, const int *nlat1, const int *nlon2, const int *nlat2,
			 const double *lon1, const double *lat1, const double *lon2, const double *lat2,
			 const double *mask1, int *grid1_i, int *grid1_j, int *grid2_i, int *grid2_j,
			 double *xgrid_area, double *grid1_di, double *grid1_dj, double *grid2_di, double *grid2_dj)
{
  int nxgrid;
  nxgrid = create_xgrid_2dx1d_order2(nlon1, nlat1, nlon2, nlat2, lon1, lat1, lon2, lat2, mask1, grid1_i,
                                     grid1_j, grid2_i, grid2_j, xgrid_area, grid1_di, grid1_dj, grid2_di, grid2_dj);
  return nxgrid;

};
#endif
int create_xgrid_2dx1d_order2(const int *nlon1, const int *nlat1, const int *nlon2, const int *nlat2,
			const double *lon1, const double *lat1, const double *lon2, const double *lat2,
			const double *mask1, int *grid1_i, int *grid1_j, int *grid2_i, int *grid2_j,
			double *xgrid_area, double *grid1_di, double *grid1_dj, double *grid2_di, double *grid2_dj)
{

  int nx1, ny1, nx2, ny2, nx1p, nx2p;
  int i1, j1, i2, j2, nxgrid;
  double ll_lon, ll_lat, ur_lon, ur_lat, lon_in[MV], lat_in[MV], lon_out[MV], lat_out[MV];
  double *ctrlon2, *ctrlat2, ctrlon1, *ctrlat1;
  
  nx1 = *nlon1;
  ny1 = *nlat1;
  nx2 = *nlon2;
  ny2 = *nlat2;

  nxgrid = 0;
  nx1p = nx1 + 1;
  nx2p = nx2 + 1;

  ctrlon2    = (double *)malloc(nx2*ny2*sizeof(double));
  ctrlat2    = (double *)malloc(nx2*ny2*sizeof(double));
  ctrlat1    = (double *)malloc(nx1*ny1*sizeof(double));
  
  for(j2=0; j2<ny2; j2++) {
    for(i2=0; i2<nx2; i2++) {
      ll_lon = lon2[i2]; ll_lat = lat2[j2];
      ur_lon = lon2[i2+1]; ur_lat = lat2[j2+1];
      ctrlon2[j2*nx2+i2] = box_ctrlon(ll_lon, ll_lat, ur_lon, ur_lat, 0.5*(ll_lon+ur_lon) );
      ctrlat2[j2*nx2+i2] = box_ctrlat(ll_lon, ll_lat, ur_lon, ur_lat );
    }
  }

  for(j1=0; j1<ny1; j1++) {
    for(i1=0; i1<nx1; i1++) {
      int n_in;
      lon_in[0] = lon1[j1*nx1p+i1];
      lon_in[1] = lon1[j1*nx1p+i1+1];
      lon_in[2] = lon1[(j1+1)*nx1p+i1+1];
      lon_in[3] = lon1[(j1+1)*nx1p+i1];
      lat_in[0] = lat1[j1*nx1p+i1];
      lat_in[1] = lat1[j1*nx1p+i1+1];
      lat_in[2] = lat1[(j1+1)*nx1p+i1+1];
      lat_in[3] = lat1[(j1+1)*nx1p+i1];
      n_in = fix_lon(lon_in, lat_in, 4, M_PI);
      ctrlat1[j1*nx1+i1]    = poly_ctrlat(lon_in, lat_in, n_in);
    }
  }

  for(j2=0; j2<ny2; j2++) for(i2=0; i2<nx2; i2++) {

    ll_lon = lon2[i2];   ll_lat = lat2[j2];
    ur_lon = lon2[i2+1]; ur_lat = lat2[j2+1];
    for(j1=0; j1<ny1; j1++) for(i1=0; i1<nx1; i1++) if( mask1[j1*nx1+i1] > MASK_THRESH ) {
      int n_in, n_out;
      double xarea, xctrlon, xctrlat;
      
      lat_in[0] = lat1[j1*nx1p+i1];
      lat_in[1] = lat1[j1*nx1p+i1+1];
      lat_in[2] = lat1[(j1+1)*nx1p+i1+1];
      lat_in[3] = lat1[(j1+1)*nx1p+i1];
      if (  (lat_in[0]<=ll_lat) && (lat_in[1]<=ll_lat)
	    && (lat_in[2]<=ll_lat) && (lat_in[3]<=ll_lat) ) continue;
      if (  (lat_in[0]>=ur_lat) && (lat_in[1]>=ur_lat)
	    && (lat_in[2]>=ur_lat) && (lat_in[3]>=ur_lat) ) continue;

      lon_in[0] = lon1[j1*nx1p+i1];
      lon_in[1] = lon1[j1*nx1p+i1+1];
      lon_in[2] = lon1[(j1+1)*nx1p+i1+1];
      lon_in[3] = lon1[(j1+1)*nx1p+i1];

      n_in = fix_lon(lon_in, lat_in, 4, (ll_lon+ur_lon)/2);
      
      if (  (n_out = clip ( lon_in, lat_in, n_in, ll_lon, ll_lat, ur_lon, ur_lat, lon_out, lat_out )) > 0 ) {
	xarea = poly_area (lon_out, lat_out, n_out ) * mask1[j1*nx1+i1];
	if(xarea > AREA_THRESH ) {	  
	  xgrid_area[nxgrid] = xarea;
	  grid1_i[nxgrid]    = i1;
	  grid1_j[nxgrid]    = j1;
	  grid2_i[nxgrid]    = i2;
	  grid2_j[nxgrid]    = j2;
	  ctrlon1 = poly_ctrlon(lon_in, lat_in, n_in, 0.5*(ll_lon+ur_lon));
	  xctrlon = poly_ctrlon ( lon_out, lat_out, n_out, 0.5*(ll_lon+ur_lon));
	  xctrlat = poly_ctrlat ( lon_out, lat_out, n_out );
	  grid1_di[nxgrid] = xctrlon - ctrlon1;
	  grid1_dj[nxgrid] = xctrlat - ctrlat1[j1*nx1+i1];
	  grid2_di[nxgrid] = xctrlon - ctrlon2[j2*nx2+i2];
	  grid2_dj[nxgrid] = xctrlat - ctrlat2[j2*nx2+i2];
	  ++nxgrid;
	  if(nxgrid > MAXXGRID) error_handler("nxgrid is greater than MAXXGRID, increase MAXXGRID");
	}
      }
    }
  }

  free(ctrlon2);
  free(ctrlat2);
  free(ctrlat1);
  
  return nxgrid;
  
}; /* create_xgrid_2dx1d_order2 */

/*******************************************************************************
  void create_xgrid_2DX2D_order1
  This routine generate exchange grids between two grids for the first order
  conservative interpolation. nlon1,nlat1,nlon2,nlat2 are the size of the grid cell
  and lon1,lat1, lon2,lat2 are geographic grid location of grid cell bounds.
  mask is on grid lon1/lat1.
*******************************************************************************/
#ifndef __AIX
int create_xgrid_2dx2d_order1_(const int *nlon1, const int *nlat1, const int *nlon2, const int *nlat2,
			       const double *lon1, const double *lat1, const double *lon2, const double *lat2,
			       const double *mask1, int *grid1_i, int *grid1_j, int *grid2_i,
			       int *grid2_j, double *xgrid_area)
{
  int nxgrid;
  
  nxgrid = create_xgrid_2dx2d_order1(nlon1, nlat1, nlon2, nlat2, lon1, lat1, lon2, lat2, mask1,
			       grid1_i, grid1_j, grid2_i, grid2_j, xgrid_area);
  return nxgrid;
    
};  
#endif
int create_xgrid_2dx2d_order1(const int *nlon1, const int *nlat1, const int *nlon2, const int *nlat2,
			      const double *lon1, const double *lat1, const double *lon2, const double *lat2,
			      const double *mask1, int *grid1_i, int *grid1_j, int *grid2_i,
			      int *grid2_j, double *xgrid_area)
{

  int nx1, ny1, nx2, ny2, nx1p, nx2p, ny1p, ny2p, nxgrid, n1_in, n2_in;
  int n0, n1, n2, n3, i1, j1, i2, j2, l;
  double lon1_in[MV], lat1_in[MV], lon2_in[MV], lat2_in[MV], lon_out[MV], lat_out[MV];
  double lon1_min, lon2_min, lon1_max, lon2_max, lat1_min, lat1_max, lat2_min, lat2_max;
  double lon1_avg;

  nx1 = *nlon1;
  ny1 = *nlat1;
  nx2 = *nlon2;
  ny2 = *nlat2;

  nxgrid = 0;
  nx1p = nx1 + 1;
  nx2p = nx2 + 1;
  
  for(j1=0; j1<ny1; j1++) for(i1=0; i1<nx1; i1++) if( mask1[j1*nx1+i1] > MASK_THRESH ) {
    n0 = j1*nx1p+i1;       n1 = j1*nx1p+i1+1;
    n2 = (j1+1)*nx1p+i1+1; n3 = (j1+1)*nx1p+i1;      
    lon1_in[0] = lon1[n0]; lat1_in[0] = lat1[n0];
    lon1_in[1] = lon1[n1]; lat1_in[1] = lat1[n1];
    lon1_in[2] = lon1[n2]; lat1_in[2] = lat1[n2];
    lon1_in[3] = lon1[n3]; lat1_in[3] = lat1[n3];
    lat1_min = minval_double(4, lat1_in);
    lat1_max = maxval_double(4, lat1_in);
    n1_in = fix_lon(lon1_in, lat1_in, 4, M_PI);
    lon1_min = minval_double(n1_in, lon1_in);
    lon1_max = maxval_double(n1_in, lon1_in);
    lon1_avg = avgval_double(n1_in, lon1_in);
      
    for(j2=0; j2<ny2; j2++) for(i2=0; i2<nx2; i2++) {
      int n_in, n_out;
      double Xarea;
      n0 = j2*nx2p+i2; n1 = j2*nx2p+i2+1;
      n2 = (j2+1)*nx2p+i2+1; n3 = (j2+1)*nx2p+i2;
      lon2_in[0] = lon2[n0]; lat2_in[0] = lat2[n0];
      lon2_in[1] = lon2[n1]; lat2_in[1] = lat2[n1];
      lon2_in[2] = lon2[n2]; lat2_in[2] = lat2[n2];
      lon2_in[3] = lon2[n3]; lat2_in[3] = lat2[n3];

      lat2_min = minval_double(4, lat2_in);
      lat2_max = maxval_double(4, lat2_in);
      if(lat2_min >= lat1_max || lat2_max <= lat1_min ) continue;
      n2_in = fix_lon(lon2_in, lat2_in, 4, lon1_avg);
      lon2_min = minval_double(n2_in, lon2_in);
      lon2_max = maxval_double(n2_in, lon2_in);    
      
      /* lon2_in should in the same range as lon1_in after lon_fix, so no need to
	 consider cyclic condition
      */
      if(lon2_min >= lon1_max || lon2_max <= lon1_min ) continue;

      if ( (n_out = clip_2dx2d ( lon1_in, lat1_in, n1_in, lon2_in, lat2_in, n2_in, lon_out, lat_out )) > 0) {
	Xarea = poly_area (lon_out, lat_out, n_out) * mask1[j1*nx1+i1];
	if( Xarea > AREA_THRESH ) {
	  xgrid_area[nxgrid] = Xarea;
	  grid1_i[nxgrid]    = i1;
	  grid1_j[nxgrid]    = j1;
	  grid2_i[nxgrid]    = i2;
	  grid2_j[nxgrid]    = j2;
	  ++nxgrid;
	  if(nxgrid > MAXXGRID) error_handler("nxgrid is greater than MAXXGRID, increase MAXXGRID");
	}
      }
    }
  }

  return nxgrid;
  
};/* get_xgrid_2Dx2D_order1 */

/********************************************************************************
  void create_xgrid_2dx1d_order2
  This routine generate exchange grids between two grids for the second order
  conservative interpolation. nlon1,nlat1,nlon2,nlat2 are the size of the grid cell
  and lon1,lat1, lon2,lat2 are geographic grid location of grid cell bounds.
  mask is on grid lon1/lat1. 
********************************************************************************/
#ifndef __AIX
int create_xgrid_2dx2d_order2_(const int *nlon1, const int *nlat1, const int *nlon2, const int *nlat2,
			 const double *lon1, const double *lat1, const double *lon2, const double *lat2,
			 const double *mask1, int *grid1_i, int *grid1_j, int *grid2_i, int *grid2_j,
			 double *xgrid_area, double *grid1_di, double *grid1_dj, double *grid2_di, double *grid2_dj)
{
  int nxgrid;
  nxgrid = create_xgrid_2dx2d_order2(nlon1, nlat1, nlon2, nlat2, lon1, lat1, lon2, lat2, mask1, grid1_i,
                                     grid1_j, grid2_i, grid2_j, xgrid_area, grid1_di, grid1_dj, grid2_di, grid2_dj);
  return nxgrid;

};
#endif
int create_xgrid_2dx2d_order2(const int *nlon1, const int *nlat1, const int *nlon2, const int *nlat2,
			const double *lon1, const double *lat1, const double *lon2, const double *lat2,
			const double *mask1, int *grid1_i, int *grid1_j, int *grid2_i, int *grid2_j,
			double *xgrid_area, double *grid1_di, double *grid1_dj, double *grid2_di, double *grid2_dj)
{

  int nx1, nx2, ny1, ny2, nx1p, nx2p, ny1p, ny2p, nxgrid, n1_in, n2_in;
  int n0, n1, n2, n3, i1, j1, i2, j2, l;
  double lon1_in[MV], lat1_in[MV], lon2_in[MV], lat2_in[MV], lon_out[MV], lat_out[MV];
  double *ctrlat1, *ctrlat2;
  double lon1_min, lon2_min, lon1_max, lon2_max, lat1_min, lat1_max, lat2_min, lat2_max;
  double lon1_avg, ctrlon1, ctrlon2, xctrlon, xctrlat;

  nx1 = *nlon1;
  ny1 = *nlat1;
  nx2 = *nlon2;
  ny2 = *nlat2;  
  nxgrid = 0;
  nx1p = nx1 + 1;
  nx2p = nx2 + 1;
  ny1p = ny1 + 1;
  ny2p = ny2 + 1;

  ctrlat1    = (double *)malloc(nx1*ny1*sizeof(double));
  ctrlat2    = (double *)malloc(nx2*ny2*sizeof(double));
  
  for(j1=0; j1<ny1; j1++) {
    for(i1=0; i1<nx1; i1++) {
      n0 = j1*nx1p+i1; n1 = j1*nx1p+i1+1;
      n2 = (j1+1)*nx1p+i1+1; n3 = (j1+1)*nx1p+i1;
      lon1_in[0] = lon1[n0]; lat1_in[0] = lat1[n0];
      lon1_in[1] = lon1[n1]; lat1_in[1] = lat1[n1];
      lon1_in[2] = lon1[n2]; lat1_in[2] = lat1[n2];
      lon1_in[3] = lon1[n3]; lat1_in[3] = lat1[n3];
      /* fix_lon should only be called for tripolar grid, may need to add
	 argument to indicate if grid1/grid2 is tripolar or not. */
      n1_in = fix_lon(lon1_in, lat1_in, 4, M_PI);
      ctrlat1[j1*nx1+i1] = poly_ctrlat(lon1_in, lat1_in, n1_in);
    }
  }

  for(j2=0; j2<ny2; j2++) {
    for(i2=0; i2<nx2; i2++) {
      n0 = j2*nx2p+i2; n1 = j2*nx2p+i2+1;
      n2 = (j2+1)*nx2p+i2+1; n3 = (j2+1)*nx2p+i2;
      lon2_in[0] = lon2[n0]; lat2_in[0] = lat2[n0];
      lon2_in[1] = lon2[n1]; lat2_in[1] = lat2[n1];
      lon2_in[2] = lon2[n2]; lat2_in[2] = lat2[n2];
      lon2_in[3] = lon2[n3]; lat2_in[3] = lat2[n3];      
      n2_in = fix_lon(lon2_in, lat2_in, 4, M_PI);
      ctrlat2[j2*nx2+i2] = poly_ctrlat(lon2_in, lat2_in, n2_in);
    }
  }

  for(j1=0; j1<ny1; j1++) for(i1=0; i1<nx1; i1++) if( mask1[j1*nx1+i1] > MASK_THRESH ) {
    n0 = j1*nx1p+i1;       n1 = j1*nx1p+i1+1;
    n2 = (j1+1)*nx1p+i1+1; n3 = (j1+1)*nx1p+i1;      
    lon1_in[0] = lon1[n0]; lat1_in[0] = lat1[n0];
    lon1_in[1] = lon1[n1]; lat1_in[1] = lat1[n1];
    lon1_in[2] = lon1[n2]; lat1_in[2] = lat1[n2];
    lon1_in[3] = lon1[n3]; lat1_in[3] = lat1[n3];
    lat1_min = minval_double(4, lat1_in);
    lat1_max = maxval_double(4, lat1_in);
    n1_in = fix_lon(lon1_in, lat1_in, 4, M_PI);
    lon1_min = minval_double(n1_in, lon1_in);
    lon1_max = maxval_double(n1_in, lon1_in);
    lon1_avg = avgval_double(n1_in, lon1_in);
      
    for(j2=0; j2<ny2; j2++) for(i2=0; i2<nx2; i2++) {
      int n_in, n_out;
      double Xarea;
      n0 = j2*nx2p+i2; n1 = j2*nx2p+i2+1;
      n2 = (j2+1)*nx2p+i2+1; n3 = (j2+1)*nx2p+i2;
      lon2_in[0] = lon2[n0]; lat2_in[0] = lat2[n0];
      lon2_in[1] = lon2[n1]; lat2_in[1] = lat2[n1];
      lon2_in[2] = lon2[n2]; lat2_in[2] = lat2[n2];
      lon2_in[3] = lon2[n3]; lat2_in[3] = lat2[n3];

      lat2_min = minval_double(4, lat2_in);
      lat2_max = maxval_double(4, lat2_in);
      if(lat2_min >= lat1_max || lat2_max <= lat1_min ) continue;
      n2_in = fix_lon(lon2_in, lat2_in, 4, lon1_avg);
      lon2_min = minval_double(n2_in, lon2_in);
      lon2_max = maxval_double(n2_in, lon2_in);    

      /* lon2_in should in the same range as lon1_in after lon_fix, so no need to
	 consider cyclic condition
      */
      if(lon2_min >= lon1_max || lon2_max <= lon1_min ) continue;

      if (  (n_out = clip_2dx2d( lon1_in, lat1_in, n1_in, lon2_in, lat2_in, n2_in, lon_out, lat_out )) > 0) {
	Xarea = poly_area(lon_out, lat_out, n_out ) * mask1[j1*nx1+i1];
	if( Xarea > AREA_THRESH ) {
	  xgrid_area[nxgrid] = Xarea;
	  ctrlon1            = poly_ctrlon(lon1_in, lat1_in, n1_in, lon1_avg);
	  ctrlon2            = poly_ctrlon(lon2_in, lat2_in, n2_in, lon1_avg);
	  xctrlon            = poly_ctrlon ( lon_out, lat_out, n_out, lon1_avg);
	  xctrlat            = poly_ctrlat ( lon_out, lat_out, n_out );	
	  grid1_di[nxgrid]   = xctrlon - ctrlon1;
	  grid1_dj[nxgrid]   = xctrlat - ctrlat1[j1*nx1+i1];
	  grid2_di[nxgrid]   = xctrlon - ctrlon2;
	  grid2_dj[nxgrid]   = xctrlat - ctrlat2[j2*nx2+i2];	
	  grid1_i[nxgrid]    = i1;
	  grid1_j[nxgrid]    = j1;
	  grid2_i[nxgrid]    = i2;
	  grid2_j[nxgrid]    = j2;
	  ++nxgrid;
	  if(nxgrid > MAXXGRID) error_handler("nxgrid is greater than MAXXGRID, increase MAXXGRID");
	}
      }
    }
  }

  free(ctrlat1);
  free(ctrlat2);
  
  return nxgrid;
  
};/* get_xgrid_2Dx2D_order2 */

/*******************************************************************************
   Sutherland-Hodgeman algorithm sequentially clips parts outside 4 boundaries
*******************************************************************************/

int clip(const double lon_in[], const double lat_in[], int n_in, double ll_lon, double ll_lat,
	 double ur_lon, double ur_lat, double lon_out[], double lat_out[])
{
  double x_tmp[MV], y_tmp[MV], x_last, y_last;
  int i_in, i_out, n_out, inside_last, inside;

  /* clip polygon with LEFT boundary - clip V_IN to V_TMP */
  x_last = lon_in[n_in-1];
  y_last = lat_in[n_in-1];
  inside_last = (x_last >= ll_lon);
  for (i_in=0,i_out=0;i_in<n_in;i_in++) {
 
    /* if crossing LEFT boundary - output intersection */
    if ((inside=(lon_in[i_in] >= ll_lon))!=inside_last) {
      x_tmp[i_out] = ll_lon;
      y_tmp[i_out++] = y_last + (ll_lon - x_last) * (lat_in[i_in] - y_last) / (lon_in[i_in] - x_last);
    }

    /* if "to" point is right of LEFT boundary, output it */
    if (inside) {
      x_tmp[i_out]   = lon_in[i_in];
      y_tmp[i_out++] = lat_in[i_in];
    }
    x_last = lon_in[i_in];
    y_last = lat_in[i_in];
    inside_last = inside;
  }
  if (!(n_out=i_out)) return(0);

  /* clip polygon with RIGHT boundary - clip V_TMP to V_OUT */
  x_last = x_tmp[n_out-1];
  y_last = y_tmp[n_out-1];
  inside_last = (x_last <= ur_lon);
  for (i_in=0,i_out=0;i_in<n_out;i_in++) {
 
    /* if crossing RIGHT boundary - output intersection */
    if ((inside=(x_tmp[i_in] <= ur_lon))!=inside_last) {
      lon_out[i_out]   = ur_lon;
      lat_out[i_out++] = y_last + (ur_lon - x_last) * (y_tmp[i_in] - y_last)
                                                 / (x_tmp[i_in] - x_last);
    }

    /* if "to" point is left of RIGHT boundary, output it */
    if (inside) {
      lon_out[i_out]   = x_tmp[i_in];
      lat_out[i_out++] = y_tmp[i_in];
    }
    
    x_last = x_tmp[i_in];
    y_last = y_tmp[i_in];
    inside_last = inside;
  }
  if (!(n_out=i_out)) return(0);

  /* clip polygon with BOTTOM boundary - clip V_OUT to V_TMP */
  x_last = lon_out[n_out-1];
  y_last = lat_out[n_out-1];
  inside_last = (y_last >= ll_lat);
  for (i_in=0,i_out=0;i_in<n_out;i_in++) {
 
    /* if crossing BOTTOM boundary - output intersection */
    if ((inside=(lat_out[i_in] >= ll_lat))!=inside_last) {
      y_tmp[i_out]   = ll_lat;
      x_tmp[i_out++] = x_last + (ll_lat - y_last) * (lon_out[i_in] - x_last) / (lat_out[i_in] - y_last);
    }

    /* if "to" point is above BOTTOM boundary, output it */
    if (inside) {
      x_tmp[i_out]   = lon_out[i_in];
      y_tmp[i_out++] = lat_out[i_in];
    }
    x_last = lon_out[i_in];
    y_last = lat_out[i_in];
    inside_last = inside;
  }
  if (!(n_out=i_out)) return(0);

  /* clip polygon with TOP boundary - clip V_TMP to V_OUT */
  x_last = x_tmp[n_out-1];
  y_last = y_tmp[n_out-1];
  inside_last = (y_last <= ur_lat);
  for (i_in=0,i_out=0;i_in<n_out;i_in++) {
 
    /* if crossing TOP boundary - output intersection */
    if ((inside=(y_tmp[i_in] <= ur_lat))!=inside_last) {
      lat_out[i_out]   = ur_lat;
      lon_out[i_out++] = x_last + (ur_lat - y_last) * (x_tmp[i_in] - x_last) / (y_tmp[i_in] - y_last);
    }

    /* if "to" point is below TOP boundary, output it */
    if (inside) {
      lon_out[i_out]   = x_tmp[i_in];
      lat_out[i_out++] = y_tmp[i_in];
    }
    x_last = x_tmp[i_in];
    y_last = y_tmp[i_in];
    inside_last = inside;
  }
  return(i_out);
}; /* clip */  


/*******************************************************************************
   Revise Sutherland-Hodgeman algorithm to find the vertices of the overlapping
   between any two grid boxes. It return the number of vertices for the exchange grid.
*******************************************************************************/

int clip_2dx2d(const double lon1_in[], const double lat1_in[], int n1_in, 
	 const double lon2_in[], const double lat2_in[], int n2_in, 
	 double lon_out[], double lat_out[])
{
  double lon_tmp[MV], lat_tmp[MV];
  double x1_0, y1_0, x1_1, y1_1, x2_0, y2_0, x2_1, y2_1;
  int i_out, n_out, inside_last, inside, i1, i2;

  /* clip polygon with each boundary of the polygon */
  /* We treat lon1_in/lat1_in as clip polygon and lon2_in/lat2_in as subject polygon */
  n_out = n1_in;
  for(i1=0; i1<n1_in; i1++) {
    lon_tmp[i1] = lon1_in[i1];
    lat_tmp[i1] = lat1_in[i1];
  }
  x2_0 = lon2_in[n2_in-1];
  y2_0 = lat2_in[n2_in-1];
  for(i2=0; i2<n2_in; i2++) {
    x2_1 = lon2_in[i2];
    y2_1 = lat2_in[i2];
    x1_0 = lon_tmp[n_out-1];
    y1_0 = lat_tmp[n_out-1];
    inside_last = inside_edge( x2_0, y2_0, x2_1, y2_1, x1_0, y1_0);
    for(i1=0, i_out=0; i1<n_out; i1++) {
      x1_1 = lon_tmp[i1];
      y1_1 = lat_tmp[i1];
      if((inside = inside_edge(x2_0, y2_0, x2_1, y2_1, x1_1, y1_1)) != inside_last ) {
	/* there is intersection, the line between <x1_0,y1_0> and  <x1_1,y1_1>
	   should not parallel to the line between <x2_0,y2_0> and  <x2_1,y2_1>
	   may need to consider truncation error
	*/
	if( fabs(x2_0 - x2_1) < EPSLN ) {  /* vertical line, x1_0 should not equal x1_1 */
	  lon_out[i_out]   = x2_0;
	  lat_out[i_out++] = y1_0 + (x2_0-x1_0)*(y1_1-y1_0)/(x1_1-x1_0);
	}
	else if( fabs(y2_0 - y2_1) < EPSLN ) {  /* horizontal line, y1_0 should not equal y1_1 */
	  lat_out[i_out] = y2_0;
	  lon_out[i_out++] = x1_0 + (y2_0 - y1_0) * (x1_1 - x1_0) / (y1_1 - y1_0);
	}
	else if( fabs(x1_0 - x1_1) < EPSLN ) {  /* vertical line, x2_0 should not equal x2_1 */
	  lon_out[i_out]   = x1_0;
	  lat_out[i_out++] = y2_0 + (x1_0-x2_0)*(y2_1-y2_0)/(x2_1-x2_0);
	}
	else if( fabs(y1_0 - y1_1) < EPSLN ) {  /* horizontal line, y1_0 should not equal y1_1 */
	  lat_out[i_out] = y1_0;
	  lon_out[i_out++] = x2_0 + (y1_0 - y2_0) * (x2_1 - x2_0) / (y2_1 - y2_0);
	}
	else {
	  double r1, r2, a1, a2;
	  r1 = (y1_1-y1_0)/(x1_1-x1_0);
	  r2 = (y2_1-y2_0)/(x2_1-x2_0);
	  if(r1 == r2) error_handler("xgrid: line between x1_0/y1_0 and x1_1/y1_1 is parallel to the line between x2_0/y2_0 and x2_1/y2_1");
	  a1 = y1_0 - r1*x1_0;
	  a2 = y2_0 - r2*x2_0;
	  lon_out[i_out]   = (a1-a2)/(r2-r1);
	  lat_out[i_out++] = r1 * lon_out[i_out] + a1;
	}
      }
      if(inside) {
	lon_out[i_out]   = x1_1;
	lat_out[i_out++] = y1_1;	
      }
      x1_0 = x1_1;
      y1_0 = y1_1;
      inside_last = inside;
    }
    if(!(n_out=i_out)) return 0;
    for(i1=0; i1<n_out; i1++) {
      lon_tmp[i1] = lon_out[i1];
      lat_tmp[i1] = lat_out[i1];
    }
    /* shift the starting point */
    x2_0 = x2_1;
    y2_0 = y2_1;
  }
  return(n_out);
}; /* clip */
    

/*------------------------------------------------------------------------------
  double poly_int3(const double x[], const double y[], int n)
  This routine is used to calculate the integral which equal the product of
   latitude of the centroid and area of any grid box 
   ---------------------------------------------------------------------------*/

double poly_int3(const double x[], const double y[], int n)
{
  double int3 = 0.0;
  int    i;

  for (i=0;i<n;i++) {
    int ip = (i+1) % n;
    double dx = (x[ip]-x[i]);
    double lat1, lat2;
    lat1 = y[ip];
    lat2 = y[i];
    if      (dx==0.0) continue;
    if(dx > M_PI)  dx = dx - 2.0*M_PI;
    if(dx < -M_PI) dx = dx + 2.0*M_PI;
    if (lat1 == lat2) /* cheap area calculation along latitude */
      int3 -= dx*(cos(lat1) + lat1*sin(lat1));
    else 
      int3 -= dx*(2*sin(lat1)-lat1*cos(lat1) - 2*sin(lat2)+lat2*cos(lat2) )/(lat1-lat2);
  }
  return int3;
}; /* poly_int3 */        

/*------------------------------------------------------------------------------
  double poly_ctrlat(const double x[], const double y[], int n)
  this routine can be used to calculate the latitude of the centroid of any grid box
  ----------------------------------------------------------------------------*/
double poly_ctrlat(const double x[], const double y[], int n)
{
  double area;
  area = poly_area_unit_radius(x, y, n);
  if(area < UNIT_AREA_THRESH)
     return (0.0);
  else  
    return (poly_int3(x, y, n)/area);
}; /* poly_ctrlat */     

/*------------------------------------------------------------------------------
  double poly_int2(const double x[], const double y[], int n, double clon)
  This routine is used to calculate the integral which equal the product of
   lontitude of the centroid and area of any grid box
   ---------------------------------------------------------------------------*/

double poly_int2(const double x[], const double y[], int n, double clon)
{
  double int2 = 0.0;
  int    i;

  clon = clon;
  for (i=0;i<n;i++) {
    int ip = (i+1) % n;
    double phi1, phi2, dphi, lat1, lat2, dphi1, dphi2;
    double f1, f2, fac, fint;
    phi1   = x[ip];
    phi2   = x[i];
    lat1 = y[ip];
    lat2 = y[i];    
    dphi   = phi1 - phi2;
    
    if      (dphi==0.0) continue;

    f1 = 0.5*(cos(lat1)*sin(lat1)+lat1);
    f2 = 0.5*(cos(lat2)*sin(lat2)+lat2);

     /* this will make sure longitude of centroid is at 
        the same interval as the center of any grid */  
    if(dphi > M_PI)  dphi = dphi - 2.0*M_PI;
    if(dphi < -M_PI) dphi = dphi + 2.0*M_PI;
    dphi1 = phi1 - clon;
    if( dphi1 > M_PI) dphi1 -= 2.0*M_PI;
    if( dphi1 <-M_PI) dphi1 += 2.0*M_PI;
    dphi2 = phi2 -clon;
    if( dphi2 > M_PI) dphi2 -= 2.0*M_PI;
    if( dphi2 <-M_PI) dphi2 += 2.0*M_PI;    

    if(abs(dphi2 -dphi1) < M_PI) {
      int2 -= dphi * (dphi1*f1+dphi2*f2)/2.0;
    }
    else {
      if(dphi1 > 0.0)
	fac = M_PI;
      else
	fac = -M_PI;
      fint = f1 + (f2-f1)*(fac-dphi1)/abs(dphi);
      int2 -= 0.5*dphi1*(dphi1-fac)*f1 - 0.5*dphi2*(dphi2+fac)*f2
	+ 0.5*fac*(dphi1+dphi2)*fint;
	}
    
  }
  return (int2);
};   /* poly_int2 */

/*------------------------------------------------------------------------------
  double poly_ctrlon(const double x[], const double y[], int n, double clon)
  this routine can be used to calculate the lontitude of the centroid of any grid box
  ----------------------------------------------------------------------------*/
double poly_ctrlon(const double x[], const double y[], int n, double clon)
{

  double area;
  area = poly_area_unit_radius(x, y, n);
  if(area < UNIT_AREA_THRESH)
     return (0.0);
  else
    return (poly_int2(x, y, n, clon)/area);
} /* poly_ctrlon */


/* -----------------------------------------------------------------------------
   double box_int3(double ll_lon, double ll_lat, double ur_lon, double ur_lat)
   This routine is used to calculate the integral which equal the product of
   latitude of the centroid and area of a uniform grid box
   ---------------------------------------------------------------------------*/
double box_int3(double ll_lon, double ll_lat, double ur_lon, double ur_lat)
{
  double dphi = ur_lon-ll_lon;
 
  if(dphi > M_PI)  dphi = dphi - 2.0*M_PI;
  if(dphi < -M_PI) dphi = dphi + 2.0*M_PI;
  return ( dphi*(cos(ur_lat) + ur_lat*sin(ur_lat)-(cos(ll_lat) +
							 ll_lat*sin(ll_lat))) );
}; /* box_int3 */

/*------------------------------------------------------------------------------
  double box_ctrlat(VTX ll, VTX ur)
  this routine can be used to calculate the latitude of the centroid of a uniform grid box
  ----------------------------------------------------------------------------*/
double box_ctrlat(double ll_lon, double ll_lat, double ur_lon, double ur_lat)
{
  double area;
  area = box_area_unit_radius(ll_lon, ll_lat, ur_lon, ur_lat);
  if(area < UNIT_AREA_THRESH)
     return (0.0);
  else    
    return ( box_int3(ll_lon, ll_lat, ur_lon, ur_lat)/area);
}; /* box_ctrlat */

/*------------------------------------------------------------------------------
  double box_int2(double ll_lon, double ll_lat, double ur_lon, double ur_lat, double clon)
  This routine is used to calculate the integral which equal the product of
   lontitude of the centroid and area of a uniform grid box
   ----------------------------------------------------------------------------*/
double box_int2(double ll_lon, double ll_lat, double ur_lon, double ur_lat, double clon)
{
  double phi1, phi2, dphi, lat1, lat2, dphi1, dphi2;
  double f1, f2, fac, fint;  
  double int2  = 0.0;
  int i;
  clon = clon;  
  for( i =0; i<2; i++) {
    if(i == 0) {
      phi1 = ur_lon;
      phi2 = ll_lon;
      lat1 = lat2 = ll_lat;
    }
    else {
      phi1 = ll_lon;
      phi2 = ur_lon;
      lat1 = lat2 = ur_lat;
    }
    dphi   = phi1 - phi2;
    f1 = 0.5*(cos(lat1)*sin(lat1)+lat1);
    f2 = 0.5*(cos(lat2)*sin(lat2)+lat2);

    if(dphi > M_PI)  dphi = dphi - 2.0*M_PI;
    if(dphi < -M_PI) dphi = dphi + 2.0*M_PI;
    /* make sure the center is in the same grid box. */
    dphi1 = phi1 - clon;
    if( dphi1 > M_PI) dphi1 -= 2.0*M_PI;
    if( dphi1 <-M_PI) dphi1 += 2.0*M_PI;
    dphi2 = phi2 -clon;
    if( dphi2 > M_PI) dphi2 -= 2.0*M_PI;
    if( dphi2 <-M_PI) dphi2 += 2.0*M_PI;    

    if(abs(dphi2 -dphi1) < M_PI) {
      int2 -= dphi * (dphi1*f1+dphi2*f2)/2.0;
    }
    else {
      if(dphi1 > 0.0)
	fac = M_PI;
      else
	fac = -M_PI;
      fint = f1 + (f2-f1)*(fac-dphi1)/abs(dphi);
      int2 -= 0.5*dphi1*(dphi1-fac)*f1 - 0.5*dphi2*(dphi2+fac)*f2
	+ 0.5*fac*(dphi1+dphi2)*fint;
    }
  }
  return (int2);    
} /* box_int2 */

/*------------------------------------------------------------------------------
  double box_ctrlon(double ll_lon, double ll_lat, double ur_lon, double ur_lat, double clon)
  this routine can be used to calculate the lontitude of the centroid of a uniform grid box
  ---------------------------------------------------------------------------*/
double box_ctrlon(double ll_lon, double ll_lat, double ur_lon, double ur_lat, double clon)
{
  double area;
  area = box_area_unit_radius(ll_lon, ll_lat, ur_lon, ur_lat);
  if(area < UNIT_AREA_THRESH)
    return (0.0);
  else      
    return(box_int2(ll_lon, ll_lat, ur_lon, ur_lat, clon)/area);
}   /* box_ctrlon */


/*******************************************************************************
  double grid_box_radius(double *x, double *y, double *z, int n);
  Find the radius of the grid box, the radius is defined the
  maximum distance between any two vertices
*******************************************************************************/ 
double grid_box_radius(const double *x, const double *y, const double *z, int n)
{
  double radius;
  int i, j;
  
  radius = 0;
  for(i=0; i<n-1; i++) {
    for(j=i+1; j<n; j++) {
      radius = max(radius, pow(x[i]-x[j],2.)+pow(y[i]-y[j],2.)+pow(z[i]-z[j],2.));
    }
  }
  
  radius = sqrt(radius);

  return (radius);
  
}; /* grid_box_radius */

/*******************************************************************************
  double dist_between_boxes(const double *x1, const double *y1, const double *z1, int n1,
			    const double *x2, const double *y2, const double *z2, int n2);
  Find the distance between any two grid boxes. The distance is defined by the maximum
  distance between any vertices of these two box
*******************************************************************************/
double dist_between_boxes(const double *x1, const double *y1, const double *z1, int n1,
			  const double *x2, const double *y2, const double *z2, int n2)
{
  double dist;
  int i, j;

  for(i=0; i<n1; i++) {
    for(j=0; j<n2; j++) {   
      dist = max(dist, pow(x1[i]-x2[j],2.)+pow(y1[i]-y2[j],2.)+pow(z1[i]-z2[j],2.));
    }
  }

  dist = sqrt(dist);
  return (dist);

}; /* dist_between_boxes */

/*******************************************************************************
 int inside_edge(double x0, double y0, double x1, double y1, double x, double y)
 determine a point(x,y) is inside or outside a given edge with vertex,
 (x0,y0) and (x1,y1). return 1 if inside and 0 if outside. <y1-y0, -(x1-x0)> is
 the outward edge normal from vertex <x0,y0> to <x1,y1>. <x-x0,y-y0> is the vector
 from <x0,y0> to <x,y>. 
 if Inner produce <x-x0,y-y0>*<y1-y0, -(x1-x0)> > 0, outside, otherwise inside.
 inner product value = 0 also treate as inside.
*******************************************************************************/
 int inside_edge(double x0, double y0, double x1, double y1, double x, double y)
 {
   double product;
   
   product = ( x-x0 )*(y1-y0) + (x0-x1)*(y-y0);

   if(product<EPSLN)
     return 1;
   else
     return 0;
 }; /* inside_edge */

