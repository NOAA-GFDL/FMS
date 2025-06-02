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
#include <string.h>
#include "grid_utils.h"
#include "tree_utils.h"
#include "constant.h"

#ifdef use_libMPI
#include <mpi.h>
#endif

/** \file
 *  \ingroup mosaic
 *  \brief Error handling and other general utilities for @ref mosaic_mod
 */

/***********************************************************
    void error_handler(char *str)
    error handler: will print out error message and then abort
***********************************************************/

void error_handler(const char *msg)
{
  fprintf(stderr, "FATAL Error: %s\n", msg );
#ifdef use_libMPI
  MPI_Abort(MPI_COMM_WORLD, -1);
#else
  exit(1);
#endif
} /* error_handler */


/*******************************************************************************
  double maxval_double(int size, double *data)
  get the maximum value of double array
*******************************************************************************/
double maxval_double(int size, const double *data)
{
  int n;
  double maxval;

  maxval = data[0];
  for(n=1; n<size; n++){
    if( data[n] > maxval ) maxval = data[n];
  }

  return maxval;

} /* maxval_double */


/*******************************************************************************
  double minval_double(int size, double *data)
  get the minimum value of double array
*******************************************************************************/
double minval_double(int size, const double *data)
{
  int n;
  double minval;

  minval = data[0];
  for(n=1; n<size; n++){
    if( data[n] < minval ) minval = data[n];
  }

  return minval;

} /* minval_double */

/*******************************************************************************
  double avgval_double(int size, double *data)
  get the average value of double array
*******************************************************************************/
double avgval_double(int size, const double *data)
{
  int n;
  double avgval;

  avgval = 0;
  for(n=0; n<size; n++) avgval += data[n];
  avgval /= size;

  return avgval;

} /* avgval_double */


/*******************************************************************************
  void latlon2xyz
  Routine to map (lon, lat) to (x,y,z)
******************************************************************************/
void latlon2xyz(int size, const double *lon, const double *lat, double *x, double *y, double *z)
{
  int n;

  for(n=0; n<size; n++) {
    x[n] = cos(lat[n])*cos(lon[n]);
    y[n] = cos(lat[n])*sin(lon[n]);
    z[n] = sin(lat[n]);
  }

} /* latlon2xyz */

/*------------------------------------------------------------
  void xyz2laton(np, p, xs, ys)
  Transfer cartesian coordinates to spherical coordinates
  ----------------------------------------------------------*/
void xyz2latlon( int np, const double *x, const double *y, const double *z, double *lon, double *lat)
{

  double xx, yy, zz;
  double dist;
  int i;

  for(i=0; i<np; i++) {
    xx = x[i];
    yy = y[i];
    zz = z[i];
    dist = sqrt(xx*xx+yy*yy+zz*zz);
    xx /= dist;
    yy /= dist;
    zz /= dist;

    if ( fabs(xx)+fabs(yy)  < EPSLN10 )
      lon[i] = 0;
    else
      lon[i] = atan2(yy, xx);
    lat[i] = asin(zz);

    if ( lon[i] < 0.) lon[i] = 2.*M_PI + lon[i];
  }

} /* xyz2latlon */

int delete_vtx(double x[], double y[], int n, int n_del)
{
  for (;n_del<n-1;n_del++) {
    x[n_del] = x[n_del+1];
    y[n_del] = y[n_del+1];
  }

  return (n-1);
} /* delete_vtx */

int insert_vtx(double x[], double y[], int n, int n_ins, double lon_in, double lat_in)
{
  int i;

  for (i=n-1;i>=n_ins;i--) {
    x[i+1] = x[i];
    y[i+1] = y[i];
  }

  x[n_ins] = lon_in;
  y[n_ins] = lat_in;
  return (n+1);
} /* insert_vtx */

void v_print(double x[], double y[], int n)
{
  int i;

  for (i=0;i<n;i++){
    printf(" %20g   %20g\n", x[i], y[i]);
  }
} /* v_print */

int fix_lon(double x[], double y[], int n, double tlon)
{
  double x_sum, dx, dx_;
  int i, nn = n, pole = 0;

  for (i=0;i<nn;i++) if (fabs(y[i])>=HPI-TOLERANCE) pole = 1;
  if (0&&pole) {
    printf("fixing pole cell\n");
    v_print(x, y, nn);
    printf("---------");
  }

  /* all pole points must be paired */
  for (i=0;i<nn;i++) if (fabs(y[i])>=HPI-TOLERANCE) {
      int im=(i+nn-1)%nn, ip=(i+1)%nn;

      if (y[im]==y[i] && y[ip]==y[i]) {
        nn = delete_vtx(x, y, nn, i);
        i--;
      } else if (y[im]!=y[i] && y[ip]!=y[i]) {
        nn = insert_vtx(x, y, nn, i, x[i], y[i]);
        i++;
      }
    }
  /* first of pole pair has longitude of previous vertex */
  /* second of pole pair has longitude of subsequent vertex */
  for (i=0;i<nn;i++) if (fabs(y[i])>=HPI-TOLERANCE) {
      int im=(i+nn-1)%nn, ip=(i+1)%nn;

      if (y[im]!=y[i]){
        x[i] = x[im];
      }
      if (y[ip]!=y[i]){
        x[i] = x[ip];
      }
    }

  if (nn){
    x_sum = x[0];
  }
  else{
    return(0);
  }
  for (i=1;i<nn;i++) {
    dx_ = x[i]-x[i-1];

    if (dx_ < -M_PI) dx_ = dx_ + TPI;
    else if (dx_ >  M_PI) dx_ = dx_ - TPI;
    x_sum += (x[i] = x[i-1] + dx_);
  }

  dx = (x_sum/nn)-tlon;
  if (dx < -M_PI){
    for (i=0;i<nn;i++){
      x[i] += TPI;
    }
  }
  else if (dx >  M_PI){
    for (i=0;i<nn;i++){
      x[i] -= TPI;
    }
  }

  if (0&&pole) {
    printf("area=%g\n", poly_area(x, y,nn));
    v_print(x, y, nn);
    printf("---------");
  }

  return (nn);
} /* fix_lon */


/*------------------------------------------------------------------------------
  double great_circle_distance()
  computes distance between two points along a great circle
  (the shortest distance between 2 points on a sphere)
  returned in units of meter
  ----------------------------------------------------------------------------*/
double great_circle_distance(double *p1, double *p2)
{
  double dist, beta;

  /* This algorithm is not accurate for small distance
     dist = RADIUS*ACOS(SIN(p1[1])*SIN(p2[1]) + COS(p1[1])*COS(p2[1])*COS(p1[0]-p2[0]));
  */
  beta = 2.*asin( sqrt( sin((p1[1]-p2[1])/2.)*sin((p1[1]-p2[1])/2.) +
                        cos(p1[1])*cos(p2[1])*(sin((p1[0]-p2[0])/2.)*sin((p1[0]-p2[0])/2.)) ) );
  dist = RADIUS*beta;
  return dist;

} /* great_circle_distance */

/*------------------------------------------------------------------------------
  double spherical_angle(const double *p1, const double *p2, const double *p3)
  p3
  /
  /
  p1 ---> angle
  \
  \
  p2
  -----------------------------------------------------------------------------*/
double spherical_angle(const double *v1, const double *v2, const double *v3)
{
  double angle;
  long double px, py, pz, qx, qy, qz, ddd;

  /* vector product between v1 and v2 */
  px = v1[1]*v2[2] - v1[2]*v2[1];
  py = v1[2]*v2[0] - v1[0]*v2[2];
  pz = v1[0]*v2[1] - v1[1]*v2[0];
  /* vector product between v1 and v3 */
  qx = v1[1]*v3[2] - v1[2]*v3[1];
  qy = v1[2]*v3[0] - v1[0]*v3[2];
  qz = v1[0]*v3[1] - v1[1]*v3[0];

  ddd = (px*px+py*py+pz*pz)*(qx*qx+qy*qy+qz*qz);
  if ( ddd <= 0.0 )
    angle = 0. ;
  else {
    ddd = (px*qx+py*qy+pz*qz) / sqrtl(ddd);
    if( fabsl(ddd-1) < EPSLN30 ) ddd = 1;
    if( fabsl(ddd+1) < EPSLN30 ) ddd = -1;
    if ( ddd>1. || ddd<-1. ) {
      /*FIX (lmh) to correctly handle co-linear points (angle near pi or 0) */
      if (ddd < 0.)
        angle = M_PI;
      else
        angle = 0.;
    }
    else
      angle = ((double)acosl( ddd ));
  }

  return angle;
} /* spherical_angle */


/*----------------------------------------------------------------------
  void vect_cross(e, p1, p2)
  Perform cross products of 3D vectors: e = P1 X P2
  -------------------------------------------------------------------*/

void vect_cross(const double *p1, const double *p2, double *e )
{

  e[0] = p1[1]*p2[2] - p1[2]*p2[1];
  e[1] = p1[2]*p2[0] - p1[0]*p2[2];
  e[2] = p1[0]*p2[1] - p1[1]*p2[0];

} /* vect_cross */


/*----------------------------------------------------------------------
  double* vect_cross(p1, p2)
  return cross products of 3D vectors: = P1 X P2
  -------------------------------------------------------------------*/

double dot(const double *p1, const double *p2)
{

  return( p1[0]*p2[0] + p1[1]*p2[1] + p1[2]*p2[2] );

}


double metric(const double *p) {
  return (sqrt(p[0]*p[0] + p[1]*p[1]+p[2]*p[2]) );
}

void normalize_vect(double *e)
{
  double pdot;
  int k;

  pdot = e[0]*e[0] + e[1] * e[1] + e[2] * e[2];
  pdot = sqrt( pdot );

  for(k=0; k<3; k++) e[k] /= pdot;
}


/*------------------------------------------------------------------
  void unit_vect_latlon(int size, lon, lat, vlon, vlat)
  calculate unit vector for latlon in cartesian coordinates
  ---------------------------------------------------------------------*/
void unit_vect_latlon(int size, const double *lon, const double *lat, double *vlon, double *vlat)
{
  double sin_lon, cos_lon, sin_lat, cos_lat;
  int n;

  for(n=0; n<size; n++) {
    sin_lon = sin(lon[n]);
    cos_lon = cos(lon[n]);
    sin_lat = sin(lat[n]);
    cos_lat = cos(lat[n]);

    vlon[3*n] = -sin_lon;
    vlon[3*n+1] =  cos_lon;
    vlon[3*n+2] =  0.;

    vlat[3*n]   = -sin_lat*cos_lon;
    vlat[3*n+1] = -sin_lat*sin_lon;
    vlat[3*n+2] =  cos_lat;
  }
} /* unit_vect_latlon */


/* Intersect a line and a plane
   Intersects between the plane ( three points ) (entries in counterclockwise order)
   and the line determined by the endpoints l1 and l2 (t=0.0 at l1 and t=1.0 at l2)
   returns true if the two intersect and the output variables are valid
   outputs p containing the coordinates in the tri and t the coordinate in the line
   of the intersection.
   NOTE: the intersection doesn't have to be inside the tri or line for this to return true
*/
int intersect_tri_with_line(const double *plane, const double *l1, const double *l2, double *p,
                            double *t) {

  long double M[3*3], inv_M[3*3];
  long double V[3];
  long double X[3];
  int is_invert=0;

  const double *pnt0=plane;
  const double *pnt1=plane+3;
  const double *pnt2=plane+6;

  /* To do intersection just solve the set of linear equations for both
     Setup M
  */
  M[0]=l1[0]-l2[0]; M[1]=pnt1[0]-pnt0[0]; M[2]=pnt2[0]-pnt0[0];
  M[3]=l1[1]-l2[1]; M[4]=pnt1[1]-pnt0[1]; M[5]=pnt2[1]-pnt0[1];
  M[6]=l1[2]-l2[2]; M[7]=pnt1[2]-pnt0[2]; M[8]=pnt2[2]-pnt0[2];


  /* Invert M */
  is_invert = invert_matrix_3x3(M,inv_M);
  if (!is_invert) return 0;

  /* Set variable holding vector */
  V[0]=l1[0]-pnt0[0];
  V[1]=l1[1]-pnt0[1];
  V[2]=l1[2]-pnt0[2];

  /* Calculate solution */
  mult(inv_M, V, X);

  /* Get answer out */
  *t=((double)X[0]);
  p[0]=((double)X[1]);
  p[1]=((double)X[2]);

  return 1;
}


void mult(long double m[], long double v[], long double out_v[]) {

  out_v[0]=m[0]*v[0]+m[1]*v[1]+m[2]*v[2];
  out_v[1]=m[3]*v[0]+m[4]*v[1]+m[5]*v[2];
  out_v[2]=m[6]*v[0]+m[7]*v[1]+m[8]*v[2];

}


/* returns 1 if matrix is inverted, 0 otherwise */
int invert_matrix_3x3(long double m[], long double m_inv[]) {


  const long double det =  m[0] * (m[4]*m[8] - m[5]*m[7])
      -m[1] * (m[3]*m[8] - m[5]*m[6])
      +m[2] * (m[3]*m[7] - m[4]*m[6]);
  if (fabsl(det) < EPSLN15 ) return 0;

  const long double deti = 1.0/det;

  m_inv[0] = (m[4]*m[8] - m[5]*m[7]) * deti;
  m_inv[1] = (m[2]*m[7] - m[1]*m[8]) * deti;
  m_inv[2] = (m[1]*m[5] - m[2]*m[4]) * deti;

  m_inv[3] = (m[5]*m[6] - m[3]*m[8]) * deti;
  m_inv[4] = (m[0]*m[8] - m[2]*m[6]) * deti;
  m_inv[5] = (m[2]*m[3] - m[0]*m[5]) * deti;

  m_inv[6] = (m[3]*m[7] - m[4]*m[6]) * deti;
  m_inv[7] = (m[1]*m[6] - m[0]*m[7]) * deti;
  m_inv[8] = (m[0]*m[4] - m[1]*m[3]) * deti;

  return 1;
}


int inside_a_polygon(double *lon1, double *lat1, int *npts, double *lon2, double *lat2)
{

  double x2[20], y2[20], z2[20];
  double x1, y1, z1;
  double min_x2, max_x2, min_y2, max_y2, min_z2, max_z2;
  int isinside, i;

  struct Node *grid1=NULL, *grid2=NULL;

  /* first convert to cartesian grid */
  latlon2xyz(*npts, lon2, lat2, x2, y2, z2);
  latlon2xyz(1, lon1, lat1, &x1, &y1, &z1);

  max_x2 = maxval_double(*npts, x2);
  if(x1 >= max_x2+RANGE_CHECK_CRITERIA) return 0;
  min_x2 = minval_double(*npts, x2);
  if(min_x2 >= x1+RANGE_CHECK_CRITERIA) return 0;

  max_y2 = maxval_double(*npts, y2);
  if(y1 >= max_y2+RANGE_CHECK_CRITERIA) return 0;
  min_y2 = minval_double(*npts, y2);
  if(min_y2 >= y1+RANGE_CHECK_CRITERIA) return 0;

  max_z2 = maxval_double(*npts, z2);
  if(z1 >= max_z2+RANGE_CHECK_CRITERIA) return 0;
  min_z2 = minval_double(*npts, z2);
  if(min_z2 >= z1+RANGE_CHECK_CRITERIA) return 0;


  /* add x2,y2,z2 to a Node */
  rewindList();
  grid1 = getNext();
  grid2 = getNext();

  addEnd(grid1, x1, y1, z1, 0, 0, 0, -1);
  for(i=0; i<*npts; i++) addEnd(grid2, x2[i], y2[i], z2[i], 0, 0, 0, -1);

  isinside = insidePolygon(grid1, grid2);

  return isinside;

}

int inside_a_polygon_(double *lon1, double *lat1, int *npts, double *lon2, double *lat2)
{

  int isinside;

  isinside = inside_a_polygon(lon1, lat1, npts, lon2, lat2);

  return isinside;

}

double get_global_area(void)
{
  double garea;
  garea = 4*M_PI*RADIUS*RADIUS;

  return garea;
}

double get_global_area_(void)
{
  double garea;
  garea = 4*M_PI*RADIUS*RADIUS;

  return garea;
}

double poly_area(const double x[], const double y[], int n)
{
  double area = 0.0;
  int    i;

  for (i=0;i<n;i++) {
    int ip = (i+1) % n;
    double dx = (x[ip]-x[i]);
    double lat1, lat2;
    double dy, dat;

    lat1 = y[ip];
    lat2 = y[i];
    if(dx > M_PI)  dx = dx - 2.0*M_PI;
    if(dx < -M_PI) dx = dx + 2.0*M_PI;
    if (dx==0.0) continue;

    if ( fabs(lat1-lat2) < SMALL_VALUE) /* cheap area calculation along latitude */
      area -= dx*sin(0.5*(lat1+lat2));
    else {
      dy = 0.5*(lat1-lat2);
      dat = sin(dy)/dy;
      area -= dx*sin(0.5*(lat1+lat2))*dat;
    }
  }
  if(area < 0)
    return -area*RADIUS*RADIUS;
  else
    return area*RADIUS*RADIUS;

} /* poly_area */

double poly_area_no_adjust(const double x[], const double y[], int n)
{
  double area = 0.0;
  int    i;

  for (i=0;i<n;i++) {
    int ip = (i+1) % n;
    double dx = (x[ip]-x[i]);
    double lat1, lat2;

    lat1 = y[ip];
    lat2 = y[i];
    if (dx==0.0) continue;

    if ( fabs(lat1-lat2) < SMALL_VALUE) /* cheap area calculation along latitude */
      area -= dx*sin(0.5*(lat1+lat2));
    else
      area += dx*(cos(lat1)-cos(lat2))/(lat1-lat2);
  }
  if(area < 0)
    return area*RADIUS*RADIUS;
  else
    return area*RADIUS*RADIUS;
} /* poly_area_no_adjust */

/*------------------------------------------------------------------------------
  double poly_area(const x[], const y[], int n)
  obtains area of input polygon by line integrating -sin(lat)d(lon)
  Vertex coordinates must be in degrees.
  Vertices must be listed counter-clockwise around polygon.
  grid is in radians.
  ----------------------------------------------------------------------------*/
double poly_area_dimensionless(const double x[], const double y[], int n)
{
  double area = 0.0;
  int    i;

  for (i=0;i<n;i++) {
    int ip = (i+1) % n;
    double dx = (x[ip]-x[i]);
    double lat1, lat2;
    double dy, dat;

    lat1 = y[ip];
    lat2 = y[i];
    if(dx > M_PI)  dx = dx - 2.0*M_PI;
    if(dx < -M_PI) dx = dx + 2.0*M_PI;
    if (dx==0.0) continue;

    if ( fabs(lat1-lat2) < SMALL_VALUE) /* cheap area calculation along latitude */
      area -= dx*sin(0.5*(lat1+lat2));
    else {
      dy = 0.5*(lat1-lat2);
      dat = sin(dy)/dy;
      area -= dx*sin(0.5*(lat1+lat2))*dat;
    }
  }
  if(area < 0)
    return (-area/(4*M_PI));
  else
    return (area/(4*M_PI));

} /* poly_area */

/* Compute the great circle area of a polygon on a sphere */
double great_circle_area(int n, const double *x, const double *y, const double *z) {
  int i;
  double pnt0[3], pnt1[3], pnt2[3];
  double sum, area;

  /* sum angles around polygon */
  sum=0.0;
  for ( i=0; i<n; i++) {
    /* points that make up a side of polygon */
    pnt0[0] = x[i];
    pnt0[1] = y[i];
    pnt0[2] = z[i];
    pnt1[0] = x[(i+1)%n];
    pnt1[1] = y[(i+1)%n];
    pnt1[2] = z[(i+1)%n];
    pnt2[0] = x[(i+2)%n];
    pnt2[1] = y[(i+2)%n];
    pnt2[2] = z[(i+2)%n];

    /* compute angle for pnt1 */
    sum += spherical_angle(pnt1, pnt2, pnt0);

  }
  area = (sum - (n-2.)*M_PI) * RADIUS* RADIUS;
  return area;
}

/*------------------------------------------------------------------------------
  double spherical_excess_area(p_lL, p_uL, p_lR, p_uR)
  get the surface area of a cell defined as a quadrilateral
  on the sphere. Area is computed as the spherical excess
  [area units are m^2]
  ----------------------------------------------------------------------------*/
double spherical_excess_area(const double* p_ll, const double* p_ul,
                             const double* p_lr, const double* p_ur, double radius)
{
  double area, ang1, ang2, ang3, ang4;
  double v1[3], v2[3], v3[3];

  /*   S-W: 1   */
  latlon2xyz(1, p_ll, p_ll+1, v1, v1+1, v1+2);
  latlon2xyz(1, p_lr, p_lr+1, v2, v2+1, v2+2);
  latlon2xyz(1, p_ul, p_ul+1, v3, v3+1, v3+2);
  ang1 = spherical_angle(v1, v2, v3);

  /*   S-E: 2   */
  latlon2xyz(1, p_lr, p_lr+1, v1, v1+1, v1+2);
  latlon2xyz(1, p_ur, p_ur+1, v2, v2+1, v2+2);
  latlon2xyz(1, p_ll, p_ll+1, v3, v3+1, v3+2);
  ang2 = spherical_angle(v1, v2, v3);

  /*   N-E: 3   */
  latlon2xyz(1, p_ur, p_ur+1, v1, v1+1, v1+2);
  latlon2xyz(1, p_ul, p_ul+1, v2, v2+1, v2+2);
  latlon2xyz(1, p_lr, p_lr+1, v3, v3+1, v3+2);
  ang3 = spherical_angle(v1, v2, v3);

  /*   N-W: 4   */
  latlon2xyz(1, p_ul, p_ul+1, v1, v1+1, v1+2);
  latlon2xyz(1, p_ur, p_ur+1, v2, v2+1, v2+2);
  latlon2xyz(1, p_ll, p_ll+1, v3, v3+1, v3+2);
  ang4 = spherical_angle(v1, v2, v3);

  area = (ang1 + ang2 + ang3 + ang4 - 2.*M_PI) * radius* radius;

  return area;

} /* spherical_excess_area */

/*******************************************************************************
void get_grid_area(const int *nlon, const int *nlat, const double *lon, const double *lat, const double *area)
  return the grid area.
*******************************************************************************/
void get_grid_area_(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area)
{
  get_grid_area(nlon, nlat, lon, lat, area);
}

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

}  /* get_grid_area */


/*******************************************************************************
void get_grid_area_ug(const int *npts, const double *lon, const double *lat, const double *area)
  return the grid area.
*******************************************************************************/
void get_grid_area_ug_(const int *npts, const double *lon, const double *lat, double *area)
{
  get_grid_area_ug(npts, lon, lat, area);
}

void get_grid_area_ug(const int *npts, const double *lon, const double *lat, double *area)
{
  int nl, l, n_in, nv;
  double x_in[20], y_in[20];

  nl = *npts;
  nv = 4;

  for(l=0; l<nl; l++) {
      x_in[0] = lon[l*nv];
      x_in[1] = lon[l*nv+1];
      x_in[2] = lon[l*nv+2];
      x_in[3] = lon[l*nv+3];
      y_in[0] = lat[l*nv];
      y_in[1] = lat[l*nv+1];
      y_in[2] = lat[l*nv+2];
      y_in[3] = lat[l*nv+3];
      n_in = fix_lon(x_in, y_in, nv, M_PI);
      area[l] = poly_area(x_in, y_in, n_in);
    }

}  /* get_grid_area_ug */


void get_grid_great_circle_area_(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area)
{
  get_grid_great_circle_area(nlon, nlat, lon, lat, area);

}

void get_grid_great_circle_area(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area)
{
  int nx, ny, nxp, nyp, i, j;
  int n0, n1, n2, n3;
  struct Node *grid=NULL;
  double *x=NULL, *y=NULL, *z=NULL;


  nx = *nlon;
  ny = *nlat;
  nxp = nx + 1;
  nyp = ny + 1;

  x = (double *)malloc(nxp*nyp*sizeof(double));
  y = (double *)malloc(nxp*nyp*sizeof(double));
  z = (double *)malloc(nxp*nyp*sizeof(double));

  latlon2xyz(nxp*nyp, lon, lat, x, y, z);

  for(j=0; j<ny; j++) for(i=0; i < nx; i++) {
      /* clockwise */
      n0 = j*nxp+i;
      n1 = (j+1)*nxp+i;
      n2 = (j+1)*nxp+i+1;
      n3 = j*nxp+i+1;
      rewindList();
      grid = getNext();
      addEnd(grid, x[n0], y[n0], z[n0], 0, 0, 0, -1);
      addEnd(grid, x[n1], y[n1], z[n1], 0, 0, 0, -1);
      addEnd(grid, x[n2], y[n2], z[n2], 0, 0, 0, -1);
      addEnd(grid, x[n3], y[n3], z[n3], 0, 0, 0, -1);
      area[j*nx+i] = gridArea(grid);
    }

  free(x);
  free(y);
  free(z);

}  /* get_grid_great_circle_area */

void get_grid_great_circle_area_ug_(const int *npts, const double *lon, const double *lat, double *area)
{
  get_grid_great_circle_area_ug(npts, lon, lat, area);

}

void get_grid_great_circle_area_ug(const int *npts, const double *lon, const double *lat, double *area)
{
  int l, nl, nv;
  int n0, n1, n2, n3;
  struct Node *grid=NULL;
  double *x=NULL, *y=NULL, *z=NULL;

  nl = *npts;
  nv = 4;

  x = (double *)malloc(nl*nv*sizeof(double));
  y = (double *)malloc(nl*nv*sizeof(double));
  z = (double *)malloc(nl*nv*sizeof(double));

  latlon2xyz(nl*nv, lon, lat, x, y, z);

  for(l=0; l<nv; l++) {
      /* clockwise */
      n0 = l*nv;
      n1 = l*nv+1;
      n2 = l*nv+2;
      n3 = l*nv+3;
      rewindList();
      grid = getNext();
      addEnd(grid, x[n0], y[n0], z[n0], 0, 0, 0, -1);
      addEnd(grid, x[n1], y[n1], z[n1], 0, 0, 0, -1);
      addEnd(grid, x[n2], y[n2], z[n2], 0, 0, 0, -1);
      addEnd(grid, x[n3], y[n3], z[n3], 0, 0, 0, -1);
      area[l] = gridArea(grid);
    }

  free(x);
  free(y);
  free(z);

}  /* get_grid_great_circle_area_ug */

void get_grid_area_dimensionless(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area)
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
      area[j*nx+i] = poly_area_dimensionless(x_in, y_in, n_in);
    }

}  /* get_grid_area */



void get_grid_area_no_adjust(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area)
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
      n_in = 4;
      area[j*nx+i] = poly_area_no_adjust(x_in, y_in, n_in);
    }

}  /* get_grid_area_no_adjust */


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
} /* clip */


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
  double dx1, dy1, dx2, dy2, determ, ds1, ds2;
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
           may need to consider truncation error */
        dy1 = y1_1-y1_0;
        dy2 = y2_1-y2_0;
        dx1 = x1_1-x1_0;
        dx2 = x2_1-x2_0;
        ds1 = y1_0*x1_1 - y1_1*x1_0;
        ds2 = y2_0*x2_1 - y2_1*x2_0;
        determ = dy2*dx1 - dy1*dx2;
        if(fabs(determ) < EPSLN30) {
          error_handler("the line between <x1_0,y1_0> and  <x1_1,y1_1> should not parallel to "
                        "the line between <x2_0,y2_0> and  <x2_1,y2_1>");
        }
        lon_out[i_out]   = (dx2*ds1 - dx1*ds2)/determ;
        lat_out[i_out++] = (dy2*ds1 - dy1*ds2)/determ;


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
} /* clip */


/*******************************************************************************
   Revise Sutherland-Hodgeman algorithm to find the vertices of the overlapping
   between any two grid boxes. It return the number of vertices for the exchange grid.
   Each edge of grid box is a part of great circle. All the points are cartesian
   coordinates. Here we are assuming each polygon is convex.
   RANGE_CHECK_CRITERIA is used to determine if the two grid boxes are possible to be
   overlap. The size should be between 0 and 0.5. The larger the range_check_criteria,
   the more expensive of the computatioin. When the value is close to 0,
   some small exchange grid might be lost. Suggest to use value 0.05 for C48.
*******************************************************************************/

int clip_2dx2d_great_circle(const double x1_in[], const double y1_in[], const double z1_in[], int n1_in,
                            const double x2_in[], const double y2_in[], const double z2_in [], int n2_in,
                            double x_out[], double y_out[], double z_out[])
{
  struct Node *grid1List=NULL;
  struct Node *grid2List=NULL;
  struct Node *intersectList=NULL;
  struct Node *polyList=NULL;
  struct Node *curList=NULL;
  struct Node *firstIntersect=NULL, *curIntersect=NULL;
  struct Node *temp1=NULL, *temp2=NULL, *temp=NULL;

  int    i1, i2, i1p, i2p, i2p2, npts1, npts2;
  int    nintersect, n_out;
  int    maxiter1, maxiter2, iter1, iter2;
  int    found1, found2, curListNum;
  int    has_inbound, inbound;
  double pt1[MV][3], pt2[MV][3];
  double *p1_0=NULL, *p1_1=NULL;
  double *p2_0=NULL, *p2_1=NULL, *p2_2=NULL;
  double intersect[3];
  double u1, u2;
  double min_x1, max_x1, min_y1, max_y1, min_z1, max_z1;
  double min_x2, max_x2, min_y2, max_y2, min_z2, max_z2;


  /* first check the min and max of (x1_in, y1_in, z1_in) with (x2_in, y2_in, z2_in) */
  min_x1 = minval_double(n1_in, x1_in);
  max_x2 = maxval_double(n2_in, x2_in);
  if(min_x1 >= max_x2+RANGE_CHECK_CRITERIA) return 0;
  max_x1 = maxval_double(n1_in, x1_in);
  min_x2 = minval_double(n2_in, x2_in);
  if(min_x2 >= max_x1+RANGE_CHECK_CRITERIA) return 0;

  min_y1 = minval_double(n1_in, y1_in);
  max_y2 = maxval_double(n2_in, y2_in);
  if(min_y1 >= max_y2+RANGE_CHECK_CRITERIA) return 0;
  max_y1 = maxval_double(n1_in, y1_in);
  min_y2 = minval_double(n2_in, y2_in);
  if(min_y2 >= max_y1+RANGE_CHECK_CRITERIA) return 0;

  min_z1 = minval_double(n1_in, z1_in);
  max_z2 = maxval_double(n2_in, z2_in);
  if(min_z1 >= max_z2+RANGE_CHECK_CRITERIA) return 0;
  max_z1 = maxval_double(n1_in, z1_in);
  min_z2 = minval_double(n2_in, z2_in);
  if(min_z2 >= max_z1+RANGE_CHECK_CRITERIA) return 0;

  rewindList();

  grid1List = getNext();
  grid2List = getNext();
  intersectList = getNext();
  polyList = getNext();

  /* insert points into SubjList and ClipList */
  for(i1=0; i1<n1_in; i1++) addEnd(grid1List, x1_in[i1], y1_in[i1], z1_in[i1], 0, 0, 0, -1);
  for(i2=0; i2<n2_in; i2++) addEnd(grid2List, x2_in[i2], y2_in[i2], z2_in[i2], 0, 0, 0, -1);
  npts1 = length(grid1List);
  npts2 = length(grid2List);

  n_out = 0;
  /* set the inside value */
  /* first check number of points in grid1 is inside grid2 */

  temp = grid1List;
  while(temp) {
    if(insidePolygon(temp, grid2List))
      temp->isInside = 1;
    else
      temp->isInside = 0;
    temp = getNextNode(temp);
  }

  /* check if grid2List is inside grid1List */
  temp = grid2List;

  while(temp) {
    if(insidePolygon(temp, grid1List))
      temp->isInside = 1;
    else
      temp->isInside = 0;
    temp = getNextNode(temp);
  }

  /* make sure the grid box is clockwise */

  /*make sure each polygon is convex, which is equivalent that the great_circle_area is positive */
  if( gridArea(grid1List) <= 0 )
    error_handler("create_xgrid.c(clip_2dx2d_great_circle): grid box 1 is not convex");
  if( gridArea(grid2List) <= 0 )
    error_handler("create_xgrid.c(clip_2dx2d_great_circle): grid box 2 is not convex");

  /* get the coordinates from grid1List and grid2List.
     Please not npts1 might not equal n1_in, npts2 might not equal n2_in because of pole
  */

  temp = grid1List;
  for(i1=0; i1<npts1; i1++) {
    getCoordinates(temp, pt1[i1]);
    temp = temp->Next;
  }
  temp = grid2List;
  for(i2=0; i2<npts2; i2++) {
    getCoordinates(temp, pt2[i2]);
    temp = temp->Next;
  }

  firstIntersect=getNext();
  curIntersect = getNext();

  /* first find all the intersection points */
  nintersect = 0;
  for(i1=0; i1<npts1; i1++) {
    i1p = (i1+1)%npts1;
    p1_0 = pt1[i1];
    p1_1 = pt1[i1p];
    for(i2=0; i2<npts2; i2++) {
      i2p = (i2+1)%npts2;
      i2p2 = (i2+2)%npts2;
      p2_0 = pt2[i2];
      p2_1 = pt2[i2p];
      p2_2 = pt2[i2p2];
      if( line_intersect_2D_3D(p1_0, p1_1, p2_0, p2_1, p2_2, intersect, &u1, &u2, &inbound) ) {
        /* from the value of u1, u2 and inbound, we can partially decide if a point is inside or outside of polygon */

        /* add the intersection into intersetList, The intersection might already be in
           intersectList and will be taken care addIntersect
        */
        if(addIntersect(intersectList, intersect[0], intersect[1], intersect[2], 1, u1, u2, inbound, i1, i1p, i2, i2p)) {
          /* add the intersection into the grid1List */

          if(u1 == 1) {
            insertIntersect(grid1List, intersect[0], intersect[1], intersect[2], 0.0, u2, inbound, p1_1[0], p1_1[1], p1_1[2]);
          }
          else
            insertIntersect(grid1List, intersect[0], intersect[1], intersect[2], u1, u2, inbound, p1_0[0], p1_0[1], p1_0[2]);
          /* when u1 == 0 or 1, need to adjust the vertice to intersect value for roundoff error */
          if(u1==1) {
            p1_1[0] = intersect[0];
            p1_1[1] = intersect[1];
            p1_1[2] = intersect[2];
          }
          else if(u1 == 0) {
            p1_0[0] = intersect[0];
            p1_0[1] = intersect[1];
            p1_0[2] = intersect[2];
          }
          /* add the intersection into the grid2List */
          if(u2==1)
            insertIntersect(grid2List, intersect[0], intersect[1], intersect[2], 0.0, u1, 0, p2_1[0], p2_1[1], p2_1[2]);
          else
            insertIntersect(grid2List, intersect[0], intersect[1], intersect[2], u2, u1, 0, p2_0[0], p2_0[1], p2_0[2]);
          /* when u2 == 0 or 1, need to adjust the vertice to intersect value for roundoff error */
          if(u2==1) {
            p2_1[0] = intersect[0];
            p2_1[1] = intersect[1];
            p2_1[2] = intersect[2];
          }
          else if(u2 == 0) {
            p2_0[0] = intersect[0];
            p2_0[1] = intersect[1];
            p2_0[2] = intersect[2];
          }
        }
      }
    }
  }

  /* set inbound value for the points in intersectList that has inbound == 0,
     this will also set some inbound value of the points in grid1List
  */

  /* get the first point in intersectList has inbound = 2, if not, set inbound value */
  has_inbound = 0;
  /* loop through intersectList to see if there is any has inbound=1 or 2 */
  temp = intersectList;
  nintersect = length(intersectList);
  if(nintersect > 1) {
    getFirstInbound(intersectList, firstIntersect);
    if(firstIntersect->initialized) {
      has_inbound = 1;
    }
  }

  /* when has_inbound == 0, get the grid1List and grid2List */
  if( !has_inbound && nintersect > 1) {
    setInbound(intersectList, grid1List);
    getFirstInbound(intersectList, firstIntersect);
    if(firstIntersect->initialized) has_inbound = 1;
  }

  /* if has_inbound = 1, find the overlapping */
  n_out = 0;

  if(has_inbound) {
    maxiter1 = nintersect;
    temp1 = getNode(grid1List, *firstIntersect);
    if( temp1 == NULL) {
      double lon[10], lat[10];
      int i;
      xyz2latlon(n1_in, x1_in, y1_in, z1_in, lon, lat);
      for(i=0; i< n1_in; i++) printf("lon1 = %g, lat1 = %g\n", lon[i]*R2D, lat[i]*R2D);
      printf("\n");
      xyz2latlon(n2_in, x2_in, y2_in, z2_in, lon, lat);
      for(i=0; i< n2_in; i++) printf("lon2 = %g, lat2 = %g\n", lon[i]*R2D, lat[i]*R2D);
      printf("\n");

      error_handler("firstIntersect is not in the grid1List");
    }
    addNode(polyList, *firstIntersect);
    nintersect--;

    /* Loop over the grid1List and grid2List to find again the firstIntersect */
    curList = grid1List;
    curListNum = 0;

    /* Loop through curList to find the next intersection, the loop will end
       when come back to firstIntersect
    */
    copyNode(curIntersect, *firstIntersect);
    iter1 = 0;
    found1 = 0;

    while( iter1 < maxiter1 ) {
      /* find the curIntersect in curList and get the next intersection points */
      temp1 =  getNode(curList, *curIntersect);
      temp2 = temp1->Next;
      if( temp2 == NULL ) temp2 = curList;

      maxiter2 = length(curList);
      found2 = 0;
      iter2  = 0;
      /* Loop until find the next intersection */
      while( iter2 < maxiter2 ) {
        int temp2IsIntersect;

        temp2IsIntersect = 0;
        if( isIntersect( *temp2 ) ) { /* copy the point and switch to the grid2List */
          struct Node *temp3;

          /* first check if temp2 is the firstIntersect */
          if( sameNode( *temp2, *firstIntersect) ) {
            found1 = 1;
            break;
          }

          temp3 = temp2->Next;
          if( temp3 == NULL) temp3 = curList;
          if( temp3 == NULL) error_handler("creat_xgrid.c: temp3 can not be NULL");
          found2 = 1;
          /* if next node is inside or an intersection,
             need to keep on curList
          */
          temp2IsIntersect = 1;
          if( isIntersect(*temp3) || (temp3->isInside == 1)  ) found2 = 0;
        }
        if(found2) {
          copyNode(curIntersect, *temp2);
          break;
        }
        else {
          addNode(polyList, *temp2);
          if(temp2IsIntersect) {
            nintersect--;
          }
        }
        temp2 = temp2->Next;
        if( temp2 == NULL ) temp2 = curList;
        iter2 ++;
      }
      if(found1) break;

      if( !found2 ) error_handler(" not found the next intersection ");

      /* if find the first intersection, the poly found */
      if( sameNode( *curIntersect, *firstIntersect) ) {
        found1 = 1;
        break;
      }

      /* add curIntersect to polyList and remove it from intersectList and curList */
      addNode(polyList, *curIntersect);
      nintersect--;


      /* switch curList */
      if( curListNum == 0) {
        curList = grid2List;
        curListNum = 1;
      }
      else {
        curList = grid1List;
        curListNum = 0;
      }
      iter1++;
    }
    if(!found1) error_handler("not return back to the first intersection");

    /* currently we are only clipping convex polygon to convex polygon */
    if( nintersect > 0) error_handler("After clipping, nintersect should be 0");

    /* copy the polygon to x_out, y_out, z_out */
    temp1 = polyList;
    while (temp1 != NULL) {
      getCoordinate(*temp1, x_out+n_out, y_out+n_out, z_out+n_out);
      temp1 = temp1->Next;
      n_out++;
    }

    /* if(n_out < 3) error_handler(" The clipped region has < 3 vertices"); */
    if( n_out < 3) n_out = 0;
  }

  /* check if grid1 is inside grid2 */
  if(n_out==0){
    /* first check number of points in grid1 is inside grid2 */
    int n, n1in2;
    /* One possible is that grid1List is inside grid2List */
    n1in2 = 0;
    temp = grid1List;
    while(temp) {
      if(temp->intersect != 1) {
        if( temp->isInside == 1) n1in2++;
      }
      temp = getNextNode(temp);
    }
    if(npts1==n1in2) { /* grid1 is inside grid2 */
      n_out = npts1;
      n = 0;
      temp = grid1List;
      while( temp ) {
        getCoordinate(*temp, &x_out[n], &y_out[n], &z_out[n]);
        n++;
        temp = getNextNode(temp);
      }
    }
    if(n_out>0) return n_out;
  }

  /* check if grid2List is inside grid1List */
  if(n_out ==0){
    int n, n2in1;

    temp = grid2List;
    n2in1 = 0;
    while(temp) {
      if(temp->intersect != 1) {
        if( temp->isInside == 1) n2in1++;
      }
      temp = getNextNode(temp);
    }

    if(npts2==n2in1) { /* grid2 is inside grid1 */
      n_out = npts2;
      n = 0;
      temp = grid2List;
      while( temp ) {
        getCoordinate(*temp, &x_out[n], &y_out[n], &z_out[n]);
        n++;
        temp = getNextNode(temp);
      }

    }
  }


  return n_out;
}


/* Intersects between the line a and the seqment s
   where both line and segment are great circle lines on the sphere represented by
   3D cartesian points.
   [sin sout] are the ends of a line segment
   returns true if the lines could be intersected, false otherwise.
   inbound means the direction of (a1,a2) go inside or outside of (q1,q2,q3)
*/

int line_intersect_2D_3D(double *a1, double *a2, double *q1, double *q2, double *q3,
                         double *intersect, double *u_a, double *u_q, int *inbound){

  /* Do this intersection by reprsenting the line a1 to a2 as a plane through the
     two line points and the origin of the sphere (0,0,0). This is the
     definition of a great circle arc.
  */
  double plane[9];
  double plane_p[2];
  double u;
  double p1[3], v1[3], v2[3];
  double c1[3], c2[3], c3[3];
  double coincident, sense, norm;
  int    i;
  int is_inter1, is_inter2;

  *inbound = 0;

  /* first check if any vertices are the same */
  if(samePoint(a1[0], a1[1], a1[2], q1[0], q1[1], q1[2])) {
    *u_a = 0;
    *u_q = 0;
    intersect[0] = a1[0];
    intersect[1] = a1[1];
    intersect[2] = a1[2];
    return 1;
  }
  else if (samePoint(a1[0], a1[1], a1[2], q2[0], q2[1], q2[2])) {
    *u_a = 0;
    *u_q = 1;
    intersect[0] = a1[0];
    intersect[1] = a1[1];
    intersect[2] = a1[2];
    return 1;
  }
  else if(samePoint(a2[0], a2[1], a2[2], q1[0], q1[1], q1[2])) {
    *u_a = 1;
    *u_q = 0;
    intersect[0] = a2[0];
    intersect[1] = a2[1];
    intersect[2] = a2[2];
    return 1;
  }
  else if (samePoint(a2[0], a2[1], a2[2], q2[0], q2[1], q2[2])) {
    *u_a = 1;
    *u_q = 1;
    intersect[0] = a2[0];
    intersect[1] = a2[1];
    intersect[2] = a2[2];
    return 1;
  }


  /* Load points defining plane into variable (these are supposed to be in counterclockwise order) */
  plane[0]=q1[0];
  plane[1]=q1[1];
  plane[2]=q1[2];
  plane[3]=q2[0];
  plane[4]=q2[1];
  plane[5]=q2[2];
  plane[6]=0.0;
  plane[7]=0.0;
  plane[8]=0.0;

  /* Intersect the segment with the plane */
  is_inter1 = intersect_tri_with_line(plane, a1, a2, plane_p, u_a);

  if(!is_inter1)
    return 0;

  if(fabs(*u_a) < EPSLN8) *u_a = 0;
  if(fabs(*u_a-1) < EPSLN8) *u_a = 1;


  if( (*u_a < 0) || (*u_a > 1) ) return 0;

  /* Load points defining plane into variable (these are supposed to be in counterclockwise order) */
  plane[0]=a1[0];
  plane[1]=a1[1];
  plane[2]=a1[2];
  plane[3]=a2[0];
  plane[4]=a2[1];
  plane[5]=a2[2];
  plane[6]=0.0;
  plane[7]=0.0;
  plane[8]=0.0;

  /* Intersect the segment with the plane */
  is_inter2 = intersect_tri_with_line(plane, q1, q2, plane_p, u_q);

  if(!is_inter2)
    return 0;

  if(fabs(*u_q) < EPSLN8) *u_q = 0;
  if(fabs(*u_q-1) < EPSLN8) *u_q = 1;


  if( (*u_q < 0) || (*u_q > 1) ) return 0;

  u =*u_a;

  /* The two planes are coincidental */
  vect_cross(a1, a2, c1);
  vect_cross(q1, q2, c2);
  vect_cross(c1, c2, c3);
  coincident = metric(c3);

  if(fabs(coincident) < EPSLN30) return 0;

  /* Calculate point of intersection */
  intersect[0]=a1[0] + u*(a2[0]-a1[0]);
  intersect[1]=a1[1] + u*(a2[1]-a1[1]);
  intersect[2]=a1[2] + u*(a2[2]-a1[2]);

  norm = metric( intersect );
  for(i = 0; i < 3; i ++) intersect[i] /= norm;

  /* when u_q =0 or u_q =1, the following could not decide the inbound value */
  if(*u_q != 0 && *u_q != 1){

    p1[0] = a2[0]-a1[0];
    p1[1] = a2[1]-a1[1];
    p1[2] = a2[2]-a1[2];
    v1[0] = q2[0]-q1[0];
    v1[1] = q2[1]-q1[1];
    v1[2] = q2[2]-q1[2];
    v2[0] = q3[0]-q2[0];
    v2[1] = q3[1]-q2[1];
    v2[2] = q3[2]-q2[2];

    vect_cross(v1, v2, c1);
    vect_cross(v1, p1, c2);

    sense = dot(c1, c2);
    *inbound = 1;
    if(sense > 0) *inbound = 2; /* v1 going into v2 in CCW sense */
  }

  return 1;
}


/*------------------------------------------------------------------------------
  double poly_ctrlat(const double x[], const double y[], int n)
  This routine is used to calculate the latitude of the centroid
  ---------------------------------------------------------------------------*/

double poly_ctrlat(const double x[], const double y[], int n)
{
  double ctrlat = 0.0;
  int    i;

  for (i=0;i<n;i++) {
    int ip = (i+1) % n;
    double dx = (x[ip]-x[i]);
    double dy, avg_y, hdy;
    double lat1, lat2;
    lat1 = y[ip];
    lat2 = y[i];
    dy = lat2 - lat1;
    hdy = dy*0.5;
    avg_y = (lat1+lat2)*0.5;
    if      (dx==0.0) continue;
    if(dx > M_PI)  dx = dx - 2.0*M_PI;
    if(dx < -M_PI) dx = dx + 2.0*M_PI;

    if ( fabs(hdy)< SMALL_VALUE ) /* cheap area calculation along latitude */
      ctrlat -= dx*(2*cos(avg_y) + lat2*sin(avg_y) - cos(lat1) );
    else
      ctrlat -= dx*( (sin(hdy)/hdy)*(2*cos(avg_y) + lat2*sin(avg_y)) - cos(lat1) );
  }
  return (ctrlat*RADIUS*RADIUS);
} /* poly_ctrlat */

/*------------------------------------------------------------------------------
  double poly_ctrlon(const double x[], const double y[], int n, double clon)
  This routine is used to calculate the lontitude of the centroid.
  ---------------------------------------------------------------------------*/
double poly_ctrlon(const double x[], const double y[], int n, double clon)
{
  double ctrlon = 0.0;
  int    i;

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

    if(fabs(dphi2 -dphi1) < M_PI) {
      ctrlon -= dphi * (dphi1*f1+dphi2*f2)/2.0;
    }
    else {
      if(dphi1 > 0.0)
        fac = M_PI;
      else
        fac = -M_PI;
      fint = f1 + (f2-f1)*(fac-dphi1)/fabs(dphi);
      ctrlon -= 0.5*dphi1*(dphi1-fac)*f1 - 0.5*dphi2*(dphi2+fac)*f2
          + 0.5*fac*(dphi1+dphi2)*fint;
    }

  }
  return (ctrlon*RADIUS*RADIUS);
}   /* poly_ctrlon */

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
  const double SMALL = 1.e-12;
  double product;

  product = ( x-x0 )*(y1-y0) + (x0-x1)*(y-y0);
  return (product<=SMALL) ? 1:0;

} /* inside_edge */
