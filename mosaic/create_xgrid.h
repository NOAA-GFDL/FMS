#ifndef CREATE_XGRID_H_
#define CREATE_XGRID_H_
#define MAXXGRID 1e7
#define MV 50
/* this value is small compare to earth area */
#define AREA_THRESH  10
#define UNIT_AREA_THRESH (1.e-12)

double poly_ctrlon(const double lon[], const double lat[], int n, double clon);
double poly_ctrlat(const double lon[], const double lat[], int n);
double box_ctrlon(double ll_lon, double ll_lat, double ur_lon, double ur_lat, double clon);
double box_ctrlat(double ll_lon, double ll_lat, double ur_lon, double ur_lat);
int get_maxxgrid(void);
void get_grid_area(const int *nlon, const int *nlat, const double *lon, const double *lat, double *area);
int clip(const double lon_in[], const double lat_in[], int n_in, double ll_lon, double ll_lat,
	 double ur_lon, double ur_lat, double lon_out[], double lat_out[]);
int clip_2dx2d(const double lon1_in[], const double lat1_in[], int n1_in, 
	       const double lon2_in[], const double lat2_in[], int n2_in, 
	       double lon_out[], double lat_out[]);
int create_xgrid_1dx2d_order1(const int *nlon1, const int *nlat1, const int *nlon2, const int *nlat2, 
                              const double *lon1, const double *lat1, const double *lon2, const double *lat2, 
                              const double *mask1, int *grid1_i, int *grid1_j, int *grid2_i, int *grid2_j, 
                              double *xgrid_area);
int create_xgrid_1dx2d_order2(const int *nlon1, const int *nlat1, const int *nlon2, const int *nlat2, 
                              const double *lon1, const double *lat1, const double *lon2, const double *lat2, 
                              const double *mask1, int *grid1_i, int *grid1_j, int *grid2_i, int *grid2_j, 
                              double *xgrid_area, double *grid1_di, double *grid1_dj,
			      double *grid2_di, double *grid2_dj);
int create_xgrid_2dx1d_order1(const int *nlon1, const int *nlat1, const int *nlon2, const int *nlat2, 
                              const double *lon1, const double *lat1, const double *lon2, const double *lat2, 
                              const double *mask1, int *grid1_i, int *grid1_j, int *grid2_i, int *grid2_j, 
                              double *xgrid_area);
int create_xgrid_2dx1d_order2(const int *nlon1, const int *nlat1, const int *nlon2, const int *nlat2, 
                              const double *lon1, const double *lat1, const double *lon2, const double *lat2, 
                              const double *mask1, int *grid1_i, int *grid1_j, int *grid2_i, int *grid2_j, 
                              double *xgrid_area, double *grid1_di, double *grid1_dj,
			      double *grid2_di, double *grid2_dj);
int create_xgrid_2dx2d_order1(const int *nlon1, const int *nlat1, const int *nlon2, const int *nlat2, 
                              const double *lon1, const double *lat1, const double *lon2, const double *lat2, 
                              const double *mask1, int *grid1_i, int *grid1_j, int *grid2_i, int *grid2_j, 
                              double *xgrid_area);
int create_xgrid_2dx2d_order2(const int *nlon1, const int *nlat1, const int *nlon2, const int *nlat2, 
                              const double *lon1, const double *lat1, const double *lon2, const double *lat2, 
                              const double *mask1, int *grid1_i, int *grid1_j, int *grid2_i, int *grid2_j, 
                              double *xgrid_area, double *grid1_di, double *grid1_dj,
			      double *grid2_di, double *grid2_dj);
#endif
