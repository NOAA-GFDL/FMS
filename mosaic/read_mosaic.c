#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "read_mosaic.h"
#include "constant.h"
#ifdef use_netCDF
#include <netcdf.h>
#endif
/*********************************************************************
    void netcdf_error( int status )
    status is the returning value of netcdf call. this routine will
    handle the error when status is not NC_NOERR.
********************************************************************/
void handle_netcdf_error(const char *msg, int status )
{
  char errmsg[512];

  sprintf( errmsg, "%s: %s", msg, nc_strerror(status) );
  error_handler(errmsg);

}; /* handle_netcdf_error */

int get_dimlen(const char* file, const char *name)
{
  int ncid, dimid, status, len;
  size_t size;
  char msg[512];
#ifdef use_netCDF  
  status = nc_open(file, NC_NOWRITE, &ncid);
  if(status != NC_NOERR) {
    sprintf(msg, "in opening file ", file);
    handle_netcdf_error(msg, status);
  }
  
  status = nc_inq_dimid(ncid, name, &dimid);
  if(status != NC_NOERR) {
    sprintf(msg, "in getting dimid of %s from file %s.", name, file);
    handle_netcdf_error(msg, status);
  }
  
  status = nc_inq_dimlen(ncid, dimid, &size);
  if(status != NC_NOERR) {
    sprintf(msg, "in getting dimension size of %s from file %s.", name, file);
    handle_netcdf_error(msg, status);
  }
  status = nc_close(ncid);
  
  
  len = size;
  if(status != NC_NOERR) {
    sprintf(msg, "in closing file ", file);
    handle_netcdf_error(msg, status);
  }
#else
  error_handler("read_mosaic: Add flag -Duse_netCDF when compiling");
#endif
  
  return len;
  
}; /* get_dimlen */

/*******************************************************************************
   void get_string_data(const char *file, const char *name, char *data)
   get string data of field with "name" from "file".
******************************************************************************/
void get_string_data(const char *file, const char *name, char *data)
{
  int ncid, varid, status;
  char msg[512];

#ifdef use_netCDF    
  status = nc_open(file, NC_NOWRITE, &ncid);
  if(status != NC_NOERR) {
    sprintf(msg, "in opening file %s", file);
    handle_netcdf_error(msg, status);
  }
  status = nc_inq_varid(ncid, name, &varid);
  if(status != NC_NOERR) {
    sprintf(msg, "in getting varid of %s from file %s.", name, file);
    handle_netcdf_error(msg, status);
  }     
  status = nc_get_var_text(ncid, varid, data);
  if(status != NC_NOERR) {
    sprintf(msg, "in getting data of %s from file %s.", name, file);
    handle_netcdf_error(msg, status);
  }    
#else
  error_handler("read_mosaic: Add flag -Duse_netCDF when compiling");
#endif
  
}; /* get_string_data */

/*******************************************************************************
   void get_string_data_level(const char *file, const char *name, const size_t *start, const size_t *nread, char *data)
   get string data of field with "name" from "file".
******************************************************************************/
void get_string_data_level(const char *file, const char *name, char *data, const int *level)
{
  int ncid, varid, status, i;
  size_t start[4], nread[4];
  char msg[512];

#ifdef use_netCDF  
  status = nc_open(file, NC_NOWRITE, &ncid);
  if(status != NC_NOERR) {
    sprintf(msg, "in opening file %s", file);
    handle_netcdf_error(msg, status);
  }
  status = nc_inq_varid(ncid, name, &varid);
  if(status != NC_NOERR) {
    sprintf(msg, "in getting varid of %s from file %s.", name, file);
    handle_netcdf_error(msg, status);
  }
  for(i=0; i<4; i++) {
    start[i] = 0; nread[i] = 1;
  }
  start[0] = *level; nread[1] = STRING;
  status = nc_get_vara_text(ncid, varid, start, nread, data);
  if(status != NC_NOERR) {
    sprintf(msg, "in getting data of %s from file %s.", name, file);
    handle_netcdf_error(msg, status);
  }    
#else
  error_handler("read_mosaic: Add flag -Duse_netCDF when compiling");
#endif
  
}; /* get_string_data_level */


/*******************************************************************************
   void get_int_data(const char *file, const char *name, int *data)
   get int data of field with "name" from "file".
******************************************************************************/
void get_int_data(const char *file, const char *name, int *data)
{
  int ncid, varid, status;
  char msg[512];

#ifdef use_netCDF    
  status = nc_open(file, NC_NOWRITE, &ncid);
  if(status != NC_NOERR) {
    sprintf(msg, "in opening file %s", file);
    handle_netcdf_error(msg, status);
  }
  status = nc_inq_varid(ncid, name, &varid);
  if(status != NC_NOERR) {
    sprintf(msg, "in getting varid of %s from file %s.", name, file);
    handle_netcdf_error(msg, status);
  }     
  status = nc_get_var_int(ncid, varid, data);
  if(status != NC_NOERR) {
    sprintf(msg, "in getting data of %s from file %s", name, file);
    handle_netcdf_error(msg, status);
  }    
#else
  error_handler("read_mosaic: Add flag -Duse_netCDF when compiling");
#endif
  
}; /* get_int_data */

/*******************************************************************************
   void get_double_data(const char *file, const char *name, double *data)
   get double data of field with "name" from "file".
******************************************************************************/
void get_double_data(const char *file, const char *name, double *data)
{

  int ncid, varid, status;  
  char msg[512];

#ifdef use_netCDF    
  status = nc_open(file, NC_NOWRITE, &ncid);
  if(status != NC_NOERR) {
    sprintf(msg, "in opening file ", file);
    handle_netcdf_error(msg, status);
  }
  status = nc_inq_varid(ncid, name, &varid);
  if(status != NC_NOERR) {
    sprintf(msg, "in getting varid of %s from file %s.", name, file);
    handle_netcdf_error(msg, status);
  }     
  status = nc_get_var_double(ncid, varid, data);
  if(status != NC_NOERR) {
    sprintf(msg, "in getting data of %s from file %s.", name, file);
    handle_netcdf_error(msg, status);
  }    
#else
  error_handler("read_mosaic: Add flag -Duse_netCDF when compiling");
#endif
  
}; /* get_double_data */

/******************************************************************************
   void get_var_text_att(const char *file, const char *name, const char *attname, char *att)
   get text attribute of field 'name' from 'file
******************************************************************************/
void get_var_text_att(const char *file, const char *name, const char *attname, char *att)
{
  int ncid, varid, status;  
  char msg[512];

#ifdef use_netCDF    
  status = nc_open(file, NC_NOWRITE, &ncid);
  if(status != NC_NOERR) {
    sprintf(msg, "in opening file ", file);
    handle_netcdf_error(msg, status);
  }
  status = nc_inq_varid(ncid, name, &varid);
  if(status != NC_NOERR) {
    sprintf(msg, "in getting varid of %s from file %s.", name, file);
    handle_netcdf_error(msg, status);
  }     
  status = nc_get_att_text(ncid, varid, attname, att);
  if(status != NC_NOERR) {
    sprintf(msg, "in getting attribute %s of %s from file %s.", attname, name, file);
    handle_netcdf_error(msg, status);
  }
#else
  error_handler("read_mosaic: Add flag -Duse_netCDF when compiling");
#endif
  
}; /* get_var_text_att */

/****************************************************************************/
#ifndef __AIX
int read_mosaic_xgrid_order1_(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2, double *area )
{
  int nxgrid;
  nxgrid = read_mosaic_xgrid_order1(xgrid_file, i1, j1, i2, j2, area);
  return nxgrid;
  
};
#endif


int read_mosaic_xgrid_order1(const char *xgrid_file, int *i1, int *j1, int *i2, int *j2, double *area )
{
  int    ncells, n;
  int    *tile1_cell, *tile2_cell;
  double garea;
  
  ncells = get_dimlen(xgrid_file, "ncells");
  if(ncells ==0) return ncells;

  tile1_cell       = (int *)malloc(ncells*2*sizeof(int));
  tile2_cell       = (int *)malloc(ncells*2*sizeof(int));
  get_int_data(xgrid_file, "tile1_cell", tile1_cell);
  get_int_data(xgrid_file, "tile2_cell", tile2_cell);
  get_double_data(xgrid_file, "xgrid_area", area);
  garea = 4*M_PI*RADIUS*RADIUS;
  
  for(n=0; n<ncells; n++) {
    i1[n] = tile1_cell[n*2] - 1;
    j1[n] = tile1_cell[n*2+1] - 1;
    i2[n] = tile2_cell[n*2] - 1;
    j2[n] = tile2_cell[n*2+1] - 1;
    area[n] /= garea; /* rescale the exchange grid area to unit earth area */
  }

  free(tile1_cell);
  free(tile2_cell);
  
  return ncells;
  
}; /* get_mosaic_xgrid */


/******************************************************************************
  int read_mosaic_ntiles(const char *mosaic_file)
  return number tiles in mosaic_file
******************************************************************************/
#ifndef __AIX
int read_mosaic_ntiles_(const char *mosaic_file)
{
  return read_mosaic_ntiles(mosaic_file);
}
#endif
int read_mosaic_ntiles(const char *mosaic_file)
{

  int ntiles;

  ntiles = get_dimlen(mosaic_file, "ntiles");

  return ntiles;
  
}; /* read_mosaic_ntiles */

/******************************************************************************
  int read_mosaic_ncontacts(const char *mosaic_file)
  return number of contacts in mosaic_file
******************************************************************************/
#ifndef __AIX
int read_mosaic_ncontacts_(const char *mosaic_file)
{
  return read_mosaic_ncontacts(mosaic_file);
}
#endif
int read_mosaic_ncontacts(const char *mosaic_file)
{

  int ncontacts;

  ncontacts = get_dimlen(mosaic_file, "ncontact");

  return ncontacts;
  
}; /* read_mosaic_ncontacts */


/*****************************************************************************
  void read_mosaic_grid_sizes(const char *mosaic_file, int *nx, int *ny)
  read mosaic grid size of each tile, currently we are assuming the refinement is 2.
*****************************************************************************/
#ifndef __AIX
void read_mosaic_grid_sizes_(const char *mosaic_file, const char *tile_dir, int *nx, int *ny)
{
  read_mosaic_grid_sizes(mosaic_file, tile_dir, nx, ny);
}
#endif
void read_mosaic_grid_sizes(const char *mosaic_file, const char *tile_dir, int *nx, int *ny)
{
  int ntiles, n;
  char gridfile[STRING], tilefile[2*STRING];
  char dir[STRING];
  const int x_refine = 2, y_refine = 2;

  if( tile_dir )
    strcpy(dir, tile_dir);
  else
    get_string_data(mosaic_file, "gridlocation", dir);
  
  ntiles = get_dimlen(mosaic_file, "ntiles");
  for(n = 0; n < ntiles; n++) {
    get_string_data_level(mosaic_file, "gridfiles", gridfile, &n);
    sprintf(tilefile, "%s/%s", dir, gridfile);
    nx[n] = get_dimlen(tilefile, "nx");
    ny[n] = get_dimlen(tilefile, "ny");
    if(nx[n]%x_refine != 0) error_handler("Error from read_mosaic_grid_sizes: nx is not divided by x_refine");
    if(ny[n]%y_refine != 0) error_handler("Error from read_mosaic_grid_sizes: ny is not divided by y_refine");
    nx[n] /= x_refine;
    ny[n] /= y_refine;
  }
  
}; /* read_mosaic_grid_sizes */
  

/******************************************************************************
  void read_mosaic_contact(const char *mosaic_file)
  read mosaic contact information
******************************************************************************/
#ifndef __AIX
void read_mosaic_contact_(const char *mosaic_file, int *tile1, int *tile2, int *istart1, int *iend1,
			 int *jstart1, int *jend1, int *istart2, int *iend2, int *jstart2, int *jend2)
{
  read_mosaic_contact(mosaic_file, tile1, tile2, istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2);
}
#endif

void read_mosaic_contact(const char *mosaic_file, int *tile1, int *tile2, int *istart1, int *iend1,
			 int *jstart1, int *jend1, int *istart2, int *iend2, int *jstart2, int *jend2)
{
  char contacts[STRING];
  char **gridtiles;
  const unsigned int MAXVAR = 40;
  char pstring[MAXVAR][STRING];
  int ntiles, ncontacts, n, m, l, nstr, found;
  const int x_refine = 2, y_refine = 2;
  
  ntiles = get_dimlen(mosaic_file, "ntiles");
  gridtiles = (char **)malloc(ntiles*sizeof(char *));
  for(n=0; n<ntiles; n++) {
    gridtiles[n] = (char *)malloc(STRING*sizeof(char));
    get_string_data_level(mosaic_file, "gridtiles", gridtiles[n], &n);
  }
    
  ncontacts = get_dimlen(mosaic_file, "ncontact"); 
  for(n = 0; n < ncontacts; n++) {
    get_string_data_level(mosaic_file, "contacts", contacts, &n);
    /* parse the string contacts to get tile number */
    tokenize( contacts, ":", STRING, MAXVAR, pstring, &nstr);
    if(nstr != 4) error_handler("Error from read_mosaic: number of elements "
				 "in contact seperated by :/:: should be 4");
    found = 0;
    for(m=0; m<ntiles; m++) {
      if(strcmp(gridtiles[m], pstring[1]) == 0) { /*found the tile name */
	found = 1;
	tile1[n] = m+1;
	break;
      }
    }
    if(!found) error_handler("error from read_mosaic: the first tile name specified "
			     "in contact is not found in tile list");
    found = 0;
    for(m=0; m<ntiles; m++) {
      if(strcmp(gridtiles[m], pstring[3]) == 0) { /*found the tile name */
	found = 1;
	tile2[n] = m+1;
	break;
      }
    }
    if(!found) error_handler("error from read_mosaic: the second tile name specified "
			     "in contact is not found in tile list");    
    get_string_data_level(mosaic_file, "contact_index", contacts, &n);
    /* parse the string to get contact index */
    tokenize( contacts, ":,", STRING, MAXVAR, pstring, &nstr);
    if(nstr != 8) error_handler("Error from read_mosaic: number of elements "
				 "in contact_index seperated by :/, should be 8");
    /* make sure the string is only composed of numbers */
    for(m=0; m<nstr; m++) for(l=0; l<strlen(pstring[m]); l++) {
      if(pstring[m][l] > '9' ||  pstring[m][l] < '0' ) {
	error_handler("Error from read_mosaic: some of the character in "
		      "contact_indices except token is not digit number");
      }
    }
    istart1[n] = atoi(pstring[0]);
    iend1[n]   = atoi(pstring[1]);
    jstart1[n] = atoi(pstring[2]);
    jend1[n]   = atoi(pstring[3]);
    istart2[n] = atoi(pstring[4]);
    iend2[n]   = atoi(pstring[5]);
    jstart2[n] = atoi(pstring[6]);
    jend2[n]   = atoi(pstring[7]);
    if(istart1[n] == iend1[n] ) {
      istart1[n] = (istart1[n]+1)/x_refine-1;
      iend1[n]   = istart1[n];
      if( jend1[n] > jstart1[n] ) {
	--(jstart1[n]);
	if(jstart1[n]%y_refine != 0)
	  error_handler("Error from read_mosaic_contact: jstart1 should be an odd number when istart1=iend1");
	if((jend1[n]-jstart1[n])%y_refine != 0)
	  error_handler("Error from read_mosaic_contact: ny1 can not be divided by y_refine when istart1=iend1");
	jstart1[n] /= y_refine;
	jend1[n] = jstart1[n] + (jend1[n]-jstart1[n])/y_refine - 1;
      }
      else if( jstart1[n] > jend1[n] ) {
	--(jend1[n]);
	if(jend1[n]%y_refine != 0)
	  error_handler("Error from read_mosaic_contact: jend1 should be an odd number when istart1=iend1");
	if((jstart1[n]-jend1[n])%y_refine != 0)
	  error_handler("Error from read_mosaic_contact: ny1 can not be divided by y_refine when istart1=iend1");
	jend1[n] /= y_refine;
	jstart1[n] = jend1[n] + (jstart1[n]-jend1[n])/y_refine - 1;
      }
      else {
	error_handler("Error from read_mosaic_contact: jstart1 and jend1 should not be equal when istart1=iend1");
      }
    }
    else if( jstart1[n] == jend1[n] ) {
      jstart1[n] = (jstart1[n]+1)/y_refine-1;
      jend1[n]   = jstart1[n];
      if(iend1[n] > istart1[n] ){
	--(istart1[n]);
	if(istart1[n]%x_refine != 0)
	  error_handler("Error from read_mosaic_contact: istart1 should be an odd number when jstart1=jend1");
	if((iend1[n]-istart1[n])%x_refine != 0)
	  error_handler("Error from read_mosaic_contact: nx1 can not be divided by y_refine when jstart1=jend1");
	istart1[n] /= x_refine;
	iend1[n] = istart1[n] + (iend1[n]-istart1[n])/x_refine - 1;
      }
      else if(istart1[n] > iend1[n] ){
	--(iend1[n]);
	if(iend1[n]%x_refine != 0)
	  error_handler("Error from read_mosaic_contact: iend1 should be an odd number when jstart1=jend1");
	if((istart1[n]-iend1[n])%x_refine != 0)
	  error_handler("Error from read_mosaic_contact: nx1 can not be divided by y_refine when jstart1=jend1");
	iend1[n] /= x_refine;
	istart1[n] = iend1[n] + (istart1[n]-iend1[n])/x_refine - 1;
      }
      else {
	error_handler("Error from read_mosaic_contact: istart1 and iend1 should not be equal when jstart1=jend1");
      }
    }
    else {
      error_handler("Error from read_mosaic_contact: only line contact is supported now, contact developer");
    }
    if(istart2[n] == iend2[n] ) {
      istart2[n] = (istart2[n]+1)/x_refine-1;
      iend2[n]   = istart2[n];
      if( jend2[n] > jstart2[n] ) {
	--(jstart2[n]);
	if(jstart2[n]%y_refine != 0)
	  error_handler("Error from read_mosaic_contact: jstart2 should be an odd number when istart2=iend2");
	if((jend2[n]-jstart2[n])%y_refine != 0)
	  error_handler("Error from read_mosaic_contact: ny2 can not be divided by y_refine when istart2=iend2");
	jstart2[n] /= y_refine;
	jend2[n] = jstart2[n] + (jend2[n]-jstart2[n])/y_refine - 1;
      }
      else if( jstart2[n] > jend2[n] ) {
	--(jend2[n]);
	if(jend2[n]%y_refine != 0)
	  error_handler("Error from read_mosaic_contact: jend2 should be an odd number when istart2=iend2");
	if((jstart2[n]-jend2[n])%y_refine != 0)
	  error_handler("Error from read_mosaic_contact: ny2 can not be divided by y_refine when istart2=iend2");
	jend2[n] /= y_refine;
	jstart2[n] = jend2[n] + (jstart2[n]-jend2[n])/y_refine - 1;
      }
      else {
	error_handler("Error from read_mosaic_contact: jstart2 and jend2 should not be equal when istart2=iend2");
      }
    }
    else if( jstart2[n] == jend2[n] ) {
      jstart2[n] = (jstart2[n]+1)/y_refine-1;
      jend2[n]   = jstart2[n];
      if(iend2[n] > istart2[n] ){
	--(istart2[n]);
	if(istart2[n]%x_refine != 0)
	  error_handler("Error from read_mosaic_contact: istart2 should be an odd number when jstart2=jend2");
	if((iend2[n]-istart2[n])%x_refine != 0)
	  error_handler("Error from read_mosaic_contact: nx2 can not be divided by y_refine when jstart2=jend2");
	istart2[n] /= x_refine;
	iend2[n] = istart2[n] + (iend2[n]-istart2[n])/x_refine - 1;
      }
      else if(istart2[n] > iend2[n] ){
	--(iend2[n]);
	if(iend2[n]%x_refine != 0)
	  error_handler("Error from read_mosaic_contact: iend2 should be an odd number when jstart2=jend2");
	if((istart2[n]-iend2[n])%x_refine != 0)
	  error_handler("Error from read_mosaic_contact: nx2 can not be divided by y_refine when jstart2=jend2");
	iend2[n] /= x_refine;
	istart2[n] = iend2[n] + (istart2[n]-iend2[n])/x_refine - 1;
      }
      else {
	error_handler("Error from read_mosaic_contact: istart2 and iend2 should not be equal when jstart2=jend2");
      }
    }
    else {
      error_handler("Error from read_mosaic_contact: only line contact is supported now, contact developer");
    }

    
  }

  

}; /* read_mosaic_contact */


