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

/** \file
 *  \ingroup tree_utils
 *  \brief utilities for create_xgrid_great_circle
 */

struct Node *nodeList=NULL;
int curListPos=0;

void rewindList(void)
{
  int n;

  curListPos = 0;
  if(!nodeList) nodeList = (struct Node *)malloc(MAXNODELIST*sizeof(struct Node));
  for(n=0; n<MAXNODELIST; n++) initNode(nodeList+n);

}

struct Node *getNext()
{
  struct Node *temp=NULL;
  int n;

  if(!nodeList) {
    nodeList = (struct Node *)malloc(MAXNODELIST*sizeof(struct Node));
    for(n=0; n<MAXNODELIST; n++) initNode(nodeList+n);
  }

  temp = nodeList+curListPos;
  curListPos++;
  if(curListPos > MAXNODELIST) error_handler("getNext: curListPos >= MAXNODELIST");

  return (temp);
}


void initNode(struct Node *node)
{
  node->x = 0;
  node->y = 0;
  node->z = 0;
  node->u = 0;
  node->intersect = 0;
  node->inbound = 0;
  node->isInside = 0;
  node->Next = NULL;
  node->initialized=0;

}

void addEnd(struct Node *list, double x, double y, double z, int intersect, double u, int inbound, int inside)
{

  struct Node *temp=NULL;

  if(list == NULL) error_handler("addEnd: list is NULL");

  if(list->initialized) {

    /* (x,y,z) might already in the list when intersect is true and u=0 or 1 */
    temp = list;
    while (temp) {
      if(samePoint(temp->x, temp->y, temp->z, x, y, z)) return;
      temp=temp->Next;
    }
    temp = list;
    while(temp->Next)
      temp=temp->Next;

    /* Append at the end of the list.  */
    temp->Next = getNext();
    temp = temp->Next;
  }
  else {
    temp = list;
  }

  temp->x = x;
  temp->y = y;
  temp->z = z;
  temp->u = u;
  temp->intersect = intersect;
  temp->inbound = inbound;
  temp->initialized=1;
  temp->isInside = inside;
}

/* return 1 if the point (x,y,z) is added in the list, return 0 if it is already in the list */

int addIntersect(struct Node *list, double x, double y, double z, int intersect, double u1, double u2, int inbound,
                 int is1, int ie1, int is2, int ie2)
{

  double u1_cur, u2_cur;
  int    i1_cur, i2_cur;
  struct Node *temp=NULL;

  if(list == NULL) error_handler("addEnd: list is NULL");

  /* first check to make sure this point is not in the list */
  u1_cur = u1;
  i1_cur = is1;
  u2_cur = u2;
  i2_cur = is2;
  if(u1_cur == 1) {
    u1_cur = 0;
    i1_cur = ie1;
  }
  if(u2_cur == 1) {
    u2_cur = 0;
    i2_cur = ie2;
  }

  if(list->initialized) {
    temp = list;
    while(temp) {
      if( temp->u == u1_cur && temp->subj_index == i1_cur) return 0;
      if( temp->u_clip == u2_cur && temp->clip_index == i2_cur) return 0;
      if( !temp->Next ) break;
      temp=temp->Next;
    }

    /* Append at the end of the list.  */
    temp->Next = getNext();
    temp = temp->Next;
  }
  else {
    temp = list;
  }

  temp->x = x;
  temp->y = y;
  temp->z = z;
  temp->intersect = intersect;
  temp->inbound = inbound;
  temp->initialized=1;
  temp->isInside = 0;
  temp->u = u1_cur;
  temp->subj_index = i1_cur;
  temp->u_clip = u2_cur;
  temp->clip_index = i2_cur;

  return 1;
}


int length(struct Node *list)
{
  struct Node *cur_ptr=NULL;
  int count=0;

  cur_ptr=list;

  while(cur_ptr)
  {
    if(cur_ptr->initialized ==0) break;
    cur_ptr=cur_ptr->Next;
    count++;
  }
  return(count);
}

/* two points are the same if there are close enough */
int samePoint(double x1, double y1, double z1, double x2, double y2, double z2)
{
  if( fabs(x1-x2) > EPSLN10 || fabs(y1-y2) > EPSLN10 || fabs(z1-z2) > EPSLN10 )
    return 0;
  else
    return 1;
}


int sameNode(struct Node node1, struct Node node2)
{
  if( node1.x == node2.x && node1.y == node2.y && node1.z==node2.z )
    return 1;
  else
    return 0;
}


void addNode(struct Node *list, struct Node inNode)
{

  addEnd(list, inNode.x, inNode.y, inNode.z, inNode.intersect, inNode.u, inNode.inbound, inNode.isInside);

}

struct Node *getNode(struct Node *list, struct Node inNode)
{
  struct Node *thisNode=NULL;
  struct Node *temp=NULL;

  temp = list;
  while( temp ) {
    if( sameNode( *temp, inNode ) ) {
      thisNode = temp;
      temp = NULL;
      break;
    }
    temp = temp->Next;
  }

  return thisNode;
}

struct Node *getNextNode(struct Node *list)
{
  return list->Next;
}

void copyNode(struct Node *node_out, struct Node node_in)
{

  node_out->x = node_in.x;
  node_out->y = node_in.y;
  node_out->z = node_in.z;
  node_out->u = node_in.u;
  node_out->intersect = node_in.intersect;
  node_out->inbound   = node_in.inbound;
  node_out->Next = NULL;
  node_out->initialized = node_in.initialized;
  node_out->isInside = node_in.isInside;
}

void printNode(struct Node *list, char *str)
{
  struct Node *temp;

  if(list == NULL) error_handler("printNode: list is NULL");
  if(str) printf("  %s \n", str);
  temp = list;
  while(temp) {
    if(temp->initialized ==0) break;
    printf(" (x, y, z, interset, inbound, isInside) = (%19.15f,%19.15f,%19.15f,%d,%d,%d)\n",
           temp->x, temp->y, temp->z, temp->intersect, temp->inbound, temp->isInside);
    temp = temp->Next;
  }
  printf("\n");
}

int intersectInList(struct Node *list, double x, double y, double z)
{
  struct Node *temp;
  int found=0;

  temp = list;
  found = 0;
  while ( temp ) {
    if( temp->x == x && temp->y == y && temp->z == z ) {
      found = 1;
      break;
    }
    temp=temp->Next;
  }
  if (!found) error_handler("intersectInList: point (x,y,z) is not found in the list");
  if( temp->intersect == 2 )
    return 1;
  else
    return 0;

}


/* The following insert a intersection after non-intersect point (x2,y2,z2), if the point
   after (x2,y2,z2) is an intersection, if u is greater than the u value of the intersection,
   insert after, otherwise insert before
*/
void insertIntersect(struct Node *list, double x, double y, double z, double u1, double u2, int inbound,
                     double x2, double y2, double z2)
{
  struct Node *temp1=NULL, *temp2=NULL;
  struct Node *temp;
  double u_cur;
  int found=0;

  temp1 = list;
  found = 0;
  while ( temp1 ) {
    if( temp1->x == x2 && temp1->y == y2 && temp1->z == z2 ) {
      found = 1;
      break;
    }
    temp1=temp1->Next;
  }
  if (!found) error_handler("inserAfter: point (x,y,z) is not found in the list");

  /* when u = 0 or u = 1, set the grid point to be the intersection point to solve truncation error isuse */
  u_cur = u1;
  if(u1 == 1) {
    u_cur = 0;
    temp1 = temp1->Next;
    if(!temp1) temp1 = list;
  }
  if(u_cur==0) {
    temp1->intersect = 2;
    temp1->isInside = 1;
    temp1->u = u_cur;
    temp1->x = x;
    temp1->y = y;
    temp1->z = z;
    return;
  }

  /* when u2 != 0 and u2 !=1, can decide if one end of the point is outside depending on inbound value */
  if(u2 != 0 && u2 != 1) {
    if(inbound == 1) { /* goes outside, then temp1->Next is an outside point */
      /* find the next non-intersect point */
      temp2 = temp1->Next;
      if(!temp2) temp2 = list;
      while(temp2->intersect) {
        temp2=temp2->Next;
        if(!temp2) temp2 = list;
      }

      temp2->isInside = 0;
    }
    else if(inbound ==2) { /* goes inside, then temp1 is an outside point */
      temp1->isInside = 0;
    }
  }

  temp2 = temp1->Next;
  while ( temp2 ) {
    if( temp2->intersect == 1 ) {
      if( temp2->u > u_cur ) {
        break;
      }
    }
    else
      break;
    temp1 = temp2;
    temp2 = temp2->Next;
  }

  /* assign value */
  temp = getNext();
  temp->x = x;
  temp->y = y;
  temp->z = z;
  temp->u = u_cur;
  temp->intersect = 1;
  temp->inbound = inbound;
  temp->isInside = 1;
  temp->initialized = 1;
  temp1->Next = temp;
  temp->Next = temp2;

}

double gridArea(struct Node *grid) {
  double x[20], y[20], z[20];
  struct Node *temp=NULL;
  double area;
  int n;

  temp = grid;
  n = 0;
  while( temp ) {
    x[n] = temp->x;
    y[n] = temp->y;
    z[n] = temp->z;
    n++;
    temp = temp->Next;
  }

  area = great_circle_area(n, x, y, z);

  return area;

}

int isIntersect(struct Node node) {

  return node.intersect;

}


int getInbound( struct Node node )
{
  return node.inbound;
}

struct Node *getLast(struct Node *list)
{
  struct Node *temp1;

  temp1 = list;
  if( temp1 ) {
    while( temp1->Next ) {
      temp1 = temp1->Next;
    }
  }

  return temp1;
}


int getFirstInbound( struct Node *list, struct Node *nodeOut)
{
  struct Node *temp=NULL;

  temp=list;

  while(temp) {
    if( temp->inbound == 2 ) {
      copyNode(nodeOut, *temp);
      return 1;
    }
    temp=temp->Next;
  }

  return 0;
}

void getCoordinate(struct Node node, double *x, double *y, double *z)
{


  *x = node.x;
  *y = node.y;
  *z = node.z;

}

void getCoordinates(struct Node *node, double *p)
{


  p[0] = node->x;
  p[1] = node->y;
  p[2] = node->z;

}

void setCoordinate(struct Node *node, double x, double y, double z)
{


  node->x = x;
  node->y = y;
  node->z = z;

}

/* set inbound value for the points in interList that has inbound =0,
   this will also set some inbound value of the points in list1
*/

void setInbound(struct Node *interList, struct Node *list)
{

  struct Node *temp1=NULL, *temp=NULL;
  struct Node *temp1_prev=NULL, *temp1_next=NULL;
  int prev_is_inside, next_is_inside;

  /* for each point in interList, search through list to decide the inbound value the interList point */
  /* For each inbound point, the prev node should be outside and the next is inside. */
  if(length(interList) == 0) return;

  temp = interList;

  while(temp) {
    if( !temp->inbound) {
      /* search in grid1 to find the prev and next point of temp, when prev point is outside and next point is inside
         inbound = 2, else inbound = 1*/
      temp1 = list;
      temp1_prev = NULL;
      temp1_next = NULL;
      while(temp1) {
        if(sameNode(*temp1, *temp)) {
          if(!temp1_prev) temp1_prev = getLast(list);
          temp1_next = temp1->Next;
          if(!temp1_next) temp1_next = list;
          break;
        }
        temp1_prev = temp1;
        temp1 = temp1->Next;
      }
      if(!temp1_next) error_handler("Error from create_xgrid.c: temp is not in list1");
      if( temp1_prev->isInside == 0 && temp1_next->isInside == 1)
        temp->inbound = 2;   /* go inside */
      else
        temp->inbound = 1;
    }
    temp=temp->Next;
  }
}

int isInside(struct Node *node) {

  if(node->isInside == -1) error_handler("Error from mosaic_util.c: node->isInside is not set");
  return(node->isInside);

}

/*  #define debug_test_create_xgrid */

/* check if node is inside polygon list or not */
int insidePolygon( struct Node *node, struct Node *list)
{
  int is_inside;
  double pnt0[3], pnt1[3], pnt2[3];
  double anglesum;
  struct Node *p1=NULL, *p2=NULL;

  anglesum = 0;

  pnt0[0] = node->x;
  pnt0[1] = node->y;
  pnt0[2] = node->z;

  p1 = list;
  p2 = list->Next;
  is_inside = 0;


  while(p1) {
    pnt1[0] = p1->x;
    pnt1[1] = p1->y;
    pnt1[2] = p1->z;
    pnt2[0] = p2->x;
    pnt2[1] = p2->y;
    pnt2[2] = p2->z;
    if( samePoint(pnt0[0], pnt0[1], pnt0[2], pnt1[0], pnt1[1], pnt1[2]) ){
      return 1;
    }
    anglesum += spherical_angle(pnt0, pnt2, pnt1);
    p1 = p1->Next;
    p2 = p2->Next;
    if(p2==NULL){
      p2 = list;
    }
  }

  if( fabs(anglesum - 2*M_PI) < EPSLN8 ){
    is_inside = 1;
  }
  else{
    is_inside = 0;
  }

  return is_inside;

}
