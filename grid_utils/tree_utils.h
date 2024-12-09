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
#ifndef TREE_UTILS_H_
#define TREE_UTILS_H_

#ifndef MAXNODELIST
#define MAXNODELIST 100
#endif

struct Node{
  double x, y, z, u, u_clip;
  int intersect; /* indicate if this point is an intersection, 0 = no, 1= yes, 2=both intersect and vertices */
  int inbound;      /* -1 uninitialized, 0 coincident, 1 outbound, 2 inbound */
  int initialized; /* = 0 means empty list */
  int isInside;   /* = 1 means one point is inside the other polygon, 0 is not, -1 undecided. */
  int subj_index; /* the index of subject point that an intersection follow. */
  int clip_index; /* the index of clip point that an intersection follow */
  struct Node *Next;
};

void rewindList(void);

struct Node *getNext();

void initNode(struct Node *node);

void addEnd(struct Node *list, double x, double y, double z, int intersect, double u, int inbound, int inside);

int addIntersect(struct Node *list, double x, double y, double z, int intersect, double u1, double u2,
                 int inbound, int is1, int ie1, int is2, int ie2);

int length(struct Node *list);

int sameNode(struct Node node1, struct Node node2);

void addNode(struct Node *list, struct Node nodeIn);

struct Node *getNode(struct Node *list, struct Node inNode);

struct Node *getNextNode(struct Node *list);

void copyNode(struct Node *node_out, struct Node node_in);

void printNode(struct Node *list, char *str);

int intersectInList(struct Node *list, double x, double y, double z);

void insertIntersect(struct Node *list, double x, double y, double z, double u1, double u2, int inbound,
                     double x2, double y2, double z2);

void insertAfter(struct Node *list, double x, double y, double z, int intersect, double u, int inbound,
                 double x2, double y2, double z2);

double gridArea(struct Node *grid);

int isIntersect(struct Node node);

int getInbound( struct Node node );

struct Node *getLast(struct Node *list);

int getFirstInbound( struct Node *list, struct Node *nodeOut);

void getCoordinate(struct Node node, double *x, double *y, double *z);

void getCoordinates(struct Node *node, double *p);

void setCoordinate(struct Node *node, double x, double y, double z);

void setInbound(struct Node *interList, struct Node *list);

int isInside(struct Node *node);

int insidePolygon( struct Node *node, struct Node *list);

#endif
