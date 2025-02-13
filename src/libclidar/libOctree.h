/*#########################*/
/*# Structures for octrees #*/
/*# in lidar programs      #*/
/*##########################*/

/*#######################################*/
/*# Copyright 2006-2017, Steven Hancock #*/
/*# The program is distributed under    #*/
/*# the terms of the GNU General Public #*/
/*# License.    svenhancock@gmail.com   #*/
/*#######################################*/


/*########################################################################*/
/*# This file is part of libCLidar.                                      #*/
/*#                                                                      #*/
/*# libCLidar is free software: you can redistribute it and/or modify    #*/
/*# it under the terms of the GNU General Public License as published by #*/
/*# the Free Software Foundation, either version 3 of the License, or    #*/
/*#  (at your option) any later version.                                 #*/
/*#                                                                      #*/
/*# libCLidar is distributed in the hope that it will be useful,         #*/
/*# but WITHOUT ANY WARRANTY; without even the implied warranty of       #*/
/*#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #*/
/*#   GNU General Public License for more details.                       #*/
/*#                                                                      #*/
/*#    You should have received a copy of the GNU General Public License #*/
/*#    along with libClidar.  If not, see <http://www.gnu.org/licenses/>.#*/
/*########################################################################*/


/*#######################################*/
/*the tree itself*/

typedef struct{
  void **tree;  /*recursive octree pointer*/
  uint32_t mapInd;   /*map index*/ 
}treeStruct;


/*#######################################*/
/*octree structure. 2D for now*/

typedef struct{
  int nLevel;   /*number of levels of octree*/
  float res;    /*resolution of top level*/
  double minX;  /*min x of octree*/
  double maxX;  /*max x of octree*/
  double minY;  /*min y of octree*/
  double maxY;  /*max y of octree*/
  int nX;       /*number of x pixels at top level*/
  int nY;       /*number of y pixels at top level*/

  treeStruct **tree;  /*the octree*/

  int **mapFile;       /*map to file indices*/
  uint32_t **mapPoint; /*map to point indices*/
  uint32_t *nIn;      /*number of points within lowest level*/
  uint32_t nMaps;      /*number of map array elements*/
}octreeStruct;


/*#########################################*/
/*structure to hold map of points to use*/

typedef struct{
  uint32_t nPoints;  /*number of points to test*/
  int *fList;        /*list of file indices*/
  uint32_t *pList;   /*list of point indices*/
}pointMapStruct;


/*#######################################################################*/
/*functions*/

octreeStruct *allocateOctree(int,int,double,double,double,double);
octreeStruct *tidyOctree(octreeStruct *);
pointMapStruct *mapFromOctree(int *,int,octreeStruct *,double,double,double,double);
void fillOctree(double,double,double,int,uint32_t,octreeStruct *);


/*the end*/
/*#######################################*/

