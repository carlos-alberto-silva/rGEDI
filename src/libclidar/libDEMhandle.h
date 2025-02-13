
/*#########################*/
/*# Functions to handle  #*/
/*# DEMs to process lidar #*/
/*#########################*/


/*#######################################*/
/*# Copyright 2006-2016, Steven Hancock #*/
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
/*structure to hold DEM*/

typedef struct{
  double minX;     /*minimum dem x*/
  double maxX;     /*maximum dem x*/
  double minY;     /*minimum dem y*/
  double maxY;     /*maximum dem y*/
  double minZ;     /*minimum dem height*/
  double maxZ;     /*maximum dem height*/
  double *z;       /*the dem*/
  int nX;          /*number of x elements*/
  int nY;          /*number of y elements*/
  float res;       /*dem bounds resolution*/
  double noData;   /*missing data flag*/
}demStruct;


/*#######################################*/
/*functions*/

demStruct *readTifDEM(char *,double,double,double,double);
demStruct *tidyDEMstruct(demStruct *);
double findDEMelev(double,double,demStruct *);

/*the end*/
/*#######################################*/

