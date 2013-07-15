#ifndef IO_H
#define IO_H

#pragma once
#include <vector>
//#include <string>
#include "point.h"
#include "ANN/ANN.h"
#include "R.h"

//void r_file_s(const char* filename, int lngCol, int latCol, int timeCol, int zCol );
//int r_file_s(const char* filename, ANNpointArray& pointArr, int longCol, int latCol, int zCol);  
//void r_file_q(const char* filename, std::vector<point_q>& list_q, int xCol, int yCol, int dayCol);
//void r_file_q(const char* filename, std::vector<point_q>& list_q, int xCol, int yCol);
void calculate(int& index, double& qlat, double& qlon, int& qtime);
void calculate(int& index, double& qlat, double& qlon);
int r_file_s(double* x, int *xn, double *y, int *yn, double *z, int *zn, ANNpointArray& pointArr);
void r_file_q(double *lon, int *lon_n, double *lat, int *lat_n, std::vector<point_q>& list_q);
#endif
