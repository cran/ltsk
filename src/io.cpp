#include "io.h"
//#include <iostream>
#include <vector>
//#include <fstream>
//#include <cassert>
//#include <string>
#include <cstdlib>
#include <cmath>

using namespace std;
int r_file_s(double *x, int *xn, double *y, int *yn, double *z, int *zn, ANNpointArray& pointArr)
{
   int i;
   int npts = *xn;
   for(i=0; i< *xn; i++)
   {
        pointArr[i][0]=x[i];
   	pointArr[i][1]=y[i];
	plist.push_back(point(x[i],y[i],z[i]));
   }

   return npts;
}

void r_file_q(double *lon, int *lon_n, double *lat, int *lat_n, std::vector<point_q>& list_q)
{
  int i;
  for(i=0; i< *lon_n; i++)
  {
    list_q.push_back(point_q(lon[i],lat[i]));
  }
}

void calculate(int& index, double& qlat, double& qlon, int& qtime)
{
  plist[index].distance = sqrt(pow((plist[index].x_coord - qlat), 2) + pow((plist[index].y_coord - qlon), 2));
  plist[index].timediff = abs(plist[index].time - qtime);
}

void calculate(int& index, double& qlat, double& qlon)
{
  plist[index].distance = sqrt(pow((plist[index].x_coord - qlat), 2) + pow((plist[index].y_coord - qlon), 2));
}
