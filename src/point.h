#ifndef POINT_H
#define POINT_H
#pragma once

#define DBL_MIN         2.2250738585072014e-308 /* min positive value */
//#include <iomanip>    
//#include <iostream>    
#include <vector>    


class point
{
 private:
 public:
  double x_coord;
  double y_coord;
  int time;
  double z_value;
  double timediff;
  double distance;		
  //		double z_2;
  //double z_3;
 public:
  //point():longitude(0),latitude(0),time(0) {}
 point(double lon, double lat, int day, double z1, int time_diff=0, double Distance_=0):x_coord(lon),y_coord(lat),time(day),z_value(z1) {}

 /* point(double lon=0, double lat=0, int day=0, double z1=0, int time_diff=0, double Distance_=0):x_coord(lon),y_coord(lat),time(day),z_value(z1) {} */

 point(double lon, double lat, double z1):x_coord(lon),y_coord(lat),z_value(z1) {}
  
  
  inline bool operator < (point  other) const
  { return this->x_coord < other.x_coord; }
  inline bool operator >= (point  other) const
  { return this->y_coord >= other.y_coord; }

  inline const point  operator=( const point&  rhs )
    {
      if( this != &rhs )
	{
	  this->x_coord=rhs.x_coord;
	  this->y_coord=rhs.y_coord;
	  this->time=rhs.time;
	}
      return *this;
    }

  inline bool empty()
  {
    return !x_coord && !y_coord && !time;
  }

  //inline void print() { std::cout << std::setprecision(9)<< x_coord << ", " << y_coord << ", " << time << std::endl; }

};


class point_q
{
 private:
  double x_coord;
  double y_coord;
  int time;
 public:
 point_q(double lon=0, double lat=0, int day=0 ):x_coord(lon),y_coord(lat),time(day) {}
  inline double& x() { return x_coord; }
  inline double& y() { return y_coord; }
  inline int& z() { return time; }
		 
};

#ifdef DEFINE_GLOBALS
#define GLOBAL
#else // !DEFINE_GLOBALS
#define GLOBAL extern
#endif

GLOBAL std::vector<point> plist;

#endif
