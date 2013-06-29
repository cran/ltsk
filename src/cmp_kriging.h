//----------------------------------------------------------------------
// File:			CMP Kriging Library.h
// Programmer:		Jin Chen
// Last modified:	Aug/03/2010 (Release 0.1)
// Description:		Basic include file for calculating kriging.
//----------------------------------------------------------------------
// Copyright (c) 2010-2020 University of Iowa and Jin Chen.
// All Rights Reserved. 
//----------------------------------------------------------------------


#ifndef CMPKRIGING_H
#define CMPKRIGING_H

/* #ifdef WIN32 */
/*   #ifdef CMP_KRIGING_LIBRARY_EXPORTS */
/* 	 #define DLL_API __declspec(dllexport) */
/*   #else */
/* 	#define DLL_API __declspec(dllimport) */
/*   #endif */
/*  #else */
/*   #define DLL_API */
/* #endif */

#define DLL_API

#include <math.h>
#include <time.h>
#include <algorithm>
#include <functional>
#include <vector>
#include "point.h"
#include "R.h"

#define CMPKrigingVersion 		"1.0"			
#define CMPKrigingVersionCmt	""

#define MINNEIGHBOR  5			// define the min number of neighbors
#define MAXNEIGHBOR  2000		// define the max number of neighbors
#define NEIGHBORCELL 10			// define the distance neighbor cells to a 10 * 10 grids
#define MINBIN	 3			// define the min number of bins on distance and time
#define DISTBIN      15		// define the number of distance bins

#define RANDOM(x) (rand()%x)		// a function which can create a random number. This number is between 0-x.


using namespace std;

typedef long double		CMPKrigingCoord;	// coordinate data type

typedef CMPKrigingCoord*	CMPKrigingPoint;	// a point

// 05/01/2013: change Pairs to double restored not fixing the problem;
typedef	struct AvDistTimeSemi_Struct
{
	double AvDistance;
	double AvTimeDiff;
	double AvSemi;
	int Pairs;
} AvDistTimeSemi;

typedef	struct AvDistSemi_Struct
{
	double AvDistance;
	double AvSemi;
	int Pairs;
} AvDistSemi;

DLL_API CMPKrigingPoint		CMPKrigingAllocPt(
		int	dim,		// dimension
		CMPKrigingCoord 	c = 0);	// coordinate value (all equal)

DLL_API void CMPKrigingDeallocPt(
		CMPKrigingPoint		&p);				// deallocate 1 point

/* class NeighborNode */
/* { */
/* public: */
/* 	NeighborNode(double x,double y,int t,double z,double d,double td) */
/* 	  :x_coord(x),y_coord(y),time(t),z_value(z),distance(d),timediff(td){} */
/* 	double	x_coord; */
/* 	double	y_coord; */
/* 	int	time; */
/* 	double	z_value; */
/* 	double	distance; */
/* 	double	timediff; */

/* 	virtual	~NeighborNode(){}; */
/* }; */

class DistTimeSemi 
{
public:
	DistTimeSemi(double a, double b,double c):distance(a), timediff(b),semi(c){}
	double distance;
	double timediff;
	double semi;
};

class NeighborVector
{
  vector<int> neighbor_vector;
 public:
  point& operator[] (const int index)
    {
      return plist[index];
    }

  int get_index(const int num)
    {
      return neighbor_vector[num];
    }
  
  void push_back (const int indx)
  { neighbor_vector.push_back(indx);}

  void clear()
  { neighbor_vector.clear();}
  
  int size() const
  { return neighbor_vector.size(); }

};	// an vector of neighbor's index

enum ModelSelection {
		Init			= 0,					// Initial value
		Spherical	    = 1,					// Calculate kriging with Spherical model
		Gaussian	    = 2,					// Calculate kriging with Gaussian model
		Exponential		= 3,					// Calculate kriging with Exponential model
		Matern			= 4,					// Calculate kriging with Matern model
		IDW			= 5};					// Calculate IDW value					

typedef struct Variogram_Parameters_Struct
{
	double Hs;
	double nuggets;
	double sills;
	ModelSelection ms;
	double Ht;
	double nuggett;
	double sillt;
	ModelSelection mt;
	double kappa;
}VariogramParameters;

enum CalculateMethod{		// when user select IDW, user must specify calculate methond
		None	= 0,	// user doesn't select IDW (initial value)
		Plus	= 1,	// use addition as calculate methond
		Multiply= 2};	// use multiply as calculate methond

enum KrigingValueQuality {
		Empty		= 0,		// empty value(initial value)
		IDWValue	= 1,		// the returned value is IDW
		Good		= 2,		// good returned kriging value 
		LessNeighbors	= 3,		// less than 5 neighbors
		LessNeighborsInSpace		= 4,		// less neighbors in space
		LessNeighborsInTime		= 5,		// less neighbors in time
		LessNeighborsInSpaceAndTime 	= 6,		// less neighbors in space and time
		NotInverseMatrix		= 7,		// the gamma matrix can't be inversed
		OutofRange			= 8,		// the claculated kriging is negtive
		WrongParameter			= 9};		// wrong input parameter

enum MethodWhenLessNeighbor{	// when less neighbors, user can specify how to calculate
		NotCalculate	= 0,	// don't calculate when less neighbors, output -99999
		Mean	        = 1,	// calculate mean value based on neighbors instead of kriging value
		Kriging		= 2	// calculate kriging by space when less neighbors in time; calculate kriging by time when less neighbors in space
};

DLL_API void CMP_Kriging_By_Space_and_Time(		// calculate kriging by space and time
					   //		CMPKrigingPoint			queryPt,
		// query point, has query point location and time information
		NeighborVector&			neighbor,
		// neighbors are stored in this vector
		int				DistBin,
		// the number of distance bins
		int				TimeBin,
		// the number of time bins
		ModelSelection			m,
		// calculation model
		double				alpha,
		// when model selection is IDW, specify two parameters, this is for distance
		double				beta,
		// when model selection is IDW, specify two parameters, this is for time
		CalculateMethod			c,
		// when model selection is IDW, specify calculation method
		double&					kriging,
		// kriging value (returned)
		double& Hs,
		double& Ht,
		double&					sigma,
		// kriging variance (returned)
		KrigingValueQuality&	kvq);			// kriging value quality (returned)


void VariogramCalculationBySpaceAndTime(		// calculate variogram by space and time function
		AvDistTimeSemi**		&BinCube,		// neighbors are stored in this two demission array
		int				DistBin,		// the number of distance bins
		int				TimeBin,		// the number of time bins
		double				d_interval,		// the interval on distance
		double				t_interval,		// the interval on time
		VariogramParameters&	parameters);	// variogram parameters (returned)

void	Kappa_Calculate(						// calculate the parameter kappa and adjust nugs, sills, nugt, sillt
		NeighborVector&			neighbor,		// neighbors are stored in this vector
		vector<int>&			sample,			// sample neighbors are stored in this vector
		VariogramParameters&	parameters);	// variogram parameters (returned)

double  InterpolateMissingValueBySpaceAndTime(	// find the nearest three bins and use this three bins' value to interpolate
		AvDistTimeSemi**		&BinCube,		// neighbors are stored in this two demission array
		int				DistBin,		// the number of distance bins
		int				TimeBin,		// the number of time bins
		int				d,				// the index on distance Bin of missing value
		int				t);				// the index on time Bin of missing value

DLL_API void Kriging_By_Space(					// calculate kriging by space
	      //CMPKrigingPoint	queryPt,// query point, has query point location information
		NeighborVector&		neighbor,	// neighbors are stored in this vector
		int			DistBin,	// the number of distance bins
		ModelSelection          m,	// calculation model
		double			alpha,	// when model selection is IDW,	specify one parameter for distance
		double&			kriging,	// kriging value (returned)
		double&			sigma,	// kriging variance (returned)
		double&			Hs,	// 
		VariogramParameters&	parameters, // kriging variogram parameter (returned)
		KrigingValueQuality&	kvq);	// kriging value quality (returned)

void    Variogram_By_Space_Calculation(			// calculate variogram by space function
		AvDistSemi*				&Bins,			// neighbors are stored in this one demission array
		int					DistBin,		// the number of distance bins
		double					d_interval,		// the interval on distance
		VariogramParameters&	parameters);	// variogram parameters (returned)

double  InterpolateMissingValueBySpace(
				       // find the nearest three bins and use this three bins' value to interpolate
		AvDistSemi*		&Bin,			// neighbors are stored in this two demission array
		int			DistBin,		// the number of distance bins
		int			d);				// the index on distance Bin of missing value


#endif
