//----------------------------------------------------------------------
// File:			CMP Kriging Library.cpp
// Programmer:		Jin Chen
// Last modified:	Aug/03/2010 (Release 0.1)
// Description:		CPP code file for calculating CMP kriging.
//----------------------------------------------------------------------
// Copyright (c) 2010-2020 University of Iowa and Jin Chen.
// All Rights Reserved. 
//----------------------------------------------------------------------
// CMP kriging library.cpp : Defines the entry point for the DLL application.
//

// #include "stdafx.h"
#include "cmp_kriging.h"
//#include <fstream>


extern   "C"  
{  
  // #include   <../include/f2c.h>  
  // #include   <../include/clapack.h>  
  // #include   <../include/blaswrap.h>  
#include   "f2c.h"  
#include   "clapack.h"  
#include   "blaswrap.h"  
} 


#ifdef _MANAGED
#pragma managed(push, off)
#endif

// bool APIENTRY DllMain( HMODULE hModule,
//                        DWORD  ul_reason_for_call,
//                        LPVOID lpReserved
// 					 )
// {
//     return TRUE;
// }

CMPKrigingPoint CMPKrigingAllocPt(int dim, CMPKrigingCoord c)		// allocate 1 point
{
  CMPKrigingPoint p = new CMPKrigingCoord[dim];
  for (int i = 0; i < dim; i++) p[i] = c;
  return p;
}

void CMPKrigingDeallocPt(CMPKrigingPoint &p)						// deallocate 1 point
{
  delete [] p;
  p = NULL;
}

// add Hs, Ht
void CMP_Kriging_By_Space_and_Time(   // calculate kriging by space and time
				   //		   CMPKrigingPoint	queryPt,
				   // query point,	has query point location and time information
				   NeighborVector&	neighbor,
				   // neighbors are stored in this vector
				   int			DistBin,
				   // the number of distance bins
				   int			TimeBin,
				   // the number of time bins
				   ModelSelection	m,
				   // calculation model
				   double		alpha,
				   // when model selection is IDW,	specify two parameters, this is for distance
				   double		beta,
				   // when model selection is IDW,	specify two parameters, this is for time
				   CalculateMethod	c,
				   // when model selection is IDW,	specify calculation method
				   double&		kriging,
				   // kriging value (returned)		   
				   double&		Hs,
				   double&		Ht,
				   double&		sigma,
				   // kriging variance (returned)
				   KrigingValueQuality&	kvq
				   // kriging value quality (returned)
				      )
{
  // declare local variables
  kriging = -99999;
  sigma = -99999;
  kvq = Empty;	

  int    npts;
  double TempX1,TempY1,TempZ1,TempT1,TempX2,TempY2,TempT2;
  double distance;
  double timediff;
  double tempAlpha;
  double tempBeta;
  int    tmax, tmin;
  double xmax, xmin, ymax, ymin;
  int    xn, yn, tn;
  int    xi, yi, ti;
  vector<int> ***neighborcube;
  double x_interval, y_interval;
  int	   sam_npts;
  vector<int> sampleNeighbor;
  vector<bool> sampling;
// int    TempR; // random number
  int    ns, nt;
  double dist_max;
  double time_max;
  double semi;
  int    di;
  double d_interval;
  double t_interval;
  vector<int> subNeighbor;
  double epsilon;
  AvDistTimeSemi** BinCube;
  VariogramParameters	parameters;
  //double krigingbar;
	
	
  npts = (int)neighbor.size();

  // if ((queryPt == NULL )||(DistBin <= 3)||(TimeBin <= 3))
  if ((DistBin <= 3)||(TimeBin <= 3))    
    {
      kriging	= -99999;
      sigma	= -99999;
      kvq	= WrongParameter;
    }
  else if ((m == IDW)&&(c != Plus)&&(c != Multiply))
    {
      kriging	= -99999;
      sigma	= -99999;
      kvq	= WrongParameter;
    }
  else if (npts == 0)
    {
      kriging	= -99999;
      sigma	= -99999;
      kvq	= LessNeighborsInSpaceAndTime;
    }
  else
    {
      if (m == IDW) // calculate IDW
	{
	  TempX2 = 0;
	  TempY2 = 0;
	  for (int i = 0; i < npts; i++)
	    {
	      //TempX1   = neighbor[i].x_coord;
	      //TempY1   = neighbor[i].y_coord;
	      TempZ1   = neighbor[i].z_value;
	      distance = neighbor[i].distance;
	      timediff = neighbor[i].timediff;

	      if (distance == 0)
		{
		  tempAlpha = 1;
		}
	      else if (distance != 0)
		{
		  tempAlpha = pow((double)distance,alpha*(-1));
		}

	      if (timediff == 0)
		{
		  tempBeta = 1;
		}
	      else if (timediff != 0)
		{
		  tempBeta = pow((double)timediff,beta*(-1));
		}

	      if (c == Plus)
		{
		  TempX2   = TempX2 + TempZ1 * (tempAlpha + tempBeta);
		  TempY2   = TempY2 + (tempAlpha + tempBeta);
		}
	      else if (c == Multiply)
		{
		  TempX2   = TempX2 + TempZ1 * tempAlpha * tempBeta;
		  TempY2   = TempY2 + tempAlpha * tempBeta;
		}
	    }
	  kriging	= -99999;
	  sigma	= -99999;
	  kvq		= IDWValue;
	  if (TempY2 != 0)
	    {
	      kriging	= TempX2 / TempY2;
	    }
	}
      else // calculate kriging
	{
	  if (npts <= MINNEIGHBOR)
	    {
	      kriging	= -99999;
	      sigma	= -99999;
	      kvq		= LessNeighborsInSpaceAndTime;
	    }
	  else // npts > MINNEIGHBOR
	    {
	      // setup random value enable 
	      // srand((unsigned)time(NULL));

	      // search the range of x, y and time
	      xmin = neighbor[0].x_coord;
	      xmax = neighbor[0].x_coord;
	      ymax = neighbor[0].y_coord;
	      ymin = neighbor[0].y_coord;
	      tmax = neighbor[0].time;
	      tmin = neighbor[0].time;
	      sampling.clear();
	      sampling.push_back(false);

	      for (int i = 1; i < npts; i++)
		{
		  if (xmin > neighbor[i].x_coord)
		    {
		      xmin = neighbor[i].x_coord;
		    }
		  if (xmax < neighbor[i].x_coord)
		    {
		      xmax = neighbor[i].x_coord;
		    }
		  if (ymin > neighbor[i].y_coord)
		    {
		      ymin = neighbor[i].y_coord;
		    }
		  if (ymax < neighbor[i].y_coord)
		    {
		      ymax = neighbor[i].y_coord;
		    }
		  if (tmin > neighbor[i].time)
		    {
		      tmin = neighbor[i].time;
		    }
		  if (tmax < neighbor[i].time)
		    {
		      tmax = neighbor[i].time;
		    }
		  sampling.push_back(false);
		}

	      // declare neighbor cubes
	      xn = NEIGHBORCELL;
	      yn = NEIGHBORCELL;
	      tn = (tmax - tmin) + 1;

	      neighborcube = new vector<int>** [xn];
	      for (int i = 0; i < xn; i++)
		{
		  neighborcube[i] = new vector<int>* [yn];
		  for (int j = 0; j < yn; j++)
		    {
		      neighborcube[i][j] = new vector<int> [tn];
		      for (int k = 0; k < tn; k++)
			{
			  neighborcube[i][j][k].clear();
			}
		    }
		}

	      // insert the neighbors into neighbor cubes
	      x_interval = abs(xmax - xmin) / xn;
	      y_interval = abs(ymax - ymin) / yn;
	      for (int i = 0; i < npts; i++)
		{
		  xi = (int)floor((neighbor[i].x_coord - xmin) / x_interval);
		  if (xi == xn) // (xi == ((neighbor[i].x_coord - xmin) / x_interval))
		    {
		      xi = xi - 1;
		    }

		  yi = (int)floor((neighbor[i].y_coord - ymin) / y_interval);
		  if (yi == yn) // (yi == ((neighbor[i].y_coord - ymin) / y_interval))
		    {
		      yi = yi - 1;
		    }

		  ti = neighbor[i].time - tmin;

		  if ((xi < xn)&&(yi < yn)&&(ti < tn))
		    {
		      neighborcube[xi][yi][ti].push_back(i);
		    }					
		}
	      // all neighbors are inserted into neighbor cubes

	      // calculate ns number of space regions, nt number of days with neighbors
	      ns = 0; 
	      nt = 0;

	      for (int i = 0; i < xn; i++)
		{
		  for (int j = 0; j < yn; j++)
		    {
		      for (int k = 0; k < tn; k++)
			{
			  if (neighborcube[i][j][k].size() > 0)
			    {
			      ns = ns + 1;
			      k = tn;
			    }
			}
		    }
		}

	      for (int k = 0; k < tn; k++)
		{
		  for (int i = 0; i < xn; i++)
		    {
		      for (int j = 0; j < yn; j++)
			{
			  if (neighborcube[i][j][k].size() > 0)
			    {
			      nt = nt + 1;
			      j = yn;
			      i = xn;
			    }
			}
		    }
		}

	      if ((ns <= MINBIN)||(nt <= MINBIN))
		{
		  kriging	= -99999;
		  sigma = -99999;
					
		  if ((ns <= MINBIN)&&(nt <= MINBIN))
		    {
		      kvq	= LessNeighborsInSpaceAndTime;
		    }
		  else if ((ns <= MINBIN)&&(nt > MINBIN))
		    {
		      kvq	= LessNeighborsInSpace;
		    }
		  else if ((ns > MINBIN)&&(nt <= MINBIN))
		    {
		      kvq	= LessNeighborsInTime;
		    }
		}
	      else if ((ns > MINBIN)&&(nt > MINBIN))
		{
		  if (npts > MAXNEIGHBOR) // the neighbors are greater than maxium, randomly select sampling neighbors to reduce calculation
		    {
		      // randomly select subneighbors from neighbor cubes
		      sampleNeighbor.clear();
		      for (int i = 0; i < xn; i++)
			{
			  for (int j = 0 ; j < yn; j++)
			    {
			      for (int k = 0; k < tn; k++)
				{
				  if (neighborcube[i][j][k].size() > 0)
				    {
				      // calculate randomly sampling number for this neighbor cube
				      sam_npts = (int)(MAXNEIGHBOR * neighborcube[i][j][k].size() / npts);

				      // if sampling number is greater than number in cube, output all neighbors in the cube
				      if (sam_npts >= (int)neighborcube[i][j][k].size())
					{
					  for (int p = 0; p < (int)neighborcube[i][j][k].size(); p++)
					    {
					      sampleNeighbor.push_back(neighborcube[i][j][k][p]);
					    }
					}
				      else
					{
					  // output one sampling neighbor at least
					  if (sam_npts == 0)
					    {
					      sam_npts = 1;
					    }

					  // output sam_npts neighbors from cube to subneighbors
					  //for (int p = 0; p < sam_npts; p++)
					  //  {
					  //    TempR = rand() % (int)neighborcube[i][j][k].size();
					  //    while (sampling[neighborcube[i][j][k][TempR]] == true)
						//{
						 // TempR = rand() % (int)neighborcube[i][j][k].size();
						//}
					     // sampling[neighborcube[i][j][k][TempR]] = true;
					      //sampleNeighbor.push_back(neighborcube[i][j][k][TempR]);
					    //}
						
					  // output first sam_npts neighbors from cube to subneighbors
					  for (int p = 0; p < sam_npts; p++)
					    {
							sampleNeighbor.push_back(neighborcube[i][j][k][p]);
					    }
					}								
				    }
				}
			    }
			}
		    }
		  else // the neighbors are less than maxium value, use all neighbors to calculate kriging value
		    {
		      sampleNeighbor.clear();
		      for (int i = 0; i < (int)neighbor.size(); i++)
			{
			  sampleNeighbor.push_back(i);
			}
		    }

		  // Bin sample neighbors in space and time
		  xmax = neighbor[sampleNeighbor[0]].x_coord;
		  xmin = neighbor[sampleNeighbor[0]].x_coord;
		  ymax = neighbor[sampleNeighbor[0]].y_coord;
		  ymin = neighbor[sampleNeighbor[0]].y_coord;
		  tmax = neighbor[sampleNeighbor[0]].time;
		  tmin = neighbor[sampleNeighbor[0]].time;
		  for (int i = 1; i < (int)sampleNeighbor.size(); i++)
		    {
		      if (xmax < neighbor[sampleNeighbor[i]].x_coord)
			{
			  xmax = neighbor[sampleNeighbor[i]].x_coord;
			}
		      if (xmin > neighbor[sampleNeighbor[i]].x_coord)
			{
			  xmin = neighbor[sampleNeighbor[i]].x_coord;
			}
		      if (ymax < neighbor[sampleNeighbor[i]].y_coord)
			{
			  ymax = neighbor[sampleNeighbor[i]].y_coord;
			}
		      if (ymin > neighbor[sampleNeighbor[i]].y_coord)
			{
			  ymin = neighbor[sampleNeighbor[i]].y_coord;
			}
		      if (tmax < neighbor[sampleNeighbor[i]].time)
			{
			  tmax = neighbor[sampleNeighbor[i]].time;
			}
		      if (tmin > neighbor[sampleNeighbor[i]].time)
			{
			  tmin = neighbor[sampleNeighbor[i]].time;
			}
		    }

		  dist_max = 0.75 * sqrt(pow((xmax - xmin), 2) + pow((ymax - ymin), 2));
		  time_max = 0.75 * (tmax - tmin); 

		  // initial Bin Cube
		  BinCube = new AvDistTimeSemi* [DistBin];
		  for (int i = 0; i < DistBin; i++)
		    {
		      BinCube[i] = new AvDistTimeSemi [TimeBin];
		      for (int j = 0; j < TimeBin; j++)
			{
			  BinCube[i][j].AvDistance = 0;
			  BinCube[i][j].AvSemi = 0;
			  BinCube[i][j].AvTimeDiff = 0;
			  BinCube[i][j].Pairs = 0;
			}
		    }

		  // ofstream outbin1("./output/bin1.csv", ofstream::out);
		  //   outbin1 << "di,ti,i,j,distance,timediff,semi" << endl;

		  d_interval = (double)dist_max / DistBin;
		  t_interval = (double)time_max / TimeBin;
		  // insert sub neighbors into distance-time bins
		  for (int i = 0; i < (int)sampleNeighbor.size(); i++)
		    {
		      for (int j = 0; j < i; j++)
			{
			  distance = sqrt(pow((neighbor[sampleNeighbor[i]].x_coord - neighbor[sampleNeighbor[j]].x_coord), 2) + pow((neighbor[sampleNeighbor[i]].y_coord - neighbor[sampleNeighbor[j]].y_coord), 2));
			  timediff = abs(neighbor[sampleNeighbor[i]].time - neighbor[sampleNeighbor[j]].time);

			  if ((distance <= dist_max)&&(timediff <= time_max))
			    {
			      semi = 0.5 * pow((neighbor[sampleNeighbor[i]].z_value - neighbor[sampleNeighbor[j]].z_value), 2);

			      di = (int)floor(distance / d_interval);
			      if (di == DistBin)
				{
				  di = di - 1;
				}

			      ti = (int)floor(timediff / t_interval);
			      if (ti == TimeBin)
				{
				  ti = ti - 1;
				}

			      //outbin1 << di << "," << ti << "," << i << "," << j << "," << distance << "," << timediff << "," << semi << endl;

			      if ((di < DistBin)&&(ti < TimeBin))
				{
				  BinCube[di][ti].AvDistance = BinCube[di][ti].AvDistance + distance;
				  BinCube[di][ti].AvSemi = BinCube[di][ti].AvSemi + semi;
				  BinCube[di][ti].AvTimeDiff = BinCube[di][ti].AvTimeDiff + timediff;
				  BinCube[di][ti].Pairs = BinCube[di][ti].Pairs + 1;
				}								
			    }
			}
		    }
		  //outbin1.close();

		  for (int i = 0; i < DistBin; i++)
		    {
		      for (int j = 0; j < TimeBin; j++)
			{
			  if (BinCube[i][j].Pairs != 0)
			    {
			      BinCube[i][j].AvDistance = BinCube[i][j].AvDistance / BinCube[i][j].Pairs;
			      BinCube[i][j].AvSemi = BinCube[i][j].AvSemi / BinCube[i][j].Pairs;
			      BinCube[i][j].AvTimeDiff = BinCube[i][j].AvTimeDiff / BinCube[i][j].Pairs;
			    }
			}
		    }
		  // Bin sample neighbors in space and time

		  // "TEST" output Bin
		  // ofstream outBinCube("./output/BinCube.csv", ofstream::out);
		  //   for (int i = 0; i < DistBin; i++)
		  //   {
		  //   for (int j = 0; j < TimeBin; j++)
		  //   {
		  //   outBinCube << "," << BinCube[i][j].Pairs;
		  //   }
		  //   outBinCube << endl;
		  //   }

		  //   outBinCube << endl;
		  //   outBinCube << endl;
		  //   outBinCube << endl;

		  //   for (int i = 0; i < DistBin; i++)
		  //   {
		  //   for (int j = 0; j < TimeBin; j++)
		  //   {
		  //   outBinCube << "," << BinCube[i][j].AvDistance;
		  //   }
		  //   outBinCube << endl;
		  //   }

		  //   outBinCube << endl;
		  //   outBinCube << endl;
		  //   outBinCube << endl;

		  //   for (int i = 0; i < DistBin; i++)
		  //   {
		  //   for (int j = 0; j < TimeBin; j++)
		  //   {
		  //   outBinCube << "," << BinCube[i][j].AvSemi;
		  //   }
		  //   outBinCube << endl;
		  //   }

		  //   outBinCube << endl;
		  //   outBinCube << endl;
		  //   outBinCube << endl;

		  //   for (int i = 0; i < DistBin; i++)
		  //   {
		  //   for (int j = 0; j < TimeBin; j++)
		  //   {
		  //   outBinCube << "," << BinCube[i][j].AvTimeDiff;
		  //   }
		  //   outBinCube << endl;
		  //   }

		  //   outBinCube.close();
		  // "TEST" output Bin

		  // calculate variogram
		  VariogramCalculationBySpaceAndTime(// calculate variogram by space and time function
						     BinCube,					// neighbors are stored in this two demission array
						     DistBin,					// the number of distance bins
						     TimeBin,					// the number of time bins
						     d_interval,					// the interval on distance
						     t_interval,					// the interval on time
						     parameters);				// variogram parameters (returned)
		  // calculate variogram

		  // estimate the interaction parameter kappa, return kappa and adjusted nugs,sills,nugt,sillt
		  Kappa_Calculate(				// calculate the parameter kappa and adjust nugs, sills, nugt, sillt
				  neighbor,					// neighbors are stored in this vector
				  sampleNeighbor,				// sample neighbors are stored in this vector
				  parameters);				// variogram parameters (returned)						

		  // Calculate Gamma matrix
		  subNeighbor.clear();
		  subNeighbor.push_back(0);
		  subNeighbor.clear();
		  for (int i = 0; i < (int)sampleNeighbor.size(); i++)
		    {
		      distance = neighbor[sampleNeighbor[i]].distance;
		      timediff = neighbor[sampleNeighbor[i]].timediff;

		      if ((distance <= parameters.Hs)||(timediff <= parameters.Ht))
			{
			  subNeighbor.push_back(sampleNeighbor[i]);
			}
		    }

		  if (subNeighbor.size() <= 1)
		    {
		      if (sampleNeighbor.size() <= 50)
			{
			  subNeighbor.clear();
			  for (int i = 0; i < (int)sampleNeighbor.size(); i++)
			    {
			      subNeighbor.push_back(sampleNeighbor[i]);
			    }
			}
		      else
			{
			  subNeighbor.clear();
			  for (int i = 0; i < (int)(0.1 * subNeighbor.size()); i++)
			    {
			      subNeighbor.push_back(sampleNeighbor[i]);
			    }
			}				
		    }

		  if ((parameters.nuggets == 0)&&(parameters.nuggett == 0))
		    {
		      epsilon = 0.00001;
		    }
		  else
		    {
		      epsilon = 0;
		    }

		  npts = (int)subNeighbor.size();
		  //cout << "npts " << npts << endl;

		  doublereal * GammaMatrix = new doublereal[(npts+1)*(npts+1)];
		  int TempK  = 0;
		  int TempP  = 0;
		  double Temp1;
		  double Temp2;

		  // TEST output neighbor
		  // ofstream outN("./output/neighbor.csv", ofstream::out);
		  //   outN << "i,j,t,z" << endl;
		  //   for (int i = 0; i < npts; i++)
		  //   {
		  //   outN << neighbor[subNeighbor[i]].x_coord << "," << neighbor[subNeighbor[i]].y_coord << "," <<  neighbor[subNeighbor[i]].time << "," <<  neighbor[subNeighbor[i]].z_value << endl;
		  //   }
		  //   outN << endl;
		  //   outN << endl;

		  //   outN << "ms,Hs,nuggets,sills,mt,Ht,nuggett,sillt,kappa,epsilon" << endl;
		  //   outN << parameters.ms << "," << parameters.Hs << "," << parameters.nuggets << "," << parameters.sills
		  // 	 << "," << parameters.mt << "," << parameters.Ht << "," << parameters.nuggett << "," << parameters.sillt
		  // 	 << "," << parameters.kappa << "," << epsilon << endl;

		  //   outN.close();
		  // TEST output neighbor
		  Hs = parameters.Hs;
		  Ht = parameters.Ht;		  

		  for (int i = 0; i < npts + 1; i++)
		    {
		      for (int j = 0; j < npts + 1; j++)
			{
			  TempK = j + (npts + 1) * i;
			  if (j == i)
			    {
			      GammaMatrix[TempK] = 0;
			    }
			  else
			    {
			      if (j > i)
				{
				  if (j == npts)
				    {
				      GammaMatrix[TempK] = 1;
				    }
				  else
				    {
				      TempX1 = neighbor[subNeighbor[i]].x_coord;
				      TempY1 = neighbor[subNeighbor[i]].y_coord;
				      TempT1 = neighbor[subNeighbor[i]].time;

				      TempX2 = neighbor[subNeighbor[j]].x_coord;
				      TempY2 = neighbor[subNeighbor[j]].y_coord;
				      TempT2 = neighbor[subNeighbor[j]].time;

				      distance = sqrt((TempX1-TempX2)*(TempX1-TempX2)+(TempY1-TempY2)*(TempY1-TempY2));
				      timediff = abs(TempT1 - TempT2);

				      switch (parameters.ms)
					{
					case 1: // Spherical
					  {
					    if (distance >= parameters.Hs)
					      {
						Temp1 = parameters.nuggets + parameters.sills;
					      }
					    else
					      {
						//Temp1 = nugs + sills * (1.5 * DistanceArray[i][j] / Hs - 0.5 * pow(DistanceArray[i][j] / Hs, 3));
						Temp1 = parameters.nuggets + parameters.sills * (1.5 * distance / parameters.Hs - 0.5 * pow(distance / parameters.Hs, 3));
					      }
					    break;
					  }
					case 2: // Gaussian
					  {
					    //Temp1 = nugs + sills * (1 - exp(-3 * pow(DistanceArray[i][j] / Hs, 2)));
					    Temp1 = parameters.nuggets + parameters.sills * (1 - exp(-3 * pow(distance / parameters.Hs, 2)));
					    break;
					  }
					case 3: // Exponential
					  {
					    //Temp1 = nugs + sills * (1 - exp(-3 * DistanceArray[i][j] / Hs));
					    Temp1 = parameters.nuggets + parameters.sills * (1 - exp(-3 * distance / parameters.Hs));
					    break;
					  }
					case 4: // Matern
					  {
					    //Temp1 = nugs + sills * (1 - (1 + 4.5 * DistanceArray[i][j] / Hs) * exp(-4.5 * DistanceArray[i][j] / Hs));
					    Temp1 = parameters.nuggets + parameters.sills * (1 - (1 + 4.5 * distance / parameters.Hs) * exp(-4.5 * distance / parameters.Hs));
					    break;
					  }
					default:
					  {
					    Temp1 = 0;
					    break;
					  }
					}
				      switch (parameters.mt)
					{
					case 1: // Spherical
					  {
					    if (timediff >= parameters.Ht)
					      {
						Temp2 = parameters.nuggett + parameters.sillt;
					      }
					    else
					      {
						//Temp2 = nugt + sillt * (1.5 * TimeDiffArray[i][j] / Ht - 0.5 * pow(TimeDiffArray[i][j] / Ht, 3));
						Temp2 = parameters.nuggett + parameters.sillt * (1.5 * timediff / parameters.Ht - 0.5 * pow(timediff / parameters.Ht, 3));
					      }
					    break;
					  }
					case 2: // Gaussian
					  {
					    //Temp2 = nugt + sillt * (1 - exp(-3 * pow(TimeDiffArray[i][j] / Ht, 2)));
					    Temp2 = parameters.nuggett + parameters.sillt * (1 - exp(-3 * pow(timediff / parameters.Ht, 2)));
					    break;
					  }
					case 3: // Exponential
					  {
					    //Temp2 = nugt + sillt * (1 - exp(-3 * TimeDiffArray[i][j] / Ht));
					    Temp2 = parameters.nuggett + parameters.sillt * (1 - exp(-3 * timediff / parameters.Ht));
					    break;
					  }
					case 4: // Matern
					  {
					    //Temp2 = nugt + sillt * (1 - (1 + 4.5 * TimeDiffArray[i][j] / Ht) * exp(-4.5 * TimeDiffArray[i][j] / Ht));
					    Temp2 = parameters.nuggett + parameters.sillt * (1 - (1 + 4.5 * timediff / parameters.Ht) * exp(-4.5 * timediff / parameters.Ht));
					    break;
					  }
					default:
					  {
					    Temp2 = 0;
					    break;
					  }
					}
				      GammaMatrix[TempK] = (doublereal)(Temp1 + Temp2 - parameters.kappa * Temp1 * Temp2 + epsilon);
				    }
				}
			      else // (p < q)
				{
				  if (i == npts)
				    {
				      GammaMatrix[TempK] = 1;
				    }
				  else
				    {
				      TempP =  i + (npts + 1) * j;
				      GammaMatrix[TempK] = GammaMatrix[TempP];
				    }
				}
			    }
			}
		    }

		  // TEST output Gamma matrix
		  // ofstream outGamma("./output/Gamma.csv", ofstream::out);
		  //   for (int i = 0; i < npts+1; i++)
		  //   {
		  //   for (int j = 0; j < npts+1; j++)
		  //   {
		  //   outGamma << "," << GammaMatrix[i*(npts+1)+j];
		  //   }
		  //   outGamma << endl;
		  //   }
		  //   outGamma.close();
		  // TEST output Gamma matrix

		  // calculate g vector by distance and time
		  doublereal *gVector;
		  doublereal *lamda;
		  doublereal dVector;
		  doublereal tVector;
		  gVector	= new doublereal [npts+1];
		  lamda = new doublereal [npts+1];

		  for (int i = 0; i < npts; i++)
		    {
		      dVector =  neighbor[subNeighbor[i]].distance;
		      tVector =  neighbor[subNeighbor[i]].timediff;

		      switch (parameters.ms)
			{
			case 1: // Spherical
			  {
			    if (dVector > parameters.Hs)
			      {
				Temp1 = parameters.nuggets + parameters.sills;
			      }
			    else
			      {
				Temp1 = parameters.nuggets + parameters.sills * (1.5 * dVector / parameters.Hs - 0.5 * pow(dVector / parameters.Hs, 3));
			      }
			    break;
			  }
			case 2: // Gaussian
			  {
			    Temp1 = parameters.nuggets + parameters.sills * (1 - exp(-3 * pow(dVector / parameters.Hs, 2)));
			    break;
			  }
			case 3: // Exponential
			  {
			    Temp1 = parameters.nuggets + parameters.sills * (1 - exp(-3 * dVector / parameters.Hs));
			    break;
			  }
			case 4: // Matern
			  {
			    Temp1 = parameters.nuggets + parameters.sills * (1 - (1 + 4.5 * dVector / parameters.Hs) * exp(-4.5 * dVector / parameters.Hs));
			    break;
			  }
			default:
			  {
			    Temp1 = 0;
			    break;
			  }
			}
		      switch (parameters.mt)
			{
			case 1: // Spherical
			  {
			    if (tVector > parameters.Ht)
			      {
				Temp2 = parameters.nuggett + parameters.sillt;
			      }
			    else
			      {
				Temp2 = parameters.nuggett + parameters.sillt * (1.5 * tVector / parameters.Ht - 0.5 * pow(tVector / parameters.Ht, 3));
			      }
			    break;
			  }
			case 2: // Gaussian
			  {
			    Temp2 = parameters.nuggett + parameters.sillt * (1 - exp(-3 * pow(tVector / parameters.Ht, 2)));
			    break;
			  }
			case 3: // Exponential
			  {
			    Temp2 = parameters.nuggett + parameters.sillt * (1 - exp(-3 * tVector / parameters.Ht));
			    break;
			  }
			case 4: // Matern
			  {
			    Temp2 = parameters.nuggett + parameters.sillt * (1 - (1 + 4.5 * tVector / parameters.Ht) * exp(-4.5 * tVector / parameters.Ht));
			    break;
			  }
			default:
			  {
			    Temp2 = 0;
			    break;
			  }
			}
		      gVector[i] = (doublereal)(Temp1 + Temp2 - parameters.kappa * Temp1 * Temp2);
		    }
		  gVector[npts] = 1;

		  // TEST output g vector
		  // ofstream outg("./output/gvector.csv", ofstream::out);
		  //   for (int i = 0; i < npts+1; i++)
		  //   {
		  //   outg << gVector[i] << endl;
		  //   }
		  //   outg.close();
		  // TEST output g vector

		  for (int i = 0; i < npts + 1; i++)
		    {
		      lamda[i] = gVector[i];
		    }

		  //METHOD 1: calculate weight vector through GammaMatrix and gVector, GammaMatrix * gVector = weights
		  /*integer M = npts;
		    integer nrhs = 1;
		    integer lda = M;
		    integer * ipiv = new integer[(npts + 1)];
		    integer ldb = M;
		    integer INFO;

		    dgesv_(&M,&nrhs,GammaMatrix,&lda,ipiv,gVector,&ldb,&INFO);

		    delete [] ipiv;
		    delete [] GammaMatrix;*/
		  // the weights value in the gVector now

		  // METHOD 2: calculate weight vector through GammaMatrix and gVector, GammaMatrix * gVector = weights
		  char uplo = 'U';
		  integer M = npts + 1;
		  integer nrhs = 1;
		  integer lda = M;
		  integer * ipiv = new integer[M];
		  integer ldb = M;

		  doublereal * c = new doublereal[(npts + 1)];
		  integer lwork = npts + 1;
		  integer INFO=0;

		  // 		  cout << INFO << endl;
		  dsysv_(&uplo,&M,&nrhs,GammaMatrix,&lda,ipiv,gVector,&ldb,c,&lwork,&INFO);//
		  // 		  cout << INFO << endl;
		  
		  delete [] c;
		  delete [] ipiv;
		  delete [] GammaMatrix;
		  // the weights value in the gVector now

		  // TEST output lamda vector
		  // ofstream outl("./output/lamda.csv", ofstream::out);
		  //   for (int i = 0; i < npts+1; i++)
		  //   {
		  //   outl << gVector[i] << endl;
		  //   }
		  //   outl.close();
		  // TEST output g vector

		  if (INFO != 0)
		    {
		      kvq	= NotInverseMatrix;

		      kriging	= -99999;
		      sigma = -99999;

		      /*kriging	= 0;
			sigma = 0;
			krigingbar = 0;

			for (int i = 0; i < npts; i++)
			{
			krigingbar = krigingbar + neighbor[i].z_value;
			}
			krigingbar = krigingbar / npts;
			kriging = krigingbar;

			for (int i = 0; i < npts; i++)
			{
			sigma = sigma + (krigingbar - neighbor[i].z_value)*(krigingbar - neighbor[i].z_value);
			}
			if (npts > 1)
			{
			sigma = sigma / (npts - 1);
			}
			sigma = sqrt(sigma);*/
		    }
		  else
		    {
		      kriging = 0;
		      sigma = gVector[npts]; // 0;
		      for (int i = 0; i < npts; i++)
			{
			  kriging = kriging + gVector[i]* neighbor[subNeighbor[i]].z_value;
			  sigma = sigma + lamda[i] * gVector[i];
			}

		      if (sigma > 0)
			{
			  sigma = sqrt(sigma);
			}

		      if (kriging >= 0)
			{
			  kvq	= Good;
			}
		      else
			{
			  kvq = OutofRange;
			}
		    }
		  delete [] gVector;
		  delete [] lamda;


		  // free BinCube memory
		  for (int i = 0; i < DistBin; i++)
		    {
		      delete [] BinCube[i];
		      BinCube[i] = NULL;
		    }
		  delete [] BinCube;
		  BinCube = NULL;


		  // free memory
		  for (int i = 0; i < xn; i++)
		    {
		      for (int j = 0; j< yn; j++)
			{
			  for (int k = 0; k < tn; k++)
			    {
			      neighborcube[i][j][k].clear();
			    }
			  delete [] neighborcube[i][j];
			  neighborcube[i][j] = NULL;
			}
		      delete [] neighborcube[i];
		      neighborcube[i] = NULL;
		    }
		  delete [] neighborcube;
		  neighborcube = NULL;
		}				
	    }
	}
    }
}

void VariogramCalculationBySpaceAndTime(			// calculate variogram by space and time function
					AvDistTimeSemi**		&BinCube,			// neighbors are stored in this two demission array
					int						DistBin,			// the number of distance bins
					int						TimeBin,			// the number of time bins
					double					d_interval,			// the interval on distance
					double					t_interval,			// the interval on time
					VariogramParameters&	parameters)			// variogram parameters (returned)
{
  double Hs = -99999;
  double nugs = -99999;
  double sills = -99999;
  ModelSelection ms = Init;
  double Ht = -99999;
  double nugt = -99999;
  double sillt = -99999;
  ModelSelection mt = Init;
  double kappa = -99999;

  double peakofsemi;
  int nn;
  double TempX1;
  double TempY1;
  double nugget_sphy;
  double sill_sphy;
  double nugget_mtey;
  double sill_mtey;
  double nugget_gausy;
  double sill_gausy;	
  double nugget_expy;
  double sill_expy;
  bool   sphy_flag;
  bool   expy_flag;
  bool   gausy_flag;
  bool   mtey_flag;
  double sse_sphy;
  double sse_expy;
  double sse_gausy;
  double sse_mtey;
  double sse = INFINITY;
  AvDistTimeSemi* Bin;

  // initial the first distance Bins
  Bin = new AvDistTimeSemi [DistBin];
  for (int i = 0; i < DistBin; i++)
    {
      Bin[i].AvDistance = (i + 1) * d_interval; // BinCube[i][0].AvDistance;
      Bin[i].AvSemi = BinCube[i][0].AvSemi;
      Bin[i].AvTimeDiff = (i + 1) * t_interval; // BinCube[i][0].AvTimeDiff;
      Bin[i].Pairs = BinCube[i][0].Pairs;
    }

  // calculate how many missing values in the first distance Bins
  nn = 0;
  for (int i = 0; i < DistBin; i++)
    {
      if (Bin[i].Pairs == 0)
	{
	  nn = nn + 1;
	}
    }

  // if the numbers of missing values in the first distance Bins are greater than 0.5*Bins, interpolate the missing values
  if (nn >= 0.5*DistBin)
    {
      for (int i = 0; i < DistBin; i++)
	{
	  if (Bin[i].Pairs == 0)
	    {
	      Bin[i].AvDistance = (i + 1) * d_interval;

	      Bin[i].AvSemi = InterpolateMissingValueBySpaceAndTime(		// find the nearest three bins and use this three bins' value to interpolate
								    BinCube,		// neighbors are stored in this two demission array
								    DistBin,		// the number of distance bins
								    TimeBin,		// the number of time bins
								    i,				// the index on distance Bin of missing value
								    0);				// the index on time Bin of missing value
				
	      Bin[i].AvTimeDiff = 0; // (i + 1) * t_interval;
	      Bin[i].Pairs = 1;
	    }		
	}		
    }

  // search the semi peak
  peakofsemi = Bin[0].AvSemi;
  for (int i = 1; i < DistBin; i++)
    {
      peakofsemi = peakofsemi > Bin[i].AvSemi ? peakofsemi : Bin[i].AvSemi;
    }
  peakofsemi = 0.8 * peakofsemi;

  // calculate Hs based on semi peak
  for (int i = 0; i < DistBin; i++)
    {
      if (Bin[i].AvSemi > peakofsemi)
	{
	  Hs = (i + 1) * d_interval; //Bin[i].AvDistance;
	  i = DistBin;
	}
    }

  // Spherical model
  // SphY[jj] = nugget + sill*((3*x)/(2*H) - 0.5*pow(x/H, 3)) if x <= Hs
  // SphY[jj] = nugget + sill if x > Hs
  double sigma_sqp_sqf = 0;
  double sigma_sqp_z = 0;
  double sigma_sqp_f = 0;
  double sigma_sqp_f_z = 0;
  double sigma_sqp = 0;
			
  for (int i = 0; i < DistBin; i++)
    {
      if (Bin[i].Pairs != 0)
	{
	  if (Bin[i].AvDistance <= Hs)
	    {
	      TempX1 = 1.5 * Bin[i].AvDistance / Hs - 0.5 * pow((Bin[i].AvDistance / Hs), 3);
	    }
	  else
	    {
	      TempX1 = 1;
	    }
	  sigma_sqp_sqf = sigma_sqp_sqf + Bin[i].Pairs * TempX1 * TempX1;
	  sigma_sqp_z = sigma_sqp_z + Bin[i].Pairs * Bin[i].AvSemi;
	  sigma_sqp_f = sigma_sqp_f + Bin[i].Pairs * TempX1;
	  sigma_sqp_f_z = sigma_sqp_f_z + Bin[i].Pairs * TempX1 * Bin[i].AvSemi;
	  sigma_sqp = sigma_sqp + Bin[i].Pairs;
	}		
    }
  if ((sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f) != 0)
    {
      nugget_sphy = (sigma_sqp_sqf * sigma_sqp_z - sigma_sqp_f * sigma_sqp_f_z) / (sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f);
      sill_sphy = (sigma_sqp * sigma_sqp_f_z - sigma_sqp_z * sigma_sqp_f) / (sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f);
    }
  else
    {
      nugget_sphy = 0;
      sill_sphy = 0;
    }
	
	
  // Exponential model
  // ExpY[jj] = nugget + sill*(1 - exp(-3*x/H));
  sigma_sqp_sqf = 0;
  sigma_sqp_z = 0;
  sigma_sqp_f = 0;
  sigma_sqp_f_z = 0;
  sigma_sqp = 0;

  for (int i = 0; i < DistBin; i++)
    {
      if (Bin[i].Pairs != 0)
	{
	  TempX1 = 1 - exp(-3 * Bin[i].AvDistance / Hs);

	  sigma_sqp_sqf = sigma_sqp_sqf + Bin[i].Pairs * TempX1 * TempX1;
	  sigma_sqp_z = sigma_sqp_z + Bin[i].Pairs * Bin[i].AvSemi;
	  sigma_sqp_f = sigma_sqp_f + Bin[i].Pairs * TempX1;
	  sigma_sqp_f_z = sigma_sqp_f_z + Bin[i].Pairs * TempX1 * Bin[i].AvSemi;
	  sigma_sqp = sigma_sqp + Bin[i].Pairs;
	}
    }
  if ((sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f) != 0)
    {
      nugget_expy = (sigma_sqp_sqf * sigma_sqp_z - sigma_sqp_f * sigma_sqp_f_z) / (sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f);
      sill_expy = (sigma_sqp * sigma_sqp_f_z - sigma_sqp_z * sigma_sqp_f) / (sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f);
    }
  else
    {
      nugget_expy = 0;
      sill_expy = 0;
    }
	

  // Gaussian model
  // GausY[jj] = nugget + sill*(1-exp(-3*(x/H)*(x/H)));
  sigma_sqp_sqf = 0;
  sigma_sqp_z = 0;
  sigma_sqp_f = 0;
  sigma_sqp_f_z = 0;
  sigma_sqp = 0;

  for (int i = 0; i < DistBin; i++)
    {
      if (Bin[i].Pairs != 0)
	{
	  TempX1 = 1 - exp(-3 * Bin[i].AvDistance * Bin[i].AvDistance / (Hs * Hs));

	  sigma_sqp_sqf = sigma_sqp_sqf + Bin[i].Pairs * TempX1 * TempX1;
	  sigma_sqp_z = sigma_sqp_z + Bin[i].Pairs * Bin[i].AvSemi;
	  sigma_sqp_f = sigma_sqp_f + Bin[i].Pairs * TempX1;
	  sigma_sqp_f_z = sigma_sqp_f_z + Bin[i].Pairs * TempX1 * Bin[i].AvSemi;
	  sigma_sqp = sigma_sqp + Bin[i].Pairs;
	}
    }
  if ((sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f) != 0)
    {
      nugget_gausy = (sigma_sqp_sqf * sigma_sqp_z - sigma_sqp_f * sigma_sqp_f_z) / (sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f);
      sill_gausy = (sigma_sqp * sigma_sqp_f_z - sigma_sqp_z * sigma_sqp_f) / (sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f);
    }
  else
    {
      nugget_gausy = 0;
      sill_gausy = 0;
    }
	

  // MteY model
  // MteY[jj] = nugget + sill*{1-[1+(4.5*x/H)]*exp[-4.5*x/H]};
  sigma_sqp_sqf = 0;
  sigma_sqp_z = 0;
  sigma_sqp_f = 0;
  sigma_sqp_f_z = 0;
  sigma_sqp = 0;

  for (int i = 0; i < DistBin; i++)
    {
      if (Bin[i].Pairs != 0)
	{
	  TempX1 = 1 - (1 + 4.5 * Bin[i].AvDistance / Hs) * exp(-4.5 * Bin[i].AvDistance / Hs);

	  sigma_sqp_sqf = sigma_sqp_sqf + Bin[i].Pairs * TempX1 * TempX1;
	  sigma_sqp_z = sigma_sqp_z + Bin[i].Pairs * Bin[i].AvSemi;
	  sigma_sqp_f = sigma_sqp_f + Bin[i].Pairs * TempX1;
	  sigma_sqp_f_z = sigma_sqp_f_z + Bin[i].Pairs * TempX1 * Bin[i].AvSemi;
	  sigma_sqp = sigma_sqp + Bin[i].Pairs;
	}
    }
  if ((sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f) != 0)
    {
      nugget_mtey = (sigma_sqp_sqf * sigma_sqp_z - sigma_sqp_f * sigma_sqp_f_z) / (sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f);
      sill_mtey = (sigma_sqp * sigma_sqp_f_z - sigma_sqp_z * sigma_sqp_f) / (sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f);
    }
  else
    {
      nugget_mtey = 0;
      sill_mtey = 0;
    }	

  sphy_flag = false;
  expy_flag = false;
  gausy_flag = false;
  mtey_flag = false;
  sse_sphy = 0;
  sse_expy = 0;
  sse_gausy = 0;
  sse_mtey = 0;

  if ((nugget_sphy > 0)&&(sill_sphy > 0))
    {
      sphy_flag = true;
      for (int i = 0; i < DistBin; i++)
	{
	  if (Bin[i].Pairs != 0)
	    {
	      if (Bin[i].AvDistance <= Hs)
		{
		  TempX1 = Bin[i].AvSemi - nugget_sphy - sill_sphy * (1.5 * Bin[i].AvDistance / Hs - 0.5 * pow((Bin[i].AvDistance / Hs), 3));
		}
	      else
		{
		  TempX1 = Bin[i].AvSemi - nugget_sphy - sill_sphy;
		}				
	      sse_sphy = sse_sphy + Bin[i].Pairs * TempX1 * TempX1;
	    }
	}
    }
  if ((nugget_expy > 0)&&(sill_expy > 0))
    {
      expy_flag = true;
      for (int i = 0; i < DistBin; i++)
	{
	  if (Bin[i].Pairs != 0)
	    {
	      TempX1 = Bin[i].AvSemi - nugget_expy - sill_expy * (1 - exp(-3 * Bin[i].AvDistance / Hs));
	      sse_expy = sse_expy + Bin[i].Pairs * TempX1 * TempX1;
	    }
	}
    }
  if ((nugget_gausy > 0)&&(sill_gausy > 0))
    {
      gausy_flag = true;
      for (int i = 0; i < DistBin; i++)
	{
	  if (Bin[i].Pairs != 0)
	    {
	      TempX1 = Bin[i].AvSemi - nugget_gausy - sill_gausy * (1 - exp(-3 * Bin[i].AvDistance * Bin[i].AvDistance / (Hs * Hs)));
	      sse_gausy = sse_gausy + Bin[i].Pairs * TempX1 * TempX1;
	    }
	}
    }
  if ((nugget_mtey > 0)&&(sill_mtey > 0))
    {
      mtey_flag = true;
      for (int i = 0; i < DistBin; i++)
	{
	  if (Bin[i].Pairs != 0)
	    {
	      TempX1 = Bin[i].AvSemi - nugget_mtey - sill_mtey * (1 - (1 + 4.5 * Bin[i].AvDistance / Hs) * exp(-4.5 * Bin[i].AvDistance / Hs));
	      sse_mtey = sse_mtey + Bin[i].Pairs * TempX1 * TempX1;
	    }
	}
    }

  if ((sphy_flag == false)&&(expy_flag == false)&&(gausy_flag == false)&&(mtey_flag == false))
    {
      nugs = 0;
      sills = 0;
      double sigma_sqp_f_z = 0;
      double sigma_sqp_sqf = 0;
		
      for (int i = 0; i < DistBin; i++)
	{
	  if (Bin[i].Pairs != 0)
	    {
	      TempX1 = (1 - exp(-3 * Bin[i].AvDistance / Hs));
	      TempY1 = Bin[i].AvSemi;

	      sigma_sqp_f_z = sigma_sqp_f_z + Bin[i].Pairs * TempX1 * TempY1;
	      sigma_sqp_sqf = sigma_sqp_sqf + Bin[i].Pairs * TempX1 * TempX1;
	    }
	}
      if (sigma_sqp_sqf != 0)
	{
	  sills = sigma_sqp_f_z / sigma_sqp_sqf;
	}
      else
	{
	  sills = 0;
	}
      ms = Exponential;
    }
  else
    {
      if (sphy_flag == true)
	{
	  nugs = nugget_sphy;
	  sills = sill_sphy;
	  ms = Spherical;
	  sse = sse_sphy;
	}
      else if (gausy_flag == true)
	{
	  nugs = nugget_gausy;
	  sills = sill_gausy;
	  ms = Gaussian;
	  sse = sse_gausy;
	}
      else if (expy_flag == true)
	{
	  nugs = nugget_expy;
	  sills = sill_expy;
	  ms = Exponential;
	  sse = sse_expy;
	}
      else if (mtey_flag == true)
	{
	  nugs = nugget_mtey;
	  sills = sill_mtey;
	  ms = Matern;
	  sse = sse_mtey;
	}

      if ((sphy_flag == true)&&(sse_sphy < sse))
	{
	  nugs = nugget_sphy;
	  sills = sill_sphy;
	  ms = Spherical;
	  sse = sse_sphy;
	}
      if ((gausy_flag == true)&&(sse_gausy < sse))
	{
	  nugs = nugget_gausy;
	  sills = sill_gausy;
	  ms = Gaussian;
	  sse = sse_gausy;
	}
      if ((expy_flag == true)&&(sse_expy < sse))
	{
	  nugs = nugget_expy;
	  sills = sill_expy;
	  ms = Exponential;
	  sse = sse_expy;
	}
      if ((mtey_flag == true)&&(sse_mtey < sse))
	{
	  nugs = nugget_mtey;
	  sills = sill_mtey;
	  ms = Matern;
	  sse = sse_mtey;
	}
    }

  // free Bin memory
  delete [] Bin;
  Bin = NULL;

  // estimate a variogram model from BinCube[0][0] to BinCube[0][DISTBIN]
  // return mt  = best fitted model
  //		  nugt = nugget effect
  //        sillt   = partial sill
  //        Ht      = range paramter

  // initial the first timr Bins
  Bin = new AvDistTimeSemi [TimeBin];
  for (int i = 0; i < TimeBin; i++)
    {
      Bin[i].AvDistance = (i + 1) * d_interval; // BinCube[0][i].AvDistance;
      Bin[i].AvSemi = BinCube[0][i].AvSemi;
      Bin[i].AvTimeDiff = (i + 1) * t_interval; // BinCube[0][i].AvTimeDiff;
      Bin[i].Pairs = BinCube[0][i].Pairs;
    }

  // calculate how many missing values in the first time Bins
  nn = 0;
  for (int i = 0; i < TimeBin; i++)
    {
      if (Bin[i].Pairs == 0)
	{
	  nn = nn + 1;
	}
    }

  // if the numbers of missing values in the first time Bins are greater than 0.5*Bins, interpolate the missing values
  if (nn >= 0.5*TimeBin)
    {
      for (int i = 0; i < TimeBin; i++)
	{
	  if (Bin[i].Pairs == 0)
	    {
	      Bin[i].AvDistance = 0; //(i + 1) * d_interval;

	      Bin[i].AvSemi = InterpolateMissingValueBySpaceAndTime(		// find the nearest three bins and use this three bins' value to interpolate
								    BinCube,		// neighbors are stored in this two demission array
								    DistBin,		// the number of distance bins
								    TimeBin,		// the number of time bins
								    0,				// the index on distance Bin of missing value
								    i);				// the index on time Bin of missing value
				
	      Bin[i].AvTimeDiff = (i + 1) * t_interval;
	      Bin[i].Pairs = 1;
	    }		
	}		
    }

  // search the semi peak
  peakofsemi = Bin[0].AvSemi;
  for (int i = 1; i < TimeBin; i++)
    {
      peakofsemi = peakofsemi > Bin[i].AvSemi ? peakofsemi : Bin[i].AvSemi;
    }
  peakofsemi = 0.8 * peakofsemi;

  // calculate Ht based on semi peak
  for (int i = 0; i < TimeBin; i++)
    {
      if (Bin[i].AvSemi > peakofsemi)
	{
	  Ht = (i + 1) * t_interval; // Bin[i].AvTimeDiff;
	  i = TimeBin;
	}
    }	

  // Spherical model
  // SphY[jj] = nugget + sill*((3*x)/(2*H) - 0.5*pow(x/H, 3)) if x <= H
  // SphY[jj] = nugget + sill if x > H
  sigma_sqp_sqf = 0;
  sigma_sqp_z = 0;
  sigma_sqp_f = 0;
  sigma_sqp_f_z = 0;
  sigma_sqp = 0;

  for (int i = 0; i < TimeBin; i++)
    {
      if (Bin[i].Pairs != 0)
	{
	  if (Bin[i].AvTimeDiff <= Ht)
	    {
	      TempX1 = 1.5 * Bin[i].AvTimeDiff / Ht - 0.5 * pow((Bin[i].AvTimeDiff / Ht), 3);
	    }
	  else
	    {
	      TempX1 = 1;
	    }
	  sigma_sqp_sqf = sigma_sqp_sqf + Bin[i].Pairs * TempX1 * TempX1;
	  sigma_sqp_z = sigma_sqp_z + Bin[i].Pairs * Bin[i].AvSemi;
	  sigma_sqp_f = sigma_sqp_f + Bin[i].Pairs * TempX1;
	  sigma_sqp_f_z = sigma_sqp_f_z + Bin[i].Pairs * TempX1 * Bin[i].AvSemi;
	  sigma_sqp = sigma_sqp + Bin[i].Pairs;
	}
    }
  if ((sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f) != 0)
    {
      nugget_sphy = (sigma_sqp_sqf * sigma_sqp_z - sigma_sqp_f * sigma_sqp_f_z) / (sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f);
      sill_sphy = (sigma_sqp * sigma_sqp_f_z - sigma_sqp_z * sigma_sqp_f) / (sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f);
    }
  else
    {
      nugget_sphy = 0;
      sill_sphy = 0;
    }

  // Exponential model
  // ExpY[jj] = nugget + sill*(1 - exp(-3*x/H));
  sigma_sqp_sqf = 0;
  sigma_sqp_z = 0;
  sigma_sqp_f = 0;
  sigma_sqp_f_z = 0;
  sigma_sqp = 0;

  for (int i = 0; i < TimeBin; i++)
    {
      if (Bin[i].Pairs != 0)
	{
	  TempX1 = 1 - exp(-3 * Bin[i].AvTimeDiff / Ht);

	  sigma_sqp_sqf = sigma_sqp_sqf + Bin[i].Pairs * TempX1 * TempX1;
	  sigma_sqp_z = sigma_sqp_z + Bin[i].Pairs * Bin[i].AvSemi;
	  sigma_sqp_f = sigma_sqp_f + Bin[i].Pairs * TempX1;
	  sigma_sqp_f_z = sigma_sqp_f_z + Bin[i].Pairs * TempX1 * Bin[i].AvSemi;
	  sigma_sqp = sigma_sqp + Bin[i].Pairs;
	}
    }
  if ((sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f) != 0)
    {
      nugget_expy = (sigma_sqp_sqf * sigma_sqp_z - sigma_sqp_f * sigma_sqp_f_z) / (sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f);
      sill_expy = (sigma_sqp * sigma_sqp_f_z - sigma_sqp_z * sigma_sqp_f) / (sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f);
    }
  else
    {
      nugget_expy = 0;
      sill_expy = 0;
    }

  // Gaussian model
  // GausY[jj] = nugget + sill*(1-exp(-3*(x/H)*(x/H)));
  sigma_sqp_sqf = 0;
  sigma_sqp_z = 0;
  sigma_sqp_f = 0;
  sigma_sqp_f_z = 0;
  sigma_sqp = 0;

  for (int i = 0; i < TimeBin; i++)
    {
      if (Bin[i].Pairs != 0)
	{
	  TempX1 = 1 - exp(-3 * Bin[i].AvTimeDiff * Bin[i].AvTimeDiff / (Ht * Ht));

	  sigma_sqp_sqf = sigma_sqp_sqf + Bin[i].Pairs * TempX1 * TempX1;
	  sigma_sqp_z = sigma_sqp_z + Bin[i].Pairs * Bin[i].AvSemi;
	  sigma_sqp_f = sigma_sqp_f + Bin[i].Pairs * TempX1;
	  sigma_sqp_f_z = sigma_sqp_f_z + Bin[i].Pairs * TempX1 * Bin[i].AvSemi;
	  sigma_sqp = sigma_sqp + Bin[i].Pairs;
	}
    }
  if ((sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f) != 0)
    {
      nugget_gausy = (sigma_sqp_sqf * sigma_sqp_z - sigma_sqp_f * sigma_sqp_f_z) / (sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f);
      sill_gausy = (sigma_sqp * sigma_sqp_f_z - sigma_sqp_z * sigma_sqp_f) / (sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f);
    }
  else
    {
      nugget_gausy = 0;
      sill_gausy = 0;
    }

  // MteY model
  // MteY[jj] = nugget + sill*[1-(1+4.5*x/H)*exp(-4.5*x/H)];
  sigma_sqp_sqf = 0;
  sigma_sqp_z = 0;
  sigma_sqp_f = 0;
  sigma_sqp_f_z = 0;
  sigma_sqp = 0;

  for (int i = 0; i < TimeBin; i++)
    {
      if (Bin[i].Pairs != 0)
	{
	  TempX1 = 1 - (1 + 4.5 * Bin[i].AvTimeDiff / Ht) * exp(-4.5 * Bin[i].AvTimeDiff / Ht);

	  sigma_sqp_sqf = sigma_sqp_sqf + Bin[i].Pairs * TempX1 * TempX1;
	  sigma_sqp_z = sigma_sqp_z + Bin[i].Pairs * Bin[i].AvSemi;
	  sigma_sqp_f = sigma_sqp_f + Bin[i].Pairs * TempX1;
	  sigma_sqp_f_z = sigma_sqp_f_z + Bin[i].Pairs * TempX1 * Bin[i].AvSemi;
	  sigma_sqp = sigma_sqp + Bin[i].Pairs;
	}
    }
  if ((sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f) != 0)
    {
      nugget_mtey = (sigma_sqp_sqf * sigma_sqp_z - sigma_sqp_f * sigma_sqp_f_z) / (sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f);
      sill_mtey = (sigma_sqp * sigma_sqp_f_z - sigma_sqp_z * sigma_sqp_f) / (sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f);
    }
  else
    {
      nugget_mtey = 0;
      sill_mtey = 0;
    }

  sphy_flag = false;
  expy_flag = false;
  gausy_flag = false;
  mtey_flag = false;
  sse_sphy = 0;
  sse_expy = 0;
  sse_gausy = 0;
  sse_mtey = 0;

  if ((nugget_sphy > 0)&&(sill_sphy > 0))
    {
      sphy_flag = true;
      for (int i = 0; i < TimeBin; i++)
	{
	  if (Bin[i].Pairs != 0)
	    {
	      if (Bin[i].AvTimeDiff <= Ht)
		{
		  TempX1 = Bin[i].AvSemi - nugget_sphy - sill_sphy * (1.5 * Bin[i].AvTimeDiff / Ht - 0.5 * pow((Bin[i].AvTimeDiff / Ht), 3));
		}
	      else
		{
		  TempX1 = Bin[i].AvSemi - nugget_sphy - sill_sphy;
		}				
	      sse_sphy = sse_sphy + Bin[i].Pairs * TempX1 * TempX1;
	    }
	}
    }
  if ((nugget_expy > 0)&&(sill_expy > 0))
    {
      expy_flag = true;
      for (int i = 0; i < TimeBin; i++)
	{
	  if (Bin[i].Pairs != 0)
	    {
	      TempX1 = Bin[i].AvSemi - nugget_expy - sill_expy * (1 - exp(-3 * Bin[i].AvTimeDiff / Ht));
	      sse_expy = sse_expy + Bin[i].Pairs * TempX1 * TempX1;
	    }
	}
    }
  if ((nugget_gausy > 0)&&(sill_gausy > 0))
    {
      gausy_flag = true;
      for (int i = 0; i < TimeBin; i++)
	{
	  if (Bin[i].Pairs != 0)
	    {
	      TempX1 = Bin[i].AvSemi - nugget_gausy - sill_gausy * (1 - exp(-3 * Bin[i].AvTimeDiff * Bin[i].AvTimeDiff / (Ht * Ht)));
	      sse_gausy = sse_gausy + Bin[i].Pairs * TempX1 * TempX1;
	    }
	}
    }
  if ((nugget_mtey > 0)&&(sill_mtey > 0))
    {
      mtey_flag = true;
      for (int i = 0; i < TimeBin; i++)
	{
	  if (Bin[i].Pairs != 0)
	    {
	      TempX1 = Bin[i].AvSemi - nugget_mtey - sill_mtey * (1 - (1 + 4.5 * Bin[i].AvTimeDiff / Ht) * exp(-4.5 * Bin[i].AvTimeDiff / Ht));
	      sse_mtey = sse_mtey + Bin[i].Pairs * TempX1 * TempX1;
	    }
	}
    }

  if ((sphy_flag == false)&&(expy_flag == false)&&(gausy_flag == false)&&(mtey_flag == false))
    {
      nugt = 0;
      sillt = 0;
      sigma_sqp_f_z = 0;
      sigma_sqp_sqf = 0;

      for (int i = 0; i < TimeBin; i++)
	{
	  if (Bin[i].Pairs != 0)
	    {
	      TempX1 = 1 - exp(-3 * Bin[i].AvTimeDiff / Ht);
	      TempY1 = Bin[i].AvSemi;

	      sigma_sqp_f_z = sigma_sqp_f_z + Bin[i].Pairs * TempX1 * TempY1;
	      sigma_sqp_sqf = sigma_sqp_sqf + Bin[i].Pairs * TempX1 * TempX1;
	    }
	}
      if (sigma_sqp_sqf != 0)
	{
	  sillt = sigma_sqp_f_z / sigma_sqp_sqf;
	}
      else
	{
	  sillt = 0;
	}
      mt = Exponential;
    }
  else
    {
      if (sphy_flag == true)
	{
	  nugt = nugget_sphy;
	  sillt = sill_sphy;
	  mt = Spherical;
	  sse = sse_sphy;
	}
      else if (gausy_flag == true)
	{
	  nugt = nugget_gausy;
	  sillt = sill_gausy;
	  mt = Gaussian;
	  sse = sse_gausy;
	}
      else if (expy_flag == true)
	{
	  nugt = nugget_expy;
	  sillt = sill_expy;
	  mt = Exponential;
	  sse = sse_expy;
	}
      else if (mtey_flag == true)
	{
	  nugt = nugget_mtey;
	  sillt = sill_mtey;
	  mt = Matern;
	  sse = sse_mtey;
	}

      if ((sphy_flag == true)&&(sse_sphy < sse))
	{
	  nugt = nugget_sphy;
	  sillt = sill_sphy;
	  mt = Spherical;
	  sse = sse_sphy;
	}
      if ((gausy_flag == true)&&(sse_gausy < sse))
	{
	  nugt = nugget_gausy;
	  sillt = sill_gausy;
	  mt = Gaussian;
	  sse = sse_gausy;
	}
      if ((expy_flag == true)&&(sse_expy < sse))
	{
	  nugt = nugget_expy;
	  sillt = sill_expy;
	  mt = Exponential;
	  sse = sse_expy;
	}
      if ((mtey_flag == true)&&(sse_mtey < sse))
	{
	  nugt = nugget_mtey;
	  sillt = sill_mtey;
	  mt = Matern;
	  sse = sse_mtey;
	}
    }

  // free Bin memory
  delete [] Bin;
  Bin = NULL;

  parameters.Hs = Hs;
  parameters.nuggets = nugs;
  parameters.sills = sills;
  parameters.ms = ms;
  parameters.Ht = Ht;
  parameters.nuggett = nugt;
  parameters.sillt = sillt;
  parameters.mt = mt;
  parameters.kappa = kappa;
}


void	Kappa_Calculate(							// calculate the parameter kappa and adjust nugs, sills, nugt, sillt
			NeighborVector&			neighbor,			// neighbors are stored in this vector
			vector<int>&			sample,				// sample neighbors are stored in this vector
			VariogramParameters&	parameters)			// variogram parameters (returned)
{
  // estimate the interaction parameter kappa
  // return kappa and adjusted nugs,sills,nugt,sillt

  double kappa = -99999;
  int npts; 
  double sumZ;
  double zbar;
  double Cst;
  double Cs;
  double Ct;
  double d1,d2,d;
  double rs,rt;
  double nugs = parameters.nuggets;
  double sills = parameters.sills;
  double nugt = parameters.nuggett;
  double sillt = parameters.sillt;


  npts = (int)sample.size(); 
  sumZ = 0;
  zbar = 0;
  if (npts > 0)
    {
      for (int i = 0; i < npts; i++)
	{
	  sumZ = sumZ + neighbor[sample[i]].z_value;
	}
      zbar = sumZ / npts;
    }
	

  Cst = 0;

  for (int i = 0; i < npts; i++)
    {
      Cst = Cst + (neighbor[sample[i]].z_value - zbar) * (neighbor[sample[i]].z_value - zbar);
    }

  Cst = Cst / (npts - 1);
  Cs = nugs + sills;
  Ct = nugt + sillt;

  d1 = Cst - Cs;
  d2 = Cst - Ct;
  d = Cst - Cs - Ct;

  if ((d1 < 0)||(d2 < 0))
    {
      Cs = min(Cst,Cs);
      Ct = min(Cst,Ct);
    } 

  if (d >= 0)
    {
      Cs = Cs + min(d1, 0.51 * (d + 0.001));
      Ct = Ct + min(d2, 0.51 * (d + 0.001));
    }

  rs = max(sills / (nugs + sills), 0.001);
  rt = max(sillt / (nugt + sillt), 0.001);

  sills  = Cs * rs;
  nugs = Cs * (1 - rs);
  sillt  = Ct * rt;
  nugt = Ct * (1 - rt);
  kappa = (Cs + Ct - Cst) / (Cs * Ct);
  // estimate the interaction parameter kappa
  // return kappa and adjusted nugs,sills,nugt,sillt

  parameters.kappa   = kappa;
  parameters.nuggets = nugs;
  parameters.sills   = sills;
  parameters.nuggett = nugt;
  parameters.sillt   = sillt;
}

double  InterpolateMissingValueBySpaceAndTime(		// find the nearest three bins and use this three bins' value to interpolate
					      AvDistTimeSemi**		&BinCube,			// neighbors are stored in this two demission array
					      int						DistBin,			// the number of distance bins
					      int						TimeBin,			// the number of time bins
					      int						d,					// the index on distance Bin of missing value
					      int						t)					// the index on time Bin of missing value
{
  double semi = 0;
  vector<double> BinCubeSemi;
  vector<int>    BinCubePairs;
  vector<double> Dist;
  double		   dx;
  double		   dy;
  double         dd;
  int n = 0;
  double smallest;
  double Temp1 = 0;
  double Temp2 = 0;

  BinCubeSemi.clear();
  BinCubePairs.clear();
  Dist.clear();

  for (int i = 0; i < DistBin; i++)
    {
      for (int j = 0; j < TimeBin; j++)
	{
	  if (BinCube[i][j].Pairs != 0)
	    {
	      dx = i - d;
	      dy = j - t;
	      dd = sqrt(dx*dx + dy*dy);

	      BinCubeSemi.push_back(BinCube[i][j].AvSemi);
	      BinCubePairs.push_back(BinCube[i][j].Pairs);
	      Dist.push_back(dd);
	    }
	}
    }

  while (n < 3)
    {
      // search the smallest value
      smallest = Dist[0];
      for (int i = 1; i < (int)Dist.size(); i++)
	{
	  smallest = smallest < Dist[i] ? smallest : Dist[i];
	}

      // calculate smallest value
      for (int i = 0; i < (int)Dist.size(); i++)
	{
	  if (Dist[i] <= smallest)
	    {
	      Temp1 = Temp1 + BinCubeSemi[i] * BinCubePairs[i] / Dist[i];
	      Temp2 = Temp2 + BinCubePairs[i] / Dist[i];
	      n = n + 1;
	      Dist[i] = 99999;
	    }
	}
    }

  if (Temp2 != 0)
    {
      semi = Temp1 / Temp2;
    }	

  return semi;
}



//TODO: also return Hs
//05/01/2013: also return variogram parameters
void Kriging_By_Space(	// calculate kriging by space
		      //CMPKrigingPoint	queryPt,	// query point,			has query point location information
		      NeighborVector&	neighbor,	// neighbors are stored in this vector
		      int		DistBin,	// the number of distance bins
		      ModelSelection    m,	// calculation model
		      double		alpha,	// when model selection is IDW, specify one parameter for distance
		      double&		kriging,	// kriging value (returned)
		      double&		sigma,	// kriging variance (returned)
		      double&           Hs, // 
			  VariogramParameters&	parameters, // variogram parameters
		      KrigingValueQuality&	kvq)	// kriging value quality (returned)
{
  kriging = -99999;
  sigma = -99999;
  kvq = Empty;

  double			TempX1,TempY1,TempX2, TempY2;
  double			distance;
  double			tempAlpha;
  //double			krigingbar;
  double			xmax,xmin,ymax,ymin;
  vector<bool>	sampling;
  double			dist_max,dbar;
  vector<int>   **neighborcell;
  double			x_interval, y_interval;
  int				xn,yn;
  int				xi,yi;
  int				npts;
  double			TempZ1;
  vector<int>		sampleNeighbor;
  int				sam_npts;
  //int				TempR; // random number
  AvDistSemi     *Bins;
  double			d_interval;
  int				di;
  double			semi;
  //VariogramParameters	parameters;
  vector<int>		subNeighbor;
  
  //ofstream dout;
  //dout.open("debug_0");	

  npts = (int)neighbor.size();

  //dout <<"neighbor input" <<endl;
  //for(int i=0; i <npts; i++){
  //  int j = neighbor.get_index(i);
	//dout << "x=" << neighbor[j].x_coord<<" y= "<< neighbor[j].y_coord << "distance= " << neighbor[j].distance << endl;
  //}
  //if ((queryPt == NULL )||(DistBin <= 3))
  // cout <<  " DistBin " << DistBin << endl;
  if (DistBin <= 3)  
    {
      kriging	= -99999;
      sigma	= -99999;
      kvq		= WrongParameter;
    }
  else if (npts == 0)
    {
      kriging	= -99999;
      sigma	= -99999;
      kvq		= LessNeighborsInSpace;
    }
  else
    {
      if (m == IDW) // calculate IDW
	{
	  TempX2 = 0;
	  TempY2 = 0;
	  for (int i = 0; i < npts; i++)
	    {
		  int j = neighbor.get_index(i);
	      TempZ1   = neighbor[j].z_value;
	      distance = neighbor[j].distance;

	      if (distance == 0)
		{
		  tempAlpha = 1;
		}
	      else if (distance != 0)
		{
		  tempAlpha = pow((double)distance,alpha*(-1));
		}

	      TempX2   = TempX2 + TempZ1 * tempAlpha;
	      TempY2   = TempY2 + tempAlpha;
	    }
	  kriging	= -99999;
	  sigma	= -99999;
	  kvq		= IDWValue;
	  if (TempY2 != 0)
	    {
	      kriging	= TempX2 / TempY2;
	    }
	}
      else // calculate kriging by space
	{
	  if (npts <= MINNEIGHBOR)
	    {
	      //Rprintf("cmp_kriging.cpp:2045 Less neighbor in space\n");
		  kvq	= LessNeighborsInSpace;
	      kriging	= -99999;
	      sigma = -99999;
	    }
	  else // npts > MINNEIGHBOR
	    {
	      // search the range of x, y and time
	      xmin = neighbor[neighbor.get_index(0)].x_coord;
	      xmax = neighbor[neighbor.get_index(0)].x_coord;
	      ymax = neighbor[neighbor.get_index(0)].y_coord;
	      ymin = neighbor[neighbor.get_index(0)].y_coord;
				
	      sampling.clear();
	      sampling.push_back(false);

	      for (int i = 1; i < npts; i++)
		{
		  int j = neighbor.get_index(i);

		  if (xmin > neighbor[j].x_coord)
		    {
		      xmin = neighbor[j].x_coord;
		    }
		  if (xmax < neighbor[j].x_coord)
		    {
		      xmax = neighbor[j].x_coord;
		    }
		  if (ymin > neighbor[j].y_coord)
		    {
		      ymin = neighbor[j].y_coord;
		    }
		  if (ymax < neighbor[j].y_coord)
		    {
		      ymax = neighbor[j].y_coord;
		    }
		  sampling.push_back(false);
		}

	      dist_max = sqrt((xmax - xmin)*(xmax - xmin)+(ymax - ymin)*(ymax - ymin));
	      // cout << "dist_max first " << dist_max << endl;
		  // Rprintf("dist_max first %f\n",dist_max);
		  //dout<< xmax <<","<<xmin<<","<<ymax<<","<<ymin<<","<<dist_max<< endl;
	      dbar = 0;
	      for (int i = 0; i < npts; i++)
		{
		  int j = neighbor.get_index(i);
		  dbar = dbar + neighbor[j].distance;
		}
	      dbar = dbar / npts;

	      if (dbar > 2 * dist_max)
		{
		  //Rprintf("cmp_kriging.cpp:2098 Sparse neighborhood.\n");
		  kriging = -99999;
		  sigma = -99999;
		  kvq = LessNeighbors;
		}
	      else
		{
		  if (npts > MAXNEIGHBOR) // the neighbors are greater than maxium, randomly select sampling neighbors to reduce calculation
		    {
			  //Rprintf("cmp_kriging.cpp:2106 %d too large. Subsampling starts\n",npts);
			  //dout<<"here 1 npts="<<npts<<endl;
		      // setup random value enable 
		      // srand((unsigned)time(NULL));

		      xn = NEIGHBORCELL;
		      yn = NEIGHBORCELL;

		      neighborcell = new vector<int>* [xn];
		      for (int i = 0; i < xn; i++)
			{
			  neighborcell[i] = new vector<int> [yn];
			  for (int j = 0; j < yn; j++)
			    {
			      neighborcell[i][j].clear();
			    }
			}

		      // insert the neighbors into neighbor cubes
		      x_interval = abs(xmax - xmin) / xn;
		      y_interval = abs(ymax - ymin) / yn;
		      for (int i = 0; i < npts; i++)
			{
			  int j = neighbor.get_index(i);
			  xi = (int)floor((neighbor[j].x_coord - xmin) / x_interval);
			  if (xi == xn)
			    {
			      xi = xi - 1;
			    }

			  yi = (int)floor((neighbor[j].y_coord - ymin) / y_interval);
			  if (yi == yn)
			    {
			      yi = yi - 1;
			    }

			  if ((xi < xn)&&(yi < yn))
			    {
				  //dout << "xi= "<< xi << "," << " yi=" << yi << " index=" <<j<< endl;
			      neighborcell[xi][yi].push_back(j);
			  	  //Rprintf("cmp_kriging.cpp:2147 Insert neighbor %d into cell (%d,%d)\n",j,xi,yi);

			    }					
			}
		      // all neighbors are inserted into neighbor cubes
			  //dout<<"here 2"<<endl;

		      // randomly select subneighbors from neighbor cubes
		      sampleNeighbor.clear();
		      for (int i = 0; i < xn; i++)
			{
			  for (int j = 0; j < yn; j++)
			    {
			      if (neighborcell[i][j].size() > 0)
				{
				  // calculate randomly sampling number for this neighbor cube
				  sam_npts = (int)(MAXNEIGHBOR * neighborcell[i][j].size() / npts);
				  //dout << "sam=" << sam_npts <<" size=" << neighborcell[i][j].size() <<endl;

				  // if sampling number is greater than number in cube, output all neighbors in the cube
				  if (sam_npts >= (int)neighborcell[i][j].size())
				    {
				      for (int k = 0; k < (int)neighborcell[i][j].size(); k++)
					{
					  sampleNeighbor.push_back(neighborcell[i][j][k]);
					  //cout << " add1 " << neighborcell[i][j][k];
					  //dout << " add1 " << neighborcell[i][j][k];					  

				    }
				}
			      else
				{
				  // output one sampling neighbor at least
				  if (sam_npts == 0)
				    {
				      sam_npts = 1;
				    }

				  // output sam_npts neighbors from cube to subneighbors
				  //for (int k = 0; k < sam_npts; k++)
				   // {
				   //   TempR = rand() % (int)neighborcell[i][j].size();
				   //   while (sampling[neighborcell[i][j][TempR]] == true)
					//{
					//  TempR = rand() % (int)neighborcell[i][j].size();
					//}
				    //  sampling[neighborcell[i][j][TempR]] = true;
				    //  sampleNeighbor.push_back(neighborcell[i][j][TempR]);
				    //  //cout << " add2 " << neighborcell[i][j][k];	
				    //  dout << " add2 " << neighborcell[i][j][TempR];
				    //}
					// Algorithm Due to 
					// http://www.cs.ucr.edu/~ciardo/teaching/CS177/section6.5.pdf;
					// c 2006 Pearson Ed., Inc. 0-13-142917-5
					// int l = 0;
					// int TempK;
					// for (unsigned int k = 0; k < neighborcell[i][j].size(); k++)
					// {
						// l = l + 1;
						// TempR = rand() % (neighborcell[i][j].size() - k) + k;
						// TempK = neighborcell[i][j][TempR];
						// neighborcell[i][j][TempR] = neighborcell[i][j][k];
						// neighborcell[i][j][k] = TempK;
						// //dout << " add2 " << TempK;
						// sampleNeighbor.push_back(TempK);
						// if (l > sam_npts){
							// break;
						// }
					// }
					for (int k = 0; k< sam_npts;k++)
					{
						sampleNeighbor.push_back(neighborcell[i][j][k]);
					}
				}								
			    }
				//Rprintf("cmp_kriging.cpp:2218 Sample %d neighbors into cell (%d,%d)\n",sam_npts,i+1,j+1);
			}
		    }
		  //dout<<"here 3"<<endl;

		  // free memory
		  for (int i = 0; i < xn; i++)
		    {
		      for (int j = 0; j< yn; j++)
			{
			  neighborcell[i][j].clear();
			}
		      delete [] neighborcell[i];
		      neighborcell[i] = NULL;
		    }
		  delete [] neighborcell;
		  neighborcell = NULL;
		}
	      else // the neighbors are less than maxium value, use all neighbors to calculate kriging value
		{
		  sampleNeighbor.clear();
		  for (int i = 0; i < (int)neighbor.size(); i++)
		    {
		      sampleNeighbor.push_back(neighbor.get_index(i));
		      // cout << " add3 " << i;			      
		    }
			//Rprintf("cmp_kriging.cpp:2244 Use %d neighbors without sampling\n",(int)neighbor.size());

		}

	      // initial Bins
	      Bins = new AvDistSemi [DistBin];
	      for (int i = 0; i < DistBin; i++)
		{
		  Bins[i].AvDistance = 0;
		  Bins[i].AvSemi = 0;
		  Bins[i].Pairs = 0;
		}

	      // cout << "sampleNeighbor size " << sampleNeighbor.size() << endl;
	      // cout << "sampleNeighbor 0 " << sampleNeighbor[0] << endl;
	      //dout << "sampleNeighbor size " << sampleNeighbor.size() << endl;
	      //dout << "sampleNeighbor 0 " << sampleNeighbor[0] << endl;		
		  
		  //dout<<"here 4"<<endl;
	  
		  //Rprintf("sample neighbor size %d \n",sampleNeighbor.size());
	      dist_max = neighbor[sampleNeighbor[0]].distance;
	      // cout << "dist_max 2nd " << dist_max << endl;
		  
	      // cout << "dist max " <<  dist_max << endl;
	      // cout << "dist max " <<  dist_max << endl;		  		  
	      for (int i = 1; i < (int)sampleNeighbor.size(); i++)
		{
		  // cout << "dist " << sampleNeighbor[i] << " " <<  neighbor[sampleNeighbor[i]].distance << endl;		  
		  //dout << "dist " << sampleNeighbor[i] << " " <<  neighbor[sampleNeighbor[i]].distance << endl;		  
		  dist_max = dist_max > neighbor[sampleNeighbor[i]].distance ? dist_max : neighbor[sampleNeighbor[i]].distance;
		}
	      // cout << "dist max " <<  dist_max << endl;
		  
	      d_interval = (double)dist_max / DistBin;
	      // cout << "dist bin " <<  DistBin << endl;
		  //dout<<"dist max " <<  dist_max << "dist bin " <<  DistBin <<endl;

	      for (int i = 0; i < (int)sampleNeighbor.size(); i++)
		{
		  for (int j = 0; j < i; j++)
		    {
		      distance = sqrt(pow((neighbor[sampleNeighbor[i]].x_coord - neighbor[sampleNeighbor[j]].x_coord), 2) + pow((neighbor[sampleNeighbor[i]].y_coord - neighbor[sampleNeighbor[j]].y_coord), 2));
							
		      semi = 0.5 * pow((neighbor[sampleNeighbor[i]].z_value - neighbor[sampleNeighbor[j]].z_value), 2);
		      // cout << "i " << i << " semi " << semi << endl;

		      // cout << distance << endl;
		      // cout << d_interval << endl;
		      //dout << "i " << i << " semi " << semi << endl;
		      //dout << distance << endl;
		      //dout << d_interval << endl;

		      if ( (distance/d_interval) < 1)
			di=0;
		      else
			di = (int)floor(distance / d_interval);
			  
		      // cout << "di " << di << endl;		
		      //dout << "di " << di << endl;			  
			  
		      if (di == DistBin)
			{
			  di = di - 1;
			}


		      if (di < DistBin)
			{
			  Bins[di].AvDistance = Bins[di].AvDistance + distance;
			  Bins[di].AvSemi = Bins[di].AvSemi + semi;
			  Bins[di].Pairs = Bins[di].Pairs + 1;
			}
		    }
		}
	      for (int i = 0; i < DistBin; i++)
		{
		  if (Bins[i].Pairs != 0)
		    {
		      Bins[i].AvDistance = Bins[i].AvDistance / Bins[i].Pairs;
		      Bins[i].AvSemi = Bins[i].AvSemi / Bins[i].Pairs;
		    }
			//Rprintf("cmp_kriging.cpp:2326 Variogram = %f npair= %d lag=%f\n",Bins[i].AvSemi,Bins[i].Pairs,Bins[i].AvDistance);

		}
		//dout<<"here 5"<<endl;

	      // calculate variogram
	      Variogram_By_Space_Calculation(			// calculate variogram by space function
					     Bins,	// neighbors are stored in this one demission array
					     DistBin,	// the number of distance bins
					     d_interval,	// the interval on distance
					     parameters);	// variogram parameters (returned)
	      // calculate variogram

		//Rprintf("cmp_kriging.cpp:2339 model = %d sill=%f nugget=%f range=%f\n",parameters.ms,parameters.sills,parameters.nuggets,parameters.Hs);
	      // Calculate Gamma matrix
		//dout<<"here 6"<<endl;
	      subNeighbor.clear();
	      for (int i = 0; i < (int)sampleNeighbor.size(); ++i)
		{
		  distance = neighbor[sampleNeighbor[i]].distance;

		  if (distance <= parameters.Hs)
		    {
		      subNeighbor.push_back(sampleNeighbor[i]);
		    }
		}

		/* when no neighbor correlated, return less neighbor instead of subsampling
			07/14/2013;
	      if (subNeighbor.size() <= 1)
		{
		  if (sampleNeighbor.size() <= 50)
		    {
		      subNeighbor.clear();
		      for (int i = 0; i < (int)sampleNeighbor.size(); i++)
			{
			  subNeighbor.push_back(sampleNeighbor[i]);
			}
		    }
		  else
		    {
		      subNeighbor.clear();
		      for (int i = 0; i < (int)(0.1 * subNeighbor.size()); i++)
			{
			  subNeighbor.push_back(sampleNeighbor[i]);
			}
		    }				
		}
		*/
		if (subNeighbor.size() <=1)
		{
		  kvq	= LessNeighbors;
		  kriging	= -99999;
		  sigma = -99999;
		  //Rprintf("cmp_kriging.cpp:2380 No neighbor within %f returning\n",parameters.Hs);
		  return;
		}


	      /*if (parameters.nuggets == 0)
		{
		epsilon = 0.00001;
		}
		else
		{
		epsilon = 0;
		}*/

	      npts = (int)subNeighbor.size();

	      doublereal * GammaMatrix = new doublereal[(npts+1)*(npts+1)];
	      int TempK  = 0;
	      int TempP  = 0;
	      double Temp1;
	      //double Temp2;

	      // TEST output neighbor
	      // ofstream outN("./output/neighbor.csv", ofstream::out);
	      //   outN << "i,j,t,z" << endl;
	      //   for (int i = 0; i < npts; i++)
	      //   {
	      //   outN << neighbor[subNeighbor[i]].x_coord << "," << neighbor[subNeighbor[i]].y_coord << "," <<  neighbor[subNeighbor[i]].time << "," <<  neighbor[subNeighbor[i]].z_value << endl;
	      //   }
	      //   outN << endl;
	      //   outN << endl;

	      //   outN << "ms,Hs,nuggets,sills,mt,Ht,nuggett,sillt,kappa,epsilon" << endl;
	      //   outN << parameters.ms << "," << parameters.Hs << "," << parameters.nuggets << "," << parameters.sills << "," << parameters.mt << "," << parameters.Ht << "," << parameters.nuggett << "," << parameters.sillt << "," << parameters.kappa << "," << endl;
	      //     //epsilon << endl;

	      //   outN.close();
	      // TEST output neighbor

	      for (int i = 0; i < npts + 1; i++)
		{
		  for (int j = 0; j < npts + 1; j++)
		    {
		      TempK = j + (npts + 1) * i;
		      if (j == i)
			{
			  GammaMatrix[TempK] = 0;
			}
		      else
			{
			  if (j > i)
			    {
			      if (j == npts)
				{
				  GammaMatrix[TempK] = 1;
				}
			      else
				{
				  TempX1 = neighbor[subNeighbor[i]].x_coord;
				  TempY1 = neighbor[subNeighbor[i]].y_coord;

				  TempX2 = neighbor[subNeighbor[j]].x_coord;
				  TempY2 = neighbor[subNeighbor[j]].y_coord;

				  distance = sqrt((TempX1-TempX2)*(TempX1-TempX2)+(TempY1-TempY2)*(TempY1-TempY2));

				  switch (parameters.ms)
				    {
				    case 1: // Spherical
				      {
					if (distance >= parameters.Hs)
					  {
					    Temp1 = parameters.nuggets + parameters.sills;
					  }
					else
					  {
					    //Temp1 = nugs + sills * (1.5 * DistanceArray[i][j] / Hs - 0.5 * pow(DistanceArray[i][j] / Hs, 3));
					    Temp1 = parameters.nuggets + parameters.sills * (1.5 * distance / parameters.Hs - 0.5 * pow(distance / parameters.Hs, 3));
					  }
					break;
				      }
				    case 2: // Gaussian
				      {
					//Temp1 = nugs + sills * (1 - exp(-3 * pow(DistanceArray[i][j] / Hs, 2)));
					Temp1 = parameters.nuggets + parameters.sills * (1 - exp(-3 * pow(distance / parameters.Hs, 2)));
					break;
				      }
				    case 3: // Exponential
				      {
					//Temp1 = nugs + sills * (1 - exp(-3 * DistanceArray[i][j] / Hs));
					Temp1 = parameters.nuggets + parameters.sills * (1 - exp(-3 * distance / parameters.Hs));
					break;
				      }
				    case 4: // Matern
				      {
					//Temp1 = nugs + sills * (1 - (1 + 4.5 * DistanceArray[i][j] / Hs) * exp(-4.5 * DistanceArray[i][j] / Hs));
					Temp1 = parameters.nuggets + parameters.sills * (1 - (1 + 4.5 * distance / parameters.Hs) * exp(-4.5 * distance / parameters.Hs));
					break;
				      }
				    default:
				      {
					Temp1 = 0;
					break;
				      }
				    }
				  GammaMatrix[TempK] = (doublereal)Temp1;
				}
			    }
			  else // (p < q)
			    {
			      if (i == npts)
				{
				  GammaMatrix[TempK] = 1;
				}
			      else
				{
				  TempP =  i + (npts + 1) * j;
				  GammaMatrix[TempK] = GammaMatrix[TempP];
				}
			    }
			}
		    }
			//if (i==0){
			//	Rprintf("cmp_kriging.cpp:2503 Gamma[%d] = %f\n",TempK,GammaMatrix[TempK]);
			//}
		}

	      // TEST output Gamma matrix
	      // ofstream outGamma("./Gamma.csv", ofstream::out);
	      //   for (int i = 0; i < npts+1; i++)
	      //   {
	      //   for (int j = 0; j < npts+1; j++)
	      //   {
	      //   outGamma << "," << GammaMatrix[i*(npts+1)+j];
	      //   }
	      //   outGamma << endl;
	      //   }
	      //   outGamma.close();
	      // TEST output Gamma matrix

	      // calculate g vector by distance and time
	      doublereal *gVector;
	      doublereal *lamda;
	      doublereal dVector;
	      //doublereal tVector;
	      gVector	= new doublereal [npts+1];
	      lamda = new doublereal [npts+1];

	      for (int i = 0; i < npts; i++)
		{
		  dVector =  neighbor[subNeighbor[i]].distance;

		  switch (parameters.ms)
		    {
		    case 1: // Spherical
		      {
			if (dVector > parameters.Hs)
			  {
			    Temp1 = parameters.nuggets + parameters.sills;
			  }
			else
			  {
			    Temp1 = parameters.nuggets + parameters.sills * (1.5 * dVector / parameters.Hs - 0.5 * pow(dVector / parameters.Hs, 3));
			  }
			break;
		      }
		    case 2: // Gaussian
		      {
			Temp1 = parameters.nuggets + parameters.sills * (1 - exp(-3 * pow(dVector / parameters.Hs, 2)));
			break;
		      }
		    case 3: // Exponential
		      {
			Temp1 = parameters.nuggets + parameters.sills * (1 - exp(-3 * dVector / parameters.Hs));
			break;
		      }
		    case 4: // Matern
		      {
			Temp1 = parameters.nuggets + parameters.sills * (1 - (1 + 4.5 * dVector / parameters.Hs) * exp(-4.5 * dVector / parameters.Hs));
			break;
		      }
		    default:
		      {
			Temp1 = 0;
			break;
		      }
		    }
						
		  gVector[i] = (doublereal)Temp1;
		 //Rprintf("cmp_kriging.cpp:2569 gVector[%d] = %f\n",i,Temp1);

		}
	      gVector[npts] = 1;

	      // TEST output g vector
	      // ofstream outg("./gvector.csv", ofstream::out);
	      //   for (int i = 0; i < npts+1; i++)
	      //   {
	      //   outg << gVector[i] << endl;
	      //   }
	      //   outg.close();
	      // TEST output g vector

	      for (int i = 0; i < npts + 1; i++)
		{
		  lamda[i] = gVector[i];
		}

	      //METHOD 1: calculate weight vector through GammaMatrix and gVector, GammaMatrix * gVector = weights
	      /*  integer M = npts;
		integer nrhs = 1;
		integer lda = M;
		integer * ipiv = new integer[(npts + 1)];
		integer ldb = M;
		integer INFO;

		dgesv_(&M,&nrhs,GammaMatrix,&lda,ipiv,gVector,&ldb,&INFO); */

		//delete [] ipiv;
		//delete [] GammaMatrix;
	      // the weights value in the gVector now

	      // METHOD 2: calculate weight vector through GammaMatrix and gVector, GammaMatrix * gVector = weights
	      char uplo = 'U';
	      integer M = npts + 1;
	      integer nrhs = 1;
	      integer lda = M;
	      integer * ipiv = new integer[M];
	      integer ldb = M;

	      doublereal * c = new doublereal[(npts + 1)];
	      integer lwork = npts + 1;
	      integer INFO=1;

	      //Rprintf("pre dsysv_: INFO = %d\n",INFO);
	      dsysv_(&uplo,&M,&nrhs,GammaMatrix,&lda,ipiv,gVector,&ldb,c,&lwork,&INFO);
	      //Rprintf("cmp_kriging.cpp:2615 dsysv_: INFO = %d\n",INFO);

	      //delete [] c;
	      //delete [] ipiv;
	      //delete [] GammaMatrix;
	      // the weights value in the gVector now

	      // TEST output lamda vector
	      // ofstream outl("./output/lamda.csv", ofstream::out);
	      //   for (int i = 0; i < npts+1; i++)
	      //   {
	      //   outl << gVector[i] << endl;
	      //   }
	      //   outl.close();
	      // TEST output g vector
	
	      // Rprintf("INFO = %d\n",INFO);
		//if (INFO > 0) {
		//	Rprintf("INFO=%d > 0\n",INFO);
		//}
		//else if (INFO <0){
		//	Rprintf("INFO=%d < 0\n",INFO);
		//}
		//else {
		//	Rprintf("INFO=%d == 0\n",INFO);
		//}
	      //if (INFO != 0)
		//{
		  //kvq	= NotInverseMatrix;
		  //kriging	= -99999;
		  //sigma = -99999;
		  //Rprintf("cmp_kriging.cpp:2637 not invertible matrix, returning\n");
		//}
	      //else
		
	      if (INFO == 0)
		{
		  kriging = 0;
		  sigma = gVector[npts]; // 0;
		  for (int i = 0; i < npts; i++)
		    {
		      kriging = kriging + gVector[i]* neighbor[subNeighbor[i]].z_value;
		      sigma = sigma + lamda[i] * gVector[i];
		    }

		  if (sigma > 0)
		    {
		      sigma = sqrt(sigma);
		    }

		  if (kriging >= 0)
		    {
		      kvq	= Good;
		    }
		  else
		    {
		      kvq = OutofRange;
		    }
		  //Rprintf("cmp_kriging.cpp:2662 Kriging=%f Sigma=%f returning\n",kriging,sigma);
		}
	      else 
	        {
		  kvq	= NotInverseMatrix;
		  kriging	= -99999;
		  sigma = -99999;
		  Rprintf("cmp_kriging.cpp:2637 INFO=%d not invertible matrix, returning\n",INFO);
		}

	      Hs = parameters.Hs;
	      delete [] c;
	      delete [] ipiv;
	      delete [] GammaMatrix;
	      delete [] gVector;
	      delete [] lamda;

	      // free Bins memory
	      delete [] Bins;
	      Bins = NULL;
	    }
	}
    }
}
//dout.close();
}

void Variogram_By_Space_Calculation(		// calculate variogram by space function
				    AvDistSemi*	&Bins,	// neighbors are stored in this one demission array
				    int		 DistBin,	// the number of distance bins
				    double	 d_interval,	// the interval on distance
				    VariogramParameters&	parameters)	// variogram parameters (returned)
{
  double Hs = -99999;
  double nugs = -99999;
  double sills = -99999;
  ModelSelection ms = Init;
  double Ht = -99999;
  double nugt = -99999;
  double sillt = -99999;
  ModelSelection mt = Init;
  double kappa = -99999;

  double peakofsemi;
  int nn;
  double TempX1;
  double TempY1;
  double nugget_sphy;
  double sill_sphy;
  double nugget_mtey;
  double sill_mtey;
  double nugget_gausy;
  double sill_gausy;	
  double nugget_expy;
  double sill_expy;
  bool   sphy_flag;
  bool   expy_flag;
  bool   gausy_flag;
  bool   mtey_flag;
  double sse_sphy;
  double sse_expy;
  double sse_gausy;
  double sse_mtey;
  double sse = INFINITY;
  //AvDistTimeSemi* Bin;

  //for (int i = 0; i < DistBin; i++)
  //{
	//Rprintf("dist=%f gamma=%f n=%d \n",Bins[i].AvDistance,Bins[i].AvSemi,Bins[i].Pairs);
  //}
  // calculate how many missing values in the distance Bins
  nn = 0;
  for (int i = 0; i < DistBin; i++)
    {
      if (Bins[i].Pairs == 0)
	{
	  nn = nn + 1;
	}
    }

  // if the numbers of missing values in the distance Bins are greater than 0.5*Bins, interpolate the missing values
  if (nn >= 0.5*DistBin)
    {
	  //Rprintf("cmp_kriging.cpp:2734 first lag missing values, interpolate\n");

      for (int i = 0; i < DistBin; i++)
	{
	  if (Bins[i].Pairs == 0)
	    {
	      Bins[i].AvDistance = (i + 1) * d_interval;

	      Bins[i].AvSemi = InterpolateMissingValueBySpace(// find the nearest three bins and use this three bins' value to interpolate
							      Bins,			// neighbors are stored in this two demission array
							      DistBin,		// the number of distance bins
							      i);				// the index on distance Bin of missing value
				
	      Bins[i].Pairs = 1;
	  	  //Rprintf("cmp_kriging.cpp:2748 interpolate bin %d = %f\n",i,Bins[i].AvSemi);
	    }		
	}		
    }

  // search the semi peak
  peakofsemi = Bins[0].AvSemi;
  for (int i = 1; i < DistBin; i++)
    {
      peakofsemi = peakofsemi > Bins[i].AvSemi ? peakofsemi : Bins[i].AvSemi;
    }
  peakofsemi = 0.8 * peakofsemi;
  //Rprintf("cmp_kriging.cpp:2760 first Peak = %f\n",peakofsemi);

  // calculate Hs based on semi peak
  for (int i = 0; i < DistBin; i++)
    {
      if (Bins[i].AvSemi > peakofsemi)
	{
	  Hs = (i + 1) * d_interval; //Bin[i].AvDistance;
	  i = DistBin;
	}
    }
  //Rprintf("cmp_kriging.cpp:2771 Hs = %f \n", Hs);
  // Spherical model
  // SphY[jj] = nugget + sill*((3*x)/(2*H) - 0.5*pow(x/H, 3)) if x <= Hs
  // SphY[jj] = nugget + sill if x > Hs
  double sigma_sqp_sqf = 0;
  double sigma_sqp_z = 0;
  double sigma_sqp_f = 0;
  double sigma_sqp_f_z = 0;
  double sigma_sqp = 0;
			
  for (int i = 0; i < DistBin; i++)
    {
      if (Bins[i].Pairs != 0)
	{
	  if (Bins[i].AvDistance <= Hs)
	    {
	      TempX1 = 1.5 * Bins[i].AvDistance / Hs - 0.5 * pow((Bins[i].AvDistance / Hs), 3);
	    }
	  else
	    {
	      TempX1 = 1;
	    }
	  sigma_sqp_sqf = sigma_sqp_sqf + Bins[i].Pairs * TempX1 * TempX1;
	  sigma_sqp_z = sigma_sqp_z + Bins[i].Pairs * Bins[i].AvSemi;
	  sigma_sqp_f = sigma_sqp_f + Bins[i].Pairs * TempX1;
	  sigma_sqp_f_z = sigma_sqp_f_z + Bins[i].Pairs * TempX1 * Bins[i].AvSemi;
	  sigma_sqp = sigma_sqp + Bins[i].Pairs;
	}		
    }
  if ((sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f) != 0)
    {
      nugget_sphy = (sigma_sqp_sqf * sigma_sqp_z - sigma_sqp_f * sigma_sqp_f_z) / (sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f);
      sill_sphy = (sigma_sqp * sigma_sqp_f_z - sigma_sqp_z * sigma_sqp_f) / (sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f);
    }
  else
    {
      nugget_sphy = 0;
      sill_sphy = 0;
    }
	
  // Exponential model
  // ExpY[jj] = nugget + sill*(1 - exp(-3*x/H));
  sigma_sqp_sqf = 0;
  sigma_sqp_z = 0;
  sigma_sqp_f = 0;
  sigma_sqp_f_z = 0;
  sigma_sqp = 0;

  for (int i = 0; i < DistBin; i++)
    {
      if (Bins[i].Pairs != 0)
	{
	  TempX1 = 1 - exp(-3 * Bins[i].AvDistance / Hs);

	  sigma_sqp_sqf = sigma_sqp_sqf + Bins[i].Pairs * TempX1 * TempX1;
	  sigma_sqp_z = sigma_sqp_z + Bins[i].Pairs * Bins[i].AvSemi;
	  sigma_sqp_f = sigma_sqp_f + Bins[i].Pairs * TempX1;
	  sigma_sqp_f_z = sigma_sqp_f_z + Bins[i].Pairs * TempX1 * Bins[i].AvSemi;
	  sigma_sqp = sigma_sqp + Bins[i].Pairs;
	}
    }
  if ((sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f) != 0)
    {
      nugget_expy = (sigma_sqp_sqf * sigma_sqp_z - sigma_sqp_f * sigma_sqp_f_z) / (sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f);
      sill_expy = (sigma_sqp * sigma_sqp_f_z - sigma_sqp_z * sigma_sqp_f) / (sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f);
    }
  else
    {
      nugget_expy = 0;
      sill_expy = 0;
    }	

  // Gaussian model
  // GausY[jj] = nugget + sill*(1-exp(-3*(x/H)*(x/H)));
  sigma_sqp_sqf = 0;
  sigma_sqp_z = 0;
  sigma_sqp_f = 0;
  sigma_sqp_f_z = 0;
  sigma_sqp = 0;

  for (int i = 0; i < DistBin; i++)
    {
      if (Bins[i].Pairs != 0)
	{
	  TempX1 = 1 - exp(-3 * Bins[i].AvDistance * Bins[i].AvDistance / (Hs * Hs));

	  sigma_sqp_sqf = sigma_sqp_sqf + Bins[i].Pairs * TempX1 * TempX1;
	  sigma_sqp_z = sigma_sqp_z + Bins[i].Pairs * Bins[i].AvSemi;
	  sigma_sqp_f = sigma_sqp_f + Bins[i].Pairs * TempX1;
	  sigma_sqp_f_z = sigma_sqp_f_z + Bins[i].Pairs * TempX1 * Bins[i].AvSemi;
	  sigma_sqp = sigma_sqp + Bins[i].Pairs;
	}
    }
  if ((sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f) != 0)
    {
      nugget_gausy = (sigma_sqp_sqf * sigma_sqp_z - sigma_sqp_f * sigma_sqp_f_z) / (sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f);
      sill_gausy = (sigma_sqp * sigma_sqp_f_z - sigma_sqp_z * sigma_sqp_f) / (sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f);
    }
  else
    {
      nugget_gausy = 0;
      sill_gausy = 0;
    }

  // MteY model
  // MteY[jj] = nugget + sill*{1-[1+(4.5*x/H)]*exp[-4.5*x/H]};
  sigma_sqp_sqf = 0;
  sigma_sqp_z = 0;
  sigma_sqp_f = 0;
  sigma_sqp_f_z = 0;
  sigma_sqp = 0;

  for (int i = 0; i < DistBin; i++)
    {
      if (Bins[i].Pairs != 0)
	{
	  TempX1 = 1 - (1 + 4.5 * Bins[i].AvDistance / Hs) * exp(-4.5 * Bins[i].AvDistance / Hs);

	  sigma_sqp_sqf = sigma_sqp_sqf + Bins[i].Pairs * TempX1 * TempX1;
	  sigma_sqp_z = sigma_sqp_z + Bins[i].Pairs * Bins[i].AvSemi;
	  sigma_sqp_f = sigma_sqp_f + Bins[i].Pairs * TempX1;
	  sigma_sqp_f_z = sigma_sqp_f_z + Bins[i].Pairs * TempX1 * Bins[i].AvSemi;
	  sigma_sqp = sigma_sqp + Bins[i].Pairs;
	}
    }
  if ((sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f) != 0)
    {
      nugget_mtey = (sigma_sqp_sqf * sigma_sqp_z - sigma_sqp_f * sigma_sqp_f_z) / (sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f);
      sill_mtey = (sigma_sqp * sigma_sqp_f_z - sigma_sqp_z * sigma_sqp_f) / (sigma_sqp * sigma_sqp_sqf - sigma_sqp_f * sigma_sqp_f);
    }
  else
    {
      nugget_mtey = 0;
      sill_mtey = 0;
    }	

  sphy_flag = false;
  expy_flag = false;
  gausy_flag = false;
  mtey_flag = false;
  sse_sphy = 0;
  sse_expy = 0;
  sse_gausy = 0;
  sse_mtey = 0;

  if ((nugget_sphy > 0)&&(sill_sphy > 0))
    {
      sphy_flag = true;
      for (int i = 0; i < DistBin; i++)
	{
	  if (Bins[i].Pairs != 0)
	    {
	      if (Bins[i].AvDistance <= Hs)
		{
		  TempX1 = Bins[i].AvSemi - nugget_sphy - sill_sphy * (1.5 * Bins[i].AvDistance / Hs - 0.5 * pow((Bins[i].AvDistance / Hs), 3));
		}
	      else
		{
		  TempX1 = Bins[i].AvSemi - nugget_sphy - sill_sphy;
		}				
	      sse_sphy = sse_sphy + Bins[i].Pairs * TempX1 * TempX1;
	    }
	}
  	//Rprintf("cmp_kriging.cpp:2934 Spherical nugget= %f sill= %f SSE= %f \n", nugget_sphy,sill_sphy,sse_sphy);

    }
  if ((nugget_expy > 0)&&(sill_expy > 0))
    {
      expy_flag = true;
      for (int i = 0; i < DistBin; i++)
	{
	  if (Bins[i].Pairs != 0)
	    {
	      TempX1 = Bins[i].AvSemi - nugget_expy - sill_expy * (1 - exp(-3 * Bins[i].AvDistance / Hs));
	      sse_expy = sse_expy + Bins[i].Pairs * TempX1 * TempX1;
	    }
	}
   	//Rprintf("cmp_kriging.cpp:2949 Exponential nugget= %f sill= %f SSE= %f \n", nugget_expy,sill_expy,sse_expy);
   }
  if ((nugget_gausy > 0)&&(sill_gausy > 0))
    {
      gausy_flag = true;
      for (int i = 0; i < DistBin; i++)
	{
	  if (Bins[i].Pairs != 0)
	    {
	      TempX1 = Bins[i].AvSemi - nugget_gausy - sill_gausy * (1 - exp(-3 * Bins[i].AvDistance * Bins[i].AvDistance / (Hs * Hs)));
	      sse_gausy = sse_gausy + Bins[i].Pairs * TempX1 * TempX1;
	    }
	}
   	//Rprintf("cmp_kriging.cpp:2961 Gaussian nugget= %f sill= %f SSE= %f \n", nugget_gausy,sill_gausy,sse_gausy);
    }
  if ((nugget_mtey > 0)&&(sill_mtey > 0))
    {
      mtey_flag = true;
      for (int i = 0; i < DistBin; i++)
	{
	  if (Bins[i].Pairs != 0)
	    {
	      TempX1 = Bins[i].AvSemi - nugget_mtey - sill_mtey * (1 - (1 + 4.5 * Bins[i].AvDistance / Hs) * exp(-4.5 * Bins[i].AvDistance / Hs));
	      sse_mtey = sse_mtey + Bins[i].Pairs * TempX1 * TempX1;
	    }
	}
   	//Rprintf("cmp_kriging.cpp:2975 Matern nugget= %f sill= %f SSE= %f \n", nugget_mtey,sill_mtey,sse_mtey);

    }

  if ((sphy_flag == false)&&(expy_flag == false)&&(gausy_flag == false)&&(mtey_flag == false))
    {
      nugs = 0;
      sills = 0;
      double sigma_sqp_f_z = 0;
      double sigma_sqp_sqf = 0;
		
      for (int i = 0; i < DistBin; i++)
	{
	  if (Bins[i].Pairs != 0)
	    {
	      TempX1 = (1 - exp(-3 * Bins[i].AvDistance / Hs));
	      TempY1 = Bins[i].AvSemi;

	      sigma_sqp_f_z = sigma_sqp_f_z + Bins[i].Pairs * TempX1 * TempY1;
	      sigma_sqp_sqf = sigma_sqp_sqf + Bins[i].Pairs * TempX1 * TempX1;
	    }
	}
      if (sigma_sqp_sqf != 0)
	{
	  sills = sigma_sqp_f_z / sigma_sqp_sqf;
	}
      else
	{
	  sills = 0;
	}
      ms = Exponential;
    }
  else
    {
      if (sphy_flag == true)
	{
	  nugs = nugget_sphy;
	  sills = sill_sphy;
	  ms = Spherical;
	  sse = sse_sphy;
	}
      else if (gausy_flag == true)
	{
	  nugs = nugget_gausy;
	  sills = sill_gausy;
	  ms = Gaussian;
	  sse = sse_gausy;
	}
      else if (expy_flag == true)
	{
	  nugs = nugget_expy;
	  sills = sill_expy;
	  ms = Exponential;
	  sse = sse_expy;
	}
      else if (mtey_flag == true)
	{
	  nugs = nugget_mtey;
	  sills = sill_mtey;
	  ms = Matern;
	  sse = sse_mtey;
	}

      if ((sphy_flag == true)&&(sse_sphy < sse))
	{
	  nugs = nugget_sphy;
	  sills = sill_sphy;
	  ms = Spherical;
	  sse = sse_sphy;
	}
      if ((gausy_flag == true)&&(sse_gausy < sse))
	{
	  nugs = nugget_gausy;
	  sills = sill_gausy;
	  ms = Gaussian;
	  sse = sse_gausy;
	}
      if ((expy_flag == true)&&(sse_expy < sse))
	{
	  nugs = nugget_expy;
	  sills = sill_expy;
	  ms = Exponential;
	  sse = sse_expy;
	}
      if ((mtey_flag == true)&&(sse_mtey < sse))
	{
	  nugs = nugget_mtey;
	  sills = sill_mtey;
	  ms = Matern;
	  sse = sse_mtey;
	}
    }

  parameters.Hs = Hs;
  parameters.nuggets = nugs;
  parameters.sills = sills;
  parameters.ms = ms;
  parameters.Ht = Ht;
  parameters.nuggett = nugt;
  parameters.sillt = sillt;
  parameters.mt = mt;
  parameters.kappa = kappa;
}

double  InterpolateMissingValueBySpace(				// find the nearest three bins and use this three bins' value to interpolate
				       AvDistSemi*				&Bins,				// neighbors are stored in this one demission array
				       int						DistBin,			// the number of distance bins
				       int						d)					// the index on distance Bin of missing value
{
  double semi = 0;
  int d1,d2;
  int n;
  int nn;	
  double Temp1 = 0;
  double Temp2 = 0;

  nn = 1;
  n = 0;
  while ((n < 3)&&(nn <= DistBin))
    {
      d1 = d - nn;
      d2 = d + nn;

      if (d1 >= 0)
	{
	  if (Bins[d1].Pairs != 0)
	    {
	      Temp1 = Temp1 + Bins[d1].AvSemi * Bins[d1].Pairs / nn;
	      Temp2 = Temp2 + Bins[d1].Pairs / nn;
	      n = n + 1;
	    }
	}
		
      if (d2 < DistBin)
	{
	  if (Bins[d2].Pairs != 0)
	    {
	      Temp1 = Temp1 + Bins[d2].AvSemi * Bins[d2].Pairs / nn;
	      Temp2 = Temp2 + Bins[d2].Pairs / nn;
	      n = n + 1;
	    }
	}
		
      nn = nn + 1;
    }

  if (Temp2 != 0)
    {
      semi = Temp1 / Temp2;
    }	

  return semi;
}

#ifdef _MANAGED
#pragma managed(pop)
#endif
