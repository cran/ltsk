#include <R_ext/Rdynload.h>

#include "ltsk.h"
#include <iomanip>    
using namespace std;					// make std:: accessible
vector<point> plist;

extern "C" { 
int lk_main(double *x, int *xn, double *y, int *yn, double *z, int *zn, double *lon, int *lon_n, double *lat, int *lat_n,double *th, int *vlen, double *krig_out, double *sigma_out, double *hs_out, double *psill_out, double *nugget_out, int *ms_out)
{
  vector<point_q> list_q;

  int distbin= *vlen;

  int MaxPt;
  MaxPt = *(xn) + 1000;

  double dist_range;
  dist_range = *th;
  //Rprintf("Distance range= %f\n",dist_range);


  plist.clear(); // clear the point list;
  
  r_file_q(lon, lon_n, lat, lat_n, list_q);

  ANNpointArray	dataPts_2;	// data points
  dataPts_2 = annAllocPts(MaxPt, 2);			// allocate data points  
  int MaxPts = r_file_s(x, xn, y, yn, z, zn, dataPts_2);

  query_instance(dataPts_2, list_q, dist_range, distbin, MaxPts, 
		&krig_out[0],&sigma_out[0], &hs_out[0], &psill_out[0], &nugget_out[0], &ms_out[0]);
  return EXIT_SUCCESS;
}
}

int query_instance(ANNpointArray& dataPts_2, vector<point_q>& list_q, double dist, int distBin, int nPts, double *krig_out, 
	double *sigma_out, double *hs_out, double *psill_out, double *nugget_out, int *ms_out)
{
  int	k      = 0;		// number of nearest neighbors
  // int	dim    = 0;		// dimension
  double	eps    = 0;		// error bound
  ANNpoint	queryPt_2;	// query point
  // int		nPts;		// actual number of data points
  ANNkd_tree*	kdTree_2;	// search structure
  
  ModelSelection m=Init;
  double alpha = 2.0;
  //double beta;
  //CalculateMethod c=None;
  double kriging;
  double sigma;
  KrigingValueQuality kvq=Empty;
  VariogramParameters	parameters;
  
  NeighborVector nv;
  nv.clear();
  
  double Hs;
  
  int npts=0;

  ANNidxArray		nnIdx;					// near neighbor indices
  ANNdistArray		dists;					// near neighbor distances
  // ANNpointArray		res_dist;

  kdTree_2= new ANNkd_tree(			// build search structure
			   dataPts_2,		// the data points
			   nPts,		// number of points
			   2);			// dimension of space

  
  for (unsigned int i=0; i<list_q.size(); i++)
    {
      // if(i!=93) continue;
      queryPt_2 = annAllocPt(2);					// allocate query point
      
      queryPt_2[0]=list_q[i].x();
      queryPt_2[1]=list_q[i].y();

      npts=kdTree_2->annkFRSearch(
				  queryPt_2, // query point
				  dist, // squared radius
				  0 , // number of near neighbors to return
				  0 , // nearest neighbor array (modified)
				  0 , // dist to near neighbors (modified)
				  0 ); // error bound)))

	// fout << npts<<endl;			  
    if(npts>0)
	{
	  k=npts;
	  nnIdx = new ANNidx[k];						// allocate near neigh indices
	  dists = new ANNdist[k];						// allocate near neighbor dists

	  kdTree_2->annkFRSearch(
				 queryPt_2, // query point
				 dist, // squared radius
				 k , // number of near neighbors to return
				 nnIdx , // nearest neighbor array (modified)
				 dists , // dist to near neighbors (modified)
				 eps ); // error bound)))

 
	  for (int j = 0; j < k; j++) 
	    {
	      nv.push_back(nnIdx[j]);
	      calculate(nnIdx[j], queryPt_2[0], queryPt_2[1]);
	    }

	  // // // cout << "kriging starting " << endl;
	  // // // cout << "nv size before passing in " << nv.size() << endl;
	  // // // cout << "nv 1st index before passing in " << nv.get_index(0) << endl;
	  
	 //Rprintf("ltsk.cpp:176 kriging starting for (%f,%f)\n",queryPt_2[0],queryPt_2[1]);
	 //Rprintf("ltsk.cpp:177 neighborhood size before passing in %d\n", nv.size());
	 
	 //Rprintf("nv 1st index before passing in \n",nv.get_index(0));
	  Kriging_By_Space( nv, distBin, m, alpha, kriging, sigma, Hs, parameters, kvq);
	  //Rprintf("query=(%f,%f), Kriging=%f sigma=%f Hs=%f \n",queryPt_2[0],queryPt_2[1],kriging,sigma,Hs);
	  krig_out[i] = kriging;
	  sigma_out[i] = sigma;
	  hs_out[i] = Hs;
	  psill_out[i] = parameters.sills;
	  nugget_out[i] = parameters.nuggets;
	  ms_out[i] = parameters.ms;
	  nv.clear();
	  delete [] nnIdx;	// clean things up
	  delete [] dists;
	  delete [] queryPt_2;
	  // kriging=0;
	}
      else
	{
	  //fout << setprecision(10) << ",null" << ",null" << ",null" <<  ",null" << ",null" << ",null" << endl; 
	  krig_out[i] = -99999;
	  sigma_out[i] = -99999;
	  hs_out[i] = -99999;
	  psill_out[i] = -99999;
	  nugget_out[i] = -99999;
	  ms_out[i] = -99999;
	}
    }
  delete kdTree_2;
  annClose();		// done with ANN

  return EXIT_SUCCESS;
}

static const R_CMethodDef cMethods[] = {
	{"lk_main", (DL_FUNC) &lk_main, 10},
		{NULL,NULL,0}
	};

void R_init_ltsk(DllInfo *info)
{
	R_registerRoutines(info, cMethods, NULL, NULL,NULL);
}


