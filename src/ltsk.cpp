#include <R_ext/Rdynload.h>

#include "ltsk.h"
#include <iomanip>    
using namespace std;					// make std:: accessible
vector<point> plist;

extern "C" { 
int lk_main(double *x, int *xn, double *y, int *yn, double *z, int *zn, double *lon, int *lon_n, double *lat, int *lat_n)
{
  ConfigFile config( "./config.txt" );

  vector<point_q> list_q;

  string RS, RT;
  int Rs, Rt;
  config.readInto(RS, "Rs");
  config.readInto(RT, "Rt");
  Rs=atoi(RS.c_str());
  Rt=atoi(RT.c_str());
  int distbin=Rs;

  string numMax;
  int MaxPt;
  config.readInto(numMax, "MAXPTS_ANN");
  MaxPt=atoi(numMax.c_str());

  string Dist_Range;
  double dist_range;
  config.readInto(Dist_Range, "dist_range");
  //Rprintf("Distance range= %s\n",Dist_Range.c_str());
  dist_range=atof(Dist_Range.c_str());
  //Rprintf("Distance range= %f\n",dist_range);

  string Dist_Step;
  double dist_step;
  config.readInto(Dist_Step, "dist_interval");
  dist_step=atof(Dist_Step.c_str());

  vector<double> distVec;

  for(int i=1; i<=dist_range/dist_step;i++) {
    distVec.push_back(dist_step*i);
  }

  plist.clear(); // clear the point list;
  
  r_file_q(lon, lon_n, lat, lat_n, list_q);

  ANNpointArray	dataPts_2;	// data points
  dataPts_2 = annAllocPts(MaxPt, 2);			// allocate data points  
  int MaxPts = r_file_s(x, xn, y, yn, z, zn, dataPts_2);

  //for(int i=0; i<5; i++){
	//Rprintf("query %d = (%f,%f)\n",i,list_q[i].x(),list_q[i].y());
  //}
  //for(int i=0; i<5; i++){
	//Rprintf("obs %d = (%f,%f)\n",i,dataPts_2[i][0],dataPts_2[i][1]);
  //}
  //Rprintf("Rs=%d\n distbin=%d MaxPts=%d",Rs,distbin,MaxPts);
   for(int i=0; i<distVec.size(); i++) {
      string s;
      stringstream out;
      //out << distVec[i];
	  out << i;
      s = out.str();
      string num="_" + s;
      char outname[100]="output";
      strcat(outname, num.c_str());
      // // cout << outname << endl;
	  
      query_instance(dataPts_2, list_q, Rs, distVec[i],  outname, distbin, MaxPts);
    }
  distVec.clear();
}

}

int query_instance(ANNpointArray& dataPts_2, vector<point_q>& list_q, int Rs, double dist, char* outname, int distBin, int nPts )
{
  int	k      = 0;		// number of nearest neighbors
  // int	dim    = 0;		// dimension
  double	eps    = 0;		// error bound
  ANNpoint	queryPt_2;	// query point
  // int		nPts;		// actual number of data points
  ANNkd_tree*	kdTree_2;	// search structure
  
  ModelSelection m=Init;
  double alpha;
  double beta;
  CalculateMethod c=None;
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

  ofstream fout;
  fout.open(outname);
  // cout << outname << endl;
  
  for (int i=0; i<list_q.size(); i++)
    {
      // if(i!=93) continue;
      queryPt_2 = annAllocPt(2);					// allocate query point
      
      queryPt_2[0]=list_q[i].x();
      queryPt_2[1]=list_q[i].y();
      // // // cout << queryPt_2[0] << " ";
      // // // cout << queryPt_2[1] << " ";
      // // // cout << dist << endl;            
      // fout << queryPt_2[0] << " ";
      // fout << queryPt_2[1] << " ";
      // fout << dist << endl;            

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

	  //// // cout << "\tNN:\tIndex\tDistance\n";

	  // vector<int> res;
	  // fout << k <<endl;			  
 
	  for (int i = 0; i < k; i++) 
	    {
	      nv.push_back(nnIdx[i]);
	      calculate(nnIdx[i], queryPt_2[0], queryPt_2[1]);
	    }
	 //ofstream dout;
	 //dout.open("debug_1");
	 //dout << "x0=" << queryPt_2[0] <<" y0=" << queryPt_2[1] <<endl;
	 //for (int i = 0; i < k; i++)
	 //{
			//dout << "x=" << nv[nnIdx[i]].x_coord<<" y= "<< nv[nnIdx[i]].y_coord << "distance= " << nv[nnIdx[i]].distance << endl;
	 //}
	 //dout.close();

	  // // // cout << "kriging starting " << endl;
	  // // // cout << "nv size before passing in " << nv.size() << endl;
	  // // // cout << "nv 1st index before passing in " << nv.get_index(0) << endl;
	 //Rprintf("kriging starting\n");
	 //Rprintf("nv size before passing in %d\n", nv.size());
	 //Rprintf("nv 1st index before passing in \n",nv.get_index(0));
	  Kriging_By_Space( nv, distBin, m, alpha, kriging, sigma, Hs, parameters, kvq);
	  // fout << setprecision(9) <<  queryPt_2[0] << ", " << queryPt_2[1] << ", " << kriging << " " << sigma << " " << Hs << " " << kvq << endl;
	  fout << setprecision(10) <<  "," << kriging <<  "," << sigma << "," << Hs << "," << parameters.nuggets << "," << parameters.sills << "," <<parameters.ms <<endl;	  
	  // // // cout << queryPt_2[0] << ", " << queryPt_2[1] << " " << kriging << " " << sigma << " " << Hs << " " << kvq << endl;
	  //Rprintf("query=(%f,%f), Kriging=%f sigma=%f Hs=%f \n",queryPt_2[0],queryPt_2[1],kriging,sigma,Hs);
	  nv.clear();
	  delete [] nnIdx;	// clean things up
	  delete [] dists;
	  // kriging=0;
	}
      else
	{
	  fout << setprecision(10) << ",null" << ",null" << ",null" <<  ",null" << ",null" << ",null" << endl; 
	  // fout << setprecision(10) <<  queryPt_2[0] << ", " << queryPt_2[1] << ", null"<< endl;
}
    }
  fout.close();
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


