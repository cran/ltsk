#include "io.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <cassert>
#include <string>
#include <cstdlib>
#include <cmath>

using namespace std;

// void r_file_s(const char* filename, int longCol, int latCol, int timeCol, int zCol)
// {
  // string s;
  // int time;
  // double lon,lat,z_one;
  // int timediff;
  // double distance;
  
  // fstream file(filename);
  // if(!file)
    // {
      // //cerr << "Unable to open file " << filename << endl;;
      // return;
    // }

  // //uid,along,alat,ayear,amonth,aday,cday,pr_pm25,pr_pm10,pr_pmc
  // //lon, lat, day, z
  // while(file.good())
    // {
      // s.clear();
      // getline(file,s);
      // if(s.empty()) break;
      // char_separator<char> sep(",;| ");
      // typedef tokenizer<char_separator<char> > tokC;
      // tokC tok(s,sep);
      // tokenizer<char_separator<char> >::iterator beg=tok.begin();
      // vector<string> vec;
      // vec.assign(tok.begin(),tok.end());
      // lon=lexical_cast<double> (vec[longCol]);
      // lat=lexical_cast<double> (vec[latCol]);
      // time=lexical_cast<int> (vec[timeCol]);
      // trim(vec[zCol]);
      // z_one=lexical_cast<double> (vec[zCol]);
      
      // plist.push_back(point(lon,lat,time, z_one, timediff, distance));

    // } 
  // file.close ();
  // // plist.pop_back();
// }

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
   //plist.pop_back();
}

// int r_file_s(const char* filename, ANNpointArray& pointArr, int longCol, int latCol,  int zCol)
// {
  // string s;
  // double lon,lat,z_one;
  // double distance;
  // int npts=0;
  
  // fstream file(filename);
  // if(!file)
    // {
      // //cerr << "Unable to open file " << filename << endl;;
      // return -1;
    // }

  // //uid,along,alat,ayear,amonth,aday,cday,pr_pm25,pr_pm10,pr_pmc
  // //lon, lat, day, z
  // while(file.good())    {
    // s.clear();
    // getline(file,s);
    // if(s.empty()) break;
    // // char_separator<char> sep(",","", keep_empty_tokens);
    // // typedef tokenizer<char_separator<char> > tokC;
    // typedef tokenizer<escaped_list_separator<char> > tokC;    
    // // tokC tok(s,sep);
    // tokC tok(s);    
    // // tokenizer<char_separator<char> >::iterator beg=tok.begin();
    // tokenizer<escaped_list_separator<char> >::iterator beg=tok.begin();    
    // vector<string> vec;
    // vec.assign(tok.begin(),tok.end());
    // // cout << vec.size() << endl;
    // lon=lexical_cast<double> (vec[longCol]);
    // pointArr[npts][0]=lon;
    // lat=lexical_cast<double> (vec[latCol]);
    // pointArr[npts][1]=lat;
    // // trim(vec[timeCol]);      
    // // time=lexical_cast<int> (vec[timeCol]);
    // trim(vec[zCol]);
    // // if(vec[zCol].empty() || vec[zCol].length()==0 )  	{
    // // cout << vec[zCol] << endl;
      // // cout << "empty " << endl;
    // //   z_one = 0;
    // // }
    // // else {
    // // cout << vec[zCol] << endl;
    // // z_one=lexical_cast<double> (vec[zCol]);
    // // if(vec.size() != 15) {
    // //  cout <<  vec.size() << " ";      
    // //  cout <<  vec[longCol] << " ";
    // //  cout <<  vec[latCol] << " ";
    // //  cout << vec[zCol] << endl;
    // // }
    
    // if ( vec[zCol] == "" )
      // // cout << " null " << endl;
      // continue;
    // else {
      // z_one=lexical_cast<double> (vec[zCol]);
      // // cout << z_one << endl;
      
      // // cout << npts << endl;      
      
      // plist.push_back(point(lon,lat,z_one));
      // npts++;
    // }

    // // if ( npts == 10 ) return -1;
    // //quit
    // // cout << plist[0].z_value << " " << z_one << endl;
    // // cout << plist.size() << endl;
    // // return 0;
    
  // } 
  // file.close ();
  // //cout << "done reading source" << endl;      
  // return npts;
  // // plist.pop_back();
// }

// void r_file_q(const char* filename, std::vector<point_q>& list_q, int xCol, int yCol, int dayCol)
// {
  // string s;
  // double x,y;
  // int z;

  // fstream file(filename);
  // if(!file)
    // {
      // //cerr << "Unable to open file " << filename << endl;;
      // return;
    // }

  // //cday,uid,id,x_coord,y_coord
  // // lon, lat, day
  // while(file.good())
    // {
      // s.clear();
      // getline(file,s);
      // if(s.empty()) break;
      // char_separator<char> sep(",;| ");
      // typedef tokenizer<char_separator<char> > tokC;
      // tokC tok(s,sep);
      // // vector<double> dataLine;
      // tokenizer<char_separator<char> >::iterator beg=tok.begin();
      // vector<string> vec;
      // vec.assign(tok.begin(),tok.end());
      // trim(vec[xCol]);
      // trim(vec[yCol]);	   
      // trim(vec[dayCol]);
      // x=lexical_cast<double> (vec[xCol]);
      // y=lexical_cast<double> (vec[yCol]);
      // z=lexical_cast<int> (vec[dayCol]);

      // // cout << x << " ";
      // // cout << y << " ";
      // // cout << z << endl;
      // // 		  dataLine.push_back(lexical_cast<double> (*beg) );
      
      // list_q.push_back(point_q(x,y,z));

    // } 
  // file.close ();
  // // plist.pop_back();
// }

void r_file_q(double *lon, int *lon_n, double *lat, int *lat_n, std::vector<point_q>& list_q)
{
  int i;
  for(i=0; i< *lon_n; i++)
  {
    list_q.push_back(point_q(lon[i],lat[i]));
  }
  //plist.pop_back();
}

// void r_file_q(const char* filename, std::vector<point_q>& list_q, int xCol, int yCol)
// {
  // string s;
  // double x,y;
  // int z;

  // fstream file(filename);
  // if(!file)
    // {
      // //cerr << "Unable to open file " << filename << endl;;
      // return;
    // }

  // //cday,uid,id,x_coord,y_coord
  // // lon, lat, day
  // // cout <<  "begin reading  "<< filename << endl;
  // // cout << xCol << " " << yCol << endl;
  // while(file.good())
    // {
      // s.clear();
      // getline(file,s);
      // if(s.empty()) break;
      
      // char_separator<char> sep(",","", keep_empty_tokens);
      // // typedef tokenizer<char_separator<char> > tokC;
      // // tokenizer<escaped_list_separator<char> > tok(s);
      // typedef tokenizer<escaped_list_separator<char> > tokC;
      
      // // tokC tok(s,sep);
      // tokC tok(s);
      // // vector<double> dataLine;
      // // tokenizer<char_separator<char> >::iterator beg=tok.begin();
      // tokenizer<escaped_list_separator<char> >::iterator beg=tok.begin();      
      // vector<string> vec;
      // vec.assign(tok.begin(),tok.end());
      // trim(vec[xCol]);
      // trim(vec[yCol]);
      // // cout <<  vec.size() << " ";      
      // // cout <<  vec[xCol] << " ";
      // // cout <<  vec[yCol] << endl;      
      // x=lexical_cast<double> (vec[xCol]);
      // y=lexical_cast<double> (vec[yCol]);

      // // cout << x << " ";
      // // cout << y << " ";
      // // cout << z << endl;
      // // 		  dataLine.push_back(lexical_cast<double> (*beg) );
      
      // list_q.push_back(point_q(x,y));

    // } 
  // file.close ();
  // // plist.pop_back();
// }

void calculate(int& index, double& qlat, double& qlon, int& qtime)
{
  plist[index].distance = sqrt(pow((plist[index].x_coord - qlat), 2) + pow((plist[index].y_coord - qlon), 2));
  plist[index].timediff = abs(plist[index].time - qtime);
};

void calculate(int& index, double& qlat, double& qlon)
{
  plist[index].distance = sqrt(pow((plist[index].x_coord - qlat), 2) + pow((plist[index].y_coord - qlon), 2));
  // Rprintf("cal dist ... %d  %f \n",index, plist[index].distance);
  // cout << "cal dist... " << index << " " << plist[index].distance << endl;
}
