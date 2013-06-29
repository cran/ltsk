#include <cstdlib>						// C standard library
#include <cstdio>						// C I/O (for sscanf)
#include <cstring>						// string manipulation
#include <fstream>						// file I/O
#include <time.h>					// ANN declarations
#include <vector>
#include "point.h"	
#include "io.h"
#include "ConfigFile.h"
#include "ANN/ANN.h"					// ANN declarations
#include "cmp_kriging.h"
#include "range_tree.h"
#include "R.h"


int query_instance(ANNpointArray& dataPts_2, vector<point_q>& list_q, int Rs, double dist, char* outname, int distBin, int nPts );

