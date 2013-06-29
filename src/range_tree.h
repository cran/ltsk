#ifndef RANGE_TREE_H
#define RANGE_TREE_H

#include <iostream>    
#include <vector>    
#include "AvlTree.h"
#include "point.h"
#include <iomanip>
using namespace std;       

inline bool compare_lat(const int idx_a, const int idx_b) 
{ return plist[idx_a].y_coord < plist[idx_b].y_coord ; }

inline bool compare_lon(const int idx_a, const int idx_b)
{ return plist[idx_a].x_coord < plist[idx_b].x_coord ; }


class spliting_node
{
 private:
 public:
  int p_idx;
  double bound[2]; //
  AvlTree* as;
  spliting_node* left_child;
  spliting_node* right_child;
  int hei_last;

 public:
  inline spliting_node():p_idx(-1), left_child(NULL), right_child(NULL) , hei_last(-1)
    { 
      as = new AvlTree(); 
    }

  inline spliting_node(spliting_node* npt):p_idx(npt->p_idx), left_child(npt->left_child), right_child(npt->right_child), hei_last(npt->hei_last)
    {
      as = new AvlTree(*npt->as);
    }

  ~spliting_node()
    {
      if(as!=NULL) delete as;
      if(left_child!=NULL) delete left_child;
      if(right_child!=NULL) delete right_child;

    }

  inline spliting_node operator=(spliting_node* n)
    {
      if(this !=n)
	{
	  this->p_idx=n->p_idx;
	  this->left_child=n->left_child;
	  this->right_child=n->right_child;
	  this->as=n->as;
	  this->bound[0]=n->bound[0];
	  this->bound[1]=n->bound[1];
	}
      return *this;
    }

 spliting_node(int po_idx):p_idx(po_idx), left_child(NULL), right_child(NULL) {}
  //spliting_node(point po_idx, spliting_node* l=NULL, spliting_node* r=NULL):p_idx(po_idx), left_child(l), right_child(r) {}
 spliting_node(int po_idx, spliting_node* l, spliting_node* r):p_idx(po_idx), left_child(l), right_child(r) {}


  void print_assoc_sturture()
  {
    if(as) { as->printTree(); }
  }

  inline void print_bound()
  {
    cout << "bound " <<bound[0] <<  " " << bound[1] <<  endl;
  }
  friend class r_tree;

};

static int height_last=0;

class r_tree   /* sorted by x_coord */
{
 public:
  //r_tree():root(NULL)
 r_tree():root(NULL)
    { };

 r_tree(vector<int>& input_vector):root(NULL)
    {
      build_tree(input_vector);
    }

  ~r_tree()
    {
      //makeEmpty(root);
    }

  void makeEmpty(spliting_node* t)
  {
    if( t != NULL )
      {
	makeEmpty( t->left_child );
	makeEmpty( t->right_child );
	delete t;
      }
    t = NULL;
  }

  // assumes input_vector is not sorted
  spliting_node* build_tree(vector<int>& input_vector)
  {
    /* cout << "first " << plist[input_vector.front()].x_coord <<endl; */
    /* cout << "last " << plist[input_vector.back()].x_coord <<endl;     */
    sort(input_vector.begin(),input_vector.end(), compare_lon);
    return build_tree(input_vector, root);
  }


  inline bool is_leaf(spliting_node*& node)
  {
    return (node!=NULL && node->left_child == NULL && node->right_child) ? true:false;

  }

  //vector<int> query(double from_x, double to_x, double from_y, double to_y)
  void query(vector<int> & query_result, double from_x, double to_x, double from_y, double to_y)
  {
    spliting_node* v=root;
    /* v->print_bound(); */
    while(v)
      {
	/* v->print_bound(); */
	if(v->left_child==NULL && v->right_child==NULL)
	  {
	    //cout << " no child, out " << endl;
	    break;
	  }
	if(v->left_child->bound[1] >= to_x )
	  {
	    //cout << " left " << endl;
	    v=v->left_child;
	  }
	else if(v->right_child->bound[0] <= from_x )
	  {
	    //cout << " right " << endl;
	    v=v->right_child;
	  }
	else
	  {
	    break;
	  }
      }

    //cout << endl;
    //cout << "after split " << endl;
    /* v->print_bound(); */
    /* v->left_child->print_bound(); */
    /* v->right_child->print_bound(); */
    
    /* v->print_assoc_sturture(); */
    /* cout << endl; */

    
    /*   if(!v) */
    /*   { */
    /*   return; */
    /*   } */
    
    /* if(is_leaf(v)) */
    /* if(v->hei_last==0) */
    if(!v->left_child && !v->right_child)
      {
	if(plist[v->p_idx].y_coord >= from_y && plist[v->p_idx].y_coord <= to_y
	   && plist[v->p_idx].x_coord >=from_x
	   && plist[v->p_idx].x_coord <=to_x
	   )
	  {
	    //assert(v->p_idx==0);
	    query_result.push_back(v->p_idx);
/* 	cout << "	find  " << setprecision (10) */
/* 	<< " " << plist[v->p_idx].x_coord */
/* 	<< " " << plist[v->p_idx].y_coord */
/* 	<< " " << plist[v->p_idx].cday */
/* 	<< endl; */
	    
	  }
      }
    else
      {
	//spliting_node* split=new spliting_node(v);
	spliting_node* tmp=v;

	// left path
	if(v->left_child)  v=v->left_child;
	
	//while(v->right_child)
	//while(v->hei_last>0)
	while(v->left_child || v->right_child)
	  {
	    if(from_x <= plist[v->p_idx].x_coord && plist[v->p_idx].x_coord <= to_x)
	      {
		//v->print_assoc_sturture();
		/* cout << "node " ; v->print_bound(); */
		/* cout << "right child " ; v->right_child->print_bound(); */
		/* cout << "begin 1D search" << endl;		 */
		query(query_result, v->right_child, from_y, to_y);
		v=v->left_child;
	      }
	    else
	      {
		v=v->right_child;
	      }
	  }
	
	if(plist[v->p_idx].y_coord >= from_y && plist[v->p_idx].y_coord <= to_y
	   && plist[v->p_idx].x_coord >=from_x
	   && plist[v->p_idx].x_coord <=to_x
	   )
	  {
	    //	    cout << "find one ";//	v->p_idx.print();
	    //assert(v->p_idx>=0);
	    query_result.push_back(v->p_idx);
/* 	    cout << "	find  " << setprecision (10) */
/* 		 << " " << plist[v->p_idx].x_coord */
/* 		 << " " << plist[v->p_idx].y_coord */
/* 		 << " " << plist[v->p_idx].cday */
/* 		 << endl; */
	    
	  }

	// right path
	v=tmp->right_child;
	//while(v->hei_last>0)
	while(v->left_child || v->right_child)
	  {
	    if(plist[v->p_idx].x_coord <= to_x && plist[v->p_idx].x_coord >= from_x)
	      {
		//v->print_assoc_sturture();
		/* cout << " node " ;v->print_bound() ;				 */
		/* cout << "left child " ; v->left_child->print_bound(); */
		/* cout << "begin 1D search" << endl; */
		query(query_result, v->left_child, from_y, to_y);
		v=v->right_child;
	      }
	    else
	      {
		v=v->left_child;
	      }
	  }
	
	if(plist[v->p_idx].y_coord >= from_y && plist[v->p_idx].y_coord <= to_y
	   && plist[v->p_idx].x_coord >=from_x   && plist[v->p_idx].x_coord <=to_x
	   )
	  {
	    //	    cout << "find one ";//	v->p_idx.print();
	    //assert(v->p_idx>=0);
	    query_result.push_back(v->p_idx);
/* 	    cout << "	find  " << setprecision (10) */
/* 	    	 << " " << plist[v->p_idx].x_coord */
/* 	    	 << " " << plist[v->p_idx].y_coord */
/* 	    	 << " " << plist[v->p_idx].cday */
/* 	    	 << endl; */
	    
	  }

      }
    //return query_result;
  }

  void query(vector<int>& res, spliting_node*& v, double from, double to)
  {
    v->as->query(res, v->as->root, from, to);
    v->print_assoc_sturture();
    /* cout << " " << plist[v->as->findMin()].y_coord */
    /* 	 <<" " << plist[v->as->findMax()].y_coord */
    /* 	 << endl; */
  }


 private:
 public:
  spliting_node* root;
		
  // build a tree rooted on v
  spliting_node* build_tree(vector<int>& input_vector, spliting_node*& v)
  {
    if(!v)
      {
	v = new spliting_node();
	v->hei_last=height_last;
      }
    AvlTree* newas= new AvlTree(input_vector);
    v->bound[0]=plist[input_vector.front()].x_coord;
    v->bound[1]=plist[input_vector.back()].x_coord;
    /* cout << "size " << input_vector.size() << " " << setprecision (10) <<v->bound[0] << "  " */
    /* 	 << v->bound[1] << endl; */
    if(input_vector.size()==1)
      {
	v->p_idx=input_vector.front();
	v->as=newas;
	v->left_child=v->right_child=NULL;
	v->hei_last=0;
      }
    else
      {
	int median = input_vector.size()/2;
	vector<int> left_subset(input_vector.begin(), input_vector.begin()+median);
	vector<int> right_subset(input_vector.begin()+median, input_vector.end());
	/* cout << "	left subtree "<<left_subset.size() << endl;	 */
	/* cout <<"	right subtree "<< right_subset.size() << endl; */

	v->p_idx=input_vector[median];

	height_last++;
	//cout << " node  " ; v->p_idx.print() ;
	v->left_child=build_tree(left_subset,v->left_child );
	v->right_child=build_tree(right_subset, v->right_child);
	v->as=newas;
      }
    return v;
  }
};

#ifdef DEFINE_GLOBALS
#define GLOBAL
#else // !DEFINE_GLOBALS
#define GLOBAL extern
#endif

GLOBAL r_tree *r;



#endif
