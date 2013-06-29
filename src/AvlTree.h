#ifndef AVL_TREE_H
#define AVL_TREE_H

#include <iostream>    // For NULL
#include <cassert>    // For NULL
#include <iomanip>
#include "dsexceptions.h"
#include "point.h"

using namespace std;

// AvlTree class
//
// CONSTRUCTION: with ITEM_NOT_FOUND object used to signal failed finds
//
// Throws UnderflowException as warranted

class AvlTree
{
 public:
  inline AvlTree( ) : root( NULL ) , size(0)
    { }

  inline AvlTree( const AvlTree & rhs ) : root( NULL ), size(rhs.size)
    {
      *this = rhs;
    }

 AvlTree( vector<int>& input_vector ) : root(NULL) ,size(0)
    {
      int size=input_vector.size();
      for(int i=0; i<size ; i++)
	{
	  insert( input_vector[i] );
	}
      //printTree();
    }

  inline ~AvlTree( )
  {
    makeEmpty( );
  }

  /**
   * Find the smallest item in the tree.
   * Throw UnderflowException if empty.
   */
  //const point & findMin( ) const
  const int & findMin( ) const
  {
    if( isEmpty( ) )
      throw UnderflowException( );
    return findMin( root )->element_idx;
  }

  /**
   * Find the largest item in the tree.
   * Throw UnderflowException if empty.
   */
  // const point & findMax( ) const
  const int & findMax( ) const
  {
    if( isEmpty( ) )
      throw UnderflowException( );
    return findMax( root )->element_idx;
  }

  /**
   * Returns true if idx is found in the tree.
   */
  bool contains( int & idx ) const
  {
    return contains( idx, root );
  }

  /**
   * Test if the tree is logically empty.
   * Return true if empty, false otherwise.
   */
  inline bool isEmpty( ) const
  {
    return root == NULL;
  }

  /**
   * Print the tree contents in sorted order.
   */
  void printTree( ) const
  {
    if( isEmpty( ) )
      cout << "Empty tree" << endl;
    else
      {
	/* cout << endl; */
	/* cout << "tree size " << size << endl; */
	/* printTree( root ); */
      }
  }

  /**
   * Make the tree logically empty.
   */
  void makeEmpty( )
  {
    makeEmpty( root );
  }

  /**
   * Insert idx into the tree; duplicates are ignored.
   */
  void insert( int idx )
  {
    insert( idx, root );
  }
     
  /**
   * Remove idx from the tree. Nothing is done if idx is not found.
   */
  void remove( int & idx )
  {
    //cout << "Sorry, remove unimplemented; " << idx <<
    //       " still present" << endl;
  }

  /**
   * Deep copy.
   */
  const AvlTree & operator=( const AvlTree & rhs )
    {
      if( this != &rhs )
        {
	  makeEmpty( );
	  root = clone( rhs.root );
        }
      return *this;
    }

 private:
 public:
  struct AvlNode
  {
    int element_idx;
    vector<int> *dup;
    bool dupFlag;
    AvlNode   *left;
    AvlNode   *right;
    int       height;

  AvlNode( int & theElement_idx, AvlNode *lt,
	   AvlNode *rt, int h = 0 ) 
  : element_idx( theElement_idx ), left( lt ), right( rt ), height( h ), dupFlag(false)
    { dup = new vector<int>; }
  };

  AvlNode *root;
  int size;


  void query(vector<int>& res, AvlNode*& node, double from_lat, double to_lat)
  {
    if(node)
      {
	query(res, node->left, from_lat, to_lat);
	query(res, node->right, from_lat, to_lat);

	if(plist[node->element_idx].y_coord>=from_lat && plist[node->element_idx].y_coord<=to_lat
	   )
	  {
	    res.push_back(node->element_idx);
	    /* cout << "	find  " << setprecision (10) */
	    /* 	 << " " << plist[node->element_idx].longitude */
	    /* 	 << endl; */
	    
	    /* if(plist[node->element_idx].longitude==-87.92355 || plist[node->element_idx].longitude== -87.87437) */
	
	    /* cout << "	find  " << setprecision (10) */
	    /* 	   << " " << plist[node->element_idx].longitude */
	    /* 	   << endl; */
	    
	    while(node->dup != NULL && node->dup->size()!=0) {
	      /* cout << "dup size " <<  node->dup->size()  << endl; */
	      int tmp=node->dup->back();
	      /* if(plist[tmp].y_coord>=from_lat && plist[tmp].y_coord<=to_lat){ */
	      res.push_back(tmp);
	      /* cout << "	find in dup " << setprecision (10) */
	      /* 	   << " " << plist[tmp].longitude */
	      /* 	   << " " << plist[tmp].y_coord */
	      /* 	   << " " << plist[tmp].cday */
	      /* 	   << endl; */
	      
	      /* cout << setprecision (10) << " " << plist[node->element_idx].longitude */
	      /* 	   << " " << plist[node->element_idx].y_coord */
	      /* 	   << " " << plist[node->element_idx].cday */
	      /* 	   << endl; */
	      node->dup->pop_back();
	    }
	  }
      }
    
    /* else */
    /*   { */
    /* 	} */
    /*   } */
  }
  /**
   * Internal method to insert into a subtree.
   * idx is the item to insert.
   * t is the node that roots the subtree.
   * Set the new root of the subtree.
   */
  void insert( int  idx, AvlNode * & t )
  {
    if( t == NULL )
      {
	t = new AvlNode( idx, NULL, NULL );
      }
    else if( plist[idx].y_coord < plist[t->element_idx].y_coord )
      {
	insert( idx, t->left );
	if( height( t->left ) - height( t->right ) == 2 )
	  {
	    if( plist[idx].y_coord < plist[t->left->element_idx].y_coord )
	      rotateWithLeftChild( t );
	    else
	      doubleWithLeftChild( t );
	  }
      }
    else if( plist[t->element_idx].y_coord < plist[idx].y_coord )
      {
	insert( idx, t->right );
	if( height( t->right ) - height( t->left ) == 2 )
	  {
	    if( plist[t->right->element_idx].y_coord < plist[idx].y_coord )
	      rotateWithRightChild( t );
	    else
	      doubleWithRightChild( t );
	  }
      }
    else
      {
	t->dupFlag=true;
	t->dup->push_back(idx);
	/* cout << t->dup->size() << endl; */
	/* cout << "over; " << endl; */
      }
    t->height = max( height( t->left ), height( t->right ) ) + 1;
  }

  /**
   * Internal method to find the smallest item in a subtree t.
   * Return node containing the smallest item.
   */
  AvlNode * findMin( AvlNode *t ) const
  {
    if( t == NULL )
      return NULL;
    if( t->left == NULL )
      return t;
    return findMin( t->left );
  }

  /**
   * Internal method to find the largest item in a subtree t.
   * Return node containing the largest item.
   */
  AvlNode * findMax( AvlNode *t ) const
  {
    if( t != NULL )
      while( t->right != NULL )
	t = t->right;
    return t;
  }


  /**
   * Internal method to test if an item is in a subtree.
   * idx is item to search for.
   * t is the node that roots the tree.
   */
  bool contains( const int & idx, AvlNode *t ) const
  {
    if( t == NULL )
      return false;
    else if( plist[idx].y_coord < plist[t->element_idx].y_coord )
      return contains( idx, t->left );
    else if( plist[t->element_idx].y_coord < plist[idx].y_coord )
      return contains( idx, t->right );
    else
      return true;    // Match
  }
  /****** NONRECURSIVE VERSION*************************
	  bool contains( const point & idx, AvlNode *t ) const
	  {
	  while( t != NULL )
	  if( idx < t->element_idx )
	  t = t->left;
	  else if( t->element_idx < idx )
	  t = t->right;
	  else
	  return true;    // Match

	  return false;   // No match
	  }
  *****************************************************/

  /**
   * Internal method to make subtree empty.
   */
  void makeEmpty( AvlNode * & t )
  {
    if( t != NULL )
      {
	makeEmpty( t->left );
	makeEmpty( t->right );
	delete t;
      }
    t = NULL;
  }

  /**
   * Internal method to print a subtree rooted at t in sorted order.
   */
  void printTree( AvlNode *t ) const
  {
    if( t != NULL )
      {
	printTree( t->left );
	/* cout << t->element_idx <<endl; */
	printTree( t->right );
      }
  }

  /**
   * Internal method to clone subtree.
   */
  AvlNode * clone( AvlNode *t ) const
  {
    if( t == NULL )
      return NULL;
    else
      return new AvlNode( t->element_idx, clone( t->left ), clone( t->right ), t->height );
  }
  // Avl manipulations
  /**
   * Return the height of node t or -1 if NULL.
   */
  int height( AvlNode *t ) const
  {
    return t == NULL ? -1 : t->height;
  }

  int max( int lhs, int rhs ) const
  {
    return lhs > rhs ? lhs : rhs;
  }

  /**
   * Rotate binary tree node with left child.
   * For AVL trees, this is a single rotation for case 1.
   * Update heights, then set new root.
   */
  void rotateWithLeftChild( AvlNode * & k2 )
  {
    AvlNode *k1 = k2->left;
    k2->left = k1->right;
    k1->right = k2;
    k2->height = max( height( k2->left ), height( k2->right ) ) + 1;
    k1->height = max( height( k1->left ), k2->height ) + 1;
    k2 = k1;
  }

  /**
   * Rotate binary tree node with right child.
   * For AVL trees, this is a single rotation for case 4.
   * Update heights, then set new root.
   */
  void rotateWithRightChild( AvlNode * & k1 )
  {
    AvlNode *k2 = k1->right;
    k1->right = k2->left;
    k2->left = k1;
    k1->height = max( height( k1->left ), height( k1->right ) ) + 1;
    k2->height = max( height( k2->right ), k1->height ) + 1;
    k1 = k2;
  }

  /**
   * Double rotate binary tree node: first left child.
   * with its right child; then node k3 with new left child.
   * For AVL trees, this is a double rotation for case 2.
   * Update heights, then set new root.
   */
  void doubleWithLeftChild( AvlNode * & k3 )
  {
    rotateWithRightChild( k3->left );
    rotateWithLeftChild( k3 );
  }

  /**
   * Double rotate binary tree node: first right child.
   * with its left child; then node k1 with new right child.
   * For AVL trees, this is a double rotation for case 3.
   * Update heights, then set new root.
   */
  void doubleWithRightChild( AvlNode * & k1 )
  {
    rotateWithLeftChild( k1->right );
    rotateWithRightChild( k1 );
  }
};

#endif
