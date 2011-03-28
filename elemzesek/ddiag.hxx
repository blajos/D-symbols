#ifndef __ddiag_h
#define __ddiag_h

#include "base.hxx"
#include "simplex.hxx"
#include "svg.hxx"
#include "param.hxx"
#include "mxfunction.hxx"
#include <iostream>
//#include <fstream>
#include <list>
#include <vector>
#include <algorithm>
//#include <set>
//#include <boost/math/common_factor.hpp>
#include <math.h>
//#include <sstream>
//#define PI 3.14159265

using namespace std;

/*
   Class: Ddiag
   The Ddiag class represents a D-diagram. A D-diagram is a <Base.dim>+1 colored
   multigraph, with <Base.car> vertices. Colors represent the different types of
   adjacency operations. Vertices represent the tiling's simplex orbits.

   FIXME Theory
   */
class Ddiag : public Base {
  private:
    /* Variables: Buffer variables
       buf_symmetries - Buffer for <Ddiag::symmetries()>.
       buf_dual_diag - Buffer for <Ddiag::dual_diag()>.
       buf_params - Buffer for <Ddiag::params()>.
       buf_cancel_operation_diag - Buffer for <Ddiag::cancel_operation_diag()>.
       buf_Rmx - Buffer for <Ddiag::Rmx()>.
       */
    int buf_symmetries;
    Ddiag* buf_dual_diag;
    list<Param*> *buf_params;
    vector<list<Ddiag*> *> buf_cancel_operation_diag;
    Mxfunction* buf_Rmx;

    /* Variable: simplex_orbits
       We use <Base.car> different Simplex* pointers for vertices/simplex
       orbits. 
       
       There are <Base.car>+1 different numberings (unsorted, and of
       revery vertex we have a numbering started at that vertex.) For every
       numbering we need a <Base.car> sized array consisting of Simplex*
       pointers.
       */
    Simplex*** simplex_orbits;

    /* Func: filter_bad_orbifolds
       Only in 3 dimensions. Filter possibly bad orbifolds introduced by
       adjacency operations 1 and 2.

	Returns:
	1 - Parameters were equalized or set to a fixed value
	0 - Nothing has happened
	-1 - Error (ex. the dimension is not 3)
       */
    int filter_bad_orbifolds(list<Param*>*);
    //int filter_bad_orbifolds03(void); Szerintem nem kell kulon...
  public:
    /*
      Constructor: Ddiag
      We define some sane defaults for a D-diagram object. Dimension and
      cardinality are parameters. We allocate the memory for the matrix
      functions, create the simplex orbits and initialize them with no real
      adjacency operations. We only use the 0 indexes, because the others depend
      on the adjacency operations.
      */
    Ddiag(int dim,int car);
    
    /*
      Constructor: Ddiag
      Copy constructor.
      */
    Ddiag(Ddiag*);
    
    /* Constructor: Ddiag
       Restore Ddiag from the machine readable format dump produces.
       */
    Ddiag(istream*);
    
    /*
      Destructor: ~Ddiag
      Reclaim allocated memory.
      */
    ~Ddiag(void);
    
    /* Func: create_edge
       Create a *color* colored edge from *from* to *to*, if possible.

       We can create the edge, if neither *from* nor *to* has a *color*
       colored edge.

       Returns:
       1 - If the edge was created
       0 - If the edge couldn't be created
       */
    int create_edge(int color,int from,int to);
    
    /* Func: remove_edge
       Remove the *color* colored edge from *from* to *to*, if there is an edge.
       */
    void remove_edge(int color,int from,int to);
    
    /* Func: create_numbering
       Create numbering that begins with *first* simplex orbit.

       Algorithm: FIXME
       */
    void create_numbering(Simplex* first);
    
    /* Func: Rmx
       This is the R matrix function, which can be generated from the D-diagram.

       FIXME? Theory
       */
    Mxfunction *Rmx(void);
    
    /* Func: dual
       Generate dual D-diagram which has the adjacency operations reversed:
       Operation number 0 becomes number <Base.dim>, operation number 1 becomes
       <Base.dim>-1 and so on.
       */
    Ddiag* dual_diag(void);

    /*
       Func: cancel_operation_diag()
      
       If we cancel the *i*th operation in the diagram, we can get some D-diagram
       components in a smaller dimension.
       */
    list<Ddiag*> *cancel_operation_diag(int i);

    /* Func: check_dual
       Check if the dual diagram is smaller.

       Returns:
       -1 - The dual is smaller
       0 - Selfdual
       1 - The dual is bigger
       */
    virtual int check_dual(void);

    /* Func: check_numberings
       Check if there are alternate good numberings which give us a smaller
       diagram. Meanwhile we fill up <Ddiag.symmetries> variable.

       Returns:
       1 - We have found the one and only best numbering
       0 - We have the best numbering, but there are multiple best numberings
       -1 - There are smaller possible numberings
       */
    virtual int check_numberings(void);

    /* Func: symmetries()
       Is the D-diagram symmetric?

       Returns:
       2+ - The D-diagram is symmetric, and there are this many starting
        vertices which give the same diagram.
       1 - The D-diagram is not symmetric.
       */
    int symmetries(void);
    
    /* Func: check_r
       Check if matrix function R is sane.

       Check if matrix function R is sane: 
       - Main diagonal has only 1's
       - Every element at least 2 away from the diagonal is not more than 2

	Returns:
	1 - R is sane
	0 - R is wrong
       */
    int check_r(void);

    /*
       Func: is_bigraph()
      
       Is the graph a bipartite graph? In 2 dimensions a bipartite graph shows
       an orientable orbifold. (FIXME We need a citeable theorem... How about
       loops?)
    */
    int is_bigraph(void);

    /*
       Func: is_smaller_diag()
       Is this D-diagram with first vertice *i* smaller than the other D-diagram
       with first vertice *j*?

       Algorithm: FIXME

       Returns:
       1 - If smaller
       0 - If equal
       -1 - If bigger

       Remarks:
       FIXME The used numberings should be generated beforehand.
       */
    int is_smaller_diag(int i,Ddiag* other,int j);

    /* Func: params()
       Generate parameters based on the simplex orbits, and filter bad orbifolds
       with making parameters equal or locking their value. It is still possible
       to get bad orbifolds in D-symbols because one can find parameters which
       don't contribute to a bad orbifold with the smallest value only with
       higher values.
       */
    list<Param*> *params(void);

    /* Func: dump
       Dump Ddiag in machine readable format. 
       
       Output looks like this: 3 6 ((1,2)(3)(4)(5)(6)) ((1,4)(2,3)(5)(6))
       ((1,2)(3)(4,5)(6)) ((1,2)(3,4)(5,6))
       The first two numbers are dimension and cardinality, the sequences are
       each adjacency operations as involutive permutations.
       */
    virtual int dump(ostream*);
    
    /* Func: print_html
       Pretty print Ddiag.
       */
    virtual int print_html(ostream*);
    friend class Dsym;
};

#endif /* __ddiag_h */
