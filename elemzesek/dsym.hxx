#ifndef __dsym_h
#define __dsym_h

#include "ddiag.hxx"
#include "mxfunction.hxx"
#include <iostream>
//#include <fstream>
#include <list>
//#include <algorithm>
//#include <set>
//#include <boost/math/common_factor.hpp>
//#include <math.h>
//#include <sstream>
//#define PI 3.14159265

using namespace std;

/*
   Class: Dsym
   The Dsym class represents a D-symbol. A D-symbol consists of a D-diagram and
   a matrix valued function, which is defined on the vertices of the D-diagram.

   FIXME Theory
   */
class Dsym : public Ddiag{
  private:
    /* Variables: Buffer variables
       buf_symmetries - Buffer for <Dsym::symmetries()>.
       buf_dual_sym - Buffer for <Dsym::dual()>.
       buf_cancel_operation_sym - Buffer for <Dsym::cancel_operation_sym()>.
       Mmx - Matrix function M.
       */
    int buf_symmetries;
    Dsym* buf_dual_sym;
    vector<list<Dsym*>*> buf_cancel_operation_sym;
    Mxfunction* Mmx;

  public:
    Dsym(int,int);
    Dsym(Dsym*);
    Dsym(Ddiag*);
    Dsym(Ddiag*, Mxfunction*);
    ~Dsym(void);

    /* Func: get_Mmx()
       This is the M matrix function, which can be generated from the parameters
       and their values, or can be preset.

       FIXME? Theory
       */
    Mxfunction *get_Mmx(void);

    /* Func: update_Mmx
       Set the values of M matrix function from the parameters and their values.
       */
    void update_Mmx(void);

    /* Func: dual
       Generate dual D-symbol which has the adjacency operations reversed:
       Operation number 0 becomes number <Base.dim>, operation number 1 becomes
       <Base.dim>-1 and so on.
       */
    Dsym* dual_sym(void);

    /* Func: check_dual
       Check if the dual symbol is smaller.

       Returns:                                                                                        
       -1 - The dual is smaller                                                                        
       0 - Selfdual
       1 - The dual is bigger
       */
    int check_dual(void);

    /* Func: check_numberings
       Check if there are alternate good numberings which give us a smaller
       symbol. Meanwhile we fill up <Dsym.symmetries> variable.

       Returns:
       1 - There are smaller possible numberings
       0 - We have the best numbering, but there are multiple best numberings
       -1 - We have found the one and only best numbering
       */
    int check_numberings(void);

    /* Func: symmetries()
       Is the D-symbol symmetric?

       Returns:
       2+ - The D-symbol is symmetric, and there are this many starting
        vertices which give the same symbol.
       1 - The D-symbol is not symmetric.
       */
    int symmetries(void);
    
    /* Func: check_m
       Check if matrix function M is sane.

       Check if matrix function M is sane: 
       - Main diagonal has only 1's
       - Every element at least 2 away from the main diagonal is exactly 2
       - Every element beside the main diagonal is at least 2
       */
    int check_m(void);

    /*
       Func: cancel_operation_sym()
      
       If we cancel the nth operation in the symbol, we can get some D-symbol
       components in a smaller dimension.
       */
    list<Dsym*> *cancel_operation_sym(int);

    /*
       Func: is_smaller()
       Is this D-symbol with first vertice *i* smaller than the other D-symbol
       with first vertice *j*?

       Algorithm: FIXME

       Returns:
       1 - If smaller
       0 - If equal
       -1 - If bigger

       Remarks:
       FIXME The used numberings should be generated beforehand.
       */
    int is_smaller(int i,Dsym* other,int j);

    /* Func: check_simplex_vertices()
       Check that tilings around barycentric-simplex vertices give us spherical
       or euclidean planes. Barycentric-simplex vertices at the vertices or the
       body centers of the fundamental domain can be real or at the edge of the
       model, so the tilings around them can be spherical or euclidean. The
       other barycentric-simplex vertices must be real, so the tilings around
       them must be spherical.

       In 2 dimensions one can check that a tiling is spherical, euclidean or
       hyperbolic using the combinatorical curvature function. In spherical case
       one has to check, there aren't any bad orbifolds.

       FIXME If the combinatorical curvature function is positive we can still
       get a torus. We should check this too...

	Returns:
        2 - Every barycentric-simplex vertice is as defined, but there are ideal
	  vertices (with euclidean tiling)
        1 - Every barycentric-simplex vertice is as defined, without ideal
	  vertices (every tiling is spherical)
	0 - There are bad barycentric-simplex vertices
       */
    int check_simplex_vertices(void);

    /* Func: combinatoric_curvature()
       Combinatoric curvature of a 2 dimensional D-symbol

       FIXME Theory
       */
    float combinatoric_curvature(void);

    /* Func: is_spherical()
       If we have a 2 dimensional D-symbol we can check if it's tiling is
       spherical.

	Returns:
	1 - We have a 2 dimensional spherical tiling.
	0 - We have a 2 dimensional non-spherical tiling.
	-1 - The dimension is not 2.
	*/
    int is_spherical(void);

    /* Func: is_euclidean()
       If we have a 2 dimensional D-symbol we can check if it's tiling is
       euclidean.

	Returns:
	1 - We have a 2 dimensional euclidean tiling.
	0 - We have a 2 dimensional non-euclidean tiling.
	-1 - The dimension is not 2.
	*/
    int is_euclidean(void);

    /* Func: is_hyperbolic()
       If we have a 2 dimensional D-symbol we can check if it's tiling is
       hyperbolic.

	Returns:
	1 - We have a 2 dimensional euclidean tiling.
	0 - We have a 2 dimensional non-euclidean tiling.
	-1 - The dimension is not 2.
	*/
    int is_hyperbolic(void);

    /* Func: dump()
       Dump Dsym in machine readable format.

       Output looks like this: 3 6 ((1,2)(3)(4)(5)(6)) ((1,4)(2,3)(5)(6))
       ((1,2)(3)(4,5)(6)) ((1,2)(3,4)(5,6)) 3 6 1 2 3 ... 6*4^2

       The first part is the D-diagram, the second is the M matrix valued
       function.
       */
    virtual int dump(ostream*);

    /* Func: print_html()
       Pretty print Dsym
       */
    virtual int print_html(ostream*);
};

#endif /* __dsym_h */
