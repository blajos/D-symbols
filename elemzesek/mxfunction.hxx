#ifndef __mxfunction_h
#define __mxfunction_h

#include "base.hxx"
#include "simplex.hxx"
#include "param.hxx"
#include <list>
#include <iostream>

/* Class: Mxfunction
   Class for a matrix valued function defined on simplex orbits.
   */
class Mxfunction : public Base {
  private:
    /* Variable: mx
       3 dimensional array. First dimension is the simplex orbit's number in the
       initial (maybe wrong) numbering. The other two dimensions represent the
       matrix's rows and columns.
       */
    int*** mx;
  public:
    /* Constructor: Mxfunction
       Create initial Mxfunction.
       */
    Mxfunction(int,int);

    /* Constructor: Mxfunction
       Copy constructor
       */
    Mxfunction(Mxfunction*);

    /* Constructor: Mxfunction
       Create Mxfunction from parameters. FIXME Kell ez?
       */
    Mxfunction(int,int,list<Param*>*);

    /* Constructor: Mxfunction
       Restore Mxfunction from <Mxfunction::dump()>
       */
    Mxfunction(istream*);

    /* Destructor: ~Mxfunction
       Reclaim memory.
       */
    ~Mxfunction(void);

    /* Funcs: Simple functions
       get() - Get the matrix element *x*,*y* from *orbit*'s matrix.
       set() - Set the matrix element *x*,*y* to *z* in *orbit*'s matrix.
       */
    int get(Simplex* orbit,int x,int y);
    int set(Simplex* orbit,int x,int y,int z);

    /* Funcs: Simple functions
       get() - Get the matrix element *x*,*y* from the *k*th matrix.
       set() - Set the matrix element *x*,*y* to *z* in the *k*th matrix.
       */
    int get(int k,int x,int y);
    int set(int k,int x,int y,int z);

    /* Func: dump
       Dump Mxfunction in machine readable format. 
       
       Output looks like this: 3 6 1 2 3 ... 6*(3+1)^2
       The first two numbers are dimension and cardinality, the rest is every
       element of the matrix function: M(D1,1,1) M(D1,1,2) M(D1,1,3) M(D1,1,4)
       M(D1,2,1) ... M(D1,4,4) M(D2,1,1) ... M(D6,4,4)
       */
    virtual int dump(ostream*);

    /* Func: print_html
       Pretty print Mxfunction.
       */
    virtual int print_html(ostream*);
};

#endif /* __mxfunction_h */
