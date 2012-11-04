#ifndef __simplex_h
#define __simplex_h

#include "base.hxx"
#include "param.hxx"
class Param;

using namespace std;

/*
Class: Simplex
This class represents one vertice of a D-diagram which represents one simplex
orbit in our imaginary tiling. 

Dimension: We have one more adjancency operations. 

Cardinality: This is the number of possible good numberings.
*/
class Simplex : public Base {
  public:
    /* Variable: szomszed
       A Simplex* pointer can be defined to represent an adjacency operation.
       There are <Base.dimension>+1 adjacency operations for every simplex orbit, so
       we use a simple array with <Base.dimension>+1 element.
       */
    Simplex** szomszed;
    /* Variable: sorszam
       We define <Base.cardinality>+1 possible numberings. The first one (with 0
       index) is not necessarily a good numbering. The others are generated good
       numberings.
       */
    int* sorszam;
    /* Variable: params
       Parameters belonging to the simplex. Params[i][j] means the parameter of
       (i,j) adjacency operation pairs.
       */
    Param*** params;
    /* Constructor: Simplex
       Initialize Simplex object with <Base.dimension>, <Base.cardinality> and the index
       in the first (not necessarily a good) numbering. Every adjacency
       operation points to the simplex itself.
       */
    Simplex(int,int,int);
    /* Destructor: ~Simplex
       Destroy Simplex, reclaim memory.
       */
    ~Simplex(void);

    /* Func: dump()
       Save simplex data.
       */
    int dump(ostream*);
};

#endif /* __simplex_h */
