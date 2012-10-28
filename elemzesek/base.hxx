#ifndef __base_h
#define __base_h

#include <iostream>

using namespace std;

/*
  Class: Base

  The base of every class. It contains some basic information: dimension and
  cardinality. And it has the basic funcions that should be defined in every
  class: <dump>, <print_html>
*/
class Base {
  protected:
    /*
      Variables: Base variables

      dim - Dimension
      car - Cardinality
    */
    int dim,car;
  public:
    /*
      Constructor: Base
      Initialize Base object without parameters. It should be filled later.
    */
    Base(void);
    /*
      Constructor: Base
      Initialize Base object's dimension, cardinality.
    */
    Base(int dim, int car);
    /* 
      Constructor: Base
      Restore data dumped with function dump from istream.
    */
    Base(istream*);
    /*
      Destructor: ~Base
      Virtual destructor of Base object, so destructor works fine if used for
      virtual objects.
    */
    virtual ~Base(void);
    /* 
      Func: dump()
      Dump data to ostream in computer readable format.

      Output format is like "3 6", dimension and cardinality.
    */
    virtual int dump(ostream*);
    /* 
      Func: print_html()
      Pretty print the class in xhtml format.
    */
    virtual int print_html(ostream*);
};

#endif /* __base_h */
