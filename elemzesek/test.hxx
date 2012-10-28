#ifndef __test_h
#define __test_h

#include <sstream>
#include "dsym.hxx"
#include "dlist.hxx"
#include "algorithms.hxx"

using namespace std;

/*
  Class: Test

  The Test class contains basic module tests. And for algorithms we can check
  for already known results.
*/
class Test{
  public:
    Test(void);
    ~Test(void);

    int all(void);

    int base(void);
    int simplex(void);
    //int fundom(void);
    int mxfunction(void);
    int param(void);
    int svg(void);
    int ddiag(void);
    int dsym(void);
    int dlist(void);
    int algorithms(void);
};

#endif /* __test_h */

