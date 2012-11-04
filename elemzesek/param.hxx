#ifndef __param_h
#define __param_h

#include "simplex.hxx"
#include <string>
#include <list>

class Simplex;

/* Class: Param
   Class for representing parameters in M matrix function.
   */
class Param {
  public:
  /* Variables: Parameter variables
     coeff - Coefficient of parameter.
     min - Smallest possible value.
     val - Current value.
     letter - Which letter represents the parameter.
     changeable - Is it changeable?
     orientable - Is it orientable (without loops in the diagram)?
     simplex_operations - Which simplex and operation pairs does this parameter
       belong to?
     */
  char letter;
  int coeff,min,val;
  bool orientable;
  bool changeable;

  struct sop {
    int op1,op2;
    Simplex* simplex;
  };
  list<sop> simplex_operations;

  /* Constructor: Param
     Create param, with values for *letter*, *coeff*, *orientable*.
   */
  Param(char letter, int coeff, bool orientable, bool changeable);
  /* Destructor: Param
     Reclaim memory.
   */
  ~Param(void);

  /* Funcs: a

     set_param() - Set parameter value to an integer
     change_param() - Change parameter value by an integer (FIXME watch out for
     smallest possible value.)
     increase_param() - Alias for change_param(1).
     decrease_param() - Alias for change_param(-1).
   */
  void set_param(int);
  void change_param(int);
  void increase_param();
  void decrease_param();
  //bool is_equal_to(Param); Nem latom, hol kell
  int check_min();
};

#endif /* __param_h */
