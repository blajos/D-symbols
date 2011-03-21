#ifndef __svg_h
#define __svg_h

#include "base.hxx"
#include <iostream>
#include <string>
//#include <fstream>
#include <list>
#include <vector>
//#include <algorithm>
//#include <set>
//#include <boost/math/common_factor.hpp>
#include <math.h>
//#include <sstream>
#define PI 3.14159265

using namespace std;

/* Class: Svg
   Svg object for writing a D-diagram's SVG graphics code some output.
   */
class Svg : public Base{
  /* Variables: Geometric variables
     koordx - Array of x coordinates.
     koordy - Array of y coordinates.
     rad - Radius of circles.
     size - Size of canvas.
     fontsize - Size of fonts.
     */
  int *koordx,*koordy,rad,size,fontsize;
  list<vector<int> > lines;

  public:
  /* Constructor: Svg
     Initialize Svg object.
     */
  Svg(int dim,int car);
  /* Destructor: ~Svg
     Destroy Svg object.
     */
  ~Svg(void);
  /* Func: create_circle()
     Create the *n*'th circle.
     */
  void create_circle(int n);
  /* Func: create_numtext()
     Create the *n*'th circle's label.
     */
  void create_numtext(int n);
  /* Func: create_line()
     Create a line between the *n0*'th and the *n1*'th circle using pattern
     *color*.
   */
  void create_line(int n0,int n1,int color);
  int add_line(int n0,int n1,int color);
  int print_html(ostream*);
};

#endif /* __svg_h */
