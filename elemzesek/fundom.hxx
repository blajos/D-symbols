#ifndef __fundom_h
#define __fundom_h

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

  /* Adatstruktura (3dim):
     - minden graf csucshoz tartozik egy tetraeder (class tetrahedron)
     - minden tetraedernek 4 lapja van (operacionkent) (class facet)
     - a lapok kozul van ami latszik es van ami menten osszeragasztottunk 2
     tetraedert
     - minden tetraedernek 4 csucsa van laponkent 3-3 (class point)
     - Uj kapcsolat letrehozasa 1 uj csucsot jelent, amit konvex modon be tudunk
     rakni
     - erdekes meg a lapok szine: (from,to,operation) harmasok
     */

/* Class: Fundom
   Fundamental domain object, for generating the fundamental domain of a
   D-diagram or a D-symbol. 
   */
class Fundom{
  /* Variables:
	dim - Dimension
	car - Number of components
	tetrahedra - Vector of tetrahedra (one tetrahedron for each diagram
	vertex)
     */
  int dim,car;
  vector<Tetrahedron*> tetrahedra;

  public:
  /* Constructor: Fundom
     Initialize Fundom object and setup first tetrahedron.
     */
  Fundom(int dim,int car);

  /* Destructor: ~Fundom
     Destroy Fundom object.
     */
  ~Fundom(void);

  /* Func: add_component()
     Glue new tetrahedron *to* to old tetrahedron *from* at face *operation*.
     */
  void add_component(int from,int to,int operation);

  /* Func: check_convexity()
     Check if the polyhedron would be convex if we add a new tetrahedron to
     *from* at face *operation*, with the fourth point *newpoint*.

     Returns:
       1 - if convex
       0 - if we can't decide due to floating point errors
       -1 - if concave
   */
  int check_convexity(int from,int operation,Point newpoint);

  /* Func: print_geomview()
     Print geomview compatible information to output.
     */
  int print_geomview(ostream*);
};

/* Class: Tetrahedron
   General n-dimensional tetrahedron class for generating fundamental domains.

     FIXME Simplex osztalyba nem lehet betenni? (Talan nem kene eroltetni, mert
     az simplex orbitokrol szol...)
   */
class Tetrahedron{
  /* Variables:
     dim - Dimension
     points - Vector of the points of the tetrahedron
     facets - List/vector of facets
  */
  int dim;
  vector<Point*> points;
  list<Facet*> facets;
};

/* Class: Facet
   *dim-1* dimensional facet of a Tetrahedron. It has points, a color
   ((from,to,operation) tuple), a normal vector and is possibly invisible.
   */
class Facet{
  /* Variables:
     dim - Dimension
     points - Vector of the points of the facet
     color - Color of the facet ((from,to,operation) tuple)
     normal - Normal vector (Point type)
     visible - Is visible?
  */
  int dim;
  vector<Point*> points;
  vector<int> color;
  Point normal;
  bool visible;
};

/* Class: Point
   *dim* dimensional vector or point.

   FIXME Kell ez? Siman vetor<int> is eleg lenne...
   */
class Point{
  /* Variables:
     dim - Dimension
     coordinates - Coordinates
     */
  int dim;
  vector<double> coordinates;
};

#endif /* __fundom_h */
