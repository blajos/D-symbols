#include "dsym.hxx"

Dsym::Dsym(int dim,int car):
  Ddiag(dim,car)
{
  buf_dual_sym=0;
  buf_Mmx=0;
  //create default Mmx
}

Dsym::Dsym(Dsym* oldsym):
  Ddiag((Ddiag*)oldsym)
{
  buf_Mmx=new Mxfunction(oldsym->Mmx());
}

Dsym::Dsym(Ddiag* olddiag):
  Ddiag(olddiag)
{
  buf_dual_sym=0;
  buf_Mmx=buf_Mmx=new Mxfunction(dim,car);
}

Dsym::Dsym(Ddiag* olddiag, Mxfunction* oldmmx):
  Ddiag(olddiag)
{
  buf_Mmx=new Mxfunction(oldmmx);
}

Dsym::~Dsym(void){
  return;
}

Mxfunction* Dsym::Mmx(void){
  return buf_Mmx;
}

void Dsym::update_Mmx(void){
  return;
}

Dsym* Dsym::dual_sym(void){
  return buf_dual_sym;
}

int Dsym::check_dual(void){
  return 0;
}

int Dsym::check_numberings(void){
  return 0;
}

int Dsym::symmetries(void){
  return 0;
}

int Dsym::check_m(void){
  return 0;
}

list<Dsym*>* Dsym::cancel_operation_sym(int i){
  return buf_cancel_operation_sym[i];
}

int Dsym::is_smaller(int i,Dsym* other,int j){
  return 0;
}

int Dsym::check_simplex_vertices(void){
  return 0;
}

double Dsym::combinatoric_curvature(void){
  return 0;
}

int Dsym::is_spherical(void){
  return 0;
}

int Dsym::is_euclidean(void){
  return 0;
}

int Dsym::is_hyperbolic(void){
  return 0;
}

int Dsym::dump(ostream*){
  return 0;
}

int Dsym::print_html(ostream*){
  return 0;
}
