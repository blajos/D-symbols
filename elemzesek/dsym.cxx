#include "dsym.hxx"

Dsym::Dsym(int dim,int car):
  Ddiag(dim,car)
{
  buf_dual_sym=0;
  Mmx=new Mxfunction(dim,car);
}

Dsym::Dsym(Dsym* oldsym):
  Ddiag((Ddiag*)oldsym)
{
  Mmx=new Mxfunction(oldsym->get_Mmx());
}

Dsym::Dsym(Ddiag* olddiag):
  Ddiag(olddiag)
{
  buf_dual_sym=0;
  Mmx=new Mxfunction(dim,car);
}

Dsym::Dsym(Ddiag* olddiag, Mxfunction* oldmmx):
  Ddiag(olddiag)
{
  Mmx=new Mxfunction(oldmmx);
}

Dsym::~Dsym(void){
  return;
}

Mxfunction* Dsym::get_Mmx(void){
  return Mmx;
}

void Dsym::update_Mmx(void){
  return;
}

Dsym* Dsym::dual_sym(void){
  if(!buf_dual_sym){
    buf_dual_sym=new Dsym(dual_diag());
    for(int i=0;i<car;i++)
      for(int j=0;j<dim+1;j++)
	for(int k=0;k<dim+1;k++)
	  buf_dual_sym->Mmx->set(i,j,k,Mmx->get(i,dim-j,dim-k));
  }
  return buf_dual_sym;
}

int Dsym::check_dual(void){
  return -is_smaller_sym(1,dual_sym(),1);
}

int Dsym::check_numberings(void){
  buf_symmetries=1;
  for (int i=2;i<car+1;i++){
    int ret=is_smaller_sym(1,this,i);
    if (ret == -1){
      return -1;
    }
    if (ret == 0){
      buf_symmetries++;
    }
  }
  if (buf_symmetries > 1)                                                                             
    return 0;                                                                                          
  else                                                                                                 
    return 1;
}

int Dsym::symmetries(void){
  if (buf_symmetries == 0)
    check_numberings();
  return buf_symmetries;
}

int Dsym::check_m(void){
  return 0;
}

list<Dsym*>* Dsym::cancel_operation_sym(int i){
  return buf_cancel_operation_sym[i];
}

int Dsym::is_smaller_sym(int thisindex,Dsym* other,int otherindex){
  if(is_smaller_diag(thisindex,(Ddiag*)other,otherindex))
    return is_smaller_diag(thisindex,(Ddiag*)other,otherindex);
  else {
    // MX fuggvenyek...
    return 0;
  }
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
