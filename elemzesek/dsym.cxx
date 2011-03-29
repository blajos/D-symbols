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
  delete buf_dual_sym;
  delete Mmx;
  for(vector<list<Dsym*>*>::iterator it1=buf_cancel_operation_sym.begin();
      it1!=buf_cancel_operation_sym.end();it1++){
    for(list<Dsym*>::iterator it2=(*it1)->begin(); it2!=(*it1)->end(); it2++)
      delete *it2;
    delete *it1;
  }
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
  for(int i=0;i<car;i++)
    for(int j=0;j<dim+1;j++)
      for(int k=0;k<dim+1;k++)
        if((j==k && Mmx->get(i,j,k)!=1) ||
	   (abs(j-k)==1 && Mmx->get(i,j,k)<2) ||
	   (abs(j-k)>1 && Mmx->get(i,j,k)!=2))
	   return 0;
  return 1;
}

list<Dsym*>* Dsym::cancel_operation_sym(int cancel_op){
  if(cancel_op<0 or cancel_op > dim+1){
    throw "Operation out of range";
  }
  buf_cancel_operation_sym.resize(dim+2,0);
  if(buf_cancel_operation_sym[cancel_op]==0){
    buf_cancel_operation_sym[cancel_op]=new list<Dsym*>;
    list<Ddiag*>* cancelled_diags=cancel_operation_diag(cancel_op);
    for(list<Ddiag*>::iterator it=cancelled_diags->begin();
	it!=cancelled_diags->end(); it++){
      Dsym* new_element=new Dsym((*it));
      buf_cancel_operation_sym[cancel_op]->push_back(new_element);
      // A cancel operation hogyan tudja megtartani a 0. szamozast?
      // Uj valtozo a Simplex osztalyban: original_label
    }
  }
  return buf_cancel_operation_sym[cancel_op];
}

int Dsym::is_smaller_sym(int thisindex,Dsym* other,int otherindex){
  if(is_smaller_diag(thisindex,(Ddiag*)other,otherindex))
    return is_smaller_diag(thisindex,(Ddiag*)other,otherindex);
  else {
    for(int i=0;i<dim+1;i++)
      for(int j=0;j<dim+1;j++)
	for(int k=0;k<car;k++)
	  if(Mmx->get(simplex_orbits[thisindex][k],i,j) <
	      other->get_Mmx()->get(other->simplex_orbits[otherindex][k],i,j))
	    return 1;
	  else if(Mmx->get(simplex_orbits[thisindex][k],i,j) >
	      other->get_Mmx()->get(other->simplex_orbits[otherindex][k],i,j))
	    return -1;
  }
  return 0;
}

int Dsym::check_simplex_vertices(void){
  if(dim!=3)
    throw "Error: dimension is not 3.";
  for(int i=0;i<dim+1;i++)
    for(list<Dsym*>::iterator it=cancel_operation_sym(i)->begin();
    it!=cancel_operation_sym(i)->end(); it++){
      if(i==0 || i==3){
	if((*it)->is_spherical()==0 && (*it)->is_euclidean()==0)
	  return 0;
      }
      else
	if((*it)->is_spherical()==0)
	  return 0;
    }
  return 1;
}

float Dsym::combinatoric_curvature(void){
  if(dim!=2)
    throw "Error: dimension is not 2.";
  float curvature=0;
  for(int i=0;i<car;i++){
    for(int j=0;j<dim;j++)
      curvature+=1.0/(float)Mmx->get(i,j,j+1);
    curvature-=0.5;
  }
  if(abs(curvature<0.0000001))
    return 0.0;
  else
    return curvature;
}

int Dsym::is_spherical(void){
  if(combinatoric_curvature()>0){
    list<Param*> real_params;
    for(list<Param*>::iterator it=params()->begin(); it!=params()->end(); it++)
      if((*it)->val > 1)
        real_params.push_back(*it);
    if(real_params.size() == 1)
      return 0;
    else if (real_params.size() == 2 && 
	(*real_params.begin())->orientable == (*real_params.begin()++)->orientable)
      return 0;
    return 1;
  }
  return 0;
}

int Dsym::is_euclidean(void){
  if(combinatoric_curvature()==0)
    return 1;
  return 0;
}

int Dsym::is_hyperbolic(void){
  if(combinatoric_curvature()<0)
    return 1;
  return 0;
}

int Dsym::dump(ostream* out){
  Ddiag::dump(out);
  *out << " ";
  Mmx->dump(out);
  return 0;
}

int Dsym::print_html(ostream* out){
  Ddiag::print_html(out);
  *out << "    <p>M matrix function:</p>" << endl;
  Mmx->print_html(out);
  return 0;
}
