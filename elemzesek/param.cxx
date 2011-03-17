#include "param.hxx"

Param::Param(char letter, int coeff, bool orientable):
  letter(letter),
  coeff(coeff),
  orientable(orientable)
{
  min=0;
  while (coeff*min<2)
    min++;
}

Param::~Param(void){
  return;
  // delete simplex_operations if needed;
}

void Param::set_param(int i){
  val=i;
  /* felesleges?
  for(list<simplex*>::iterator szim=currparam->szek.begin();
      szim!=currparam->szek.end();szim++){
    int op0=currparam->op;
    int op1=op0+1;
    (*szim)->mx[op0][op1]=abs(currparam->eh)*currparam->ertek;
    (*szim)->mx[op1][op0]=abs(currparam->eh)*currparam->ertek;
  }
  */
}

void Param::change_param(int i){
  set_param(val+i);
}

void Param::increase_param(){
  change_param(+1);
}

void Param::decrease_param(){
  change_param(-1);
}
