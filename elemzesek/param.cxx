#include "param.hxx"

Param::Param(char letter, int coeff, bool orientable, bool changeable):
  letter(letter),
  coeff(coeff),
  orientable(orientable),
  changeable(changeable)
{
  check_min();
}

Param::~Param(void){
  return;
  // delete simplex_operations if needed;
}

void Param::set_param(int i){
  if (changeable){
    val=i;
    if (! check_min())
      val=min;
  }
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

int Param::check_min(){
  min=0;
  while (coeff*min<2)
    min++;
  if (val >= min)
    return 1;
  else
    return 0;
}
