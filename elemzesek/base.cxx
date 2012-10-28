#include "base.hxx"

Base::Base(void):
  dim(0),
  car(0)
{
  return;
}

Base::Base(int dimin, int carin):
  dim(dimin),
  car(carin)
{
  return;
}

Base::Base(istream* in){
  *in >> dim >> car;
}

Base::~Base(){
  return;
}

int Base::dump(ostream* out){
  *out << dim << " " << car;
  return 0;
}

int Base::print_html(ostream*){
  return 0;
}
