#include "mxfunction.hxx"

Mxfunction::Mxfunction(int dim,int car):
  Base(dim,car)
{
  mx=new int**[car];
  for(int i=0;i<car;i++){
    mx[i]=new int*[dim+1];
    for(int j=0;j<dim+1;j++){
      mx[i][j]=new int[dim+1];
      for(int k=0;k<dim+1;k++)
	mx[i][j][k]=0;
    }
  }
}

Mxfunction::Mxfunction(Mxfunction* oldmx):
  Mxfunction(oldmx->dim,oldmx->car)
{
  //Mxfunction(dim,car);
  for(int i=0;i<car;i++)
    for(int j=0;j<dim+1;j++)
      for(int k=0;k<dim+1;k++)
	mx[i][j][k]=oldmx->get(i,j,k);
}

// FIXME Kell ez?
Mxfunction::Mxfunction(int,int,list<Param*>*):
  Base(dim,car)
{
  return;
}

Mxfunction::Mxfunction(istream* in){
  *in >> dim >> car;
  Mxfunction(dim,car);
  for(int i=0;i<car;i++){
    for(int j=0;j<dim+1;j++){
      for(int k=0;k<dim+1;k++)
	*in >> mx[i][j][k];
    }
  }
}

Mxfunction::~Mxfunction(void){
  for(int i=0;i<car;i++){
    for(int j=0;j<dim+1;j++)
      delete[] mx[i][j];
    delete[] mx[i];
  } 
  delete[] mx;
}

int Mxfunction::dump(ostream *out){
  *out << dim << " " << car;
  for(int i=0;i<car;i++){
    for(int j=0;j<dim+1;j++){
      for(int k=0;k<dim+1;k++)
	*out << " " << mx[i][j][k];
    }
  }
  return 0;
}

int Mxfunction::get(Simplex* orbit,int x,int y){
  if ( x > dim || y > dim)
    throw 42;
  return mx[orbit->sorszam[0]][x][y];
}

int Mxfunction::set(Simplex* orbit,int x,int y,int z){
  if ( x > dim || y > dim)
    throw 42;
  mx[orbit->sorszam[0]][x][y]=z;
  return 0;
}

int Mxfunction::get(int k,int x,int y){
  if ( k > car-1 || x > dim || y > dim)
    throw 42;
  return mx[k][x][y];
}

int Mxfunction::set(int k,int x,int y,int z){
  if ( k > car-1 || x > dim || y > dim)
    throw 42;
  mx[k][x][y]=z;
  return 0;
}

int Mxfunction::print_html(ostream* out){
  *out << "<table border=\"1\"><tr>";
  for(int i=0;i<car;i++){
    *out << "<td><table cellpadding=\"3\">";
    for(int j=0;j<dim+1;j++){
      *out << "<tr>";
      for(int k=0;k<dim+1;k++)
	*out << "<td align=\"center\">" << mx[i][j][k] << "</td>";
      *out << "</tr>";
    }
    *out << "</table></td>";
  }
  *out << "</tr></table>";
  return 0;
}

