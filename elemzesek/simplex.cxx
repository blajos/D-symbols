#include "simplex.hxx"

//simplex konstruktor: megadjuk a dimenziot, az elemszamot, lefoglaljuk a
//szukseges memoriatartomanyokat a szomszed es a sorszam tomboknek.
//Beallitjuk a 0. variacio sorszamat az adottra. (Ezt a sorszamot a felsorolo
//algoritmus tolti ki, nem feltetlenul ugyanaz, mint az 1-es elemmel valo
//kezdeskor a sorszam.
//A dim+1 darab szomszedsagi operacio eredmenye kezdetben maga a szimplex.
//A parameterek listajat inicializaljuk.
Simplex::Simplex(int dimin,int carin,int ssz):
  Base(dimin,carin)
{
  params=new Param**[dim+1];
  for (int i=0;i<dim+1;i++){
    params[i]=new Param*[dim+1];
    for (int j=0;j<dim+1;j++)
      params[i][j]=0;
  }
  sorszam=new int[car+1];
  sorszam[0]=ssz;
  szomszed=new Simplex*[dim+1];
  for(int i=0;i<dim+1;i++) szomszed[i]=this;
}

//simplex destruktor: a lefoglalt memoriakat felszabaditjuk
Simplex::~Simplex(void) {
  if (params){
    for (int i=0;i<dim+1;i++)
      if (params[i]){
	for (int j=0;j<dim+1;j++)
	  if (params[i][j])
	    delete params[i][j];
	delete params[i];
      }
    delete[] params;
  }
  delete[] sorszam;
  delete[] szomszed;
}

int Simplex::dump(ostream* out) {
  Base::dump(out);
  *out << " " << sorszam[0];
  return 0;
}
