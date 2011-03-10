#include "d-diag.hxx"

Ddiag::Ddiag(int dim,int car):
  Base(dim,car)
{
  simplex_orbits=new simplex**[car+1];
  for (int k=0;k<car+1;k++) simplex_orbits[k]=new simplex*[car];
  for (int i=0;i<car;i++) simplex_orbits[0][i]=new simplex(dim,car,i);
}

Ddiag::Ddiag(Ddiag* olddiag){
  dim=olddiag->dim;
  car=olddiag->car;
  for(int i=0;i<car;i++) 
    for(int j=0;j<dim+1;j++) {
      // We set helper variable to k where D_k=\sigma_j(D_i)
      int helper=olddiag->simplex_orbits[0][i]->szomszed[j]->sorszam[0];
      // Set adjacency operations
      simplex_orbits[0][i]->szomszed[j]=simplex_orbits[0][helper];
      // Set base numbering
      simplex_orbits[0][i]->sorszam[0]=olddiag->simplex_orbits[0][i]->sorszam[0];
    }
  for(int k=1;k<car+1;k++)
    for(int i=0;i<car;i++) {
      // Set helper to the number D_i gets, when numbering starts with k.
      int helper=olddiag->simplex_orbits[k][i]->sorszam[0];
      simplex_orbits[k][i]=simplex_orbits[0][helper];
      simplex_orbits[k][i]->sorszam[k]=i;
    }
}

Ddiag::Ddiag(istream*){
  *in >> dim >> car;
  Ddiag(dim,car);
  char level_control;
  for(int i=0;i<dim+1;i++){
    *in >> level_control;
    if (level_control == '('){
      *in >> level_control;
      while (level_control == '('){
	int first,second;
        *in >> first;
	*in >> level_control;
	if (level_control == ','){
	  *in >> second;
	  *in >> level_control;
	  // Set adjacency: i is the operation, first and second are the orbits
	  create_edge(i,first,second);
	}
	if (level_control != ')')
	  throw 1;  //Something's wrong
	*in >> level_control;
      }
    }
    else {
      throw 2; //Something's wrong
    }
  }
}

Ddiag::~Ddiag(void){
  for (int i=0;i<car;i++) delete simplex_orbits[0][i];
  for (int i=0;i<car;i++) delete[] simplex_orbits[i];
  delete[] simplex_orbits;
  delete buf_dual;
  delete buf_Rmx;

  for(list<Param*>::iterator it=buf_params->begin();it1!=buf_params->end();
      ++it1){
    // it1 is a pointer to Param* elements, so *it1 is a Param* element
    delete *it1;
  }
  delete buf_params;

  for(list<Ddiag*>::iterator it=buf_cancel_operation->begin();
      it1!=buf_cancel_operation->end(); ++it1){
    // it1 is a pointer to Ddiag* elements, so *it1 is a Ddiag* element
    delete *it1;
  }
  delete buf_cancel_operation;
}

int Ddiag::create_edge(int color,int from,int to) {
  if (color<0 || color>dim || from<0 || from >=car || to<0 || to>=car )
    return 0;
  if (simplex_orbits[0][from]->szomszed[color] != simplex_orbits[0][from] ||
      simplex_orbits[0][to]->szomszed[color] != simplex_orbits[0][to]) {
    // Can't add edge, because there is another edge blocking it
    return 0;
  }
  simplex_orbits[0][from]->szomszed[color]=simplex_orbits[0][to];
  simplex_orbits[0][to]->szomszed[color]=simplex_orbits[0][from];
  return 1;
}

void Ddiag::remove_edge(int color,int from,int to){
  if (color<0 || color>dim || from<0 || from >=car || to<0 || to>=car )
    return;
  if (simplex_orbits[0][from]->szomszed[color] != simplex_orbits[0][to]) {
    // There's no such edge
    return;
  }
  simplex_orbits[0][from]->szomszed[color]=simplex_orbits[0][from];
  simplex_orbits[0][to]->szomszed[color]=simplex_orbits[0][to];
}

void Ddiag::create_numbering(Simplex* first){
  vector<Simplex*> D;
  D.push_back(first);
  int r=0;
  // The generated numberings are indexed from 1, but the first simplex orbit is
  // still the 0th.
  int i=first->sorszam[0]+1; // Number of first in the base numbering
  first->sorszam[i]=r++;
  while (D.size() < car){
    vector<Simplex*>::iterator it;
    if (find(D.begin(),D.end(),*(D.rbegin())->szomszed[0]) != D.end()) {
      D.push_back(*(D.rbegin())->szomszed[0]);
    }
    else if (find(D.begin(),D.end(),*(D.rbegin())->szomszed[1]) != D.end()) {
      D.push_back(*(D.rbegin())->szomszed[1]);
    }
    else
      try {
	for (int j=0;j<dim+1;j++) {
	  for (vector<Simplex*>::iterator it=D.rbegin(),it!=D.rend(),it++){
	    if (find(D.begin(),D.end(),*it->szomszed[j]) != D.end()) {
	      D.push_back(*it->szomszed[j]);
	      throw 0;
	    }
	  }
	}
      }
      catch (int a){
	void;
      }
    *(D.rbegin())->sorszam[i]=r++;
  }

  // nem csak sorszam szerint lehet vegigmenni...
  for (vector<Simplex*>::iterator it=D.begin(),it!=D.end(),it++){
    int a=*it->sorszam[i];
    simplex_orbits[i][a]=*it;
  }
}

Mxfunction *Ddiag::Rmx(void) {
  if (!buf_Rmx){
    buf_Rmx=new Mxfunction(dim,car);
    for (int r=0;r<car;r++)
      for (int i=0;i<dim+1;i++)
	for (int j=0;j<dim+1;j++) {
	  Simplex* start=simplex_orbits[0][r];
	  Simplex* next=start->szomszed[i]->szomszed[j];
	  int steps=1;
	  while (next!=start){
	    next=next->szomszed[i]->szomszed[j];
	    steps++;
	  }
	  buf_Rmx->set(r,i,j,steps);
	}
  }
  return buf_Rmx;
}

Ddiag* Ddiag::dual(void) {
  if (!buf_dual){
    buf_dual=new Ddiag(dim,car);
    for (int r=0;r<car;r++)
      for (int i=0;i<dim+1;i++){
	int icsucsjszomszedja=simplex_orbits[0][r]->szomszed[i]->sorszam[0];
        buf_dual->simplex_orbits[0][r]->szomszed[dim-i]=buf_dual->csucsok[0][icsucsjszomszedja];
	buf_dual->simplex_orbits[0][r]->sorszam[0]=r;
      }
  }
  return buf_dual;
}

list<Ddiag*> *Ddiag::cancel_operation(int i);
int Ddiag::filter_bad_orbifolds(list<Param*>*);
int Ddiag::check_dual(void);
int Ddiag::check_numberings(void);
int Ddiag::symmetries(void);
int Ddiag::check_r(void);
int Ddiag::check_all(void);
int Ddiag::is_bigraph(void);
int Ddiag::is_smaller(int i,Ddiag* other,int j);
list<Param*> *Ddiag::params(void);
int Ddiag::dump(ostream*);
int Ddiag::print_html(ostream*);
