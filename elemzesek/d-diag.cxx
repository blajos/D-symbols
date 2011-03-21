#include "d-diag.hxx"

Ddiag::Ddiag(int dim,int car):
  Base(dim,car),
  buf_symmetries(0)
{
  buf_dual=0;
  buf_Rmx=0;
  simplex_orbits=new Simplex**[car+1];
  for (int k=0;k<car+1;k++) simplex_orbits[k]=new Simplex*[car];
  for (int i=0;i<car;i++) simplex_orbits[0][i]=new Simplex(dim,car,i);
}

Ddiag::Ddiag(Ddiag* olddiag){
  dim=olddiag->dim;
  car=olddiag->car;
  Ddiag(dim,car);
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

Ddiag::Ddiag(istream* in){
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

  for(list<Param*>::iterator it1=buf_params->begin();it1!=buf_params->end();
      ++it1){
    // it1 is a pointer to Param* elements, so *it1 is a Param* element
    delete *it1;
  }
  delete buf_params;

  for(vector<list<Ddiag*> *>::iterator it1=buf_cancel_operation.begin();
      it1!=buf_cancel_operation.end(); ++it1){
    // it1 is a pointer to list<Ddiag*> * elements
    for(list<Ddiag*>::iterator it2=(*it1)->begin(); it2!=(*it1)->end(); ++it2){
      // it2 is a pointer to Ddiag* elements, so *it1 is a Ddiag* element
      delete *it2;
    }
    delete *it1;
  }
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
  list<Simplex*> D;
  D.push_back(first);
  int r=0;
  // The generated numberings are indexed from 1, but the first simplex orbit is
  // still the 0th.
  int i=first->sorszam[0]+1; // Number of first in the base numbering
  first->sorszam[i]=r++;
  while (D.size() < car){
    // (*D.rbegin()) is a Simplex*
    if (find(D.begin(),D.end(),(*D.rbegin())->szomszed[0]) == D.end()) {
      D.push_back((*D.rbegin())->szomszed[0]);
    }
    else if (find(D.begin(),D.end(),(*D.rbegin())->szomszed[1]) == D.end()) {
      D.push_back((*D.rbegin())->szomszed[1]);
    }
    else
      try {
	for (int j=0;j<dim+1;j++) {
	  for (list<Simplex*>::reverse_iterator it=D.rbegin();it!=D.rend();it++){
	    if (find(D.begin(),D.end(),(*it)->szomszed[j]) == D.end()) {
	      D.push_back((*it)->szomszed[j]);
	      throw 0;
	    }
	  }
	}
      }
      catch (int a){
	a=0;
      }
    (*D.rbegin())->sorszam[i]=r++;
  }

  // nem csak sorszam szerint lehet vegigmenni...
  for (list<Simplex*>::iterator it=D.begin();it!=D.end();it++){
    int a=(*it)->sorszam[i];
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
	  buf_Rmx->set(simplex_orbits[0][r],i,j,steps);
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
        buf_dual->simplex_orbits[0][r]->szomszed[dim-i]=buf_dual->simplex_orbits[0][icsucsjszomszedja];
	buf_dual->simplex_orbits[0][r]->sorszam[0]=r;
      }
  }
  return buf_dual;
}

list<Ddiag*> *Ddiag::cancel_operation(int cancel_op) {
  if(cancel_op<0 or cancel_op > dim+1){
    throw "Operation out of range";
  }
  buf_cancel_operation.resize(dim+2,0);
  if(buf_cancel_operation[cancel_op]==0){
    buf_cancel_operation[cancel_op]=new list<Ddiag*>;

    list<int> unreachable;
   
    for(int i=0;i<car;i++) unreachable.push_back(i);    //fill unreachable
    while (!unreachable.empty()) {
      list<int> current_component;
      list<int> new_elements;
      new_elements.push_back(*unreachable.begin());
      while (!new_elements.empty()) {
	list<int> previous_elements;
        for ( list<int>::iterator p = new_elements.begin(); p !=
	    new_elements.end(); ++p ) {
          previous_elements.push_back(*p);
	  current_component.push_back(*p);
	}
        new_elements.clear();
        for ( list<int>::iterator r = previous_elements.begin();
	    r!=previous_elements.end(); ++r ){
          //utolso szomszedai metszet nemelerheto->uj
          for (int j=0;j<dim+1;j++) {                                                                    
	    if ( j!=cancel_op ){
	      list<int>::iterator new_unreachable=find(unreachable.begin(),
		  unreachable.end(),simplex_orbits[0][*r]->szomszed[j]->sorszam[0]);
	      if (new_unreachable != unreachable.end() ) {
		new_elements.push_back(*new_unreachable);
		unreachable.erase(new_unreachable);    //talalt nemelerhetoek torlese
	      }
	    }
	  }
	}
      }

      Ddiag* curr=new Ddiag(dim-1,current_component.size());
      int num=0;

      //Copy simplexes...
      // indexes[n]= Index of the simplex in the original diagram, whose index in the component
      // is n.
      vector<int> indexes;
      for ( list<int>::iterator r = current_component.begin();
	  r!=current_component.end();r++)
	indexes.push_back(simplex_orbits[0][*r]->sorszam[0]);

      for ( list<int>::iterator r = current_component.begin();
	  r!=current_component.end();r++){
	for ( int j=0;j<dim+1;j++){
	  // helper=index of the j-th adjacent simplex in the original diagram.
	  int helper = simplex_orbits[0][indexes[num]]->szomszed[j]->sorszam[0];

	  // helper1= index of the j-th adjacent simplex in the component
	  vector<int>::iterator helper1 = find(indexes.begin(),indexes.end(),helper);
	  curr->simplex_orbits[0][num]->szomszed[j] = curr->simplex_orbits[0][*helper1];
	}
	num++;
      }
      buf_cancel_operation[cancel_op]->push_back(curr);
    }
  }
  return buf_cancel_operation[cancel_op];
}

int Ddiag::check_dual(void){
  return -is_smaller(1,dual(),1);
}

int Ddiag::check_numberings(void){
  buf_symmetries=1;
  for (int i=2;i<car+1;i++){
    int ret=is_smaller(1,this,i);
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

int Ddiag::symmetries(void){
  if (buf_symmetries == 0)
    check_numberings();
  return buf_symmetries;
}

int Ddiag::check_r(void) {
  Mxfunction* R=Rmx();
  for(int i=0;i<car;i++)
    for (int u=0;u<dim+1;u++)
      for (int v=u;v<dim+1;v++){
	int val=R->get(simplex_orbits[0][i],u,v);
	if (val < 1)
	  return 0;
        if (val > 2 and abs(u-v)>=2)
	  return 0;
	if (val > 1 and u==v)
	  return 0;
      }
  return 1;
}

int Ddiag::is_bigraph(void) {
  list<Simplex*> component1,component2,*curr_comp,*other_comp,unreached;
  curr_comp=&component1;
  other_comp=&component2;
  curr_comp->push_back(simplex_orbits[0][0]);
  for(int i=1;i<car;i++)
    unreached.push_back(simplex_orbits[0][i]);
  
  while (!unreached.empty()){
    for(list<Simplex*>::iterator curr_sim=curr_comp->begin();
	curr_sim!=curr_comp->end(); curr_sim++)
      for(int i=0;i<dim+1;i++)
	// For every simplex in the component we check, if any of the operations
	// points to an other simplex in the same component; if so we don't have
	// a bipartite graph.
	if((*curr_sim)->szomszed[i]!=*curr_sim){
	  if(find(curr_comp->begin(),curr_comp->end(),(*curr_sim)->szomszed[i])
	      != curr_comp->end())
	    return 0;
	  else {
	    // Else we add only the new simplexes to the other component, and
	    // remove them from the unreached list.
	    list<Simplex*>::iterator unreached_it=find(unreached.begin(),
		unreached.end(), (*curr_sim)->szomszed[i]);
	    if(unreached_it != unreached.end()){
	      other_comp->push_back((*curr_sim)->szomszed[i]);
	      unreached.erase(unreached_it);
	    }
	  }
	}
    // Swap current and other components
    list<Simplex*> *temp=curr_comp;
    curr_comp=other_comp;
    other_comp=temp;
  }
}

int Ddiag::is_smaller(int thisindex,Ddiag* other,int otherindex){
  if (dim > other->dim)
    return -1;
  if (dim < other->dim)
    return 1;

  if (car > other->car)
    return -1;
  if (car < other->car)
    return 1;

  Simplex** thissimplex=simplex_orbits[thisindex];
  Simplex** othersimplex=other->simplex_orbits[otherindex];
  for (int j=dim;j>=0;j--)
    for (int i=0;i<car-1;i++){
      if (thissimplex[i]->szomszed[j]->sorszam[thisindex] >
          othersimplex[i]->szomszed[j]->sorszam[otherindex])
	return -1;
      if (thissimplex[i]->szomszed[j]->sorszam[thisindex] <
          othersimplex[i]->szomszed[j]->sorszam[otherindex])
	return 1;
    }
  return 0;
}

list<Param*> *Ddiag::params(void){
  if (!buf_params){
    char my_character='m';
    for(int op0=0;op0<dim;op0++){
      int op1=op0+1;
      list<Simplex*> unreached;
      for(int i=0;i<car;i++) unreached.push_back(simplex_orbits[0][i]);
      while (!unreached.empty()){
	int coeff;
	bool orientable;
	list<Simplex*> reached;
	list<Simplex*>::iterator unreached_it;
	reached.push_back(*unreached.begin());
	Simplex* starting_point=*unreached.begin();
	unreached.erase(unreached.begin());
	while((*reached.rbegin())->szomszed[op0]->szomszed[op1] !=
	    starting_point) {
	  // We don't reach the end in 1 pair of steps
	  // Add every simplex on the route to reached
	  unreached_it=find(unreached.begin(), unreached.end(),
	      (*reached.rbegin())->szomszed[op0]);
	  if(unreached_it!=unreached.end()){
	    reached.push_back(*unreached_it);
	    unreached.erase(unreached_it);
	  }
	  else 
	    orientable=false;

	  unreached_it=find(unreached.begin(), unreached.end(),
	      (*reached.rbegin())->szomszed[op0]->szomszed[op1]);
	  if(unreached_it!=unreached.end()){
	    reached.push_back(*unreached_it);
	    unreached.erase(unreached_it);
	  }
	  else 
	    orientable=false;
	}
	// There can be one more step with operation op0
	unreached_it=find(unreached.begin(), unreached.end(),
	    (*reached.rbegin())->szomszed[op0]);
	if(unreached_it!=unreached.end()){
	  reached.push_back(*unreached_it);
	  unreached.erase(unreached_it);
	}
	else 
	  orientable=false;

	coeff=Rmx()->get(starting_point,op0,op1);

	Param* curr_param=new Param(my_character++,coeff,orientable);
	for(list<Simplex*>::iterator reached_it=reached.begin();
	    reached_it!=reached.end(); reached_it++)
	  curr_param->simplex_operations.push_back(pair<Simplex*,int>(*reached_it,op0));
	buf_params->push_back(curr_param);
      }
    }
  }
  filter_bad_orbifolds(buf_params);
  return buf_params;
}

int Ddiag::dump(ostream* out){
  *out << dim << car;
  for(int i=0;i<dim+1;i++){
    *out << '(';
    list<int> simplex_indexes;
    for(int j=0;j<car;j++)
      simplex_indexes.push_back(j);
    while(! simplex_indexes.empty()){
      int first=*simplex_indexes.begin();
      simplex_indexes.erase(simplex_indexes.begin());
      int second=simplex_orbits[0][first]->szomszed[i]->sorszam[0];
      *out << '(' << first;
      if (second!=first){
	*out << ',' << second;
	simplex_indexes.erase(find(simplex_indexes.begin(),
	      simplex_indexes.end(), second));
      }
      *out << ')';
    }
    *out << ')';
  }
  return 0;
}

int Ddiag::print_html(ostream* out){
  *out << "<html>" <<endl << "  <body>" << endl << "    <table border=\"2\">" <<
    endl << "      <tr><td>";

  Svg my_svg(dim,car);
  for(int i=0;i<car;i++)
    for(int j=0;j<dim+1;j++){
      int other=simplex_orbits[0][i]->szomszed[j]->sorszam[0];
      if(other!=i)
        my_svg.add_line(i,other,j);
    }
  my_svg.print_html(out);

  *out << "</td>" << endl;
  *out << "      <td><table>" << endl;
  *out << "        <tr><td>";
  Rmx()->print_html(out);
  *out << "</td></tr>" << endl;
  *out << "        <tr><td colspan=\"" << car << "\">Number of " << dim-1 << 
    " dimensional components: (" << cancel_operation(0)->size();
  for (int j=1;j<dim+1;j++)
    *out << "," << cancel_operation(j)->size();
  *out << ")</td></tr>" << endl;
  *out << "        <tr><td colspan=\"" << car << "\">Parameters:";
  for (list<Param*>::iterator it=params()->begin();it!=params()->end();it++)
    *out << " " << (*it)->coeff << (*it)->letter << ((*it)->orientable ? "+" : "");
  *out << "</td></tr>" << endl;
  *out << "      </table></td></tr>" << endl;
  *out << "    </table>" << endl;
  *out << "  </body>" << endl;
  *out << "</html>" << endl;
  return 0;
}

int Ddiag::filter_bad_orbifolds(list<Param*>*){
  // Let's ignore this for now, because it's possible to have a Euclidean
  // tiling, which doesn't suffer from the bad orbifold problem.
  return 0;
}

