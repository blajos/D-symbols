#include "backtrack.hxx"
#include <boost/math/common_factor.hpp>
#include <iostream>
#include <fstream>
#include <list>
#include <algorithm>
#include <set>
#include <math.h>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#define PI 3.14159265

//simplex konstruktor: megadjuk a dimenziot, az elemszamot, lefoglaljuk a
//szukseges memoriatartomanyokat a szomszed es a sorszam tomboknek.
//Beallitjuk a 0. variacio sorszamat az adottra. (Ezt a sorszamot a felsorolo
//algoritmus tolti ki, nem feltetlenul ugyanaz, mint az 1-es elemmel valo
//kezdeskor a sorszam.
//A dim+1 darab szomszedsagi operacio eredmenye kezdetben maga a szimplex.
//A parameterek listajat inicializaljuk.
simplex::simplex(int dimin,int carin,int ssz):
  dim(dimin),
  car(carin)
{
  sorszam=new int[car+1];
  sorszam[0]=ssz;
  szomszed=new simplex*[dim+1];
  for(int i=0;i<dim+1;i++) szomszed[i]=this;
  mx=new int*[dim+1];
  for (int i=0;i<dim+1;i++) mx[i]=new int[dim+1];
}

//simplex destruktor: a lefoglalt memoriakat felszabaditjuk
simplex::~simplex(void) {
  delete[] sorszam;
  delete[] szomszed;
  for (int i=0;i<dim+1;i++) delete[] mx[i];
  delete[] mx;
}

//Dsym konstruktor: dimenzio, elemszam, matrixok szamanak kitoltese, involucio
//matrix 0-ra definialasa, csucsok 2 dimenzios simplex-pointerekbol allo tomb
//feltoltese. (Az elso "koord" a szamlalas kezdopontjanak valasztott szimplex
//sorszama+1, mert a 0-adikat szinte veletlenszeruen toltjuk ki; a masodik koord
//az adott kezdopont szerinti sorszama a megfelelo szimplexnek.)
//Megfelelo memoriacimek lefoglalasa.
void Dsym::create(void){
  involucio=0;
  mxnum=0;
  pmaxertek=4;		//kisebbre nem választhatjuk, es majd megnő
  //az eredetivel egyutt az atsorszamozasok szama
  csucsok=new simplex**[car+1];	
  for (int k=0;k<car+1;k++) csucsok[k]=new simplex*[car];
  for (int i=0;i<car;i++) csucsok[0][i]=new simplex(dim,car,i);
  params=new std::list<param>::iterator*[car];
  for (int i=0;i<car;i++) params[i]=new std::list<param>::iterator[dim];
  dual=0;
  klist=new std::list<kisebbdim>[dim+1];
}

Dsym::Dsym(int dimin, int carin):
  dim(dimin),
  car(carin)
{
  create();
}

// Restore Dsym from a dump
Dsym::Dsym(std::istream* in){
  *in >> dim >> car;
  create();
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
          elhozzaad(i,first-1,second-1);
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

int Dsym::dump(std::ostream *out){
  *out<<dim<<" "<<car;
  for(int j=0;j<dim+1;j++){
    *out<<" (";
    for(int i=0;i<car;i++){
      int icsucsjszomszedja=csucsok[0][i]->szomszed[j]->sorszam[0];
      if (icsucsjszomszedja == i)
	*out<<"("<<i+1<<")";
      else if (i < icsucsjszomszedja)
	*out<<"("<<i+1<<","<<icsucsjszomszedja+1<<")";
    }
    *out<<")";
  }
  return 0;
}

int Dsym::dump(std::ostream *out, int first){
  *out<<dim<<" "<<car;
  for(int j=0;j<dim+1;j++){
    *out<<" (";
    for(int i=0;i<car;i++){
      int icsucsjszomszedja=csucsok[first][i]->szomszed[j]->sorszam[first];
      if (icsucsjszomszedja == i)
	*out<<"("<<i+1<<")";
      else if (i < icsucsjszomszedja)
	*out<<"("<<i+1<<","<<icsucsjszomszedja+1<<")";
    }
    *out<<")";
  }
  return 0;
}

int Dsym::dump(char *output){
  sprintf(output,"%d %d ",dim,car);
  for(int j=0;j<dim+1;j++){
    strcat(output,"(");
    for(int i=0;i<car;i++){
      int icsucsjszomszedja=csucsok[0][i]->szomszed[j]->sorszam[0];
      char temp[20];
      if (icsucsjszomszedja == i)
	sprintf(temp,"(%d)",i+1);
      else if (i < icsucsjszomszedja)
	sprintf(temp,"(%d,%d)",i+1,icsucsjszomszedja+1);
      strcat(output,temp);
    }
    strcat(output,")");
  }
  return 0;
}

int Dsym::dump(char *output,int first){
  sprintf(output,"%d %d ",dim,car);
  for(int j=0;j<dim+1;j++){
    strcat(output,"(");
    for(int i=0;i<car;i++){
      int icsucsjszomszedja=csucsok[first][i]->szomszed[j]->sorszam[first];
      char temp[20];
      if (icsucsjszomszedja == i)
	sprintf(temp,"(%d)",i+1);
      else if (i < icsucsjszomszedja)
	sprintf(temp,"(%d,%d)",i+1,icsucsjszomszedja+1);
      strcat(output,temp);
    }
    strcat(output,")");
  }
  return 0;
}

//Dsym destruktor: memoria felszabaditasa
Dsym::~Dsym(void) {
  if(dual!=0) delete dual;
  if(involucio) {
    for(int j=0;j<dim+1;j++){
      for(int k=0;k<*involucio[j][car];k++)
	delete[] involucio[j][k];
      delete involucio[j][car];
      delete[] involucio[j];
    }
    delete[] involucio;
  }
  if(params){
    for (int i=0;i<car;i++) delete[] params[i];
    delete[] params;
  }
  delete[] klist;

  // Ez a vegere kell...
  for (int i=0;i<car;i++) delete csucsok[0][i];
  for (int i=0;i<car+1;i++) delete[] csucsok[i];
  delete[] csucsok;
}

//Dsym::save: letrehoz egy uj Dsym objektumra mutato pointert, es egy uj
//objektumot is. 
//Atmasoljuk bele az osszes adatot az eredetibol, figyelve arra, hogy veletlenul
//se az eredeti egy simplexere mutato pointert mentsunk, hanem a sajaton belul a
//megfelelo sorszamut talaljuk meg.
Dsym* Dsym::save(void) {
  Dsym* saved=new Dsym(dim,car);
  for(int i=0;i<car;i++) for(int j=0;j<dim+1;j++) {
    //Ehhez kell tudni azonositani a csucsokat:
    int icsucsjszomszedja=csucsok[0][i]->szomszed[j]->sorszam[0];  
    //szomszedsagok beallitasa
    saved->csucsok[0][i]->szomszed[j]=saved->csucsok[0][icsucsjszomszedja];	
    //sorszamok
    saved->csucsok[0][i]->sorszam[0]=csucsok[0][i]->sorszam[0];	
  }
  for(int k=1;k<car+1;k++)
    for(int i=0;i<car;i++) {
      int icsucsksorszama=csucsok[k][i]->sorszam[0];
      saved->csucsok[k][i]=saved->csucsok[0][icsucsksorszama];
      saved->csucsok[k][i]->sorszam[k]=i;
    }
  //saved->invol_create(0);
  //saved->invol_create(1);
  //mx mentese
  for (int r=0;r<car;r++){
    for (int i=0;i<dim+1;i++)
      for (int j=0;j<dim+1;j++)
	saved->csucsok[0][r]->mx[i][j]=csucsok[0][r]->mx[i][j];
  }
  return saved;
}

//Dsym::save_with_start: letrehoz egy uj Dsym objektumra mutato pointert, es egy uj
//objektumot is. 
//Atmasoljuk bele az osszes adatot az eredetibol, figyelve arra, hogy veletlenul
//se az eredeti egy simplexere mutato pointert mentsunk, hanem a sajaton belul a
//megfelelo sorszamut talaljuk meg. A 0. simplex mas lesz.
Dsym* Dsym::save_with_start(int start) {
  Dsym* saved=new Dsym(dim,car);
  for(int i=0;i<car;i++) for(int j=0;j<dim+1;j++) {
    //Ehhez kell tudni azonositani a csucsokat:
    int icsucsjszomszedja=csucsok[start][i]->szomszed[j]->sorszam[start];  
    //szomszedsagok beallitasa
    saved->csucsok[0][i]->szomszed[j]=saved->csucsok[0][icsucsjszomszedja];	
    //sorszamok
    saved->csucsok[0][i]->sorszam[0]=csucsok[start][i]->sorszam[start];	
  }
  for(int i=0;i<car;i++) saved->atsorszamoz(i+1);
  //saved->invol_create(0);
  //saved->invol_create(1);
  //mx mentese
  for (int r=0;r<car;r++){
    for (int i=0;i<dim+1;i++)
      for (int j=0;j<dim+1;j++)
	saved->csucsok[0][r]->mx[i][j]=csucsok[start][r]->mx[i][j];
  }
  return saved;
}

//Dsym::atsorszamoz(i): i=1..car (0-val nem hivjuk meg, mert azt tetszolegesen
//toltjuk fel.)
//Az i. csucsbol(szimplexbol) indulva alkalmazzuk a sorszamozo algoritmust, es a
//csucsok[i] listat feltoltjuk az algoritmus szerint.
//Az algoritmus roviden: D1,D2...Dr mar meg van szamozva. Dr-re alkalmazzuk a 0.
//operaciot, ha olyan csucsot(szimplexet) kapunk, ami meg nem szerepel, az a
//kovetkezo; aztan Dr-re alkalmazzuk az 1. operaciot (szinten, ha ujat kapunk,
//az lesz az r+1-edik); 2,3,...dim operacio Dr-re; 0,1,...dim operacio
//D(r-1)-re; 0..dim D(r-2)-re;...;0..dim D1-re
void Dsym::atsorszamoz(int i) {
  int r=0;
  simplex** D;
  D=new simplex*[car];
  //azert i-1, mert a masodik koord 0..car-1, az elso pedig 0..car:
  D[0]=csucsok[0][i-1];
  D[0]->sorszam[i]=0;
  while (r<car-1){
    r++;
    D[r]=NULL;
    if (my_find(D[r-1]->szomszed[0],r,D) == 0) {
      D[r]=D[r-1]->szomszed[0];
    }                
    else if (my_find(D[r-1]->szomszed[1],r,D) == 0) {
      D[r]=D[r-1]->szomszed[1];
    }                
    else             
      for (int j=0;j<dim+1;j++) {
	for (int k=r-1;k>=0;k--) {
	  if (my_find(D[k]->szomszed[j],r,D) == 0) {
	    D[r]=D[k]->szomszed[j];
	    break;    
	  }           
	}             
	if (D[r] != NULL) {
	  break;      
	}             
      }              
    D[r]->sorszam[i]=r;
  }
  for (r=0;r<car;r++) {
    csucsok[i][r]=D[r];
  }
  delete[] D;
}

//Trukk: A felsorolo algoritmusunk szerint, ha az i. az elso olyan szimplex
//osztaly, akinek nincs ele, akkor a nala kisebbek osszefuggoek, a nagyobbaknak
//pedig szinten nincs ele.
void Dsym::atsorszamoz(int i, int max) {
  int r=0;
  simplex** D;
  D=new simplex*[car];
  //azert i-1, mert a masodik koord 0..car-1, az elso pedig 0..car:
  D[0]=csucsok[0][i-1];
  D[0]->sorszam[i]=0;
  while (r<max){
    r++;
    D[r]=NULL;
    if (my_find(D[r-1]->szomszed[0],r,D) == 0) {
      D[r]=D[r-1]->szomszed[0];
    }                
    else if (my_find(D[r-1]->szomszed[1],r,D) == 0) {
      D[r]=D[r-1]->szomszed[1];
    }                
    else             
      for (int j=0;j<dim+1;j++) {
	for (int k=r-1;k>=0;k--) {
	  if (my_find(D[k]->szomszed[j],r,D) == 0) {
	    D[r]=D[k]->szomszed[j];
	    break;    
	  }           
	}             
	if (D[r] != NULL) {
	  break;      
	}             
      }              
    D[r]->sorszam[i]=r;
  }
  for (r=0;r<=max;r++) {
    csucsok[i][r]=D[r];
  }
  for(r=max+1;r<car;r++){
    csucsok[i][r]=csucsok[0][r];
    csucsok[i][r]->sorszam[i]=r;
  }
  delete[] D;
}

//my_find(...): mit a keresendo dolog, hossz a lista hossza, hol a lista elejere
//mutato pointer.
//Csak egy egyszeru segedfv keresesre.
int my_find(simplex* mit, int hossz, simplex** hol){
  for(int i=0;i<hossz;i++)
    if(mit==hol[i])
      return 1;
  return 0;
}

//kisebb(...): Adott egyik es masik D-szimbolum, ev es mv kezdo
//csucsokkal(szimplexekkel) Erre alkalmazzuk a rendezesi algoritmusunkat, aminek
//alapja a tavolsag ket csucs kozott: az egyik csucsbol vett szamozas szerinti
//sorszama a masik csucsnak. Ezutan ugy rendezunk, hogy vesszuk az operaciokat
//csokkeno sorrendben majd a csucsokat novekvoben es megvizsgaljuk, hogy az
//adott operaciot alkalmazva a csucsra, milyen tavolsagra kerulunk a csucstol.
//Az elso hely, ahol a ket szimbolum kulonbozik, dont: amelyik kisebb, az a
//szimbolun a kisebb (azaz lexikografikusan rendezunk.)
int kisebb(Dsym* egyik,int ev,Dsym* masik,int mv){
  if (egyik->dim != masik->dim) return (masik->dim > egyik->dim) ? 1: -1;
  if (egyik->car != masik->car) return (masik->car > egyik->car) ? 1: -1;
  int dim=egyik->dim;
  int car=egyik->car;
  simplex** egyiks=egyik->csucsok[ev];
  simplex** masiks=masik->csucsok[mv];
  for (int j=dim;j>=0;j--)
    for (int i=0;i<car-1;i++){
      if (egyiks[i]->szomszed[j]->sorszam[ev] != 
	  masiks[i]->szomszed[j]->sorszam[mv])
      {
	return (masiks[i]->szomszed[j]->sorszam[mv] > 
	    egyiks[i]->szomszed[j]->sorszam[ev]) ? 1 : -1;
      }
    }
  return 0;
}

//Dsym::dualis: dual valtozo feltoltese. Felcsereljuk az i. es a (dim-i).
//operaciokat minden i-re. Visszateresi ertek: -1, ha a dualis kisebb, 0 ha
//egyenloek, 1 ha a dualis nagyobb.
int Dsym::dualis(void) {
  if(!dual) dual=new Dsym(dim,car);
  for(int i=0;i<car;i++) for(int j=0;j<dim+1;j++) {
    int icsucsjszomszedja=csucsok[0][i]->szomszed[j]->sorszam[0];
    //felcsereljuk az elek szineit... (dim-j)
    dual->csucsok[0][i]->szomszed[dim-j]=dual->csucsok[0][icsucsjszomszedja];
    dual->csucsok[0][i]->sorszam[0]=i; 
  }
  atsorszamoz(1);
  dual->atsorszamoz(1);
  return kisebb(dual,1,this,1);
}

//Dsym::sorszamozas: az osszes csucsbol indulva megnezzuk, hogy talalunk-e
//kisebbet a D-szimbolumunknal. Az egyformakat szamoljuk.
int Dsym::sorszamozas(void) {
  int k;
  int egyformak=1;
  atsorszamoz(1);
  if(kisebb(this,0,this,1)!=0) return 1;
  for (int i=2;i<car+1;i++) {
    atsorszamoz(i);
    k=kisebb(this,i,this,1);
    if (k==1) {
      return k;
    }
    else if (k==0) egyformak++;
  }
  if (egyformak>1 && egyformak<car) return 0;
  else return -1;
}

//Dsym::osszefuggo: Egymas utan alkalmazott szelessegi bejarasokkal
//megszamoljuk, hogy az elhagy valtozo altal jelzett operaciot elhagyva hany
//reszre esik szet a multigraf. (elhagy=-1 eseten pl nem hagyunk el valtozot.)
//A nemelerheto listat toltjuk fel az elejen, es torolgetjuk belole az elsot,
//annak szomszedait, majd annak szomszedait... ameddig lehet, ha elfogyott,
//veszunk egy uj elso elemet, noveljuk eggyel a komponensek szamat, es kezdjuk
//elolrol, amig nem ures a nemelerheto lista.
int Dsym::osszefuggo(int elhagy) {
  std::list<int> nemelerheto,utolso,uj;
  std::list<int>::iterator hanyadik;
  int hanydarab=0;

  for(int i=0;i<car;i++) nemelerheto.push_back(i);	//nemelerheto feltoltese
  while (!nemelerheto.empty()) {
    uj.clear();
    uj.push_back(*nemelerheto.begin());
    hanydarab++;
    while (!uj.empty()) {
      utolso.clear();
      for ( std::list<int>::iterator p = uj.begin(); p != uj.end(); ++p ) 
	utolso.push_back(*p);	//uj->utolso
      uj.clear();
      for (std::list<int>::iterator r = utolso.begin(  ); r != utolso.end(  ); ++r ){
	//utolso szomszedai metszet nemelerheto->uj
	for (int j=0;j<dim+1;j++) {
	  if ( j!=elhagy ){
	    hanyadik=find(nemelerheto.begin(),nemelerheto.end(),
		csucsok[0][*r]->szomszed[j]->sorszam[0]);
	    if (hanyadik != nemelerheto.end() ) {
	      uj.push_back(*hanyadik);
	      //talalt nemelerhetoek torlese
	      nemelerheto.erase(hanyadik);
	    }
	  }
	}
      }
    }
  }
  return hanydarab;
}

//Dsym::uvw: Ellenorizzuk igaz-e, hogy barmely ket nem szomszedos (es nem
//azonos) operaciot 4-szer egymas utan alkalmazva visszajutunk a kiindulasi
//pontba. Es feltoltjuk a matrixokban az informaciot (a szabad parameterekkel
//meg nem foglalkozunk, csak beirjuk oket a matrixba.) 0-t adunk vissza, ha nem
//teljesul a fenti feltetel, 1-et, ha teljesul.
int Dsym::uvw(void) {
  for (int r=0;r<car;r++) {				//atlotol tavolabbi
    for (int i=0;i<dim-1;i++)
      for (int i1=i+2;i1<dim+1;i1++){
	if(csucsok[1][r]->szomszed[i]->szomszed[i1]->szomszed[i]->szomszed[i1]
	    != csucsok[1][r]){
	  return 0;
	}
	else {
	  csucsok[1][r]->mx[i][i1]=2;
	  csucsok[1][r]->mx[i1][i]=2;
	}
      }
    for (int i=0;i<dim+1;i++) csucsok[1][r]->mx[i][i]=1;     //atlo
    for (int i=0;i<dim;i++) {				       //atlo szomszedok
      int u=1;
      simplex* csucs=csucsok[1][r];
      while (csucs->szomszed[i]->szomszed[i+1] != csucsok[1][r]) {
	u++;
	csucs=csucs->szomszed[i]->szomszed[i+1];
      }
      csucsok[1][r]->mx[i][i+1]=u;
      csucsok[1][r]->mx[i+1][i]=u;
    }
  }
  return 1;
}


//Dsym::ellenoriz: Az elobb definialt ellenorzesek lefuttatasa.
int Dsym::ellenoriz(void) {
  if (osszefuggo(-1) > 1){
    //std::cout << "Nem osszefuggo" <<std::endl;
    return 1; //csak az osszefuggoek erdekesek
  }
  sorszamozas();
  //if (sor==1) return 1; 
  if (!uvw()){
    //std::cout << "Nem uvw0" <<std::endl;
    //print(0);
    return 1;
  }
  return -1;
}

//Dsym::elhozzaad(...): A szin, honnan, hova valtozoharmas altal leirt elt adja
//hozza a multigrafhoz, ha ezt meg lehet tenni, azaz egyik csucsanal sincs meg
//definialva a szin altal leirt szomszedsagi operacio.
int Dsym::elhozzaad(int szin, int honnan, int hova) {
  if (szin<0 || szin>dim || honnan<0 || honnan >=car || hova<0 || hova>=car ) 
    return 0;
  if (csucsok[0][honnan]->szomszed[szin] != csucsok[0][honnan] ||
      csucsok[0][hova]->szomszed[szin] != csucsok[0][hova]) {
    //std::cout << "Nem tudok elt hozzaadni, mert nincs torolve az elozo" << std::endl;
    return 0;
  }
  csucsok[0][honnan]->szomszed[szin]=csucsok[0][hova];
  csucsok[0][hova]->szomszed[szin]=csucsok[0][honnan];
  return 1;
}

//Dsym::eltorol: Torli a megfelelo elt, ha ez ertelmesen megteheto.
void Dsym::eltorol(int szin, int honnan, int hova) {
  if (szin<0 || szin>dim || honnan<0 || honnan >=car || hova<0 || hova>=car ) 
    return;
  if (csucsok[0][honnan]->szomszed[szin] != csucsok[0][hova]) {
    std::cerr << "Nincs el (szin, honnan, hova): " <<szin<<" "<<honnan<<" "<<hova
      << std::endl;
    return;
  }
  csucsok[0][honnan]->szomszed[szin]=csucsok[0][honnan];
  csucsok[0][hova]->szomszed[szin]=csucsok[0][hova];
}


//Dsym::invol_create: involuciokent egyszeruen abrazolhato a multigraf konzolon.
//Ezt hozzuk itt letre. Csak az elore elek kellenek nekunk 
//(kisebb sorszam->nagyobb sorszam)
void Dsym::invol_create(int var) {
  if(involucio)
    for (int j=0;j<dim+1;j++) {
      for (int d=0;d<*involucio[j][car];d++)
	delete[] involucio[j][d];
      delete involucio[j][car];
      delete[] involucio[j];
    }
  involucio=new int**[dim+1];
  for (int j=0;j<dim+1;j++) {
    involucio[j]=new int*[car+1];
    involucio[j][car]=new int;	//*involucio[j][car] az aktualis hossz
    *involucio[j][car]=0;
    for (int i=0;i<car;i++) {
      int kell=1;
      int masik=csucsok[var][i]->szomszed[j]->sorszam[var];
      if(masik==i) continue;
      for(int k=0;k<*involucio[j][car];k++){
	if (involucio[j][k][1]==i ) {
	  kell=0;
	  break;
	}
      }
      if (kell==1) {
	involucio[j][*involucio[j][car]]=new int[2];
	involucio[j][*involucio[j][car]][0]=i;
	involucio[j][*involucio[j][car]][1]=masik;
	++(*involucio[j][car]);
      }
    }
  }
}


//Dsym::print: A fenti involucio konzolra irasa, illetve az egyes operaciok 
//elhagyasaval hany darabra esik a multigraf.
void Dsym::print(int var) {
  invol_create(var);
  std::cout << "=========================" << std::endl;
  for (int j=0;j<dim+1;j++) {
    for(int k=0;k<*involucio[j][car];k++){
      std::cout << "[" << involucio[j][k][0]+1 << "," << involucio[j][k][1]+1 <<"] ";
    }
    std::cout << std::endl;
  }
  std::cout << "Number of " << dim-1 << " dimensional components: ";
  std::cout << "("<<osszefuggo(0);
  for (int i=1;i<dim+1;i++) std::cout <<","<<osszefuggo(i);
  std::cout << ")"<<std::endl;
}

//Dsym::print_params: Kiirja az op0 es az op0+1 operaciokhoz tartozo
//parametereket, letter1 betuvel kezdve az out kimenetre.
void Dsym::print_params(int op0,std::ostream* out) {
  for (std::list<param>::iterator currparam=plist.begin();currparam!=plist.end();
      currparam++)
    if(op0==currparam->op){
      *out << abs(currparam->eh) << currparam->kar;
      if ( currparam->eh>0 ) *out<<"+";
      *out<<" ";
    }
}

//Dsym::write_xfig: file-ba letrehozza a megfelelo tartalmat, amit a fig2dev
//programmal tetszoleges formatumra lehet konvertalni. (Csak az eleket kell
//berajzolni, mert minden mas megegyezik az azonos dimenzioju es elemszamu
//multigrafok kozott, es minden elt csak egy iranyba szabad megrajzolni,
//kulonben a mintak fedesbe kerulhetnek.)
void Dsym::write_xfig(std::string file) {
  xfig out(file,dim,car);
  for(int i=0;i<car-1;i++) 
    for(int j=0;j<dim+1;j++)
      if(csucsok[1][i]->szomszed[j]->sorszam[1]>csucsok[1][i]->sorszam[1])
	out.create_line(i,csucsok[1][i]->szomszed[j]->sorszam[1],j);
}


//Dsym::create_kdim: klist[i] a D-szimbolum 0. operaciojanak elhagyasaval
//keletkezo komponensek listaja, egy-egy komponensben a benne levo szimplexek
//listajaval. 
//A komponenseket, mint korabban az osszefuggoseg kerdesenel is, tobbszori
//szelessegi bejarassal hatarozzuk meg.
//Vegul megnezzuk, hogy a komponens iranyitott vagy iranyitatlan feluletet ad
//meg.
void Dsym::create_kdim(void) {
  for (int elhagy=0;elhagy<dim+1;elhagy++){
    int currkis=elhagy;

    std::list<int> nemelerheto,utolso,uj;
    std::list<int>::iterator hanyadik;

    for(int i=0;i<car;i++) nemelerheto.push_back(i);	//nemelerheto feltoltese
    while (!nemelerheto.empty()) {
      kisebbdim curr;
      uj.clear();
      uj.push_back(*nemelerheto.begin());
      while (!uj.empty()) {
	utolso.clear();
	for ( std::list<int>::iterator p = uj.begin(); p != uj.end(); ++p ) 
	  utolso.push_back(*p);//uj->utolso
	uj.clear();
	for ( std::list<int>::iterator r = utolso.begin(); r!=utolso.end(); ++r ){
	  //utolso szomszedai metszet nemelerheto->uj
	  for (int j=0;j<dim+1;j++) {
	    if ( j!=elhagy ){
	      hanyadik=find(nemelerheto.begin(),nemelerheto.end(),
		  csucsok[0][*r]->szomszed[j]->sorszam[0]);
	      if (hanyadik != nemelerheto.end() ) {
		uj.push_back(*hanyadik);
		curr.szek.push_back(csucsok[0][*hanyadik]);
		nemelerheto.erase(hanyadik);	//talalt nemelerhetoek torlese
	      }
	    }
	  }
	}
      }

      //Megnezzuk, hogy iranyitott-e a 2 dimenzios felulet, azaz nincs 
      //keresztsapka. Ha a graf paros, akkor iranyitott, ha nem, akkor nem.
      std::list<simplex*> egyik, masik, *akt,*nemakt,nemfelsorolt;
      for(std::list<simplex*>::iterator szit=curr.szek.begin();
	  szit!=curr.szek.end();szit++) 
	if(szit==curr.szek.begin())
	  egyik.push_back(*szit);
	else
	  nemfelsorolt.push_back(*szit);

      akt=&egyik;
      nemakt=&masik;
      curr.iranyitott=true;
      while( ! nemfelsorolt.empty()){
	for(std::list<simplex*>::iterator szit=akt->begin();
	    szit!=akt->end();szit++)
	  for(int i=0;i<dim+1;i++)
	    if(i!=elhagy && (*szit)->szomszed[i]!=(*szit)){
	      if(find(akt->begin(),akt->end(),(*szit)->szomszed[i])!=
		  akt->end())
		curr.iranyitott=false;
	      std::list<simplex*>::iterator hanyadik=find(nemfelsorolt.begin(),
		  nemfelsorolt.end(),(*szit)->szomszed[i]);
	      if(hanyadik!=nemfelsorolt.end()){
		nemakt->push_back(*hanyadik);
		nemfelsorolt.erase(hanyadik);
	      }
	    }
	if(! curr.iranyitott){
	  break;
	}
	//most felcsereljuk aktot, nemakttal
	std::list<simplex*> *temp;
	temp=akt;
	akt=nemakt;
	nemakt=temp;
      }
      //Vegul ujra atfutjuk, mert az utolso lepesnel keletkezhettek elek
      //nemakt-ban, ami most akt.
      for(std::list<simplex*>::iterator szit=akt->begin();
	  szit!=akt->end();szit++)
	for(int i=0;i<dim+1;i++)
	  if(i!=elhagy && (*szit)->szomszed[i]!=(*szit)){
	    if(find(akt->begin(),akt->end(),(*szit)->szomszed[i])!=
		akt->end())
	      curr.iranyitott=false;
	  }

      klist[currkis].push_back(curr);
    }
  }
}

//Dsym::kdimsf: Igaz-e, hogy ha elhagyjuk a 0. vagy az n. operaciot, akkor a
//kisebb dimenzios szimbolum szferikus vagy euklideszi sikon megvalosulhat.
//Komponensenkent szamolunk: osszeadjuk az 1/m01+1/m12-1/2 ertekeket, ha nem
//negativ, akkor igaz az allitas. Mivel lebegopontos szamitasokat vegzunk, nem
//vizsgalhatunk ==0-t.
//Kisebb dimenzios reszek szferikusak: 1, van Euk: 0, van hip:-1
int Dsym::kdimsf(std::list<param>::iterator inf_param){
  if (dim!=3) return -100;
  int ret=1;
  for (int elhagy=0;elhagy<dim+1;elhagy+=dim){
    int currkis=elhagy;
    for(std::list<kisebbdim>::iterator curr=klist[currkis].begin();
	curr!=klist[currkis].end();curr++){
      float sum=0;
      for(std::list<simplex*>::iterator currszim=curr->szek.begin();
	  currszim!=curr->szek.end();currszim++){
	for(int j=0;j<dim;j++){
	  if(!((elhagy==0 && j==0) || (elhagy==dim && j==dim-1))){
	    if(params[(*currszim)->sorszam[1]][j]!=inf_param)
	      sum+=1.0/float((*currszim)->mx[j][j+1]);
	    else
	      sum+=10*THRESH; // Ezzel jelezzuk, hogy a vegtelen lancnak az elemei erdekesek, nem a hatarerteke
	  }
	}
	sum-=0.5;
      }
      if(sum<-THRESH) {
	ret=-1;
	return ret;
      }
      if(sum<THRESH && sum>-THRESH) 
	ret=0;
    }
  }
  return ret;
}

//Dsym::min: Megnezzuk, hogy a multigraf szimmetriain tul van-e a
//matrix-rendszernek is szimmetriaja, illetve ha nincs akkor az aktualis
//sorszamozas-e a legkisebb a matrix-rendszert is figyelembe veve. Azaz
//egyuttveve minimalis-e a D-szimbolum.
int Dsym::min(void){
  osszevlist.clear();
  int i=1;  //csak az elso csucs kell sorszamozashoz, nem max eseten is eleg.
  int nemmax=1;
  //Ha van a matrixok kozott kulonbozo, akkor lehet csak max:
  for(int op=0;op<dim;op++)
    for(int k=0;k<car-1;k++)
      if(csucsok[1][k]->mx[op][op+1]!=csucsok[1][k+1]->mx[op][op+1])
	nemmax=0;
  if (nemmax==1) 	//azaz minden matrix egyenlo
    for(int j=i+1;j<car+1;j++)
      osszevlist.push_back(j);
  else
    for(int j=i+1;j<car+1;j++)
      if(kisebb(this,i,this,j)==0){
	//felig szimmetrikus, azaz a matrixot meg nem vizsgaltuk
	int max=0;
	for(int op=0;op<dim;op++){
	  for(int k=0;k<car;k++)
	    if(csucsok[i][k]->mx[op][op+1]<csucsok[j][k]->mx[op][op+1]){
	      max=1;
	      break;
	    }
	    else if(csucsok[i][k]->mx[op][op+1]>csucsok[j][k]->mx[op][op+1]){
	      return -1;
	    }
	  if(max==1)
	    break;     //Most az aktualis parban az 1 kisebb
	}
	if(max!=1){
	  nemmax=1;
	  osszevlist.push_back(j);
	}
      }
  if(nemmax==1)
    return -2;	//nem atsorszamozhato, de nem is max
  else
    return 1;
}

//Dsym::create_params: plist parameter lista feltoltese az egyutthatoval, a
//megfelelo operacio-parral, es a szimplexekkel, amiket a parameter erint;
//illetve a parameter ertek beallitasa arra a minimalis ertekre, ahol az
//egyutthato es az ertek szorzata legalabb 3. Majd a mx-ok feltoltese minden
//szimplexben. Egy parameter egyutthatoja negativ, ha iranyitatlan, pozitiv, ha
//iranyitott. (A mx-okban csak az abszolut erteket tuntejuk fel.)
void Dsym::create_params(void){
  char karakter='a';
  for(int op0=0;op0<dim;op0++){ 
    int op1=op0+1;
    std::list<int> nemelerheto;
    std::list<int>::iterator hanyadik;

    for(int i=0;i<car;i++) nemelerheto.push_back(i);	//nemelerheto feltoltese
    while (!nemelerheto.empty()) {
      param currparam;
      currparam.eh=0;	//a param egyutthatoja
      currparam.op=op0;
      currparam.kar=karakter;
      karakter++;
      int akt=-1; //az aktualis hely a koron
      int ir=1; //1 iranyitott, 0 iranyitatlan
      int kezdet=*nemelerheto.begin();
      while ( kezdet != akt ) {
	if (akt==-1) akt=kezdet;

	hanyadik=find(nemelerheto.begin(),nemelerheto.end(),akt);
	if (hanyadik != nemelerheto.end() ) {
	  currparam.szek.push_back(csucsok[1][*hanyadik]);
	  nemelerheto.erase(hanyadik);
	}
	if (akt==csucsok[1][akt]->szomszed[op0]->sorszam[1]) ir=0;
	akt=csucsok[1][akt]->szomszed[op0]->sorszam[1];

	hanyadik=find(nemelerheto.begin(),nemelerheto.end(),akt);
	if (hanyadik != nemelerheto.end() ) {
	  currparam.szek.push_back(csucsok[1][*hanyadik]);
	  nemelerheto.erase(hanyadik);
	}
	if (akt==csucsok[1][akt]->szomszed[op1]->sorszam[1]) ir=0;
	akt=csucsok[1][akt]->szomszed[op1]->sorszam[1];

	currparam.eh++;
      }
      for(currparam.min_ertek=0;currparam.eh*currparam.min_ertek < DEGENERATION_LIMIT;
	  currparam.min_ertek++);
      if ( ir!=1 ) currparam.eh=-currparam.eh;
      plist.push_back(currparam);
      std::list<param>::iterator lastp=plist.end();
      lastp--;
      for(std::list<simplex*>::iterator szit=currparam.szek.begin();
	  szit!=currparam.szek.end();szit++)
	params[(*szit)->sorszam[1]][op0]=lastp;
    }
  }
  for(std::list<param>::iterator currparam=plist.begin();currparam!=plist.end();
      currparam++){
    currparam->ertek=currparam->min_ertek;
    change_param(currparam,0); //mx-fuggveny frissitese
  }
}

//Dsym::filter_bad_orbifolds: Az elobb eloallitott parameteres ertekek kozul
//egyesitjuk azokat, amiknek a kulonbozosege eseten rossz orbifoldokat kapnank.
//Az M matrix-fuggvenybeli ertekek kell megegyezzenek, ezert az egyutthatok
//legkisebb kozos tobbszorose lesz az uj parameter egyutthatoja.
void Dsym::filter_bad_orbifolds(void){
  if(dim!=3) return;

  for(int elhagy=1;elhagy<3;elhagy++)
    for(std::list<kisebbdim>::iterator komp=klist[elhagy].begin();
	komp!=klist[elhagy].end();komp++){
      int op;
      switch (elhagy){
	case 1: op=2;break;
	case 2: op=0;break;
	default: std::cerr<<"Nem jo"<<std::endl;
      }
      std::list<param>::iterator egyik,masik,harmadik; //ha mar harman vannak az jo
      egyik=plist.end();
      masik=plist.end();
      harmadik=plist.end();
      for(std::list<simplex*>::iterator szit=komp->szek.begin();
	  szit!=komp->szek.end();szit++)
	if( egyik==plist.end() ) 
	  egyik=params[(*szit)->sorszam[1]][op];
	else if (masik==plist.end()){ 
	  if( egyik!=params[(*szit)->sorszam[1]][op])
	    masik=params[(*szit)->sorszam[1]][op];
	}
	else if (egyik!=params[(*szit)->sorszam[1]][op] && 
	    masik!=params[(*szit)->sorszam[1]][op]){
	  harmadik=params[(*szit)->sorszam[1]][op];
	  break;
	}
      if(masik!=plist.end() && harmadik==plist.end()){
	if(egyik->eh != masik->eh)
	  std::cerr<<"Hm";
	egyik->eh=(egyik->eh<0 ? -1 : 1 ) *
	  boost::math::lcm(egyik->eh,masik->eh);

	for(egyik->min_ertek=0;abs(egyik->eh)*egyik->min_ertek < DEGENERATION_LIMIT;
	    egyik->min_ertek++);

	for(std::list<simplex*>::iterator szit=masik->szek.begin();
	    szit!=masik->szek.end();szit++)
	  params[(*szit)->sorszam[1]][op]=egyik;

	egyik->szek.splice(egyik->szek.end(),masik->szek);

	plist.erase(masik);
      }
    } 
}

//Dsym::change_param: A parameterlista megfelelo elemere mutato pointert adjuk
//at, hogy megsporoljuk a sorszam alapjan torteno azonositas plusz koltseget. A
//pointer altal mutatott parametert i-vel noveljuk (vagy csokkentjuk, ha i
//negativ,) es azon szimplexek matrixait megvaltoztatjuk, amelyekre hatassal van
//a parameter.
void Dsym::change_param(std::list<param>::iterator currparam,int i){
  currparam->ertek+=i;
  for(std::list<simplex*>::iterator szim=currparam->szek.begin();
      szim!=currparam->szek.end();szim++){
    int op0=currparam->op;
    int op1=op0+1;
    (*szim)->mx[op0][op1]=abs(currparam->eh)*currparam->ertek;
    (*szim)->mx[op1][op0]=abs(currparam->eh)*currparam->ertek;
  }
}

//A kovetkezo ket fuggveny csak az elnevezest segiti.
void Dsym::increase_param(std::list<param>::iterator currparam){ 
  if (currparam->ertek>100000)
    change_param(currparam,-currparam->ertek+currparam->min_ertek);
  else
    change_param(currparam,1);
}

void Dsym::decrease_param(std::list<param>::iterator currparam){
  if (currparam->ertek==currparam->min_ertek)
    change_param(currparam,inf++);
  else
    change_param(currparam,-1);
}

//Dsym::print_possible_mxes: Backtrack algoritmus segitsegevel felsoroljuk az
//osszes lehetseges parameter-beallitast, amire meg teljesul a 0. es a 3.
//operaciok elhagyasa utan a ket dimenzios szferikussag, vagy euklidesziseg.
//De csak azokat a lehetosegeket irjuk ki, amik minimalisak, azaz
//atsorszamozassal nem szeretnenk tobbszor megkapni ugyanazt a matrix-rendszert,
//illetve az elemeket nem szeretnenk, ha osszevonhatoak lennenek.
//Mindig csak a backtrack-fa legmelyen, azaz a lista vegen irjuk ki a
//matrixokat, mert oda minden matrix-rendszer eljut, mivel eloszor csak a
//vizsgalt sorszamot noveljuk, majd a parameter erteket.
//Mivel semmi sem garantalja, hogy nem letezhetnek vegtelen sorok, ezert a
//parametereknek adunk egy ertelmes felso korlatot (ez most 10.)
void Dsym::print_possible_mxes(std::list<param>::iterator currparam,std::ostream *out){
  if(currparam->ertek==10) return;	//A vegtelen sorokat nem szeretjuk...
  if(kdimsf(plist.end())<0) return;
  else{
    if(currparam==plist.end() && min()>0){
      mxnum++;
      *out<<"<tr>";
      for(int j=0;j<car;j++){
	*out<<"<td><table cellpadding=\"3\">";
	for(int i=0;i<dim+1;i++){
	  *out<<"<tr>";
	  for(int k=0;k<dim+1;k++){
	    *out<<"<td align=\"center\">";
	    *out<<csucsok[1][j]->mx[i][k];
	    *out<<"</td>";
	  }
	  *out<<"</tr>";
	}
	*out<<"</table></td>";
      }
      *out<<"</tr>"<<std::endl;
    }
  }
  if(currparam!=plist.end()){
    std::list<param>::iterator curr1=currparam;	//curr1=currparam+1
    curr1++;
    print_possible_mxes(curr1,out);
    increase_param(currparam);
    print_possible_mxes(currparam,out);
    decrease_param(currparam);
    return;
  }
  else return;
}

//Dsym::print_possible_params: Backtrack algoritmus segitsegevel felsoroljuk az
//osszes lehetseges parameter-beallitast, amire meg teljesul a 0. es a 3.
//operaciok elhagyasa utan a ket dimenzios szferikussag, vagy euklidesziseg.
//Mindig csak a backtrack-fa legmelyen, azaz a lista vegen irjuk ki a
//parametereket, mert oda minden allapot eljut, mivel eloszor csak a
//vizsgalt sorszamot noveljuk, majd a parameter erteket.
//Mivel semmi sem garantalja, hogy nem letezhetnek vegtelen sorok, ezert a
//parametereknek adunk egy ertelmes felso korlatot (ezt dinamikusan noveljuk,
//ahogy szukseges)
int Dsym::print_possible_params(std::list<param>::iterator currparam,std::ostream *out){
  //A vegtelen sorokat nem szeretjuk...  
  if(currparam->ertek==pmaxertek && currparam!=plist.end() ) return 0;
  if(kdimsf(plist.end())<0) return 0;
  else{
    if(currparam==plist.end()) 
      return print_possible_infs(1,plist.begin(),out);
    else {
      std::list<param>::iterator curr1=currparam;	//curr1=currparam+1
      curr1++;
      int curr_mxnum=print_possible_params(curr1,out);
      increase_param(currparam);
      curr_mxnum+=print_possible_params(currparam,out);
      decrease_param(currparam);
      return curr_mxnum;
    }
  }
}

//Dsym::print_possible_infs: Az elobb kiszamolt lehetseges parameterekben
//kiszurjuk a vegtelen lancokat a kovetkezo algoritmussal: eloszor minden 
//lehetseges modon behelyettesitunk a minimalis erteku parameterekbe vegtelent.
//Masodik korben ellenorizzuk, hogy lehetne-e meg vegtelent behelyettesiteni, ha
//nem akkor kiirjuk.
int Dsym::print_possible_infs(int pass,std::list<param>::iterator currparam,
    std::ostream* out) {
  if(pass==1){
    if(currparam==plist.end()){
      if(kdimsf(plist.end())==0)
	idcs=1;
      else
	idcs=0;
      return print_possible_infs(2,plist.begin(),out);
    }
    else{
      int ret=0;
      std::list<param>::iterator curr1=currparam;
      curr1++;
      ret+=print_possible_infs(1,curr1,out);
      if(kdimsf(currparam)>=0 && currparam->ertek==currparam->min_ertek){
	change_param(currparam,-currparam->ertek+(inf++));
	ret+=print_possible_infs(1,curr1,out);
	currparam->ertek=currparam->min_ertek;
	change_param(currparam,0);
      }
      return ret;
    }
  }
  else if(pass==2){
    std::ostringstream badorb;
    int m=min();
    if(currparam==plist.end())
      if(kdimsf(plist.end())<0 || m==-1)
	return 0;
      else{
	*out<<"<tr>";
	for(std::list<param>::iterator pit=plist.begin();pit!=plist.end();pit++){
	  if(pit->ertek>=100000)
	    *out<<"<td>&infin;</td>";
	  else{
	    if(pit->ertek==pmaxertek-1)
	      pmaxertek++;
	    *out<<"<td>"<<pit->ertek<<"</td>";
	  }
	}
	*out<<"<td>";
	if(m==-2){
	  *out<<"No (1";
	  for(std::list<int>::iterator it=osszevlist.begin();it!=osszevlist.end();
	      it++)
	    *out<<"="<<*it;
	  *out<<")";
	}
	*out<<"</td>";
	/*else if(min()==-1)
	 *out<<"<td>Atsorszamozhato</td>";*/
	*out<<"<td>";
	if(idcs==1)
	  *out<<"Exists";
	*out<<"</td>";
	*out<<"<td>";
	if(filter_bad_orbifolds03(&badorb)==0)
	  *out<<"Impossible";
	else
	  *out<<badorb.str();
	*out<<"</td>";
	//backup info
	*out<<"<td>";
	*out<<dim<<" "<<car;
	for(int j=0;j<dim+1;j++){
	  *out<<" (";
	  for(int i=0;i<car;i++){
	    int icsucsjszomszedja=csucsok[0][i]->szomszed[j]->sorszam[0];
	    if (icsucsjszomszedja == i)
	      *out<<i+1;
	    else if (i < icsucsjszomszedja)
	      *out<<"("<<i+1<<icsucsjszomszedja+1<<")";
	  }
	  *out<<")";
	}
	for(std::list<param>::iterator pit=plist.begin();pit!=plist.end();pit++){
	  *out<<" "<<pit->ertek;
	}
	*out<<"</td>";
	*out<<"</tr>"<<std::endl;
	return 1;
      }
    else{
      std::list<param>::iterator curr1=currparam;
      curr1++;
      //Ha meg egy valtozo helyere tudnank vegtelent irni, visszaterunk
      if(kdimsf(currparam)>=0 && currparam->ertek<=100000){
	return 0;
      }
      return print_possible_infs(2,curr1,out);
    }
  }
  else
    return 0;
}

//Dsym::filter_bad_orbifolds03: A 0. es a 3. operaciot elhagyva is
//keletkezhetnek rossz orbifoldok. Eloszor is ha euklideszi esetben vagyunk,
//akkor persze nem. Kulonben megnezzuk, hogy a nem parameteres matrixfuggveny
//ertekek hany peremet adnak (az ilyen peremek korul pontosan 1 fele szimplex
//lesz, iranyitott esetben biztos, hogy nem valodi forgascentrum lesz,
//iranyitatlan esetben, ha mar ket szimplex szerepel, akkor szinten kisimul.)
//Vegul a parameteres reszeket elemezzuk, hogy hany valodi forgascentrum es hany
//valodi perem lesz.
//Ha 1 vagy 2 perem es nulla forgas, illetve 1 vagy 2 forgas es nulla perem van,
//akkor megnezzuk, hogy a (mindket) centrum pontosan egyszeresen szerepel-e, ha
//igen es csak 1 van, akkor rossz orbifold, kulonben egyenlo kell legyen a 
//rendjuk.
int Dsym::filter_bad_orbifolds03(std::ostringstream *badorb){
  if(dim!=3) return -100;

  int ret=1;
  for(int op=0;op<dim+1;op+=dim)
    for(std::list<kisebbdim>::iterator komp=klist[op].begin();komp!=klist[op].end();
	komp++){
      if(! komp->iranyitott) continue; //mert a gomb iranyitott
      //Euklideszi esetekre figyelni kell:
      float sum=0;
      for(std::list<simplex*>::iterator currszim=komp->szek.begin();
	  currszim!=komp->szek.end();currszim++){
	for(int j=0;j<dim;j++)
	  if(!((op==0 && j==0) || (op==dim && j==dim-1)))
	    sum+=1.0/float((*currszim)->mx[j][j+1]);
	sum-=0.5;
      }
      //#define THRESH 0.0000001 korabban tortent
      if(sum<THRESH && sum>-THRESH) 
	continue;

      int peremek=0;
      int forgatasok=0;
      std::list<param> parperemek,parforgatasok;
      //nem parameteres polusok (itt iranyitott 2-es: semmi nincs, iranyitatlan
      //2-es: siktukrozes; iranyitatlan 1-es: 2-odrendu perem,iranyitott 1-es:
      //masodrendu forgas, amik erdekesek):
      for(std::list<simplex*>::iterator szit=komp->szek.begin();
	  szit!=komp->szek.end();szit++)
	for(int i=0;i<dim-1;i++)
	  if(i!=op)
	    for(int j=i+2;j<dim+1;j++)
	      if(j!=op){
		if((*szit)->szomszed[i]->sorszam[0]==(*szit)->sorszam[0] &&
		    (*szit)->szomszed[j]->sorszam[0]==(*szit)->sorszam[0])
		  peremek++;
		else if((*szit)->szomszed[i]->szomszed[j]->sorszam[0]==
		    (*szit)->sorszam[0])
		  forgatasok++;
	      }
      forgatasok=forgatasok/2; //mert mindig ketszer annyi van belole, mint kell
      if(peremek+forgatasok>=2) continue;

      int noop=0;
      if(op==dim)
	noop=dim-1;

      for(std::list<simplex*>::iterator szit=komp->szek.begin();
	  szit!=komp->szek.end();szit++)
	for(int i=0;i<dim;i++)
	  if(noop!=i && params[(*szit)->sorszam[1]][i]->ertek>1){
	    if(params[(*szit)->sorszam[1]][i]->eh>0){
	      if(find(parforgatasok.begin(),parforgatasok.end(), 
		    *params[(*szit)->sorszam[1]][i])==parforgatasok.end())
		parforgatasok.push_back(*params[(*szit)->sorszam[1]][i]);
	    }
	    else{
	      if(find(parperemek.begin(),parperemek.end(), 
		    *params[(*szit)->sorszam[1]][i])==parperemek.end())
		parperemek.push_back(*params[(*szit)->sorszam[1]][i]);
	    }
	  }

      if(
	  ( (peremek+parperemek.size()<=2 && 
	     forgatasok+parforgatasok.size()==0) ||
	    (peremek+parperemek.size()==0 && 
	     forgatasok+parforgatasok.size()<=2) ) &&
	  parperemek.size()+parforgatasok.size()>0)
      {
	std::list<param> parek;	//parameterek
	parek.splice(parek.end(),parperemek);
	parek.splice(parek.end(),parforgatasok);

	//Kell meg, hogy a parameterek reciprok osszegeben 1/p es 1/q legyen,
	//azaz az aktiv komponensben ugyanannyi szimplexhez tartozzon a
	//parameter, mint az egyutthatoja siktukrozes esetben, illetve ketszer
	//annyihoz forgatas esetben.
	int jo_orb=0;
	for(std::list<param>::iterator pit=parek.begin();pit!=parek.end();pit++){
	  int count=0;
	  for(std::list<simplex*>::iterator szit=pit->szek.begin();
	      szit!=pit->szek.end();szit++)
	    if(find(komp->szek.begin(),komp->szek.end(),*szit)!=
		komp->szek.end())
	      count++;
	  if(pit->eh<0 && count>abs(pit->eh))
	    jo_orb=1;
	  if(pit->eh>0 && count>2*pit->eh)
	    jo_orb=1;
	}


	if(jo_orb==0){
	  if(parek.size()==1){
	    if(peremek==1 || forgatasok==1)
	      *badorb<<parek.begin()->kar<<"=2";
	    else
	      return 0;
	  }
	  else
	    if(peremek!=1 && forgatasok!=1){
	      std::list<param>::iterator egyik,masik;
	      egyik=parek.begin();
	      masik=parek.begin();
	      masik++;
	      if(egyik->ertek>=100000 && masik->ertek>=100000)
		*badorb<<" "<<egyik->kar<<"="<<masik->kar;
	      else if(egyik->ertek>=100000)
		*badorb<<" "<<egyik->kar<<"="<<masik->ertek;
	      else if(masik->ertek>=100000)
		*badorb<<" "<<masik->kar<<"="<<egyik->ertek;
	      else if(egyik->ertek!=masik->ertek)
		return 0;
	    }
	}
      }
    }
  return ret;
}

//Dsym::print_param_mx: Kiiratjuk a parameteres matrix-fuggvenyt (html-kodkent.)
void Dsym::print_param_mx(std::ostream *out){
  *out<<"<tr>";
  for(int j=0;j<car;j++){
    *out<<"<td><table cellpadding=\"3\">";
    for(int i=0;i<dim+1;i++){
      *out<<"<tr>";
      for(int k=0;k<dim+1;k++){
	*out<<"<td align=\"center\">";
	if (abs(i-k)==1) {
	  std::list<param>::iterator cp=params[j][(i<k) ? i : k];
	  if(abs(cp->eh)>1) *out<<abs(cp->eh);
	  *out<<cp->kar;
	}
	else
	  *out<<csucsok[1][j]->mx[i][k];
	*out<<"</td>";
      }
      *out<<"</tr>";
    }
    *out<<"</table></td>";
  }
  *out<<"</tr>"<<std::endl;
}

bool operator == (Dsym::param a,Dsym::param b) {
  return a.kar==b.kar;
}

//backtrack: Felsoroljuk az osszes lehetseges el-kombinaciot, el-szam szerint
//rendezve
Dsymlista* saved;
int bt1;
void backtrack(int dim, int car) {
  int number_of_edges=0;
  bool found_new=true;
  stringdb* seen_prev=0;
  char temp[50];
  sprintf(temp,"d%dc%d_",dim,car);
  stringdb* seen_now=new stringdb(number_of_edges,temp);
  sprintf(temp,"d%dc%d_edge_",dim,car);
  stringdb* seen_edge_new=new stringdb(0,temp);
  Dsym D(dim,car);
  char dumpstr[500];
  D.atsorszamoz(1,0);
  D.dump(dumpstr,1);
  seen_now->append(dumpstr,500);
  while (found_new and ++number_of_edges) {
    found_new=false;
    if (seen_prev)
      delete seen_prev;
    seen_prev=seen_now;
    sprintf(temp,"d%dc%d_",dim,car);
    seen_now=new stringdb(number_of_edges,temp);

    //Osszes elozoleg megtalalt Dsym
    Dbc* prev_cursor=seen_prev->get_cursor();

    while (1){
      char* str=new char[2*seen_prev->keylength+1];
      int max;

      Dbt key;
      key.set_data(str);
      key.set_ulen(seen_prev->keylength);
      key.set_flags(DB_DBT_USERMEM);                                                                                     
      Dbt data(&max, sizeof(max));

      int ret = prev_cursor->get(&key, &data, DB_NEXT);
      if (ret == DB_NOTFOUND){
	delete[] str;
	break;
      }
      bt1++;
      std::istringstream dumpsstr;
      dumpsstr.str(std::string(str));
      delete[] str;
      Dsym D(&dumpsstr);
      max=*(int*)data.get_data();

      //Az osszes lehetseges elet hozzaadjuk
      for(int szin=0; szin<dim+1; szin++){
	if (szin==1) //Eltranzitiv eset
	  continue;
	for(int honnan=0; honnan<max+1 && honnan<car-1; honnan++)
	  for(int hova=honnan+1; hova<max+2 && hova<car; hova++)
	    if (D.elhozzaad(szin,honnan,hova)){
	      int newmax=max;
	      if (hova==max+1)
		newmax++;
	      bool voltmar=false;

	      //1. ellenorzes: seen_now-hoz erdemes-e hozzaadni
	      D.atsorszamoz(1,newmax);
	      char dumpstr[500];
	      D.dump(dumpstr,1);
	      if(seen_now->check(dumpstr))
		voltmar=true;
	      for (int i=2;i<newmax+2;i++){
		D.atsorszamoz(i,newmax);
		char dumpstr[500];
		D.dump(dumpstr,i);
		if(seen_now->check(dumpstr))
		  voltmar=true;
	      }
	      if (not voltmar) {
		found_new=true;
		//if ( not backtrack_breaks_uvw(&D,newmax))
		  seen_now->append(dumpstr,500,newmax);
	      }

	      //2. ellenorzes: seen_edge_new-hoz erdemes-e hozzaadni
	      if (newmax == car-1 && not voltmar && D.ellenoriz() == -1){
		seen_edge_new->append(dumpstr,500,0);
		int start=1;
		for(int i=1;i<car;i++){
		  D.atsorszamoz(i+1);
		  if(kisebb(&D,i+1,&D,start)==1)
		    start=i+1;
		}
		Dsym* ujD=D.save_with_start(start);
		if (saved->check(ujD)==0){
		  saved->append(ujD);
		}
		delete ujD;
	      }

	      //Vegul toroljuk az elet
	      D.eltorol(szin,honnan,hova);
	    }
      }
    }
  }
  delete seen_now;
  backtrack_edges(dim,car,seen_edge_new);
}

// backtrack_edges: In which combinations can we add the operations for
// edge-center-adjacencies
// szin is not needed (it's always 1)
int bt2;
void backtrack_edges(int dim,int car,stringdb* seen_edge_new) {
  int szin=1;
  int number_of_edges=0;
  bool found_new=true;
  stringdb* seen_prev=0;
  stringdb* seen_now=seen_edge_new;
  while (found_new and ++number_of_edges) {
    found_new=false;
    if (seen_prev)
      delete seen_prev;
    seen_prev=seen_now;
    char temp[50];
    sprintf(temp,"d%dc%d_edge_",dim,car);
    seen_now=new stringdb(number_of_edges,temp);
    sprintf(temp,"d%dc%d_edge_permut",dim,car);
    stringdb seen_now_permut(number_of_edges,temp);

    //Osszes elozoleg megtalalt Dsym
    Dbc* prev_cursor=seen_prev->get_cursor();

    while (1){
      char* str=new char[seen_prev->keylength];
      int max;

      Dbt key;
      key.set_data(str);
      key.set_ulen(seen_prev->keylength);
      key.set_flags(DB_DBT_USERMEM);                                                                                     
      Dbt data(&max, sizeof(max));

      int ret = prev_cursor->get(&key, &data, DB_NEXT);
      if (ret == DB_NOTFOUND){
	delete[] str;
	break;
      }
      bt2++;
      std::istringstream dumpsstr;
      dumpsstr.str(std::string(str));
      delete[] str;
      Dsym D(&dumpsstr);

      //Az osszes lehetseges elet hozzaadjuk
      for(int honnan=0; honnan < car-1; honnan++)
	for(int hova=honnan+1; hova<car; hova++)
	  if (D.elhozzaad(szin,honnan,hova)){
	    bool voltmar=false;

	    //1. ellenorzes: seen_new-hoz erdemes-e hozzaadni
	    D.atsorszamoz(1);
	    char dumpstr[500];
	    D.dump(dumpstr,1);
	    if (seen_now_permut.check(dumpstr))
	      voltmar=true;
	    else{
	      found_new=true;
	      //if (not backtrack_breaks_uvw(&D,car-1))
		seen_now->append(dumpstr,500,0);
	      seen_now_permut.append(dumpstr,500,0);
	      for (int i=2;i<=car;i++){
		D.atsorszamoz(i);
		char dumpstr[500];
		D.dump(dumpstr,i);
		seen_now_permut.append(dumpstr,500,0);
	      }
	    }

	    //2. ellenorzes: seen_edge_new-hoz erdemes-e hozzaadni
	    if (not voltmar && D.ellenoriz() == -1){
	      int start=1;
	      for(int i=1;i<car;i++){
		D.atsorszamoz(i+1);
		if(kisebb(&D,i+1,&D,start)==1)
		  start=i+1;
	      }
	      Dsym* ujD=D.save_with_start(start);
	      if (saved->check(ujD)==0){
		saved->append(ujD);
	      }
	      delete ujD;
	    }

	    //Vegul toroljuk az elet
	    D.eltorol(szin,honnan,hova);
	  }
    }
  }
  delete seen_now;
}

/* Megvizsgaljuk, hogy az aktualis diagramban van nem m_ij=2 erteku
 * operacio-par, kezdopont harmas.
 */
bool backtrack_breaks_uvw(Dsym* D,int max){
  // Nem mukodik
  int dim=D->dim;
  simplex*** csucsok=D->csucsok;
  for (int r=0;r<=max;r++)
    for (int i=0;i<dim-1;i++)
      for (int i1=i+2;i1<dim+1;i1++){
	/* Ha 4 lepesbol egyszer sem fordultunk vissza es nem is jutottunk a
	 * kiindulasi pontra, akkor nem jo. "|_|-"-szeru minta.
	 */
	simplex* csucs=csucsok[0][r];
	int j,j1;
	if ( csucs != csucs->szomszed[i] ){
	  j=i;
	  j1=i1;
	}
	else {
	  j=i1;
	  j1=i;
	}
	int steps=0;
	while (steps++ < 4 && 
	    csucs != csucs->szomszed[j] &&
	    csucs != csucsok[0][r]){
	  csucs=csucs->szomszed[j];
	  int temp=j;
	  j=j1;
	  j1=temp;
	}
	if (steps == 4)
	  return true;
      }
  return false;
}

//Dsymlista::check: Az adott elem szerepel-e mar a listainkban?
int Dsymlista::check(Dsym* element){
  std::ostringstream dumpsstr;
  element->dump(&dumpsstr);
  std::string dumpstr=dumpsstr.str();
  int a;

  Dbt key((void*)dumpstr.c_str(), dumpstr.size() + 1);
  Dbt data(&a, sizeof(a));

  int ret=fastdb.get(NULL, &key, &data, 0);
  if (ret == DB_NOTFOUND)
    return 0;
  else
    return 1;
}

//Dsymlista::check_with_start: Az adott elem szerepel-e mar a listainkban?
int Dsymlista::check_with_start(Dsym* element,int start){
  Dsym* ujD=element->save_with_start(start);
  int ret=check(ujD);
  delete ujD;
  return ret;
}

//Dsymlista::check_sorted: Az adott elem hanyadik a rendezett listaban
int Dsymlista::check_sorted(Dsym* element,int start){
  std::ostringstream dumpsstr;
  Dsym* ujD=element->save_with_start(start);
  int ret=-3;
  if (check(ujD)){
    ujD->dump(&dumpsstr);
    std::string dumpstr=dumpsstr.str();
    int a;

    Dbt key((void*)dumpstr.c_str(), dumpstr.size() + 1);
    Dbt data(&a, sizeof(a));

    int ret=sorteddb.get(NULL, &key, &data, 0);
    if (ret == DB_NOTFOUND)
      ret=-1;
    else
      ret=a;
  }
  else{
    // Elvileg nem erhetunk ide
    throw -2;
    ret=-2;
  }
  delete ujD;
  return ret;
}


//Dsymlista::append: az uj elemet beszurjuk az adatbazisokba. (Ha a fastdb-ben mar szerepel, akkor nem.)
void Dsymlista::append(Dsym* new_element){
  std::ostringstream dumpsstr;
  new_element->dump(&dumpsstr);
  std::string dumpstr=dumpsstr.str();
  int a=0;

  Dbt key((void*)dumpstr.c_str(), dumpstr.size() + 1);
  Dbt data(&a, sizeof(a));

  int ret = fastdb.put(NULL, &key, &data, DB_NOOVERWRITE);
  if (ret == 0) {
    std::cerr<<"\r"<<count++;
    sorteddb.put(NULL, &key, &data, DB_NOOVERWRITE);

    if ((int)dumpstr.size() + 1 > keylength)
      keylength=2*dumpstr.size();
  }
  // Else: do nothing...
  //else if (ret == DB_KEYEXIST) {
  //}
}                                                                                                                 

Dsymlista::Dsymlista(int dimin,int carin):
  dim(dimin),
  car(carin),
  fastdb(NULL,0),
  sorteddb(NULL,0),
  current(0),
  keylength(2048),
  count(0)
{
  Dsymlista(dimin,carin,std::string("example"));
}

Dsymlista::Dsymlista(int dimin,int carin,std::string fn):
  dim(dimin),
  car(carin),
  filename_base(fn),
  fastdb(NULL,0),
  sorteddb(NULL,0),
  current(0),
  keylength(2048),
  count(0)
{ 
  try {
    Db a(NULL,0);
    Db b(NULL,0);
    a.remove((filename_base + "_fast.db").c_str(), NULL, 0);
    b.remove((filename_base + "_sort.db").c_str(), NULL, 0);
  }
  catch (DbException& e){
    ;
  }

  fastdb.set_cachesize(1,500*1024*1024,1);
  sorteddb.set_cachesize(1,500*1024*1024,1);
  sorteddb.set_bt_compare(&compare_d);

  fastdb.open(NULL, ("data/" + filename_base + "_fast.db").c_str(), NULL, DB_HASH, DB_CREATE, 0);
  sorteddb.open(NULL, ("data/" + filename_base + "_sort.db").c_str(), NULL, DB_BTREE, DB_CREATE, 0);
}

Dsymlista::~Dsymlista(void){
  if (current!=0)
    current->close();
  fastdb.sync(0);
  fastdb.close(0);
  sorteddb.sync(0);
  sorteddb.close(0);
}

int compare_d(Db *dbp, const Dbt *a, const Dbt *b){
  char* dstr1=new char[a->get_size()];
  char* dstr2=new char[b->get_size()];

  memcpy(dstr1, a->get_data(), a->get_size());
  memcpy(dstr2, b->get_data(), b->get_size());                                                                    

  std::istringstream dsstr1, dsstr2;
  dsstr1.str(dstr1);
  dsstr2.str(dstr2);

  Dsym* dobj1=new Dsym(&dsstr1);
  Dsym* dobj2=new Dsym(&dsstr2);

  int ret=kisebb(dobj2,0,dobj1,0);
  delete dobj1;
  delete dobj2;
  delete[] dstr1;
  delete[] dstr2;
  return ret;
}

std::pair<Dsym*,int> Dsymlista::getnextsorted(void){
  if (current == 0){
    reset_cursor();
    generate_ordered_numbering();
  }

  int ssz=0;
  char* str;
  int ret;

  str=new char[keylength];

  Dbt key;
  Dbt data(&ssz, sizeof(ssz));

  key.set_data(str);
  key.set_ulen(keylength);
  key.set_flags(DB_DBT_USERMEM);

  ret = current->get(&key, &data, DB_NEXT);
  if (ret == DB_NOTFOUND){
    delete[] str;
    return std::pair<Dsym*,int>(NULL,0);
  }

  std::istringstream dumpsstr;
  dumpsstr.str(std::string(str));
  Dsym* output=new Dsym(&dumpsstr);
  ssz=*(int*)data.get_data();

  delete[] str;
  return std::pair<Dsym*,int>(output,ssz);
}

void Dsymlista::reset_cursor(void){
  sorteddb.cursor(0,&current,DB_CURSOR_BULK);
}

// Create ordered numbering and set keylength if needed
void Dsymlista::generate_ordered_numbering(void){
  Dbc *cursor;
  sorteddb.cursor(0,&cursor,DB_CURSOR_BULK);

  int ssz=1;
  int temp;
  char* str;
  int ret;

  str=new char[keylength];

  Dbt key;
  Dbt data(&temp, sizeof(temp));

  key.set_data(str);
  key.set_ulen(keylength);
  key.set_flags(DB_DBT_USERMEM);

  try {
    ret = cursor->get(&key, &data, DB_NEXT);
    while (ret != DB_NOTFOUND){
      temp=ssz++;
      data.set_data(&temp);
      cursor->put(&key, &data, DB_CURRENT);
      fastdb.put(NULL, &key, &data, 0);
      ret = cursor->get(&key, &data, DB_NEXT);
    }
  }
  catch (DbMemoryException& e){
    if (e.get_errno() == DB_BUFFER_SMALL) {
      keylength=2*key.get_size();
      generate_ordered_numbering();
    }
    else {
      throw 1;
    }
  }

  delete[] str;
}

//Dsymlista::print_html: Hasonloan mint az elobb, kiirjuk html-be az adatokat,
//annyi pluszt teszunk hozza, hogy kulon fajlba (currD) kiirjuk az osszes
//lehetseges matrix-rendszert is.
void Dsymlista::print_html(void){
  std::ostringstream filename;
  char dir[20];
  sprintf(dir,"d%dc%d",dim,car);
  mkdir(dir,0755);
  filename<<"d"<<dim<<"c"<<car<<"/index.html";
  std::ofstream html_file;
  html_file.open(filename.str().c_str());

  //Header
  html_file<<"<html>"<<std::endl<<"<body>"<<std::endl;

  while (true) {
    std::pair<Dsym*,int> it=getnextsorted();
    Dsym* curr=it.first;
    if(curr==NULL)
      break;
    curr->ellenoriz();
    int ssz=it.second;
    std::cerr<<"\r"<<ssz;
    curr->dual=0;
    curr->dualis();
    inf=1000000;
    curr->create_params();
    curr->create_kdim();
    curr->filter_bad_orbifolds();

    std::ostringstream currDname;
    currDname<<"d"<<dim<<"c"<<car<<"/"<<ssz<<".html";
    std::ofstream currD;
    currD.open(currDname.str().c_str());

    currD<<"<html>"<<std::endl<<"<body>"<<std::endl<<"<table border=\"2\">"<<std::endl;
    currD<<"<caption>"<<ssz<<"</caption>"<<std::endl;

    std::ostringstream xfigfile;
    xfigfile<<"d"<<dim<<"c"<<car<<"/"<<ssz<<".fig";
    curr->write_xfig(xfigfile.str());
    std::ostringstream figtojpg;
    figtojpg<<"fig2dev -L jpeg "<<"d"<<dim<<"c"<<car<<"_"<<ssz<<".fig "
      <<"d"<<dim<<"c"<<car<<"_"<<ssz<<".jpg";
    //system(figtojpg.str().c_str());
    //remove(xfigfile.str().c_str());
    currD<<"<tr><td><img src=\""<<ssz<<".jpg\"/></td>"
      <<std::endl;

    currD<<"<td><table border=\"1\">"<<std::endl;
    curr->print_param_mx(&currD);
    currD<< "<tr><td colspan=\""<<car<<"\">Number of " << dim-1 << 
      " dimensional components: "<< "(" <<curr->osszefuggo(0);
    for (int i=1;i<dim+1;i++) currD<<","<<curr->osszefuggo(i);
    currD<< ")</td></tr>"<<std::endl;
    currD<<"<tr><td colspan=\""<<car<<"\">"<<"Dual: ";
    int dualchk;
    for (int i=1;i<car+1;i++) {
      curr->dual->atsorszamoz(i);
      dualchk=check_with_start(curr->dual,i);
      if (dualchk!=0) {
	break;
      }
    }
    if (dualchk==0)
      currD<<"Not found";
    else
      if (dualchk==ssz)
	currD<<"Selfdual";
      else
	currD<<"D."<<dualchk;
    currD<<"</td></tr>"<<std::endl;

    currD<<"<tr><td colspan=\""<<car<<"\">"<<"Parameters: <br>"<<std::endl;
    for(int i=0;i<dim;i++) {
      currD << "(" <<i<<","<<i+1<<") ";
      curr->print_params(i,&currD);
      if (i<dim-1) currD <<"<br>"<<std::endl;
      else currD<<std::endl<<"</td></tr>"<<std::endl;
    }
    currD<<"</table></tr></table><br><table border=\"1\" cellpadding=\"3\">"
      <<std::endl<<"<caption>Possible parameter values:</caption>"<<std::endl
      <<"<thead><tr>";
    for(std::list<Dsym::param>::iterator pit=curr->plist.begin();
	pit!=curr->plist.end();pit++)
      currD<<"<td align=center>"<<pit->kar<<"</td>";
    currD<<"<td>Maximal</td><td>Ideal vertex</td>"  //Plusz infok
      <<"<td>Good orbifold criteria</td>";
    currD<<"<td>Backup info</td>";
    currD<<"</tr></thead><tbody>"<<std::endl;
    curr->mxnum=curr->print_possible_params(curr->plist.begin(),
	&currD);
    currD<<std::endl<<"</tbody></table></body></html>"<<std::endl;

    //html file:
    html_file<<"<table border=\"2\">"<<std::endl;
    html_file<<"<caption>"<<ssz<<"</caption>"<<std::endl;
    html_file<<"<tr><td><img src=\""<<ssz
      <<".jpg\"/></td>"<<std::endl;
    html_file<<"<td><table border=\"1\">"<<std::endl;
    html_file<<"<tr><td colspan=\""<<car<<"\"><a href=\""<<ssz<<".html"
      <<"\">Number of matrices: "<<curr->mxnum<<"</a></td></tr>"<<std::endl;
    curr->print_param_mx(&html_file);
    html_file << "<tr><td colspan=\""<<car<<"\">Number of " << dim-1 << 
      " dimensional components: "<< "(" <<curr->osszefuggo(0);
    for (int i=1;i<dim+1;i++) html_file <<","<<curr->osszefuggo(i);
    html_file << ")</td></tr>"<<std::endl;
    html_file<<"<tr><td colspan=\""<<car<<"\">"<<"Dual: ";
    if (dualchk==0)
      html_file<<"Not "<<dim<<"-connected";
    else
      if (dualchk==ssz)
	html_file<<"Selfdual";
      else
	html_file<<"D."<<dualchk;
    html_file<<"</td></tr>"<<std::endl;
    html_file<<"<tr><td colspan=\""<<car<<"\">"<<"Parameters: <br>"<<std::endl;
    for(int i=0;i<dim;i++) {
      html_file << "(" <<i<<","<<i+1<<") ";
      curr->print_params(i,&html_file);
      if (i<dim-1) html_file <<"<br>"<<std::endl;
      else html_file<<std::endl<<"</td></tr>"<<std::endl;
    }
    html_file<<std::endl<<"</table></tr></table><br>"<<std::endl;
    delete curr;
  }

  //Footer
  html_file<<"</body>"<<std::endl<<"</html>"<<std::endl;

}

//Simple string hash database
void stringdb::create(char* dbname)
{
  //  filename_base=filename;
  filename_base=std::string(dbname);
  keylength=2048;
  db.set_cachesize(1,500*1024*1024,1);

  char temp[50];
  sprintf(temp,"data/%s",dbname);
  db.open(NULL, temp, NULL, DB_HASH, DB_CREATE, 0);
  db.cursor(0,&cursor,DB_CURSOR_BULK);
}

stringdb::stringdb(char* dbname):
  db(NULL,0)
{
  create(dbname);
}

stringdb::stringdb(int num, char* dbname):
  db(NULL,0)
{
  char fn[50];
  sprintf(fn,"%s%d",dbname,num);
  create(fn);
}

stringdb::~stringdb(void){
  cursor->close();
  db.sync(0);
  db.close(0);
}

Dbc* stringdb::get_cursor(void){
  cursor->close();
  db.cursor(0,&cursor,DB_CURSOR_BULK);
  return cursor;
}

void stringdb::append(std::string what, int a){
  Dbt key((void*)what.c_str(), what.size() + 1);
  Dbt data(&a, sizeof(a));

  if (keylength <= (int)what.size() + 1)
    keylength = 2*what.size();

  db.put(NULL, &key, &data, DB_NOOVERWRITE);
}

void stringdb::append(std::string what){
  append(what,0);
}

void stringdb::append(char* what, int size, int a){
  Dbt key((void*)what, size);
  Dbt data(&a, sizeof(a));

  if (keylength <= size)
    keylength = 2*size;

  db.put(NULL, &key, &data, DB_NOOVERWRITE);
}

void stringdb::append(char* what, int size){
  append(what,size,0);
}

bool stringdb::check(std::string what){
  int a;

  Dbt key((void*)what.c_str(), what.size() + 1);
  Dbt data(&a, sizeof(a));

  int ret=db.get(NULL, &key, &data, 0);
  if (ret == 0)
    return true;
  else
    return false;
};


//xfig konstruktor: feltoltjuk a valtozokat, es rajzolunk egy el nelkuli
//multigrafot. rad a korok sugara, koordx illetve koordy pedig a korok x es y
//koordinatainak listaja.
xfig::xfig(std::string filename,int dim1,int car1):dim(dim1),car(car1) {
  rad=200;
  koordx=new int[car];
  koordy=new int[car];
  for (int i=0;i<car;i++) {
    koordx[i]=int(round(4000+1500*cos(PI-i*2*PI/car)));
    koordy[i]=int(round(4000-1500*sin(PI-i*2*PI/car)));
  }
  outfile.open(filename.c_str());
  outfile << "#FIG 3.2  Produced by xfig version 3.2.5-alpha5" <<std::endl;
  outfile << "Landscape" <<std::endl;
  outfile << "Center" <<std::endl;
  outfile << "Metric" <<std::endl;
  outfile << "A4" <<std::endl;
  outfile << "100.00" <<std::endl;
  outfile << "Single" <<std::endl;
  outfile << "-2" <<std::endl;
  outfile << "1200 2" <<std::endl;
  for(int i=0;i<car;i++){
    create_circle(i);
    create_numtext(i);
  }
}

//xfig::create_circle: az n-edik csucs korenek megrajzolasa.
void xfig::create_circle(int n) {
  outfile << "1 3 0 1 0 7 49 -1 20 0.000 1 0.0000 "<<koordx[n]<<" "<<koordy[n]
    <<" "<<rad<<" "<<rad<<" "<<koordx[n]<<" "<<koordy[n]<<" "
    <<koordx[n]+rad<<" "<<koordy[n]<<std::endl;
}

//xfig::create_numtext: az n-edik sorszam beleirasa.
void xfig::create_numtext(int n) {
  outfile << "4 0 0 48 -1 0 12 0.0000 4 150 105 "<<koordx[n]-52<<" "
    <<koordy[n]+75<< " "<<n+1<<"\\001"<<std::endl;
}

//xfig::create_line: n0 es n1 koze szin "szinu" (sokkal inkabb mintaju) el
//berajzolasa. Fontos, hogy ha ket pont koze tobbfele elet rajzolunk, azok ne
//mosodjanak ossze, ezert az elek iranyara meroleges iranyban fix ertekkel
//vannak eltolva a kulonbozo szinu elek.
void xfig::create_line(int n0,int n1,int szin) {
  int diff=(szin*2*rad/dim-rad)*3/4;
  float length=sqrt(pow(koordx[n1]-koordx[n0],2)+pow(koordy[n1]-koordy[n0],2));
  int diffy=(koordx[n1]-koordx[n0])*diff/(int)length;
  int diffx=-(koordy[n1]-koordy[n0])*diff/(int)length;
  int minta,mintaszelesseg;
  switch (szin){
    case 0:
      minta=2;
      mintaszelesseg=10;
      break;
    case 1:
      minta=1;
      mintaszelesseg=12;
      break;
    case 2:
      minta=0;
      mintaszelesseg=0;
      break;
    case 3:
      minta=3;
      mintaszelesseg=20;
      break;
    default:
      minta=szin;
      mintaszelesseg=20;
  }
  outfile << "2 1 "<<minta<<" 4 0 7 50 -1 -1 "<<mintaszelesseg
    <<".000 0 0 -1 0 0 2\n\t "<<koordx[n0]+diffx<<" "<<koordy[n0]+diffy
    <<" "<<koordx[n1]+diffx<<" "<<koordy[n1]+diffy<<std::endl;
}

//xfig destruktor: listak torlese, fajl bezarasa.
xfig::~xfig(void) {
  delete[] koordx;
  delete[] koordy;
  outfile.close();
}


//main: Megkerdezi a dimenziot es a bonyolultsagot, majd "y" hatasara felsorolja
//a felteteleknek megfelelo multigrafokat, es html fajlokba irja a kinyert
//informaciot. Ha nem felsorolast szeretnenk, akkor lehet probalkozni elek
//hozzaadasaval, es az informaciok kinyeresevel (ez a resz kisse elmaradt a 
//fejlesztes soran, miutan debugolasra minimalis szukseg volt a vegehez 
//kozeledve.)
int main(int,char**,char**){
  int dim;
  int car;
  dim=3;
  std::cout << "Cardinality: ";std::cin >> car;std::cout<<car<<std::endl;
  char buffer [50];
  sprintf (buffer, "d%dc%d", dim, car);
  std::string fn(buffer);
  saved=new Dsymlista(dim,car,fn);
  bt1=0;
  bt2=0;
  backtrack(dim,car);
  //for (Dsymlinklist* it=saved->first;it!=NULL;it=it->next)
  //it->ssz=ujssz++;
  std::cout<<"saved"<<std::endl;
  std::cout<<saved->count<<std::endl;
  std::cout<<"backtrack_steps: "<< bt1<<", backtrack_edges_steps: "<<bt2<<std::endl;
  saved->print_html();
  delete saved;
  return 0;
}
