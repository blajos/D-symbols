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
  pmaxertek=4;		//kisebbre nem választhatjuk, es majd megnõ
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

int Dsym::fok(int csucs,int sorszamozas){
  // fokszam: i-edik tipusu el eseten 2^i-t hozzaadunk az eddigiekhez (0-val
  // kezdunk)
  int fok=0;
  for(int j=0;j<dim+1;j++)
    if(csucsok[sorszamozas][csucs]->szomszed[j] != csucsok[sorszamozas][csucs]){
      fok+=1 << j;
    }
  return fok;
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

void Dsym::uvw1(void) {
  for (int r=0;r<car;r++) {				//atlotol tavolabbi
    for (int i=0;i<dim-1;i++)
      for (int i1=i+2;i1<dim+1;i1++){
	csucsok[0][r]->mx[i][i1]=2;
	csucsok[0][r]->mx[i1][i]=2;
      }
    for (int i=0;i<dim+1;i++) csucsok[0][r]->mx[i][i]=1;     //atlo
    for (int i=0;i<dim;i++) {				       //atlo szomszedok
      int u=1;
      simplex* csucs=csucsok[0][r];
      while (csucs->szomszed[i]->szomszed[i+1] != csucsok[0][r]) {
	u++;
	csucs=csucs->szomszed[i]->szomszed[i+1];
      }
      csucsok[0][r]->mx[i][i+1]=u;
      csucsok[0][r]->mx[i+1][i]=u;
    }
  }
}

bool Dsym::lehet_eltranzitiv(void){
  if (dim!=3)
    throw "Nope.";
  for (int r=0;r<car;r++) {
    if (csucsok[1][r]->mx[2][3] > 6)
      return false;
  }
  return true;
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
  if (!lehet_eltranzitiv()){
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

//Dsym::write_svg: out streambe letrehozza a megfelelo tartalmat. (Csak az eleket kell
//berajzolni, mert minden mas megegyezik az azonos dimenzioju es elemszamu
//multigrafok kozott, es minden elt csak egy iranyba szabad megrajzolni,
//kulonben a mintak fedesbe kerulhetnek.)
void Dsym::write_svg(std::ostream* out) {
  Svg rajz(dim,car);
  for(int i=0;i<car-1;i++) 
    for(int j=0;j<dim+1;j++)
      if(csucsok[1][i]->szomszed[j]->sorszam[1]>csucsok[1][i]->sorszam[1])
	rajz.add_line(i,csucsok[1][i]->szomszed[j]->sorszam[1],j);
  rajz.print_html(out);
}

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

Base::Base(std::istream* in){
  *in >> dim >> car;
}

Base::~Base(){
  return;
}

int Base::dump(std::ostream* out){
  *out << dim << " " << car;
  return 0;
}

int Base::print_html(std::ostream*){
  return 0;
}

Svg::Svg(int dimin,int carin):
  Base(dimin,carin)
{
  rad=20;
  fontsize=28;
  size=car*50;
  if (size < 300)
    size=300;
  koordx=new int[car];
  koordy=new int[car];
  for (int i=0;i<car;i++) {
    koordx[i]=int(round(size/2+(size/2-rad)*cos(PI-i*2*PI/car)));
    koordy[i]=int(round(size/2-(size/2-rad)*sin(PI-i*2*PI/car)));
  }
}

void Svg::create_circle(int n,std::ostream* out) {
  *out << "<circle cx=\"" << koordx[n] << "\" cy=\"" << koordy[n] <<
    "\" r=\"" << rad << "\" fill=\"white\" stroke=\"black\" stroke-width=\"1\"/>";
}

void Svg::create_numtext(int n,std::ostream* out) {
  *out << "<text x=\"" << koordx[n] << "\" y=\"" <<
    (koordy[n]+round(fontsize/3)) << "\" font-size=\""<< fontsize <<
    "\" text-anchor=\"middle\" dominant-baseline=\"mathematical\">" << n+1 <<
    "</text>";
}

void Svg::create_line(int n0,int n1,int szin,std::ostream* out) {
  int diff=(szin*2*rad/dim-rad)*3/4;
  float length=sqrt(pow(koordx[n1]-koordx[n0],2)+pow(koordy[n1]-koordy[n0],2));
  int diffy=(koordx[n1]-koordx[n0])*diff/(int)length;
  int diffx=-(koordy[n1]-koordy[n0])*diff/(int)length;
  std::string style;
  switch (szin){
    case 0:
      style="stroke-dasharray:2,10";
      break;
    case 1:
      style="stroke-dasharray:8,8";
      break;
    case 2:
      style="";
      break;
    case 3:
      style="stroke-dasharray:8,8,2,10";
      break;
    default:
      style="stroke-dasharray:8,8,2,10";
      for (int i=3;i<szin;i++) style+=",2,10";
  }
  *out << "<line style=\"" << style <<
    "\" x1=\"" << koordx[n0]+diffx <<
    "\" y1=\"" << koordy[n0]+diffy <<
    "\" x2=\"" << koordx[n1]+diffx <<
    "\" y2=\"" << koordy[n1]+diffy <<
    "\" stroke=\"black\" stroke-width=\"1\"/>";
}

Svg::~Svg(void) {
  delete[] koordx;
  delete[] koordy;
}

int Svg::print_html(std::ostream* out){
  *out << "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 "<< size<< " "<< size<<"\" width=\"" <<
    size << "px\" height=\""<<size<<"px\" version=\"1.1\">";

  for(std::list<std::vector<int> >::iterator it=lines.begin();it!=lines.end();it++)
    create_line((*it)[0],(*it)[1],(*it)[2],out);

  for(int i=0;i<car;i++){
    create_circle(i,out);
    create_numtext(i,out);
  }
  *out << "</svg>";
  return 0;
}

int Svg::add_line(int n0,int n1,int szin){
  std::vector<int> *a=new std::vector<int>;
  a->push_back(n0);
  a->push_back(n1);
  a->push_back(szin);
  lines.push_back(*a);
  return 0;
}

AnimSvg::AnimSvg(int dimin,int carin):
  Svg(dimin,carin)
{
  show=new std::list<int>**[dim+1];
  hide=new std::list<int>**[dim+1];
  for (int d=0;d<dim+1;d++){
    show[d]=new std::list<int>*[car];
    hide[d]=new std::list<int>*[car];
    for (int c1=0;c1<car;c1++){
      show[d][c1]=new std::list<int>[car];
      hide[d][c1]=new std::list<int>[car];
    }
  }
}

AnimSvg::~AnimSvg(){
  for (int d=0;d<dim+1;d++){
    for (int c1=0;c1<car;c1++){
      delete[] show[d][c1];
      delete[] hide[d][c1];
    }
    delete[] show[d];
    delete[] hide[d];
  }
  delete[] show;
  delete[] hide;
}

void AnimSvg::add_line(int n0,int n1,int color,int timetick){
  show[color][n0][n1].push_back(timetick);
}

void AnimSvg::remove_line(int n0,int n1,int color,int timetick){
  hide[color][n0][n1].push_back(timetick);
}

void AnimSvg::add_sleep(int timetick){
  sleep.push_back(timetick);
}

int AnimSvg::print_html(std::ostream* out){
  *out << "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 "<< size<< " "<< size<<"\" width=\"" <<
    size << "px\" height=\""<<size<<"px\" version=\"1.1\">"<<std::endl;

  print_sleep(out);
  for(int szin=0;szin<dim+1;szin++)
    for(int n0=0;n0<car;n0++)
      for(int n1=n0+1;n1<car;n1++)
	print_line(n0,n1,szin,out);

  for(int i=0;i<car;i++){
    create_circle(i,out);
    create_numtext(i,out);
  }
  *out << "</svg>";
  return 0;
}

void AnimSvg::create_circle(int n,std::ostream* out) {
  *out << "  <circle cx=\"" << koordx[n] << "\" cy=\"" << koordy[n] <<
    "\" r=\"" << rad << "\" fill=\"white\" stroke=\"black\" stroke-width=\"1\"/>"<<std::endl;
}

void AnimSvg::create_numtext(int n,std::ostream* out) {
  *out << "  <text x=\"" << koordx[n] << "\" y=\"" <<
    (koordy[n]+round(fontsize/3)) << "\" font-size=\""<< fontsize <<
    "\" text-anchor=\"middle\" dominant-baseline=\"mathematical\">" << n+1 <<
    "</text>"<<std::endl;
}

void AnimSvg::print_sleep(std::ostream* out) {
  *out << "  <rect x=\"0\" y=\"0\" width=\"" << size << "px\" height=\""<<size<<"px\" fill=\"green\" opacity=\"0\">" << std::endl;
  for(std::list<int>::iterator it=sleep.begin();
      it!=sleep.end(); it++){
    std::string prev;
    if (*it == 0)
      prev="0";
    else {
      std::ostringstream temp;
      temp << "A" << (*it)-1 << ".end";
      prev=temp.str();
    }
    *out << "    <animate id=\"A"<< *it+9 <<"\" attributeName=\"opacity\" values=\"0;1;0\" dur=\"2s\" begin=\"" << prev << "\" fill=\"freeze\"/>" << std::endl;
  }
  *out << "  </rect>" << std::endl;
}

void AnimSvg::print_line(int n0,int n1,int szin,std::ostream* out) {
  int diff=(szin*2*rad/dim-rad)*3/4;
  float length=sqrt(pow(koordx[n1]-koordx[n0],2)+pow(koordy[n1]-koordy[n0],2));
  int diffy=(koordx[n1]-koordx[n0])*diff/(int)length;
  int diffx=-(koordy[n1]-koordy[n0])*diff/(int)length;
  std::string style;
  switch (szin){
    case 0:
      style="stroke-dasharray:2,10";
      break;
    case 1:
      style="stroke-dasharray:8,8";
      break;
    case 2:
      style="";
      break;
    case 3:
      style="stroke-dasharray:8,8,2,10";
      break;
    default:
      style="stroke-dasharray:8,8,2,10";
      for (int i=3;i<szin;i++) style+=",2,10";
  }
  *out << "  <line style=\"" << style <<
    "\" x1=\"" << koordx[n0]+diffx <<
    "\" y1=\"" << koordy[n0]+diffy <<
    "\" x2=\"" << koordx[n1]+diffx <<
    "\" y2=\"" << koordy[n1]+diffy <<
    "\" stroke=\"black\" stroke-width=\"1\" opacity=\"0\">"<<std::endl;
  for(std::list<int>::iterator it=show[szin][n0][n1].begin();
      it!=show[szin][n0][n1].end(); it++){
    std::string prev;
    if (*it == 0)
      prev="0";
    else {
      std::ostringstream temp;
      temp << "A" << (*it)-1 << ".end";
      prev=temp.str();
    }
    *out << "    <animate id=\"A"<< *it <<"\" attributeName=\"opacity\" values=\"0;1\" dur=\"1s\" begin=\"" << prev << "\" fill=\"freeze\"/>" << std::endl;
  }
  for(std::list<int>::iterator it=hide[szin][n0][n1].begin();
      it!=hide[szin][n0][n1].end(); it++){
    std::string prev;
    if (*it == 0)
      prev="0";
    else {
      std::ostringstream temp;
      temp << "A" << (*it)-1 << ".end";
      prev=temp.str();
    }
    *out << "    <animate id=\"A"<< *it <<"\" attributeName=\"opacity\" values=\"1;0\" dur=\"0.5s\" begin=\"" << prev << "\" fill=\"freeze\"/>" << std::endl;
  }
  *out << "  </line>" << std::endl;
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
    klist[currkis].clear();

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

//Dsym::kdimsf: Igaz-e, hogy ha elhagyjuk a 0. vagy az 3. operaciot, akkor a
//kisebb dimenzios szimbolum szferikus vagy euklideszi sikon megvalosulhat.
//Komponensenkent szamolunk: osszeadjuk az 1/m01+1/m12-1/2 ertekeket, ha nem
//negativ, akkor igaz az allitas. Mivel lebegopontos szamitasokat vegzunk, nem
//vizsgalhatunk ==0-t.
//Kisebb dimenzios reszek szferikusak: 1, van Euk: 0, van hip:-1
//Miert nem nezzuk meg a 1-2 operaciokat: ekkor pl. 1/m23+1/2+1/2 kell nagyobb
//legyen 1-nel, ami mindig igaz.
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

//Dsym::kdimsf i<max eseten igaz-e, hogy a 0-3 operaciok elhagyasaval kapott
//"ki nem logo" resz D-szimbolumok megfelelnek a felteteleinknek (0,3 eseten
//euklideszi vagy szferikus; 1,2 eseten szferikus)
//1-2 esetet nem is kene nezni, mert az uvw feltetel miatt mindig szferikus lesz
bool Dsym::kdimsf(int max){
  if (dim!=3) return -100;
  uvw1();
  create_kdim();
  for (int elhagy=0;elhagy<dim+1;elhagy++){
    int currkis=elhagy;
    for(std::list<kisebbdim>::iterator curr=klist[currkis].begin();
	curr!=klist[currkis].end();curr++){
      float sum=0;
      bool skip=false;
      for(std::list<simplex*>::iterator currszim=curr->szek.begin();
	  currszim!=curr->szek.end();currszim++)
	if((*currszim)->sorszam[0] < max){
	  for(int j=0;j<dim;j++)
	    for(int j1=j+1;j1<dim+1;j1++){
	      if(j!=elhagy and j1!=elhagy)
		sum+=1.0/float((*currszim)->mx[j][j1]);
	    }
	  sum-=1;
	}
	else{
	  skip=true;
	  break;
	}
      if(not skip){
	if(sum<-THRESH){
	  if (elhagy > 0 and elhagy < dim){
	    std::cout << "Nem lehet";
	  }
	  return false;
	}
	if(sum<THRESH && sum>-THRESH and elhagy > 0 and elhagy < dim){
	  std::cout << "Nem lehet";
	  return false;
	}
      }
    }
  }
  return true;
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

//backtrack: Felsoroljuk az osszes lehetseges car elemszamu, dim+1 elszinu
//multigrafot, ahol teljesul, hogy egy csucsban nincs 2 azonos szinu, de
//kulonbozo el. A lehetseges eleket fel tudjuk sorolni x y z szamharmasokkal,
//ahol: x a szin, y a kezdo, z a vegpont. (x,y,z>=0; y<z; x<dim+1; y,z<car)
//Ezeket soroljuk fel (akar tetszoleges sorrendben is lehetne, de az alabbi tunt
//a legegyszerubbnek) eloszor el nelkul, majd az elt berajzolva. Ezzel megkapjuk
//a backtrack-fat. Nem optimalis: minden multigrafot leirunk, igy ha az adott
//kezdoponthoz beallitjuk a csucsok sorszamat, akkor (n-1)! azonos grafot
//kapunk. A rendezesi algoritmust sokkal jobban ki kene hasznalni.
//
// Modositas az eredetihez kepest: Nem adunk hozza "1" szinu elet, azokat csak a
// levelekben fogjuk tovabb csocsalni.
int bt;
int bt1;
int bt2;
long long bt0;
int timetick;
AnimSvg* animation;
void backtrack(Dsym* D,Dsymlista* saved,int szin,int honnan,int hova) {
  int car=D->car;
  int dim=D->dim;
  bt0++;
  if (honnan==car-1){	//a vegen megallunk
    //ellenorzesek
    bt++;
    int joe=D->ellenoriz();
    //mentes, ha kell
    if (joe==-1) {
      bt1++;
      int start=1;
      for(int i=1;i<car;i++){
	D->atsorszamoz(i+1);
	if(kisebb(D,i+1,D,start)==1)
	  start=i+1;
      }
      //std::cout << " " << bt1 << " " << start << std::endl;
      Dsym* ujD=D->save_with_start(start);
      if (saved->check(ujD)==0){
	//saved->append(ujD);
	bt2++;
	backtrack_edges(D,saved,0,1);
      }
      delete ujD;
    }
    return;
  }

  // Nem a legkisebb 0-ad foku csucsot nem akarjuk hozzaadni, a nala 1-el kisebb
  // csucs kell nem ures legyen, kiveve 0-1 kozti elso elet.
  bool erdemes=false;
  if (hova == 1)
    erdemes = true;
  if ( hova >= 2 ){
    for(int j=0;j<dim+1;j++)
      if (D->csucsok[0][hova-1]->szomszed[j] != D->csucsok[0][hova-1]){
	erdemes=true;
	break;
      }
  }

  //ha erdemes elt hozzaadni, ujra meghivjuk onmagunkat
  if (erdemes && D->elhozzaad(szin,honnan,hova)){
    //animation->add_line(honnan, hova, szin, timetick++);
    backtrack(D,saved,szin,honnan,hova);
    D->eltorol(szin,honnan,hova);
    //animation->remove_line(honnan, hova, szin, timetick++);
  }

  if(hova+1<car) backtrack(D,saved,szin,honnan,hova+1);
  else if(szin+1<dim+1 && szin==0) backtrack(D,saved,szin+2,honnan,honnan+1); //eltranzitiv esetben kihagyjuk az 1-es szint
  else if(szin+1<dim+1) backtrack(D,saved,szin+1,honnan,honnan+1);
  else {
    if (honnan > 0){
      //Heurisztika: 
      // mivel honnan-nal nem nagyobb vegu elet mar biztos nem adunk hozza, es az
      // osszes permutaciot megnezzuk, ezert a 0..honnan rendszer lehet fok
      // szerint csokkeno sorrendben
      // Es ehhez eleg annyi, hogy (honnan-1)-nel nem nagyobb honnan foka
      //
      // Update: Ez egy szep elkepzeles, de sajnos nem igaz, hogy minden
      // permutaciot megnezunk: ha ket egyforman legnagyobb foku csucs nem
      // szomszedos, akkor nem fogjuk oket egymas utan bevalasztani. 

      int szukseges_fok=D->fok(0,0);
      int fok=D->fok(honnan,0);
      if ( fok <= szukseges_fok && not backtrack_breaks_uvw(D,honnan) && not
	  backtrack_breaks_eltranzitiv(D,honnan) and D->kdimsf(honnan+1))
	backtrack(D,saved,0,honnan+1,honnan+2);
    }
    else{
      backtrack(D,saved,0,honnan+1,honnan+2);
    }
  }
}

// backtrack_edges: In which combinations can we add the operations for
// edge-center-adjacencies
// szin is not needed (it's always 1)
int bte;
int bte1;
int bte2;
long long bte0;
void backtrack_edges(Dsym* D,Dsymlista* saved,int honnan,int hova) {
  int car=D->car;
  int szin=1;
  bte0++;
  if (honnan==car-1){	//a vegen megallunk
    //ellenorzesek
    bte++;
    int joe=D->ellenoriz();
    //mentes, ha kell
    if (joe==-1) {
      bte1++;
      int start=1;
      for(int i=1;i<car;i++){
	D->atsorszamoz(i+1);
	if(kisebb(D,i+1,D,start)==1)
	  start=i+1;
      }
      //std::cout << " " << bt1 << " " << start << std::endl;
      Dsym* ujD=D->save_with_start(start);
      if (saved->check(ujD)==0){
	bte2++;
	//animation->add_sleep(timetick);
	timetick+=10;
	saved->append(ujD);
      }
      delete ujD;
    }
    return;
  }

  //ha tudunk elt hozzaadni, ujra meghivjuk onmagunkat
  if (D->elhozzaad(szin,honnan,hova)){
    //animation->add_line(honnan, hova, szin, timetick++);
    backtrack_edges(D,saved,honnan,hova);
    D->eltorol(szin,honnan,hova);
    //animation->remove_line(honnan, hova, szin, timetick++);
  }

  if(hova+1<car) backtrack_edges(D,saved,honnan,hova+1);
  else 
    if (not backtrack_breaks_uvw(D,honnan) and D->kdimsf(honnan+1)){
      backtrack_edges(D,saved,honnan+1,honnan+2);
    }
}

/* Megvizsgaljuk, hogy az aktualis diagramban van nem m_ij=2 erteku
 * operacio-par, kezdopont harmas.
 */
bool backtrack_breaks_uvw(Dsym* D,int honnan){
  int dim=D->dim;
  int max=honnan+1;
  simplex*** csucsok=D->csucsok;
  /* A kovetkezo mintak felismerese: |_ |_|
     Ha talalunk nem szomszedos operacio parokkal legalabb 2 hosszu lancot
     (mindket vege hurok es mindket vege kisebb, mint max)
     Vagy legalabb 5 hosszu barmit; az rossz. */
  for (int r=0;r<=max;r++)
    for (int i=0;i<dim-1;i++)
      for (int i1=i+2;i1<dim+1;i1++){
	simplex* csucs=csucsok[0][r];
	int j,j1;
	bool chain=true;
	if ( csucs != csucs->szomszed[i] ){
	  if ( csucs != csucs->szomszed[i1])
	    chain=false;
	  j=i;
	  j1=i1;
	}
	else if ( csucs != csucs->szomszed[i1]){
	  j=i1;                                                                                                      
	  j1=i;                                                                                                      
	}
	else { //i, i1 nem vezet ki csucsbol
	  continue;
	}
	int steps=0;
	do {
	  csucs=csucs->szomszed[j];
	  int temp=j;
	  j=j1;
	  j1=temp;
	} while (++steps <= 5 && csucs != csucs->szomszed[j] && csucs != csucsok[0][r]);
	if (steps == 5 or (steps >= 2 and chain and csucs->sorszam[0]<=max))
	  return true;
      }
  return false;
}

// Legfeljebb 6 alakzat talalkozhat egy elnel (mert a dualis a laptranzitiv,
// ahol legfeljebb 6-szog lapok lehetnek)
bool backtrack_breaks_eltranzitiv(Dsym* D,int honnan){
  int dim=D->dim;
  if (dim!=3)
    throw "That doesn't work.";
  int max=honnan+1;
  simplex*** csucsok=D->csucsok;
  /* Ha talalunk (2,3) operacio parokkal legalabb 7 hosszu lancot, vagy 13
   * hosszu barmit; az rossz.*/
  for (int r=0;r<=max;r++){
    int i=2;
    int i1=3;
    simplex* csucs=csucsok[0][r];
    int j,j1;
    bool chain=true;
    if ( csucs != csucs->szomszed[i] ){
      if ( csucs != csucs->szomszed[i1])
	chain=false;
      j=i;
      j1=i1;
    }
    else if ( csucs != csucs->szomszed[i1]){
      j=i1;                                                                                                      
      j1=i;                                                                                                      
    }
    else { //i, i1 nem vezet ki csucsbol
      continue;
    }
    int steps=0;
    do {
      csucs=csucs->szomszed[j];
      int temp=j;
      j=j1;
      j1=temp;
    } while (++steps <= 5 && csucs != csucs->szomszed[j] && csucs != csucsok[0][r]);
    if (steps == 13 or (steps >= 7 and chain and csucs->sorszam[0]<=max))
      return true;
  }
  return false;
}

//Dsymlista::sync: Flush sorted db to disk
void Dsymlista::sync(void){
  sorteddb.sync(0);
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
  output_limit(100),
  verify_ok(false),
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
  output_limit(100),
  verify_ok(false),
  count(0)
{ 
  try {
    Db a(NULL,0);
    a.set_bt_compare(&compare_d);
    if (a.verify((filename_base + "_sort.db").c_str(), NULL, NULL, 0) == 0)
      verify_ok=true;
    a.close(0);
  }
  catch (DbException& e){
    ;
  }

  if ( not verify_ok ){
    try {
      Db a(NULL,0);
      Db b(NULL,0);
      a.remove((filename_base + "_fast.db").c_str(), NULL, 0);
      b.remove((filename_base + "_sort.db").c_str(), NULL, 0);
    }
    catch (DbException& e){
      ;
    }
  }

  fastdb.set_cachesize(0,200*1024*1024,2);
  sorteddb.set_cachesize(0,200*1024*1024,2);
  sorteddb.set_bt_compare(&compare_d);

  fastdb.open(NULL, (filename_base + "_fast.db").c_str(), NULL, DB_BTREE, DB_CREATE, 0);
  sorteddb.open(NULL, (filename_base + "_sort.db").c_str(), NULL, DB_BTREE, DB_CREATE, 0);
}

Dsymlista::~Dsymlista(void){
  if (current!=0)
    current->close();
  fastdb.close(0);
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

void Dsymlista::create_directories(std::string fbase, int count, int prefix=0 ){
  //A kovetkezo jellegu mappa strukturat keszitjuk el:
  //fbase/0/0/234.html
  //fbase/0/2000/2347.html
  //fbase/1000000/12000/10124321.html
  // Igy minden mappaban max 1000 fajl lesz, ezt meg kenyelmesen elviselik a
  // file rendszerek
  // log10((float)count)/3 1 millio eseten 2-t ad vissza
  if ( count == 0 )
    return;
  int levels=(int)(log((float)count)/log((float)output_limit));
  if (levels == 0)
    return;
  else{
    for(int i=0;i<=count;i+=pow(output_limit,levels)){
      std::ostringstream dirst;
      dirst << fbase << prefix+i << "/";
      std::string dir=dirst.str();
      mkdir(dir.c_str(),0755);
      int max=pow(output_limit,levels)-1;
      if (max+i > count)
	max=count-i;
      create_directories(dir,max,prefix+i);
    }
  }
}

std::string Dsymlista::get_filename(std::string fbase, int num){
  std::ostringstream out;
  out << get_path(fbase,num) << num;
  return out.str();
}

std::string Dsymlista::get_path(std::string fbase, int num){
  //A kovetkezo jellegu mappa strukturat keszitjuk el:
  //fbase/0/0/234.html
  //fbase/0/2000/2347.html
  //fbase/1000000/12000/10124321.html
  // Igy minden mappaban max 1000 fajl lesz, ezt meg kenyelmesen elviselik a
  // file rendszerek
  // log10((float)count)/3 1 millio eseten 2-t ad vissza
  std::ostringstream out;
  out << fbase;
  int levels=(int)(log((float)count)/log((float)output_limit));

  for (int level=levels;level > 0;level--){
    out << (int)(((int)(num/pow(output_limit,level))) * pow(output_limit,level)) << "/";
  }
  return out.str();
}

//Dsymlista::print_html: Hasonloan mint az elobb, kiirjuk html-be az adatokat,
//annyi pluszt teszunk hozza, hogy kulon fajlba (currD) kiirjuk az osszes
//lehetseges matrix-rendszert is.
void Dsymlista::print_html(std::string filebase){
  std::string path="-1";
  std::string prev_path="-1";
  std::string fbase;
  if(filebase.length() != 0)
    fbase=filebase+"/";
  else
    fbase=std::string(".")+"/";
  mkdir(fbase.c_str(),0755);
  std::ostringstream filename;
  filename<<fbase<<"index.html";
  std::ofstream index_file;
  std::ofstream list_file;
  index_file.open(filename.str().c_str());
  create_directories(fbase,count);

  //Header
  index_file<<"<html>"<<std::endl<<"<body>"<<std::endl;
  index_file<<"<h1>Dimension: "<<dim<<"</h1>"<<std::endl;
  index_file<<"<h2>Cardinality: "<<car<<"</h2>"<<std::endl;

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

    std::string fn=get_filename("",ssz);
    path=get_path("",ssz);
    if ( path != prev_path ){
      prev_path=path;
      if(list_file.is_open()){
	list_file<<"</body>"<<std::endl<<"</html>"<<std::endl;
	list_file.close();
      }
      list_file.open((fbase+path+"list.html").c_str());
      list_file<<"<html>"<<std::endl<<"<body>"<<std::endl;
      index_file << "<a href=\"" << (path+"list.html") << "\">From " << ssz << " to " << ssz+output_limit-1 << "<br>" << std::endl;
    }
    std::ofstream currD;
    currD.open((fbase+fn+".html").c_str());

    currD<<"<html>"<<std::endl<<"<body>"<<std::endl;
    currD<<"<h1>"<<ssz<<"</h1>"<<std::endl;

    curr->write_svg(&currD);

    currD<<"<table border=\"1\">"<<std::endl;
    curr->print_param_mx(&currD);
    currD<< "<tr><td colspan=\""<<car<<"\">Number of " << dim-1 << 
      " dimensional components: "<< "(" <<curr->osszefuggo(0);
    for (int i=1;i<dim+1;i++) currD<<","<<curr->osszefuggo(i);
    currD<< ")</td></tr>"<<std::endl;
    currD<<"<tr><td colspan=\""<<car<<"\">"<<"Dual: ";
    int dualchk=0;
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
    currD<<"</table><br><table border=\"1\" cellpadding=\"3\">"
      <<std::endl<<"<caption>Possible parameter values:</caption>"<<std::endl
      <<"<thead><tr>";
    for(std::list<Dsym::param>::iterator pit=curr->plist.begin();
	pit!=curr->plist.end();pit++)
      currD<<"<td align=center>"<<pit->kar<<"</td>";
    currD<<"<td>Maximal</td><td>Ideal vertex</td>"  //Plusz infok
      <<"<td>Good orbifold criteria</td>";
    currD<<"</tr></thead><tbody>"<<std::endl;
    curr->mxnum=curr->print_possible_params(curr->plist.begin(),
	&currD);
    currD<<std::endl<<"</tbody></table></body></html>"<<std::endl;

    //html file:
    list_file<<"<h1>"<<ssz<<"</h1>"<<std::endl;
    curr->write_svg(&list_file);
    list_file<<"<table border=\"1\">"<<std::endl;
    list_file<<"<tr><td colspan=\""<<car<<"\"><a href=\""<<ssz<<".html"
      <<"\">Number of matrices: "<<curr->mxnum<<"</a></td></tr>"<<std::endl;
    curr->print_param_mx(&list_file);
    list_file << "<tr><td colspan=\""<<car<<"\">Number of " << dim-1 << 
      " dimensional components: "<< "(" <<curr->osszefuggo(0);
    for (int i=1;i<dim+1;i++) list_file <<","<<curr->osszefuggo(i);
    list_file << ")</td></tr>"<<std::endl;
    list_file<<"<tr><td colspan=\""<<car<<"\">"<<"Dual: ";
    if (dualchk==0)
      list_file<<"Not "<<dim<<"-connected";
    else
      if (dualchk==ssz)
	list_file<<"Selfdual";
      else
	list_file<<"D."<<dualchk;
    list_file<<"</td></tr>"<<std::endl;
    list_file<<"<tr><td colspan=\""<<car<<"\">"<<"Parameters: <br>"<<std::endl;
    for(int i=0;i<dim;i++) {
      list_file << "(" <<i<<","<<i+1<<") ";
      curr->print_params(i,&list_file);
      if (i<dim-1) list_file <<"<br>"<<std::endl;
      else list_file<<std::endl<<"</td></tr>"<<std::endl;
    }
    list_file<<std::endl<<"</table><br><hr>"<<std::endl;
    delete curr;
  }

  //Footer
  list_file<<"</body>"<<std::endl<<"</html>"<<std::endl;
  index_file<<"</body>"<<std::endl<<"</html>"<<std::endl;
  std::cerr<<std::endl;

}

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
  Dsym* D=new Dsym(dim,car);
  char buffer [50];
  sprintf (buffer, "d%dc%d", dim, car);
  std::string fn(buffer);
  Dsymlista* saved=new Dsymlista(dim,car,fn);
  if( not saved->verify_ok ){
    bt=0;
    bt1=0;
    bt2=0;
    bt0=0;
    bte=0;
    bte1=0;
    bte2=0;
    bte0=0;
    timetick=0;
    animation=new AnimSvg(dim,car);
    backtrack(D,saved,0,0,1);
    saved->sync();
    std::cout<<std::endl<<"saved: "<<saved->count<<std::endl;
    std::cout<<"Statistics:"<<std::endl;
    std::cout<<"Backtrack steps: "<<bt0<<", leaves: "<<bt<<", valid: "<<bt1<<", saved:"<<bt2<<std::endl;
    std::cout<<"    edges steps: "<<bte0<<", leaves: "<<bte<<", valid: "<<bte1<<", saved:"<<bte2<<std::endl;
    std::ofstream animf;
    //animf.open(fn+"_anim.svg");
    //animation->print_html(&animf);
    delete animation;
  }
  saved->print_html(fn);
  delete D;
  delete saved;
  return 0;
}
