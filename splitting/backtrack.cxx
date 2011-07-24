#include "backtrack.hxx"
#include <boost/math/common_factor.hpp>
#include <iostream>
#include <fstream>
#include <list>
#include <algorithm>
#include <set>
#include <math.h>
#include <sstream>
using namespace std;
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
Dsym::Dsym(int dimin, int carin):
  dim(dimin),
  car(carin),
  involucio(0),
  mxnum(0),
  pmaxertek(4)		//kisebbre nem választhatjuk, es majd megnõ
{
  //az eredetivel egyutt az atsorszamozasok szama
  csucsok=new simplex**[car+1];	
  for (int k=0;k<car+1;k++) csucsok[k]=new simplex*[car];
  for (int i=0;i<car;i++) csucsok[0][i]=new simplex(dim,car,i);
  params=new list<param>::iterator*[car];
  for (int i=0;i<car;i++) params[i]=new list<param>::iterator[dim];
  dual=0;
  klist=new list<kisebbdim>[dim+1];
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
  list<int> nemelerheto,utolso,uj;
  list<int>::iterator hanyadik;
  int hanydarab=0;

  for(int i=0;i<car;i++) nemelerheto.push_back(i);	//nemelerheto feltoltese
  while (!nemelerheto.empty()) {
    uj.clear();
    uj.push_back(*nemelerheto.begin());
    hanydarab++;
    while (!uj.empty()) {
      utolso.clear();
      for ( list<int>::iterator p = uj.begin(); p != uj.end(); ++p ) 
	utolso.push_back(*p);	//uj->utolso
      uj.clear();
      for (list<int>::iterator r = utolso.begin(  ); r != utolso.end(  ); ++r ){
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
    //cout << "Nem osszefuggo" <<endl;
    return 1; //csak az osszefuggoek erdekesek
  }
  int sor=sorszamozas();
  //if (sor==1) return 1; 
  if (!uvw()){
    //cout << "Nem uvw0" <<endl;
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
    //cout << "Nem tudok elt hozzaadni, mert nincs torolve az elozo" << endl;
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
    cerr << "Nincs el (szin, honnan, hova): " <<szin<<" "<<honnan<<" "<<hova
						<< endl;
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
  cout << "=========================" << endl;
  for (int j=0;j<dim+1;j++) {
    for(int k=0;k<*involucio[j][car];k++){
      cout << "[" << involucio[j][k][0]+1 << "," << involucio[j][k][1]+1 <<"] ";
    }
    cout << endl;
  }
  cout << "Number of " << dim-1 << " dimensional components: ";
  cout << "("<<osszefuggo(0);
  for (int i=1;i<dim+1;i++) cout <<","<<osszefuggo(i);
  cout << ")"<<endl;
}

//Dsym::print_params: Kiirja az op0 es az op0+1 operaciokhoz tartozo
//parametereket, letter1 betuvel kezdve az out kimenetre.
void Dsym::print_params(int op0,ostream* out) {
  for (list<param>::iterator currparam=plist.begin();currparam!=plist.end();
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
void Dsym::write_xfig(string file) {
  xfig out(file,dim,car);
  for(int i=0;i<car-1;i++) 
    for(int j=0;j<dim+1;j++)
      if(csucsok[1][i]->szomszed[j]->sorszam[1]>csucsok[1][i]->sorszam[1])
	out.create_line(i,csucsok[1][i]->szomszed[j]->sorszam[1],j);
}


//Dsym::create_kdim: klist[i] a D-szimbolum i. operaciojanak elhagyasaval
//keletkezo komponensek listaja, egy-egy komponensben a benne levo szimplexek
//listajaval. 
//A komponenseket, mint korabban az osszefuggoseg kerdesenel is, tobbszori
//szelessegi bejarassal hatarozzuk meg.
//Vegul megnezzuk, hogy a komponens iranyitott vagy iranyitatlan feluletet ad
//meg.
void Dsym::create_kdim(void) {
  for (int elhagy=0;elhagy<dim+1;elhagy++){
    list<int> nemelerheto,utolso,uj;
    list<int>::iterator hanyadik;

    for(int i=0;i<car;i++) nemelerheto.push_back(i);	//nemelerheto feltoltese
    while (!nemelerheto.empty()) {
      kisebbdim curr;
      curr.op=elhagy;
      uj.clear();
      uj.push_back(*nemelerheto.begin());
      while (!uj.empty()) {
	utolso.clear();
	for ( list<int>::iterator p = uj.begin(); p != uj.end(); ++p ) 
	  utolso.push_back(*p);//uj->utolso
	uj.clear();
	for ( list<int>::iterator r = utolso.begin(); r!=utolso.end(); ++r ){
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
      list<simplex*> egyik, masik, *akt,*nemakt,nemfelsorolt;
      for(list<simplex*>::iterator szit=curr.szek.begin();
	  szit!=curr.szek.end();szit++) 
	if(szit==curr.szek.begin())
	  egyik.push_back(*szit);
	else
	  nemfelsorolt.push_back(*szit);

      akt=&egyik;
      nemakt=&masik;
      curr.iranyitott=true;
      while( ! nemfelsorolt.empty()){
	for(list<simplex*>::iterator szit=akt->begin();
	    szit!=akt->end();szit++)
	  for(int i=0;i<dim+1;i++)
	    if(i!=elhagy && (*szit)->szomszed[i]!=(*szit)){
	      if(find(akt->begin(),akt->end(),(*szit)->szomszed[i])!=
		  akt->end())
		curr.iranyitott=false;
	      list<simplex*>::iterator hanyadik=find(nemfelsorolt.begin(),
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
	list<simplex*> *temp;
	temp=akt;
	akt=nemakt;
	nemakt=temp;
      }
      //Vegul ujra atfutjuk, mert az utolso lepesnel keletkezhettek elek
      //nemakt-ban, ami most akt.
      for(list<simplex*>::iterator szit=akt->begin();
	  szit!=akt->end();szit++)
	for(int i=0;i<dim+1;i++)
	  if(i!=elhagy && (*szit)->szomszed[i]!=(*szit)){
	    if(find(akt->begin(),akt->end(),(*szit)->szomszed[i])!=
		akt->end())
	      curr.iranyitott=false;
	  }

      klist[elhagy].push_back(curr);
    }
  }
}

//Dsym::kdimsf: Igaz-e, hogy ha elhagyjuk a 0. vagy az n. operaciot, akkor a
//kisebb dimenzios szimbolum szferikus vagy euklideszi sikon megvalosulhat.
//Komponensenkent szamolunk: osszeadjuk az 1/m01+1/m12-1/2 ertekeket, ha nem
//negativ, akkor igaz az allitas. Mivel lebegopontos szamitasokat vegzunk, nem
//vizsgalhatunk ==0-t.
//Kisebb dimenzios reszek szferikusak: 1, van Euk: 0, van hip:-1
int Dsym::kdimsf(list<param>::iterator inf_param){
  if (dim!=3) return -100;
  int ret=1;
  for (int elhagy=0;elhagy<dim+1;elhagy+=dim){
    int currkis=elhagy;
    for(list<kisebbdim>::iterator curr=klist[currkis].begin();
	curr!=klist[currkis].end();curr++){
      float sum=0;
      for(list<simplex*>::iterator currszim=curr->szek.begin();
	  currszim!=curr->szek.end();currszim++){
	for(int j=0;j<dim;j++)
	  if(!((elhagy==0 && j==0) || (elhagy==dim && j==dim-1)))
	    if(params[(*currszim)->sorszam[1]][j]!=inf_param)
	      sum+=1.0/float((*currszim)->mx[j][j+1]);
	    else
	      sum+=10*THRESH; // Ezzel jelezzuk, hogy a vegtelen lancnak az elemei erdekesek, nem a hatarerteke
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
  int atsorsz=0;
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
    list<int> nemelerheto;
    list<int>::iterator hanyadik;

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
      list<param>::iterator lastp=plist.end();
      lastp--;
      for(list<simplex*>::iterator szit=currparam.szek.begin();
	  szit!=currparam.szek.end();szit++)
	params[(*szit)->sorszam[1]][op0]=lastp;
    }
  }
  for(list<param>::iterator currparam=plist.begin();currparam!=plist.end();
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
    for(list<kisebbdim>::iterator komp=klist[elhagy].begin();
	komp!=klist[elhagy].end();komp++){
      int op;
      switch (elhagy){
	case 1: op=2;break;
	case 2: op=0;break;
	default: cerr<<"Nem jo"<<endl;
      }
      list<param>::iterator egyik,masik,harmadik; //ha mar harman vannak az jo
      egyik=plist.end();
      masik=plist.end();
      harmadik=plist.end();
      for(list<simplex*>::iterator szit=komp->szek.begin();
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
	  cerr<<"Hm";
	egyik->eh=(egyik->eh<0 ? -1 : 1 ) *
	  boost::math::lcm(egyik->eh,masik->eh);

	for(egyik->min_ertek=0;abs(egyik->eh)*egyik->min_ertek < DEGENERATION_LIMIT;
	    egyik->min_ertek++);

	for(list<simplex*>::iterator szit=masik->szek.begin();
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
void Dsym::change_param(list<param>::iterator currparam,int i){
  currparam->ertek+=i;
  for(list<simplex*>::iterator szim=currparam->szek.begin();
      szim!=currparam->szek.end();szim++){
    int op0=currparam->op;
    int op1=op0+1;
    (*szim)->mx[op0][op1]=abs(currparam->eh)*currparam->ertek;
    (*szim)->mx[op1][op0]=abs(currparam->eh)*currparam->ertek;
  }
}

//A kovetkezo ket fuggveny csak az elnevezest segiti.
void Dsym::increase_param(list<param>::iterator currparam){ 
  if (currparam->ertek>100000)
    throw "Csapda1";
  //change_param(currparam,-currparam->ertek+currparam->min_ertek);
  else
    change_param(currparam,1);
}

void Dsym::decrease_param(list<param>::iterator currparam){
  if (currparam->ertek==currparam->min_ertek)
    throw "Csapda2";
  //change_param(currparam,inf++);
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
void Dsym::print_possible_mxes(list<param>::iterator currparam,ostream *out){
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
      *out<<"</tr>"<<endl;
    }
  }
  if(currparam!=plist.end()){
    list<param>::iterator curr1=currparam;	//curr1=currparam+1
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
int Dsym::print_possible_params(list<param>::iterator currparam,ostream *out){
  //A vegtelen sorokat nem szeretjuk...  
  if(currparam->ertek==pmaxertek && currparam!=plist.end() ) return 0;
  if(kdimsf(plist.end())<0) return 0;
  else{
    if(currparam==plist.end()) 
      return print_possible_infs(1,plist.begin(),out);
    else {
      list<param>::iterator curr1=currparam;	//curr1=currparam+1
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
int Dsym::print_possible_infs(int pass,list<param>::iterator currparam,
    ostream* out) {
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
      list<param>::iterator curr1=currparam;
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
    ostringstream badorb;
    int m=min();
    if(currparam==plist.end())
      if(kdimsf(plist.end())<0 || m==-1)
	return 0;
      else{
	*out<<"<tr>";
	for(list<param>::iterator pit=plist.begin();pit!=plist.end();pit++){
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
	  for(list<int>::iterator it=osszevlist.begin();it!=osszevlist.end();
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
	for(list<param>::iterator pit=plist.begin();pit!=plist.end();pit++){
	  *out<<" "<<pit->ertek;
	}
	*out<<"</td>";
	*out<<"</tr>"<<endl;
	return 1;
      }
    else{
      int ret=0;
      list<param>::iterator curr1=currparam;
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
int Dsym::filter_bad_orbifolds03(ostringstream *badorb){
  if(dim!=3) return -100;

  int ret=1;
  for(int op=0;op<dim+1;op+=dim)
    for(list<kisebbdim>::iterator komp=klist[op].begin();komp!=klist[op].end();
	komp++){
      if(! komp->iranyitott) continue; //mert a gomb iranyitott
      //Euklideszi esetekre figyelni kell:
      float sum=0;
      for(list<simplex*>::iterator currszim=komp->szek.begin();
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
      list<param> parperemek,parforgatasok;
      //nem parameteres polusok (itt iranyitott 2-es: semmi nincs, iranyitatlan
      //2-es: siktukrozes; iranyitatlan 1-es: 2-odrendu perem,iranyitott 1-es:
      //masodrendu forgas, amik erdekesek):
      for(list<simplex*>::iterator szit=komp->szek.begin();
	  szit!=komp->szek.end();szit++)
	for(int i=0;i<dim-1;i++)
	  if(i!=op)
	    for(int j=i+2;j<dim+1;j++)
	      if(j!=op)
		if((*szit)->szomszed[i]->sorszam[0]==(*szit)->sorszam[0] &&
		    (*szit)->szomszed[j]->sorszam[0]==(*szit)->sorszam[0])
		  peremek++;
		else if((*szit)->szomszed[i]->szomszed[j]->sorszam[0]==
		    (*szit)->sorszam[0])
		  forgatasok++;
      forgatasok=forgatasok/2; //mert mindig ketszer annyi van belole, mint kell
      if(peremek+forgatasok>=2) continue;

      int noop=0;
      if(op==dim)
	noop=dim-1;

      for(list<simplex*>::iterator szit=komp->szek.begin();
	  szit!=komp->szek.end();szit++)
	for(int i=0;i<dim;i++)
	  if(noop!=i && params[(*szit)->sorszam[1]][i]->ertek>1)
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

      if(
	  ( (peremek+parperemek.size()<=2 && 
	     forgatasok+parforgatasok.size()==0) ||
	    (peremek+parperemek.size()==0 && 
	     forgatasok+parforgatasok.size()<=2) ) &&
	  parperemek.size()+parforgatasok.size()>0)
      {
	list<param> parek;	//parameterek
	parek.splice(parek.end(),parperemek);
	parek.splice(parek.end(),parforgatasok);

	//Kell meg, hogy a parameterek reciprok osszegeben 1/p es 1/q legyen,
	//azaz az aktiv komponensben ugyanannyi szimplexhez tartozzon a
	//parameter, mint az egyutthatoja siktukrozes esetben, illetve ketszer
	//annyihoz forgatas esetben.
	int jo_orb=0;
	for(list<param>::iterator pit=parek.begin();pit!=parek.end();pit++){
	  int count=0;
	  for(list<simplex*>::iterator szit=pit->szek.begin();
	      szit!=pit->szek.end();szit++)
	    if(find(komp->szek.begin(),komp->szek.end(),*szit)!=
		komp->szek.end())
	      count++;
	  if(pit->eh<0 && count>abs(pit->eh))
	    jo_orb=1;
	  if(pit->eh>0 && count>2*pit->eh)
	    jo_orb=1;
	}


	if(jo_orb==0)
	  if(parek.size()==1){
	    if(peremek==1 || forgatasok==1)
	      *badorb<<parek.begin()->kar<<"=2";
	    else
	      return 0;
	  }
	  else
	    if(peremek!=1 && forgatasok!=1){
	      list<param>::iterator egyik,masik;
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
  return ret;
}

//Dsym::print_param_mx: Kiiratjuk a parameteres matrix-fuggvenyt (html-kodkent.)
void Dsym::print_param_mx(ostream *out){
  *out<<"<tr>";
  for(int j=0;j<car;j++){
    *out<<"<td><table cellpadding=\"3\">";
    for(int i=0;i<dim+1;i++){
      *out<<"<tr>";
      for(int k=0;k<dim+1;k++){
	*out<<"<td align=\"center\">";
	if (abs(i-k)==1) {
	  list<param>::iterator cp=params[j][(i<k) ? i : k];
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
  *out<<"</tr>"<<endl;
}

void Dsym::print_possible_splittings(ostream *out){
  /*
     list of points=...
     part1 points=first point
     notpart1 points={}
     list of edges going out=edges of first point
     recursion(out,list of points,part1 points,notpart1 points,list of edges)
     */
  list<kisebbdim*> pointlist;
  for(int op=0;op<dim+1;op++)
    for(list<kisebbdim>::iterator komp=klist[op].begin();komp!=klist[op].end(); komp++){
      pointlist.push_back(&(*komp));
    }

  list<kisebbdim*> part1list;
  part1list.push_back(*(pointlist.begin()));

  list<kisebbdim*> notpart1list;
  
  list<pair<kisebbdim*,kisebbdim*> > outbound_edges;
  for(int op=0;op<dim+1;op++)
    if (op!=(*part1list.begin())->op)
      for(list<kisebbdim>::iterator komp=klist[op].begin(); komp!=klist[op].end(); komp++){
	bool vankozos=false;
	// Gyorsitas: try catch
	for(list<simplex*>::iterator szimp1=komp->szek.begin(); szimp1!=komp->szek.end(); szimp1++)
	  for(list<simplex*>::iterator szimp2=(*part1list.begin())->szek.begin(); szimp2!=(*part1list.begin())->szek.end(); szimp2++)
	    if (*szimp1 == *szimp2)
	      vankozos=true;
	if (vankozos)
	  outbound_edges.push_back(pair<kisebbdim*,kisebbdim*>(*part1list.begin(),&(*komp)));
      }

  print_possible_splittings_recursion(out, pointlist, part1list, notpart1list, outbound_edges);
}

void Dsym::print_possible_splittings_recursion(ostream *out, list<kisebbdim*> pointlist, list<kisebbdim*> part1list, list<kisebbdim*> notpart1list, list<pair<kisebbdim*,kisebbdim*> > outbound_edges){
  /* 
     make copy of list of edges going out
     make copy of notpart1 points
     while |list of edges going out| > 0;do
     make copy of part1 points
     current edge=pop(list of edges going out)
     unless  current edge's other end in notpart1 points:
     part1 points+=current edge's other end
     make copy of list of edges going out
     list of edges going out+=new points edges going out
     list of edges going out-=new points edges going in
     print_splitting(out, list of edges going out)
     recurse(list of points,part1 points,notpart1 points,list of edges going
     out,out)
     notpart1 points+=current edge's other end
  //remove every instance of current edge's other end from list of edges
  //going out(double copy solves this problem)
  done

  outbound_edges' elements' first half is in part1; second part is not in
  part1.
   */

  list<kisebbdim*> notpart1list1 = list<kisebbdim*>(notpart1list);
  for (list<pair<kisebbdim*,kisebbdim*> >::iterator edge_it=outbound_edges.begin(); edge_it!=outbound_edges.end(); edge_it++){
    list<kisebbdim*> part1list1 = list<kisebbdim*>(part1list);
    pair<kisebbdim*,kisebbdim*> current_edge = *edge_it;
    kisebbdim* new_point = current_edge.second;
    bool do_we_care_about_the_edge=true;
    for(list<kisebbdim*>::iterator it=notpart1list1.begin(); it!=notpart1list1.end(); it++)
      if(new_point == *it){
	do_we_care_about_the_edge=false;
	break;
      }
    if (do_we_care_about_the_edge){
      part1list1.push_back(new_point);
      list<pair<kisebbdim*,kisebbdim*> > outbound_edges1 = list<pair<kisebbdim*,kisebbdim*> >(outbound_edges);
      outbound_edges1.remove(current_edge);

      //list of edges going out+=new point's edges going out
      for(int op=0;op<dim+1;op++)
	if (op!=new_point->op)
	  for(list<kisebbdim>::iterator komp=klist[op].begin(); komp!=klist[op].end(); komp++){
	    bool vankozos=false;
	    // Gyorsitas: try catch
	    for(list<simplex*>::iterator szimp1=komp->szek.begin(); szimp1!=komp->szek.end(); szimp1++)
	      for(list<simplex*>::iterator szimp2=new_point->szek.begin(); szimp2!=new_point->szek.end(); szimp2++)
		if (*szimp1 == *szimp2)
		  vankozos=true;
	    if (vankozos and find(part1list.begin(),part1list.end(),&(*komp))==part1list.end()){
	      outbound_edges1.push_back(pair<kisebbdim*,kisebbdim*>(new_point,&(*komp)));
	    }
	  }

      //list of edges going out-=new point's edges going in
      for(list<pair<kisebbdim*,kisebbdim*> >::iterator edge=outbound_edges1.begin(); edge!=outbound_edges1.end();edge++)
	while (edge!=outbound_edges1.end() and edge->second == new_point)
	  edge=outbound_edges1.erase(edge);

      if (part1list1.size() > 1 and pointlist.size()-part1list1.size() > 1 )
	print_splitting(out, outbound_edges1, part1list1);
      print_possible_splittings_recursion(out, pointlist, part1list1, notpart1list1, outbound_edges1);
      notpart1list1.push_back(new_point);
    }
  }
}

void Dsym::print_splitting(ostream *out, list<pair<kisebbdim*,kisebbdim*> > outbound_edges, list<kisebbdim*> part1list){
  /*
     sum=0;
     set of simplexes={};
     simplex_sum=0;
     foreach edge in outbound_edges:
     foreach simplex having edge:
     sum+=1/m_ij
     simplex_sum+=1
     push(set of simplexes,simplex)
     Combinatorial curvature=simplex_sum-2*|set of simplexes|-sum
     If cc<0:
     Check good orbifold crit.-> S2
     elsif cc=0:
     E2
     else
     NOP (We don't care about H2)
   */
  float sum=0;
  list<simplex*> simpleces;
  float simplex_sum=0;
  list<param*> splitting_params;
  for (list<pair<kisebbdim*,kisebbdim*> >::iterator edge=outbound_edges.begin(); edge!=outbound_edges.end(); edge++){
    list<simplex*> simpleces_of_edge;
    //Egy elnek van elso es masodik csucsa, ezekbol meg tudom mondani a
    //vonatkozo operacio parokat. Majd meg kell talalni egy kozos szimplexet,
    //ebbol elindulva az operacioparokkal megtalaljuk a tobbit.
    int i,j;
    list<int> ints;
    for(int i=0;i<dim+1;i++)
      if (i!=edge->first->op and i!=edge->second->op)
	ints.push_back(i);
    if (dim == 3){
      list<int>::iterator it=ints.begin();
      i=*it;
      it++;
      j=*it;
    }
    else
      throw "We are not in 3 dimensions";

    for(list<simplex*>::iterator szimp1=edge->first->szek.begin(); szimp1!=edge->first->szek.end();szimp1++)
      for(list<simplex*>::iterator szimp2=edge->second->szek.begin(); szimp2!=edge->second->szek.end();szimp2++)
	if (*szimp1 == *szimp2)
	  simpleces_of_edge.push_back(*szimp1);

    if ( j-i == 1 ){
      param* param_of_edge=&(*params[(*simpleces_of_edge.begin())->sorszam[0]][i]);
      if ( find(splitting_params.begin(),splitting_params.end(),param_of_edge) == splitting_params.end()) {
	splitting_params.push_back(param_of_edge);
      }
    }

    for (list<simplex*>::iterator szimp=simpleces_of_edge.begin(); szimp!=simpleces_of_edge.end(); szimp++){
      //FIXME mx[i][j] helyett a parameter minimalis erteket kene venni...
      //Es tarolni, hogy egy-egy parameterrel hanyszor foglalkoztunk
      sum+=1.0/(*szimp)->mx[i][j];
      simplex_sum+=1;
      if (find(simpleces.begin(),simpleces.end(),*szimp) == simpleces.end())
	simpleces.push_back(*szimp);
    }
  }
  float cc=simplex_sum-2*simpleces.size()-sum;
  if ( cc <= THRESH ){
  //<tr><td>Vertices of one part</td><td>Type of splitting</td><td>Essential
  //parameters</td></tr>
  *out << "<tr>";
  *out << "<td>";
  /*for(int op=0;op<dim+1;op++){
    if ( op != 0 )
      *out << "; ";
    bool colon=false;
    for ( list<kisebbdim*>::iterator currpoint=part1list.begin(); currpoint!=part1list.end(); currpoint++ ){
      if ((*currpoint)->op == op){
	if (colon)
	  *out << ", ";
	else
	  colon=true;
	int min_sorsz=car+1;
	for ( list<simplex*>::iterator szimp=(*currpoint)->szek.begin(); szimp!=(*currpoint)->szek.end(); szimp++ )
	  if ( (*szimp)->sorszam[0] < min_sorsz )
	    min_sorsz=(*szimp)->sorszam[0];
	*out << min_sorsz+1;
      }
    }
  }*/
  for ( list<kisebbdim*>::iterator currpoint=part1list.begin(); currpoint!=part1list.end(); currpoint++ ){
    *out << (*currpoint)->op << ": ";
    for ( list<simplex*>::iterator szimp=(*currpoint)->szek.begin(); szimp!=(*currpoint)->szek.end(); szimp++ )
      *out << (*szimp)->sorszam[0]+1 << " ";
    *out << "<br>";
  }
  *out << "</td>";

  *out << "<td>";
  if ( cc < -THRESH )
    *out << "S2" <<endl;
  else if (cc <= THRESH)
    *out << "E2" <<endl;
  else 
    *out << "H2" <<endl;
  *out << simplex_sum <<"-2*"<< simpleces.size() <<"-"<< sum <<"="<< cc <<"</td>";

  *out << "<td>";
  for (list<param*>::iterator param=splitting_params.begin(); param!=splitting_params.end(); param++){
    if (param!=splitting_params.begin())
      *out << ", ";
    //FIXME Hanyszorosan szamit egy-egy parameter?
    *out << (*param)->kar;
  }
  *out << "</td>";

  //debug
  *out << "<td>";
  for (list<pair<kisebbdim*,kisebbdim*> >::iterator edge=outbound_edges.begin(); edge!=outbound_edges.end(); edge++){
    *out << edge->first->op << ": ";
    for ( list<simplex*>::iterator szimp=edge->first->szek.begin(); szimp!=edge->first->szek.end(); szimp++ )
      *out << (*szimp)->sorszam[0]+1 << " ";
    *out << " &rarr; ";

    *out << edge->second->op << ": ";
    for ( list<simplex*>::iterator szimp=edge->second->szek.begin(); szimp!=edge->second->szek.end(); szimp++ )
      *out << (*szimp)->sorszam[0]+1 << " ";
    *out << "<br>";
  }
  *out << "</td>";

  *out << "</tr>"<<endl;
  }
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
int bt;
int bt1;
void backtrack(Dsym* D,Dsymlista* saved,int szin,int honnan,int hova,int current_largest) {
  int car=D->car;
  int dim=D->dim;
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
      //cout << " " << bt1 << " " << start << endl;
      Dsym* ujD=D->save_with_start(start);
      if (saved->check(ujD,1)==0) saved->append(ujD);
    }
    return;
  }

  // Nem a legkisebb 0-ad foku csucsot nem akarjuk hozzaadni
  bool ures=true;
  for(int j=0;j<dim+1;j++)
    if (D->csucsok[0][hova]->szomszed[j] != D->csucsok[0][hova])
      ures=false;

  bool erdemes=true;
  if ( ures )
    for (int i=1;i<hova;i++){ // A 0-adik csucs uressege nem szamit
      bool ures1=true;
      for(int j=0;j<dim+1;j++)
	if (D->csucsok[0][i]->szomszed[j] != D->csucsok[0][i]){
	  ures1=false;
	  break;
	}
      if (ures1){
	erdemes=false;
	break;
      }
    }

  //ha erdemes elt hozzaadni, ujra meghivjuk onmagunkat
  if (erdemes && D->elhozzaad(szin,honnan,hova)){
    int largest=current_largest;
    if (hova > largest)
      largest=hova;
    backtrack(D,saved,szin,honnan,hova,largest);

    D->eltorol(szin,honnan,hova);
  }

  if(hova+1<car) backtrack(D,saved,szin,honnan,hova+1,current_largest);
  else if(szin+1<dim+1) backtrack(D,saved,szin+1,honnan,honnan+1,current_largest);
  else {
    int szukseges_fok=0;
    int fok=0;
    //Heurisztika: az elso csucs fokszamanal nem veszunk nagyobb fokszamu esetet
    //(lexikografikusan: szinek szerint)

    for(int j=0;j<dim+1;j++) 
      if(D->csucsok[0][honnan]->szomszed[j] != D->csucsok[0][honnan])
	fok+=1 << j;
    if (honnan > 0){
      for(int j=0;j<dim+1;j++)
	if (D->csucsok[0][0]->szomszed[j] != D->csucsok[0][0])
	  szukseges_fok+=1 << j;
      if(fok <= szukseges_fok &&
	  fok >= 1 &&
	  current_largest >= honnan+1 &&
	  not backtrack_breaks_uvw(D,honnan))
	backtrack(D,saved,0,honnan+1,honnan+2,current_largest);
    }
    else if(fok>=1)  //current largest ekkor egyertelmu
      backtrack(D,saved,0,honnan+1,honnan+2,current_largest);
  }
}

bool backtrack_breaks_uvw(Dsym* D,int honnan){
  int dim=D->dim;
  int max=honnan+1;
  simplex*** csucsok=D->csucsok;
  for (int r=0;r<=max;r++)
    for (int i=0;i<dim-1;i++)
      for (int i1=i+2;i1<dim+1;i1++){
	/* A kovetkezo minta felismerese: |_
	   Ha talalunk nem szomszedos operacio parokkal 2 hosszu lancot, ahol mindket
	   tag nem nagyobb, mint max; az rossz. */
	if(csucsok[0][r]->szomszed[i] == csucsok[0][r] &&
	    csucsok[0][r]->szomszed[i1] != csucsok[0][r] &&
	    csucsok[0][r]->szomszed[i1]->sorszam[0] <= max &&
	    csucsok[0][r]->szomszed[i1]->szomszed[i] != csucsok[0][r]->szomszed[i1] &&
	    csucsok[0][r]->szomszed[i1]->szomszed[i]->sorszam[0] <= max &&
	    csucsok[0][r]->szomszed[i1]->szomszed[i]->szomszed[i1] == csucsok[0][r]->szomszed[i1]->szomszed[i]){
	  //cout << "Nem uvw1" <<endl;
	  //D->print(0);
	  return true;
	}
	/* A kovetkezo minta felismerese: |_|
	   Ha talalunk nem szomszedos operacio parokkal 3 hosszu valamit, ahol mindket
	   tag nem nagyobb, mint max, es nem zarodik vissza; az is rossz. */
	if(csucsok[0][r]->szomszed[i] == csucsok[0][r] &&
	    csucsok[0][r]->szomszed[i1] != csucsok[0][r] &&
	    csucsok[0][r]->szomszed[i1]->sorszam[0] <= max &&
	    csucsok[0][r]->szomszed[i1]->szomszed[i] != csucsok[0][r]->szomszed[i1] &&
	    csucsok[0][r]->szomszed[i1]->szomszed[i]->sorszam[0] <= max &&
	    csucsok[0][r]->szomszed[i1]->szomszed[i]->szomszed[i1] != csucsok[0][r]->szomszed[i1]->szomszed[i] &&
	    csucsok[0][r]->szomszed[i1]->szomszed[i]->szomszed[i1]->sorszam[0] <= max &&
	    csucsok[0][r]->szomszed[i1]->szomszed[i]->szomszed[i1]->szomszed[i] != csucsok[0][r]){
	  //cout << "Nem uvw2" <<endl;
	  //D->print(0);
	  return true;
	}
      }
  return false;
}

//Dsymlista::check: Hatulrol indulva megnezzuk, hogy hanyadik helyre kene tenni
//az uj elemet (a megfelelo csucsot kijelolve elsonek.) Azert hatulrol, mert
//olyan felsorolo algoritmust szeretnenk, ami sorba rakva sorolja fel oket. it
//kulso valtozot feltoltjuk.
int Dsymlista::check(Dsym* uj_elem,int var_uj) {
  it=last;
  while(it!=NULL) {
    int kis=kisebb(uj_elem,var_uj,it->curr,1);
    if (kis==0) {
      return it->ssz;
    }
    else if (kis==1) it=it->prev;
    else break;		//kis==-1
  }
  return 0;
}

//Dsymlista::append: az uj elemet a linkelt lista it-vel mutatott eleme moge
//szurja be. Fontos, hogy check(uj_elem,1)-et meg kell hivni elotte.
void Dsymlista::append(Dsym* uj_elem) {
  cerr<<"\r"<<count++;
  Dsymlinklist* uj=new Dsymlinklist;
  uj->curr=uj_elem;
  uj->ssz=1;
  uj->prev=it;
  if(uj->prev==NULL){
    uj->next=first;
    first=uj;
  }
  else {
    uj->next=it->next;			//it.....uj.....it->next
    it->next=uj;
  }
  if(uj->next==NULL) {
    last=uj;
  }
  else {
    uj->next->prev=uj;
  }
}

//Dsymlista konstruktor: ures lista alapertelemezesei.
Dsymlista::Dsymlista(void):first(0),last(0),count(0) {}

//Dsymlista destruktor: 
Dsymlista::~Dsymlista(void){
  for (Dsymlinklist* it=first;it!=NULL;it=it->next) {
    delete it->curr;
    delete it->prev;
  }
  delete last;
}

//Dsymlista::print: kiirja a konzolra a multigrafot (involucio alakban) egy
//sorszammal; kiszamolja, hogy a listaban hanyadik elem a dulisa, majd kiirja;
//operacio-paronkent kiirja egy sorba az oda tartozo parametereket (+-szal
//jelolve az iranyitott korbejarasokat;) illetve egy fajlba irja a multigraf
//rajzat (fig formatumban, amit kesobb jpg-ge alakithatunk.)
void Dsymlista::print(void){
  for (Dsymlinklist* it=first;it!=NULL;it=it->next) {
    int dim=it->curr->dim;
    int car=it->curr->car;
    cout<<endl<<it->ssz;
    it->curr->print(1);
    it->curr->dual=0;
    it->curr->dualis();

    cout<<"Dual: D.";
    for (int i=1;i<car+1;i++) {
      it->curr->dual->atsorszamoz(i);
      int c=check(it->curr->dual,i);
      if (c!=0) {
	cout<<c;
	break;
      }
    }
    cout<<endl;

    cout<<"Parameters: "<<endl;
    for(int i=0;i<dim;i++) {
      cout << "(" <<i<<","<<i+1<<") ";
      it->curr->print_params(i,&cout);
      cout<<endl;
    }
    ostringstream xfigfile;
    xfigfile<<"d"<<dim<<"c"<<car<<"_"<<it->ssz<<".fig";
    it->curr->write_xfig(xfigfile.str());
  }
  cout<<endl;
  cout << count<<endl;
}

//Dsymlista::print_html: Hasonloan mint az elobb, kiirjuk html-be az adatokat,
//annyi pluszt teszunk hozza, hogy kulon fajlba (currD) kiirjuk az osszes
//lehetseges matrix-rendszert is.
void Dsymlista::print_html(void){
  if (first==NULL)
    return;
  int dim=first->curr->dim;
  int car=first->curr->car;
  ostringstream filename;
  filename<<"d"<<dim<<"c"<<car<<".html";
  ofstream html_file;
  html_file.open(filename.str().c_str());

  //Header
  html_file<<"<html>"<<endl<<"<body>"<<endl;

  for (Dsymlinklist* it=first;it!=NULL;it=it->next) {
    cerr<<"\r"<<it->ssz;
    it->curr->dual=0;
    it->curr->dualis();
    inf=1000000;
    it->curr->create_params();
    it->curr->create_kdim();
    it->curr->filter_bad_orbifolds();

    ostringstream currDname;
    currDname<<"d"<<dim<<"c"<<car<<"_"<<it->ssz<<".html";
    ofstream currD;
    currD.open(currDname.str().c_str());

    currD<<"<html>"<<endl<<"<body>"<<endl<<"<table border=\"2\">"<<endl;
    dim=it->curr->dim;
    car=it->curr->car;
    currD<<"<caption>"<<it->ssz<<"</caption>"<<endl;

    ostringstream xfigfile;
    xfigfile<<"d"<<dim<<"c"<<car<<"_"<<it->ssz<<".fig";
    it->curr->write_xfig(xfigfile.str());
    ostringstream figtojpg;
    figtojpg<<"fig2dev -L jpeg "<<"d"<<dim<<"c"<<car<<"_"<<it->ssz<<".fig "
      <<"d"<<dim<<"c"<<car<<"_"<<it->ssz<<".jpg";
    //system(figtojpg.str().c_str());
    //remove(xfigfile.str().c_str());
    currD<<"<tr><td><img src=\"d"<<dim<<"c"<<car<<"_"<<it->ssz<<".jpg\"/></td>"
      <<endl;

    currD<<"<td><table border=\"1\">"<<endl;
    it->curr->print_param_mx(&currD);
    currD<< "<tr><td colspan=\""<<car<<"\">Number of " << dim-1 << 
      " dimensional components: "<< "(" <<it->curr->osszefuggo(0);
    for (int i=1;i<dim+1;i++) currD<<","<<it->curr->osszefuggo(i);
    currD<< ")</td></tr>"<<endl;
    currD<<"<tr><td colspan=\""<<car<<"\">"<<"Dual: ";
    int dualchk;
    for (int i=1;i<car+1;i++) {
      it->curr->dual->atsorszamoz(i);
      dualchk=check(it->curr->dual,i);
      if (dualchk!=0) {
	break;
      }
    }
    if (dualchk==0)
      currD<<"Not found";
    else
      if (dualchk==it->ssz)
	currD<<"Selfdual";
      else
	currD<<"D."<<dualchk;
    currD<<"</td></tr>"<<endl;

    currD<<"<tr><td colspan=\""<<car<<"\">"<<"Parameters: <br>"<<endl;
    for(int i=0;i<dim;i++) {
      currD << "(" <<i<<","<<i+1<<") ";
      it->curr->print_params(i,&currD);
      if (i<dim-1) currD <<"<br>"<<endl;
      else currD<<endl<<"</td></tr>"<<endl;
    }
    currD<<"</table></tr></table><br><table border=\"1\" cellpadding=\"3\">"
      <<endl<<"<caption>Possible parameter values:</caption>"<<endl
      <<"<thead><tr>";
    for(list<Dsym::param>::iterator pit=it->curr->plist.begin();
	pit!=it->curr->plist.end();pit++)
      currD<<"<td align=center>"<<pit->kar<<"</td>";
    currD<<"<td>Maximal</td><td>Ideal vertex</td>"  //Plusz infok
      <<"<td>Good orbifold criteria</td>";
    currD<<"<td>Backup info</td>";
    currD<<"</tr></thead><tbody>"<<endl;
    it->curr->mxnum=it->curr->print_possible_params(it->curr->plist.begin(), &currD);
    currD<<endl<<"</tbody></table>";
    // Possible splittings
    currD<<"<br><table border=\"1\" cellpadding=\"3\">"
      <<endl<<"<caption>Possible splittings:</caption>"<<endl
      <<"<thead><tr>";
    currD<<endl<<"<td>Vertices of one part</td><td>Type of splitting</td><td>Essential parameters</td><td>Edges split</td></tr>";
    it->curr->print_possible_splittings(&currD);
    currD<<endl<<"</table>";
    currD<<endl<<"</body></html>"<<endl;

    //html file:
    html_file<<"<table border=\"2\">"<<endl;
    html_file<<"<caption>"<<it->ssz<<"</caption>"<<endl;
    html_file<<"<tr><td><img src=\"d"<<dim<<"c"<<car<<"_"<<it->ssz
      <<".jpg\"/></td>"<<endl;
    html_file<<"<td><table border=\"1\">"<<endl;
    html_file<<"<tr><td colspan=\""<<car<<"\"><a href=\""<<currDname.str()
      <<"\">Number of matrices: "<<it->curr->mxnum<<"</a></td></tr>"<<endl;
    it->curr->print_param_mx(&html_file);
    html_file << "<tr><td colspan=\""<<car<<"\">Number of " << dim-1 << 
      " dimensional components: "<< "(" <<it->curr->osszefuggo(0);
    for (int i=1;i<dim+1;i++) html_file <<","<<it->curr->osszefuggo(i);
    html_file << ")</td></tr>"<<endl;
    html_file<<"<tr><td colspan=\""<<car<<"\">"<<"Dual: ";
    if (dualchk==0)
      html_file<<"Not "<<dim<<"-connected";
    else
      if (dualchk==it->ssz)
	html_file<<"Selfdual";
      else
	html_file<<"D."<<dualchk;
    html_file<<"</td></tr>"<<endl;
    html_file<<"<tr><td colspan=\""<<car<<"\">"<<"Parameters: <br>"<<endl;
    for(int i=0;i<dim;i++) {
      html_file << "(" <<i<<","<<i+1<<") ";
      it->curr->print_params(i,&html_file);
      if (i<dim-1) html_file <<"<br>"<<endl;
      else html_file<<endl<<"</td></tr>"<<endl;
    }
    html_file<<endl<<"</table></tr></table><br>"<<endl;
  }

  //Footer
  html_file<<"</body>"<<endl<<"</html>"<<endl;

}

//xfig konstruktor: feltoltjuk a valtozokat, es rajzolunk egy el nelkuli
//multigrafot. rad a korok sugara, koordx illetve koordy pedig a korok x es y
//koordinatainak listaja.
xfig::xfig(string filename,int dim1,int car1):dim(dim1),car(car1) {
  rad=200;
  koordx=new int[car];
  koordy=new int[car];
  for (int i=0;i<car;i++) {
    koordx[i]=int(round(4000+1500*cos(PI-i*2*PI/car)));
    koordy[i]=int(round(4000-1500*sin(PI-i*2*PI/car)));
  }
  outfile.open(filename.c_str());
  outfile << "#FIG 3.2  Produced by xfig version 3.2.5-alpha5" <<endl;
  outfile << "Landscape" <<endl;
  outfile << "Center" <<endl;
  outfile << "Metric" <<endl;
  outfile << "A4" <<endl;
  outfile << "100.00" <<endl;
  outfile << "Single" <<endl;
  outfile << "-2" <<endl;
  outfile << "1200 2" <<endl;
  for(int i=0;i<car;i++){
    create_circle(i);
    create_numtext(i);
  }
}

//xfig::create_circle: az n-edik csucs korenek megrajzolasa.
void xfig::create_circle(int n) {
  outfile << "1 3 0 1 0 7 49 -1 20 0.000 1 0.0000 "<<koordx[n]<<" "<<koordy[n]
    <<" "<<rad<<" "<<rad<<" "<<koordx[n]<<" "<<koordy[n]<<" "
    <<koordx[n]+rad<<" "<<koordy[n]<<endl;
}

//xfig::create_numtext: az n-edik sorszam beleirasa.
void xfig::create_numtext(int n) {
  outfile << "4 0 0 48 -1 0 12 0.0000 4 150 105 "<<koordx[n]-52<<" "
    <<koordy[n]+75<< " "<<n+1<<"\\001"<<endl;
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
    <<" "<<koordx[n1]+diffx<<" "<<koordy[n1]+diffy<<endl;
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
  char fels;
  cout << "Dimension: "; cin >> dim;cout<<dim<<endl;
  cout << "Cardinality: ";cin >> car;cout<<car<<endl;
  //cout << "Felsorolás? ";cin >> fels;cout<<fels<<endl;
  Dsym* D=new Dsym(dim,car);
  //if (fels=='y') {
  Dsymlista* saved=new Dsymlista;
  bt=0;
  bt1=0;
  backtrack(D,saved,0,0,1,0);
  int ujssz=1;
  for (Dsymlinklist* it=saved->first;it!=NULL;it=it->next)
    it->ssz=ujssz++;
  cout<<"saved"<<endl;
  cout<<bt<<endl;
  cout<<saved->count<<endl;
  cout<<bt1<<endl;
  cout<<float(bt1)/float(saved->count)<<endl;
  saved->print_html();
  delete D;
  delete saved;
  return 0;
  //}
  /*else{
    int szin,honnan,hova;
    char utolso;
    do{
    cout << "Kérem az él paramétereit (szin honnan hova) ";
    cin >>szin>>honnan>>hova;
    cout << szin<<" "<<honnan<<" "<<hova<<endl;
    cout << "Utolsó? ";cin >>utolso;cout <<utolso<<endl;
    D->elhozzaad(szin,honnan,hova);
    }while (utolso != 'y');

    if (D->osszefuggo(-1) == 1){
    D->atsorszamoz(1);
    D->print(1);
    ostringstream xfigfile;
    xfigfile<<"d"<<dim<<"c"<<car<<".fig";
    D->write_xfig(xfigfile.str());
    }
    else D->print(0);
    if (D->uvw()) {
    cout<<"Paraméterek: "<<endl;
    for(int i=0;i<dim;i++) {
    cout << "(" <<i<<","<<i+1<<") ";
    D->print_params(i,&cout);
    cout<<endl;
    }
    }
    else cout<<"Nem megfelelõ D-szimbólum."<<endl;
    return 0;
    }*/
}
