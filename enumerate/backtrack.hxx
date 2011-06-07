#include <iostream>
#include <fstream>
#include <list>
#include <algorithm>
#include <set>

#define THRESH 0.0000001
#define DEGENERATION_LIMIT 3

using namespace std;

//simplex osztaly: A D-szimbolumok szimplexeinek avagy csucsainak tarolasara.
//Fontos informacio a dimenzio (ennel eggyel tobb szomszedsagi operacio van) 
//es az elemszam (cardinality, amitol a D-szimbolum lehetseges sorszamozasai es
//igy a csucs lehetseges sorszamai fuggenek.) Operacionkent egy, szomszedra
//mutato pointert tarolunk. Ezenkivul meg a csucshoz tartozo M-matrix is itt
//van. A konstruktor a dimenziot, az elemszamot es az atsorszamozatlan eset
//sorszamat kapja meg.
//
//params: az adott szimplexre vonatkozo paramterek.
class simplex {
  public:
    simplex** szomszed;
    int* sorszam;
    int** mx;
    int dim,car;
    simplex(int,int,int);
    ~simplex(void);
};

//Dsym osztaly: 
// adatok: dimenzio es elemszam, szimplexre mutato pointerek matrixa (osszesen
//  elemszam+1 darab kulonbozo, minden rogzitett sorszamozasra szerepel 
//  mindegyik, illetve a kivulrol feltoltott sorszamozast, ami gyakran nem felel
//  meg szuksegleteinknek,) es a dualisra mutato pointer.
// Alapvetoen valtoztato fuggvenyek:
//  konstruktor: dimenzio es elemszam megadasa utan ertelmesen feltolti a
//   szimbolum ertekeit
//  destruktor: memoria torles
//  save: keszit egy masolatot, es visszaadja a ra mutato pointert
//  elhozzaad: ha nem sert meg senkit, hozzadja a (szin, honnan, hova) elt es
//   1-at ad vissza, kulonben 0-t.
//  eltorol: ha van ilyen el (szin, honnan, hova) akkor torli, azaz onmagukra
//   iranyitja a szomszedsagokat
//  atsorszamoz: kitolti a megfelelo elemmel kezdodo sorszamozasokat
//   (simplexekben is.) 0-t nem szerencses beirni, mert az a nem valtozo resz
//   (kivulrol lett feltoltve, sorszamozas fuggetlenul.)
//
// Alapvetoen ellenorzesek:
//  osszefuggo: az argumentumkent kapott operaciot (nem feltetlenul letezot...)
//   elhagyva hany reszre esik szet a szimbolum.
//  dualis: letrehozza a dualist, es megnezi, hogy melyikuk kisebb (1 ha a
//   dualis, 0, ha egyenloek, -1 ha az eredeti kisebb)
//  sorszamozas: Ha van olyan atsorszamozas, ami kisebb, akkor 1, kulonben ha
//   van szimmetria, akkor 0, kulonben 1.
//  uvw: Megnezi minden nem szomszedos operacio-parra, hogy ketszer alkalmazva a
//   part, visszajutunk-e az eredeti csucsba. Kozben a matrixok nem valtozo
//   ertekeit feltolti.
//  ellenoriz: Az ellenorzeseket fogja ossze, amikben meg nem szerepel a
//   matrixok vizsgalata. 1 ha biztosan nincs ra szukseg (atsorszamozva
//   kisebb, nem osszefuggo, vagy a nem szomszedos operaciok elrontjak.)
// 
// Kiiratasok:
//  write_xfig: A string-kent megadott nevu fajlba irja multigrafot kirajzolo
//   fig-kodot.
//
//  involucio: A grafot involuciokent egy 3 melysegu szam listakent
//   legegyszerubb leirni. (Ez a 3 melysegu szamlista.)
//  invol_create: Az adott sorszamhoz kesziti el a sorszamot, mint 1-est
//   valasztva az involuciot.
//  print: A megfelelo sorszamozas szerint kiirja az involuciot a kepernyore.
//
//
// Parameterek es matrix-fuggvenyekek kezelese:
//  param tipus:
//   Egy-egy parametert ir le, az egyutthatojat, a lehetseges minimalis erteket,
//   a pillanatnyi erteket, melyik
//   operacio-parra vonatkozik (op,op+1) es hogy melyik szimplexeket
//   befolyasolja. Kesobb egy egyenloseget is definialunk a karakter alapjan.
//  plist: Az adott multigraf (a matrix-fuggveny fugg a parameterektol, es nem
//   forditva) parametereinek listaja.
//  kisebbdim tipus:
//   Olyan szimplexek listajat tartalmazza, amik egy bizonyos operacio elhagyasa
//   utan egy komponensbe esnek. Egy ilyen komponenst ir le.
//  klist: dim+1 elemu tomb, aminek elemei kisebbdim elemu listak. A i. elemeben
//   az i. operacio elhagyasaval keletkezett  komponensek listaja van.
//  mxnum: A kiszamolt matrix-fuggvenyek szama.
//
// Valtoztatasok: 
//  create_kdim: klist letrehozasa
//  create_params: plist letrehozasa, es az ertekek beallitasa egy minimalis
//   ertekre, hogy eh*ertek>=3
//  filter_bad_orbifolds: plist megfelelo parametereinek egyesitese, a rossz
//   orbifoldok szuresere
//  change_param: Az iterator altal mutatott parameter ertekenek novelese az
//   egeszkent megadott argumentummal (ami lehet negativ is, ekkor csokkentes.)
//   Es a matrixrendszer aktualizalasa a valtoztatassal.
//  increase_param: change_param(it,+1)
//  decrease_param: change_param(it,-1), mindig csak noveles utan csokkentunk
//
// Ellenorzesek:
//  kdimsf: Igaz-e, hogy minden vizsgalt kisebb dimenzios resszimbolum szferikus
//   vagy Euklideszi terben valosul meg (az argumentumkent megadott parametert
//   vegtelennek valasztva.) (Ezt csak a 2 dimenziosakra tudjuk
//   eldonteni.) Minden szferikus:1, van Euklideszi (es nincs hip.): 0, van 
//   hiperbolikus: -1.
//  min: A matrix-fuggvenyt is figyelembe veve, van-e szimmetriaja, illetve
//   van-e kisebb szamozasa a szimbolumnak. Ha van, akkor nem minimalis (0-t ad
//   vissza,) kulonben 1-et.
//  filter_bad_orbifolds03: A 0. es a 3. operacio elhagyaskor esetleg keletkezo
//   rossz orbifold esetek kezelese.
//
// Kiiratasok:
//  print_possible_mxes: backtrack algoritmus segitsegevel vegigmegyunk a
//   parameterek lehetseges ertekein es kiiratjuk a matrixokat (ha minimalis a
//   szimbolum,) a parametereket egyesevel novelgetve. (Az esetleges vegtelen
//   ciklusra figyelni kell.) Minden lehetseges parameter-allast max egyszer
//   szeretnenk latni. (Nem csak kiiratunk, hanem valtoztatunk is, de a
//   kiiratason van a hangsuly.)
//  print_param_mx: A matrix-fuggvenyt kiirjuk parameteresen
//  print_possible_params: backtrack algoritmus segitsegevel vegigmegyunk a
//   parameterek lehetseges ertekein es kiiratjuk oket (ha nem minimalis a
//   szimbolum jelezzuk,) a parametereket egyesevel novelgetve. (Az esetleges 
//   vegtelen
//   ciklusra figyelni kell.) Minden lehetseges parameter-allast max egyszer
//   szeretnenk latni. (Nem csak kiiratunk, hanem valtoztatunk is, de a
//   kiiratason van a hangsuly.)
//  print_possible_infs: a vegtelen lancokat irja ki a possible params altal
//   bealitott parameterek alapjan.
//  print_params: (int,int+1) operacio-parhoz tartozo parameterek es
//   egyutthatoik kiirasa ("+"-t irunk moge, ha iranyitott.) A karakter ahhoz
//   kell, hogy atlathatobb legyen a kimenet.
class Dsym {
  public:
    int dim,car;
    simplex*** csucsok;
    Dsym*  dual;

    Dsym(int,int);
    ~Dsym(void);
    Dsym* save(void);
    Dsym* save_with_start(int);
    int elhozzaad(int,int,int);
    void eltorol(int,int,int);
    void atsorszamoz(int);

    int osszefuggo(int);
    int dualis(void);
    int sorszamozas(void);
    int uvw(void);
    int ellenoriz(void);

    void write_xfig(string);

    int*** involucio;
    void invol_create(int);
    void print(int);

    
    struct param {
      int eh,min_ertek,ertek,op;
      char kar;
      list<simplex*> szek;
    };
    list<param> plist;
    struct kisebbdim {
      list<simplex*> szek;
      bool iranyitott;
    };
    list<kisebbdim> *klist;
    int mxnum;
    list<param>::iterator** params;

    void create_kdim(void);
    void create_params(void);
    void filter_bad_orbifolds(void);
    void change_param(list<param>::iterator,int);
    void increase_param(list<param>::iterator);
    void decrease_param(list<param>::iterator);

    int kdimsf(list<param>::iterator);
    int min(void);
    list<int> osszevlist;

    void print_param_mx(ostream*);
    int print_possible_params(list<param>::iterator,ostream*);
    int pmaxertek;  //a maximalis veges parametererteknel eggyel nagyobb
    int print_possible_infs(int,list<param>::iterator,ostream*);
    int filter_bad_orbifolds03(ostringstream*); //0 es 3 esetben rossz orb
    int idcs;	//segedvaltozo idealis csucsok detektalasahoz
    void print_possible_mxes(list<param>::iterator,ostream*);
    void print_params(int,ostream*);
};

bool operator == (Dsym::param,Dsym::param);

//inf: "vegtelen" nagy szam, amit folyamatosan valtoztatunk, nehogy szimmetria
//legyen.
int inf;

//Dsymlinklist: oda-vissza lepkedos linkelt lista, D-szimbolumokra mutato
//pointerekkel, mint adattal.
struct Dsymlinklist {
  struct Dsymlinklist *next,*prev;
  int ssz;
  Dsym *curr;
};

//Dsymlista osztaly:
// first,last,it: elso, utolso, iterator a szamunkra hasznos D-szimbolumok
//  linkelt listajan. A listat rendezetten tartjuk (ezert lenne jo a set...)
// count: a lista hossza (csak hogy tudjuk...)
//
// konstruktor: alapertelmezesek: mindent 0-ra allitunk
// append: A beszurasi pontba szurja az uj elemet (azaz a megfelelo pointereket
//  atallitgatja.) Fontos, hogy elotte check(uj_elem,1)==0 kell legyen.
//
// check: szerepel-e az adott szimbolum az adott kezdoelemmel a listaban: Ha
//  igen hanyadik (a dualis megtalalasahoz kell,) ha nem, akkor 0. 
//  (Beallitja az it valtozot a beszurasi pontra.)
//
// print_html: letrehoz egy html-fajlt es bele irja a szimbolumokat minden
//  kideritett informacioval egyutt (a matrixokat kulon fajlban reszletezve,) es
//  a multigrafokat leiro kep-fajlokat bele linkeli.
// print: A kepernyore irja az involuciokat es az informaciok egy reszet.
class Dsymlista {
  public:
    Dsymlinklist *first,*last,*it;	//csinaljuk meg avl-lel
    int count;

    Dsymlista(void);
    ~Dsymlista(void);
    void append(Dsym* uj_elem);

    int check(Dsym*,int);

    void print_html(void);
    void print(void);		//mindet kiirja
};

//backtrack: felsorolja az osszes lehetseges (szamunkra esetleg erdekes)
// multigrafot (barmely elszint veve igaz, hogy max 2 elemu komponensei vannak
// az igy kapott grafnak.) Lefuttatja ra az ellenorzeseket, es ha valoban
// erdekes, akkor hozzad egy masolatot a listahoz.
// argumentumok: Dsym*: benne dolgozunk (hozzaadunk es torlunk eleket)
//               Dsymlista*: az erdekesek listaja
//               int,int,int: vizsgalt el: (szin,honnan,hova,current_largest)
void backtrack(Dsym*,Dsymlista*,int,int,int,int);

//Lehet-e meg jo a diagram?
bool backtrack_breaks_uvw(Dsym*,int);

//Rendezes matrixok nelkul: elso multigraf es kezdopontja, masodik multigraf es
// kezdopontja. 
// Elso kisebb:1	egyenlok:0	masodik kisebb:-1
//
// A megfelelo atsorszamozast elotte meg kell csinalni.
int kisebb(Dsym*,int,Dsym*,int);

//my_find: Atsorszamozas megkonnyitesere segedfv: mit keresunk a hol kezdopontu,
//hossz hosszu tombben
int my_find(simplex* mit,int hossz,simplex** hol);

//xfig osztaly:
//multigraf fig kodjanak fajlba irasara
//
// konstruktor: file megnyitasa, dimenzio, elemszam kitoltese, korok
//  kozeppontjainak kiszamolasa, sugaranak beallitasa, korok es szamok rajzolasa.
// create_circle: n=0..car-1 a megfelelo sorszamu szimplexhez tartozo kor
//  rajzolasa
// create_numtext: n=0..car-1 a megfelelo sorszamu sorszam beleirasa
// create_line: ket kor osszekotese szinnel (szinek szerint az elre meroleges 
//  iranyba eltolva, az atfedesek megakadalyozasa miatt.)
// destruktor: fajl bezarasa, memoria felszabaditasa.
class xfig {
  ofstream outfile;
  int *koordx,*koordy,rad,dim,car;
  public:
  xfig(string filename,int dim,int car);
  void create_circle(int n);		//az n-edik kor megrajzolasa
  void create_numtext(int n);		//az n-edik kor szama
  void create_line(int n0,int n1,int szin);	//osszekotjuk a ket kort
  ~xfig(void);
};

//TODO: 
//  - Dsym::print_html-be atteni a Dsymlist::print_html-beli, nem odatartozo
//    dolgokat.
//  - Dsym adatait nem kozvetlenul elerni, hanem private-ba tenni
//  - backtrack fuggvenyt a rendezo algoritmusra optimalizalni
//  - Dsymlinklist helyett a beepitett set osztalyt hasznalni (ez sem trivialis)
//
//Ujrastrukturalas:
//  - Angol valtozok, angol kommentekkel
//  - Ahol a sima pointer nem egyszerubb minden masnal, meg kene nezni, milyen
//  beepitett szerkezet hasznalhato
//  - 0-rol ujrairni, sallangok nelkul
//  - Vegtelen paramterek korrekt kezelese (nem holmi 100000-es hack-kel)
//  - szem elott tartani a fundamentalis tartomany es a splitting kereses
//  lehetoseget is
//  - Az osszegyujtott infokat nem irjuk ki egybol, hogy kesobbi ponton
//  feldolgozhatoak/finomithatoak legyenek
//  - Dokumentacio Doxygen hasznalataval (angolul...)
//  - Osztalystruktura:
//    * simplex (vagy graf csucs)
//    * D-graf
//    * D-symbolum
//    * Dsym* save helyett konstruktor kene
//    * parameter, parameterlista (set)
//    * D-symbolum lista (set)
//    * <, >, == operatorok definialasa (set-ekhez)
//    * html fajl
//    * osszegzo oldal
//    * reszletezo oldal
//    * svg objektum
//    * backtrack algo altalanosan? (nem hiszem, hogy kene) (osztalykent, amibol
//    orokitheto a param backtrack es a Dsym backtrack?)
