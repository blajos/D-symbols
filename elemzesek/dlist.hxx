#ifndef __dlist_h
#define __dlist_h

#include<db_cxx.h>
#include<sstream>

using namespace std;

/*
   Class: Dlist

   The Dlist class represents a list of D-diagrams or D-symbols. We use Berkeley
   DB for storing the information, because it's fast and has an internal caching
   routine.

   <Dlist.fastdb> is a hash storing the strings dumped from D-objects.
   <Dlist.sorteddb> is a btree storing the same strings but utilizing our custom
   comparison functions.
   */
template <class D> class Dlist {
  private:
    Db fastdb,sorteddb;
    Dbc *current;
    string filename_base;
    int count;
    int keylength;
  public:
    int count(void); //Db::stat()
    Dlist(void);
    ~Dlist(void);
    int append(D* new_element);
    int check(D* element);
    D* getnextsorted(void);
    void reset_cursor(void);
    int compare_d(Db *, const Dbt *, const Dbt *);
};


int dbmain(int argc, char *argv[]){
  Db db(0,0);
  db.open(NULL, "./dbfile.db", NULL, DB_HASH, DB_CREATE, 0664);

  std::string *teszt = new std::string("Grocery bill.");
  bool a=true;

  Dbt key((void*)teszt->c_str(), teszt->size() + 1);
  Dbt data(&a, sizeof(a));

  int ret = db.put(NULL, &key, &data, DB_NOOVERWRITE);
  if (ret == DB_KEYEXIST) {
        db.err(ret, "Put failed because key %s already exists", teszt->c_str());
  }
  return 0;
}

