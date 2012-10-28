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
