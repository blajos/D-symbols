#include "dlist.hxx"

template <class D> Dlist::Dlist(string fn):
  filename_base(fn),
  fastdb(NULL,0),
  sorteddb(NULL,0),
  count(0),
  current(0),
  keylength(2048)
{
  sorteddb.set_bt_compare(compare_d);

  fastdb.open(NULL, filename_base + "_fast.db", NULL, DB_HASH, DB_CREATE, 0);
  sorteddb.open(NULL, filename_base + "_sort.db", NULL, DB_BTREE, DB_CREATE, 0);
}

template <class D> Dlist::~Dlist(){
  current->close();
  fastdb.close(0);
  sorteddb.close(0);
}

template <class D> int Dlist::compare_d(Db *dbp, const Dbt *a, const Dbt *b){
  D* dobj1,dobj2;
  char* dstr1=new char[a->get_size()];
  char* dstr2=new char[b->get_size()];

  memcpy(&dstr1, a->get_data(), a->get_size()); 
  memcpy(&dstr2, b->get_data(), b->get_size()); 

  sstream dsstr1,dsstr2;
  dsstr1 << dstr1;
  dsstr2 << dstr2;

  dobj1=new D(dsstr1);
  dobj2=new D(dsstr2);

  return dobj2->is_smaller(0,dobj1,0);
}

template <class D> int Dlist::count(void){
  return count;
}

template <class D> void Dlist::append(D* new_element){
  sstream dumpstr;
  new_element->dump(&dumpstr);
  bool a=true;

  Dbt key((void*)&dumpstr.str().c_str(), dumpstr.str().size() + 1);
  Dbt data(&a, sizeof(a));

  int ret = fastdb.put(NULL, &key, &data, DB_NOOVERWRITE);
  if (ret == 0) {
    count++;
    sorteddb.put(NULL, &key, &data, DB_NOOVERWRITE);
  }
  // Else: do nothing...
  //else if (ret == DB_KEYEXIST) {
  //}
}

template <class D> int Dlist::check(D* element){
  sstream dumpstr;
  new_element->dump(&dumpstr);
  bool a;

  Dbt key((void*)&dumpstr.str().c_str(), dumpstr.str().size() + 1);
  Dbt data(&a, sizeof(a));

  int ret=fastdb.get(NULL, &key, &data, 0);
  if (ret == DB_NOTFOUND)
    return 0;
  else
    return 1;
}

template <class D> D* Dlist::getnextsorted(void){
  if (current == 0)
    reset_cursor();

  sstream dumpstr;
  bool a;
  char* str;
  int ret;

  str=new char[keylength];
  
  Dbt key;
  Dbt data(&a, sizeof(a));

  key.set_data(str);
  key.set_ulen(keylength);
  key.set_flags(DB_DBT_USERMEM);

  try {
    ret = current->get(&key, &data, DB_NEXT);
  }
  catch (DbMemoryException& e){
    if (e.get_errno() == DB_BUFFER_SMALL) {
      keylength=2*key.get_size();
      delete[] str;

      str=new char[keylength];
      key.set_data(str);
      key.set_ulen(keylength);
      key.set_flags(DB_DBT_USERMEM);
      ret = current->get(&key, &data, DB_NEXT);
    }
    else {
      throw 1;
    }
  }
  if (ret == DB_NOTFOUND)
    return NULL;

  dumpstr << *str;
  D* output(dumpstr);

  delete[] str;

  return output;
}

template <class D> void Dlist::reset_cursor(void){
  sorteddb.cursor(0,&current,DB_CURSOR_BULK);
}
