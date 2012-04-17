#include<db_cxx.h>

int main(int argc, char *argv[]){
  Db db(0,0);
  db.open(NULL, "./dbfile.db", NULL, DB_BTREE, DB_CREATE, 0664);

  std::string *teszt = new std::string("Grocery bill.");
  int a=0;

  Dbt key((void*)teszt->c_str(), teszt->size() + 1);
  Dbt data(&a, sizeof(a));

  int ret = db.put(NULL, &key, &data, DB_NOOVERWRITE);
  if (ret == DB_KEYEXIST) {
        db.err(ret, "Put failed because key %s already exists", teszt->c_str());
  }
  return 0;
}
