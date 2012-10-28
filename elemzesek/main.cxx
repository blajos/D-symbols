#include "test.hxx"

int main(){
  Test tester;
  try{
    return tester.all();
  }
  catch (int a){
    return a;
  }
}
