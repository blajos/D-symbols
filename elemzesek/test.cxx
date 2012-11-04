#include "test.hxx"

Test::Test(void){
  return;
}

Test::~Test(void){
  return;
}

int Test::all(void){
  return base()+simplex()+mxfunction()+param()+svg()+ddiag()+dsym()+dlist()+\
    algorithms();
}

int Test::base(void){
  cerr << "Testing Base module: ";
  Base tester(3,6);
  stringstream str1,str2;
  tester.dump(&str1);
  Base tester1(&str1);
  tester1.dump(&str2);
  if (str2.str() != "3 6"){
    cerr << str2.str() << endl;
    return 1;
  }
  cerr << "Ok" << endl;
  return 0;
}

int Test::simplex(void){
  cerr << "Testing Simplex module: ";
  Simplex tester(3,6,1);
  stringstream str;
  tester.dump(&str);
  if (str.str() != "3 6 1"){
    cerr << str.str() << endl;
    return 1;
  }
  cerr << "Ok" << endl;
  return 0;
}

int Test::mxfunction(void){
  cerr << "Testing Mxfunction module: ";
  Mxfunction tester(3,6);
  tester.set(2,2,3,15);
  tester.set(1,1,1,5);
  int a=tester.get(2,2,3);
  if (a != 15){
    cerr << " " << a;
    return 1;
  }
  Mxfunction tester1(&tester);
  stringstream str1,str2;
  tester1.dump(&str1);
  if (str1.str() != "3 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 15 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"){
    cerr << " " << str1.str() << endl;
    return 1;
  }
  tester1.print_html(&str2);
  if (str2.str() != "<table border=\"1\"><tr><td><table cellpadding=\"3\"><tr><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td></tr><tr><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td></tr><tr><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td></tr><tr><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td></tr></table></td><td><table cellpadding=\"3\"><tr><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td></tr><tr><td align=\"center\">0</td><td align=\"center\">5</td><td align=\"center\">0</td><td align=\"center\">0</td></tr><tr><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td></tr><tr><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td></tr></table></td><td><table cellpadding=\"3\"><tr><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td></tr><tr><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td></tr><tr><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">15</td></tr><tr><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td></tr></table></td><td><table cellpadding=\"3\"><tr><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td></tr><tr><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td></tr><tr><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td></tr><tr><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td></tr></table></td><td><table cellpadding=\"3\"><tr><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td></tr><tr><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td></tr><tr><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td></tr><tr><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td></tr></table></td><td><table cellpadding=\"3\"><tr><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td></tr><tr><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td></tr><tr><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td></tr><tr><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td><td align=\"center\">0</td></tr></table></td></tr></table>"){
    cerr << str2.str() << endl;
    return 1;
  }
  cerr << "Ok" << endl;
  return 0;
}

int Test::param(void){
  //FIXME Ezen nincs igazan mit tesztelni
  return 0;
}

int Test::svg(void){
  cerr << "Testing Svg module: ";
  stringstream str;
  Svg my_svg(4,6);
  my_svg.add_line(0,1,0);
  my_svg.add_line(1,3,3);
  my_svg.add_line(2,4,2);
  my_svg.add_line(3,4,1);
  my_svg.add_line(4,5,4);
  my_svg.add_line(0,5,2);
  my_svg.print_html(&str);
  if (str.str() != "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 300 300\" width=\"300px\" height=\"300px\" version=\"1.1\"><line style=\"stroke-dasharray:2,10\" x1=\"7\" y1=\"143\" x2=\"72\" y2=\"30\" stroke=\"black\" stroke-width=\"1\"/><line style=\"stroke-dasharray:8,8,2,10\" x1=\"82\" y1=\"43\" x2=\"277\" y2=\"156\" stroke=\"black\" stroke-width=\"1\"/><line style=\"\" x1=\"215\" y1=\"37\" x2=\"215\" y2=\"263\" stroke=\"black\" stroke-width=\"1\"/><line style=\"stroke-dasharray:8,8\" x1=\"286\" y1=\"153\" x2=\"221\" y2=\"266\" stroke=\"black\" stroke-width=\"1\"/><line style=\"stroke-dasharray:8,8,2,10,2,10\" x1=\"215\" y1=\"248\" x2=\"85\" y2=\"248\" stroke=\"black\" stroke-width=\"1\"/><line style=\"\" x1=\"20\" y1=\"150\" x2=\"85\" y2=\"263\" stroke=\"black\" stroke-width=\"1\"/><circle cx=\"20\" cy=\"150\" r=\"20\" fill=\"white\" stroke=\"black\" stroke-width=\"1\"/><text x=\"20\" y=\"159\" font-size=\"28\" text-anchor=\"middle\" dominant-baseline=\"mathematical\">1</text><circle cx=\"85\" cy=\"37\" r=\"20\" fill=\"white\" stroke=\"black\" stroke-width=\"1\"/><text x=\"85\" y=\"46\" font-size=\"28\" text-anchor=\"middle\" dominant-baseline=\"mathematical\">2</text><circle cx=\"215\" cy=\"37\" r=\"20\" fill=\"white\" stroke=\"black\" stroke-width=\"1\"/><text x=\"215\" y=\"46\" font-size=\"28\" text-anchor=\"middle\" dominant-baseline=\"mathematical\">3</text><circle cx=\"280\" cy=\"150\" r=\"20\" fill=\"white\" stroke=\"black\" stroke-width=\"1\"/><text x=\"280\" y=\"159\" font-size=\"28\" text-anchor=\"middle\" dominant-baseline=\"mathematical\">4</text><circle cx=\"215\" cy=\"263\" r=\"20\" fill=\"white\" stroke=\"black\" stroke-width=\"1\"/><text x=\"215\" y=\"272\" font-size=\"28\" text-anchor=\"middle\" dominant-baseline=\"mathematical\">5</text><circle cx=\"85\" cy=\"263\" r=\"20\" fill=\"white\" stroke=\"black\" stroke-width=\"1\"/><text x=\"85\" y=\"272\" font-size=\"28\" text-anchor=\"middle\" dominant-baseline=\"mathematical\">6</text></svg>"){
    cerr << str.str() << endl;
    return 1;
  }
  cerr << "Ok" << endl;
  return 0;
}

int Test::ddiag(void){
  stringstream str1,str2,str3;
  cerr << "Testing Ddiag module: ";
  Ddiag tester(3,6);
  tester.create_edge(0,0,1);
  tester.create_edge(1,1,2);
  tester.create_edge(0,2,3);
  tester.create_edge(1,3,4);
  tester.create_edge(0,4,5);

  tester.dump(&str1);
  if (str1.str() != "3 6 ((0,1)(2,3)(4,5))((0)(1,2)(3,4)(5))((0)(1)(2)(3)(4)(5))((0)(1)(2)(3)(4)(5))"){
    cerr << str1.str() << endl;
    return 1;
  }

  tester.print_html(&str2);
  if (str2.str() != ""){
    cerr << " " << str2.str() << endl;
    return 1;
  }
  return 0;
}

int Test::dsym(void){
  //FIXME
  return 0;
}

int Test::dlist(void){
  //FIXME
  return 0;
}

int Test::algorithms(void){
  //FIXME
  return 0;
}

