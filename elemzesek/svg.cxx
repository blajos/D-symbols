#include "svg.hxx"

Svg::Svg(int dimin,int carin):
  Base(dimin,carin)
{
  rad=20;
  fontsize=28;
  size=max(300, car*50);
  koordx=new int[car];
  koordy=new int[car];
  for (int i=0;i<car;i++) {
    koordx[i]=int(round(size/2+(size/2-rad)*cos(PI-i*2*PI/car)));
    koordy[i]=int(round(size/2-(size/2-rad)*sin(PI-i*2*PI/car)));
  }
}

void Svg::create_circle(int n,ostream* out) {
  *out << "<svg:circle cx=\"" << koordx[n] << "\" cy=\"" << koordy[n] <<
    "\" r=\"" << rad << "\" fill=\"white\" stroke=\"black\" stroke-width=\"1\"/>"<<endl;
}

void Svg::create_numtext(int n,ostream* out) {
  *out << "<svg:text x=\"" << koordx[n] << "\" y=\"" <<
    koordy[n]+round(2/3*fontsize/2) << "\" font-size=\""<< fontsize <<
    "\" text-anchor=\"middle\" dominant-baseline=\"mathematical\">" << n+1 <<
    "</svg:text>"<<endl;
}

void Svg::create_line(int n0,int n1,int szin,ostream* out) {
  int diff=(szin*2*rad/dim-rad)*3/4;
  float length=sqrt(pow(koordx[n1]-koordx[n0],2)+pow(koordy[n1]-koordy[n0],2));
  int diffy=(koordx[n1]-koordx[n0])*diff/(int)length;
  int diffx=-(koordy[n1]-koordy[n0])*diff/(int)length;
  string style;
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
  *out << "<svg:line style=\"" << style << 
    "\" x1=\"" << koordx[n0]+diffx << 
    "\" y1=\"" << koordy[n0]+diffy << 
    "\" x2=\"" << koordx[n1]+diffx << 
    "\" y2=\"" << koordy[n1]+diffy << 
    "\" stroke=\"black\" stroke-width=\"1\"/>" << endl;
}

Svg::~Svg(void) {
  delete[] koordx;
  delete[] koordy;
}

int Svg::print_html(ostream* out){
  *out << "<svg:svg viewBox=\"0 0 "<< size<< " "<< size<<"\" width=\"" << 
    size << "px\" height=\""<<size<<"px\" version=\"1.1\">" <<endl;

  for(list<vector<int> >::iterator it=lines.begin();it!=lines.end();it++)
    create_line((*it)[0],(*it)[1],(*it)[2],out);

  for(int i=0;i<car;i++){
    create_circle(i,out);
    create_numtext(i,out);
  }
  *out << "</svg:svg>" <<endl;
  return 0;
}

int Svg::add_line(int n0,int n1,int szin){
  vector<int> *a=new vector<int>;
  a->push_back(n0);
  a->push_back(n1);
  a->push_back(szin);
  lines.push_back(*a);
  return 0;
}
