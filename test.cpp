#include <propagator.hpp>
using namespace std;

int main(int argc,char* argv[])
{
  //TEST OF MAS
  /*
  double tini=0.0;
  double X0[]={1.0,0.0};
  double h=0.1;
  int npoints=100;
  double duration=10.0;
  int nsys=2;
  double w=1.0;
  double params[]={nsys,w};
  
  double *ts=newVector(npoints);
  double **X=newMatrix(npoints,nsys);
  
  //CALL THE INTEGRATION ROUTINE
  integrateEoM(tini,X0,h,npoints,duration,nsys,EoM_MAS,params,
	       ts,X);

  //SAVE OUTPUT
  FILE* f=fopen("solucion.dat","w");
  for(int i=0;i<npoints;i++){
    fprintf(f,"%.17e %.17e %.17e\n",ts[i],X[i][0],X[i][1]);
  }
  //*/

  //TEST OF TWO BODY PROBLEM
  //*
  double RE=6.371e6;//m
  double ME=5.98e24;//kg
  double H=300e3;//m
  double ro=RE+H;
  double vcirc=sqrt(GCONST*ME/ro);//m/s

  //CANONIC UNITS
  double ul=RE;
  double um=ME;
  double ut=sqrt(ul*ul*ul/(GCONST*um));
  double uv=ul/ut;

  double tini=0.0;
  double X0[]={+ro/ul,0.0,0.0,0.0,0.5*vcirc/uv,0.0};
  double h=0.1;
  int npoints=100;
  double duration=90*60.0/ut;
  int nsys=6;
  double params[]={nsys,1.0,ME/um};
  
  double *ts=newVector(npoints);
  double **X=newMatrix(npoints,nsys);
  
  //CALL THE INTEGRATION ROUTINE
  integrateEoM(tini,X0,h,npoints,duration,nsys,EoM_2B,params,
	       ts,X);

  //SAVE OUTPUT
  FILE* f=fopen("solucion.dat","w");
  for(int i=0;i<npoints;i++){
    fprintf(f,"%.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
	    ts[i],
	    X[i][0]*ul/1e3,X[i][1]*ul/1e3,X[i][2]*ul/1e3,
	    X[i][3]*uv,X[i][4]*uv,X[i][5]*uv
	    );
  }
  //*/

  return 0;
}
