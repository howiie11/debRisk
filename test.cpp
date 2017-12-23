#include <propagator.hpp>
using namespace std;

int main(int argc,char* argv[])
{
  //COMMON
  double RE=6.371e6;//m
  double ME=5.98e24;//kg

  //CANONIC UNITS
  double ul=RE;
  double um=ME;
  double ut=sqrt(ul*ul*ul/(GCONST*um));
  double uv=ul/ut;

  //WHICH MODULE DO YOU WANT TO TEST
  goto geopot;//GEOPOTENTIAL
  goto coord;//COORDINATE TRANSFORMATION
  goto integ_2B;//INTEGRATOR TWO-BODY PROBLEM
  goto atmos;//ATMOSPHERIC MODEL
  goto integ_mas;//INTEGRATOR MAS

 coord:
  {
    //FIRST OCTANT
    //          x   y   z   vx  vy  vz
    double c[]={1.0,0.0,0.0,0.0,+1.0,+1.0};
    double s[6];
    cart2sph(s,c);
    fprintf(stdout,"Cartesian: %s\n",vec2strn(c,6,"%f "));
    fprintf(stdout,"Spherical: %s\n",vec2strn(s,6,"%f "));
    exit(0);
  }


 geopot:
  {
    double H=300e3;//m
    double ro=RE+H;
    double vcirc=sqrt(GCONST*ME/ro);//m/s

    //INITIAL CONDITION IN CARTESIAN COORDINATES
    double X0[6];

    double Xc0[]={+ro/ul,0.0,0.0,0.0,+vcirc/uv,0.0};
    fprintf(stdout,"Initial conditions (cartesian): %s\n",vec2strn(Xc0,6,"%f "));

    cart2sph(X0,Xc0);
    fprintf(stdout,"Initial conditions (spherical): %s\n",vec2strn(X0,6,"%f "));

    double dydt[6];
    double params[]={6.0,1.0};

    //Value of the equation
    EoM_Full(0,X0,dydt,params);

    //Integration
    double tini=0.0;
    double h=0.1;
    int npoints=100;
    double duration=90*60.0/ut;
    int nsys=6;
    double params[]={nsys,1.0};
    double *ts=newVector(npoints);
    double **X=newMatrix(npoints,nsys);
    integrateEoM(tini,y,h,npoints,duration,nsys,EoM_Full,params,
		 ts,X);
    
  }

 atmos:
  {
    //*MSIS90E ATMOSPHERIC MODEL
    int i,n=100,qx;
    double x,xo,xe,dx;
    double rho,T;
    double KAP[]=KAP_REF;
    FILE* f=fopen("atmosphere.dat","w");

    int iyd=00001;
    double sec=0.0,alt=1.4,lat=+36.2,lon=-75.34;

    /*Altitude*/xo=0.0;xe=100.0;qx=3;
    /*Latitude*/xo=0.0;xe=90.0;qx=4;
    /*Longitude*/xo=0.0;xe=180.0;qx=5;
    /*Day of year*/xo=1.0;xe=120.0;qx=1;

    /*Second of the day*/xo=0.0;xe=86400;qx=2;
    fprintf(f,"%-25s %-25s %-25s %-25s %-25s\n",
	    "#1:x","2:rho(90)","3:T(90)","4:rho(00)","5:T(00)");
    for(i=0;i<n;i++){
      x=xo+i*(xe-xo)/n;
      iyd=qx==1?(int)x:iyd;
      sec=qx==2?x:sec;
      alt=qx==3?x:alt;
      lat=qx==4?x:lat;
      lon=qx==5?x:lon;
      //printf("Variable : %e\n",x);
      fprintf(f,"%-25.17e ",x);
      MSIS90E(iyd,sec,alt,lat,lon,KAP,&rho,&T);
      fprintf(f,"%-25.17e %-25.17e ",rho,T);
      NRLMSISE(iyd,sec,alt,lat,lon,KAP,&rho,&T);
      fprintf(f,"%-25.17e %-25.17e ",rho,T);
      fprintf(f,"\n");
    }
    fclose(f);
    exit(0);
  }
  //TEST OF MAS
  //*
 integ_mas:
  {
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
    exit(0);
  }
  //*/

  //TEST OF TWO BODY PROBLEM
  //*
 integ_2B:
  {
    double H=300e3;//m
    double ro=RE+H;
    double vcirc=sqrt(GCONST*ME/ro);//m/s

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
    exit(0);
  }
  //*/

  return 0;
}
