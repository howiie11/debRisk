#include <propagator.hpp>
using namespace std;

/*
  Rate of variation of longitude of the nodes and periapse argument:
  http://www.braeunig.us/space/orbmech.htm
  http://www.braeunig.us/space/atmos.htm

 */

int main(int argc,char* argv[])
{
  //COMMON
  double RE=REARTH;//m
  double ME=MUEARTH/GCONST;//kg

  //INITIALIZE PROPAGATOR
  initPropagator();

  //WHICH MODULE DO YOU WANT TO TEST
  goto elem;//ELEMENTS
  goto integ_2B;//INTEGRATOR TWO-BODY PROBLEM
  goto geopot;//GEOPOTENTIAL
  goto coord;//COORDINATE TRANSFORMATION
  goto atmos;//ATMOSPHERIC MODEL
  goto integ_mas;//INTEGRATOR MAS

 elem:
  {
    /*
      You can see elements variations in:
      http://www.braeunig.us/space/orbmech.htm
     */
    
    //INITIAL CONDITIONS
    double H=3200e3;//m
    double ro=REARTH+H;

    //INITIAL ORBITAL ELEMENTS
    double ao=ro/UL;
    double eo=0.2;
    double io=30.0*DEG;
    double Wo=30.0*DEG;
    double wo=60.0*DEG;
    double Mo=0.0*DEG;
    double to=0.0;
    double mu=1.0;
    double P=2*M_PI/sqrt(mu/(ao*ao*ao));
    double h=0.01;
    double duration=10.0*P;
    int npoints=(int)(50*duration/P);//50 is the number of points per orbit
    fprintf(stdout,"P = %e hours\n",P*UT/3600);
    double elts[]={ao*(1-eo),eo,io,Wo,wo,Mo,to,mu};
    fprintf(stdout,"Initial elements: %s\n",vec2strn(elts,8,"%f "));

    //PREDICTED RATE OF VARIATION
    //See: http://www.braeunig.us/space/orbmech.htm and Danby, Pag. 347
    double n=(2*M_PI/(P*UT));
    fprintf(stdout,"n = %e deg/day\n",n*RAD*DAY);
    double J2=J[2][0];
    double J3=J[3][0];
    double face=(1-eo*eo)*(1-eo*eo);
    fprintf(stdout,"J_2 = %e, J_3 = %e\n",J2,J3);
    double dWdt=-1.5*n*J2*pow((REARTH/(ao*UL)),2.0)*cos(io)/face;
    fprintf(stdout,"dW/dt (expected) = %e deg/h\n",dWdt*RAD*DAY);
    double dwdt=+3.0*n*J2*pow((REARTH/(ao*UL)),2.0)*(1-5./4*sin(io)*sin(io))/face;
    fprintf(stdout,"dw/dt (expected) = %e deg/h\n",dwdt*RAD*DAY);
    
    //STATE VECTOR
    double X0[6],Xc0[6];
    double params[]={6.0,mu};
    conics_c(elts,to,Xc0);
    fprintf(stdout,"Initial conditions (cartesian): %s\n",vec2strn(Xc0,6,"%f "));
    car2sph(X0,Xc0);
    fprintf(stdout,"Initial conditions (spherical): %s\n",vec2strn(X0,6,"%f "));

    double Eo=totalEnergy(to,X0,params);
    fprintf(stdout,"Initial energy: %e\n",Eo);

    //Integration
    double tini=0.0;
    fprintf(stdout,"Duration = %e\n",duration);
    double *ts=newVector(npoints);
    double *Es=newVector(npoints);
    double **X=newMatrix(npoints,6);
    double **Xc=newMatrix(npoints,6);
    
    try{
      integrateEoM(tini,X0,h,npoints,duration,6.0,EoM_Fulla,params,
		   ts,X);
    }catch(int e){
      fprintf(stderr,"An exception ocurred: %d (NPS=%d)\n",e,NPS);
      npoints=NPS;
    }

    FILE* fc=fopen("conservation.dat","w");
    //double Eo=totalEnergy(ts[0],X[0],params);
    for(int i=0;i<npoints;i++){
      //COMPUTE TOTAL ENERGY
      Es[i]=totalEnergy(ts[i],X[i],params);
      fprintf(fc,"%.3e\n",fabs((Es[i]-Eo)/Eo));
      //CONVERT TO CARTESIAN
      sph2car(Xc[i],X[i]);
    }
      fclose(fc);
    savetxt("solution-elements.dat",Xc,npoints,6,ts);

    double **E=newMatrix(npoints,8);
    //ELEMENTS ARE: q,e,i,W,w,M,t,mu
    for(int i=0;i<npoints;i++){
      oscelt_c(Xc[i],ts[i],mu,E[i]);
      ts[i]*=UT/3600;//TRANSFORM TO HOURS
    }
    savetxt("elements.dat",E,npoints,8,ts);
    
    //CALCULATE RATE OF VARIATION
    dWdt=(E[npoints-1][3]-E[0][3])*RAD/(duration*UT/DAY);
    fprintf(stdout,"dWdt (observed) = %e deg/day\n",dWdt);
    dwdt=(E[npoints-1][4]-E[0][4])*RAD/(duration*UT/DAY);
    fprintf(stdout,"dWdt (observed) = %e deg/day\n",dwdt);
    exit(0);
  }

 coord:
  {
    {
      //CARTESIAN TO SPHERICAL
      //          x   y   z   vx  vy  vz
      double c[]={1.0,0.0,0.0,0.0,+1.0,+1.0};
      double s[6];
      car2sph(s,c);
      fprintf(stdout,"Cartesian to Spherical:\n");
      fprintf(stdout,"\tCartesian: %s\n",vec2strn(c,6,"%f "));
      fprintf(stdout,"\tSpherical: %s\n",vec2strn(s,6,"%f "));
    }

    {
      //SPHERICAL TO CARTESIAN
      //          r    f    q    vr   vf   vq
      double s[]={+1.0,+1.57,+3.14,+0.0,+1.0,+0.0};
      double c[6];
      sph2car(c,s);
      fprintf(stdout,"Spherical to Cartesian:\n");
      fprintf(stdout,"\tSpherical: %s\n",vec2strn(s,6,"%f "));
      fprintf(stdout,"\tCartesian: %s\n",vec2strn(c,6,"%f "));
    }
    exit(0);
  }

 geopot:
  {
    double H=300e3;//m
    double ro=REARTH+H;
    double vcirc=sqrt(GCONST*MEARTH/ro);//m/s

    //INITIAL CONDITION IN CARTESIAN COORDINATES
    double X0[6];

    double Xc0[]={+ro/UL,0.0,0.0,0.0,0.0,vcirc/UV};
    fprintf(stdout,"Initial conditions (cartesian): %s\n",vec2strn(Xc0,6,"%f "));

    car2sph(X0,Xc0);
    fprintf(stdout,"Initial conditions (spherical): %s\n",vec2strn(X0,6,"%f "));

    double dydt[6];
    double params[]={6.0,1.0};

    //Value of the equation
    EoM_Full(0,X0,dydt,params);

    //Integration
    double tini=0.0;
    double h=0.01;
    int npoints=100;
    double duration=90*60.0/UT;
    double *ts=newVector(npoints);
    double **X=newMatrix(npoints,6);
    double **Xc=newMatrix(npoints,6);
    integrateEoM(tini,X0,h,npoints,duration,6.0,EoM_Fulla,params,
		 ts,X);
    fprintf(stdout,"Duration = %e\n",duration);
    
    //TRANSFORM TO CARTESIAN
    //*
    for(int i=0;i<npoints;i++) sph2car(Xc[i],X[i]);
    savetxt("solution-spherical.dat",Xc,npoints,6,ts);
    //*/
    //savetxt("solution-spherical.dat",X,npoints,6,ts);

    //goto integ_2B;
    exit(0);
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
    double ro=REARTH+H;
    double vcirc=sqrt(GCONST*MEARTH/ro);//m/s

    double tini=0.0;
    double X0[]={+ro/UL,0.0,0.0,0.0,1.0*vcirc/UV,0.0};
    fprintf(stdout,"Initial conditions (cartesian): %s\n",vec2strn(X0,6,"%f "));

    double h=0.01;
    int npoints=1000;
    double ao=ro/UL;
    double mu=1.0;
    double P=2*M_PI/sqrt(mu/(ao*ao*ao));
    fprintf(stdout,"P = %e\n",P);
    double duration=P/2;
    fprintf(stdout,"Duration = %e\n",duration);
    int nsys=6;
    double params[]={nsys,1.0,MEARTH/UM};
  
    double *ts=newVector(npoints);
    double **X=newMatrix(npoints,nsys);
    double **Xs=newMatrix(npoints,nsys);
  
    //CALL THE INTEGRATION ROUTINE
    fprintf(stdout,"Integrating cartesian coordinates\n");
    integrateEoM(tini,X0,h,npoints,duration,nsys,EoM_2B,params,
		 ts,X);

    //TRANSFORM TO SPHERICAL
    /*
    for(int i=0;i<npoints;i++) car2sph(Xs[i],X[i]);
    savetxt("solution-cartesian.dat",Xs,npoints,6,ts);    
    //*/
    savetxt("solution-cartesian.dat",X,npoints,6,ts);

    difftxt("solution-cartesian.dat","solution-spherical.dat",npoints,7);
    /*
    //SAVE OUTPUT
    FILE* f=fopen("solution-cartesian.dat","w");
    for(int i=0;i<npoints;i++){
      fprintf(f,"%.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
	      ts[i],
	      X[i][0]*UL/1e3,X[i][1]*UL/1e3,X[i][2]*UL/1e3,
	      X[i][3]*UV,X[i][4]*UV,X[i][5]*UV
	      );
    }
    */
    exit(0);
  }
  //*/

  return 0;
}
