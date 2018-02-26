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

  //VVVV--- PUT HERE THE MODULE YOU WANT TO TEST
  goto satint;//INTEGRATE WITH ALL THE EFFECTS
  goto sat;//SATELLITE PACKAGE
  
  //THIS LINES WILL NOT BE EXECUTED
  goto elem;//ELEMENTS
  goto atmos;//ATMOSPHERIC MODEL
  goto integ_2B;//INTEGRATOR TWO-BODY PROBLEM
  goto geopot;//GEOPOTENTIAL
  goto coord;//COORDINATE TRANSFORMATION
  goto integ_mas;//INTEGRATOR MAS

 satint:
  {
    //*
    //All effects
    struct satparams params={
      .nsys=6,
      .Area=5.0,
      .mass=1000.0,
      .CR=1.3,
      .CD=2.3,
      .n=3,
      .m=3,
      .qSun  = true,
      .qMoon = true,
      .qSRad = true,
      .qDrag = true,
    };
    //*/
    /*
    //Partial effects
    struct satparams params={
      .nsys=6,
      .Area=5.0,
      .mass=1000.0,
      .CR=1.3,
      .CD=2.3,
      .n=3,
      .m=3,
      .qSun  = false,
      .qMoon = false,
      .qSRad = false,
      .qDrag = false,
    };
    //*/
    double c1,c2;

    //INITIAL CONDITIONS
    double H=350e3;//m
    double ro=REARTH+H;

    //INITIAL TIME
    double Mjd0_UTC = Mjd(1999,03,01,00,00,0.0);
    double Mjd0_TT = Mjd0_UTC + IERS::TT_UTC(Mjd0_UTC)/86400.0;
    fprintf(stdout,"Mjd0(UTC) = %.9e\n",Mjd0_UTC);
    fprintf(stdout,"Mjd0(TT) = %.9e\n",Mjd0_TT);

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
    double h=0.1;
    double duration=10.0*P;
    
    int npoints=(int)(50*duration/P);//50 is the number of points per orbit
    fprintf(stdout,"Npoints: %d\n",npoints);
    fprintf(stdout,"P = %e hours\n",P*UT/3600);

    //BUILD ELEMENTS VECTOR
    double elts[]={ao*(1-eo),eo,io,Wo,wo,Mo,to,mu};
    fprintf(stdout,"Initial elements: %s\n",vec2strn(elts,8,"%f "));
    
    //INITIAL STATE VECTOR IN CARTESIAN COORDINATES
    double Xc0[6],dXc0dt[]={0,0,0,0,0,0};
    conics_c(elts,to,Xc0);
    fprintf(stdout,"Initial conditions (cartesian): %s\n",vec2strn(Xc0,6,"%f "));

    //Integration
    double tini=0.0;
    fprintf(stdout,"Duration = %e\n",duration);

    double *ts=newVector(npoints);
    double *Es=newVector(npoints);
    double **X=newMatrix(npoints,6);
    double **Xc=newMatrix(npoints,6);

    //EoM_Satellite(0,Xc0,dXc0dt,&params);
    //INTEGRATE
    try{
      integrateEoM(tini,Xc0,h,npoints,duration,6.0,EoM_Satellite,&params,ts,Xc);
    }catch(int e){
      fprintf(stderr,"An exception ocurred: %d (NPS=%d)\n",e,NPS);
      npoints=NPS;
    }

    //Save solution
    savetxt("solution-coordinates.dat",Xc,npoints,6,ts);
    exit(0);
  }

 sat:
  {
    //Testing the vector components
    Vector x=Vector(1,0,1);
    Vector y=Vector(2,1,3);
    Vector z=x+y;
    fprintf(stdout,"%e\n",z(0));
    
    //Calculating the acceleration
    //Compute acceleration
    double Mjd0_UTC = Mjd(1999,03,01,00,00,0.0);
    double Mjd0_TT = Mjd0_UTC + IERS::TT_UTC(Mjd0_UTC)/86400.0;
    fprintf(stdout,"Mjd0(UTC) = %.9e\n",Mjd0_UTC);
    fprintf(stdout,"Mjd0(TT) = %.9e\n",Mjd0_TT);

    double Area    = 5.0;     // [m^2]  Remote sensing satellite
    double mass    = 1000.0;  // [kg]
    double CR      = 1.3;     
    double CD      = 2.3;
    int n       = 20;
    int m       = 20;
    double qSun     = true;
    double qMoon    = true;
    double qSRad    = true;
    double qDrag    = true;

    Vector r(3),v(3);
    double dens;

    double h=1000000.0;//m
    r=Vector(1.0,0.0,0.0)*(Grav.R_ref + h);
    fprintf(stdout,"Position = %e\n",r(0));

    //Model integrated with sat only compute density above 100 km
    dens=Density_HP(Mjd0_TT,r);
    fprintf(stdout,"Density = %e\n",dens);

    Vector Kep(6),Y0(6),a(3);
    //Kep = Vector ( 7178.0e3, 0.0010, 98.57*Rad, 0.0, 0.0, 0.0 ); 
    //             a              e    i         W         w         M
    Kep = Vector ( REARTH+3200e3, 0.2, 30.0*Rad, 30.0*Rad, 60.0*Rad, 0.0*Rad ); 
    Y0 = State ( GM_Earth, Kep, 0.0 );
    r = Y0.slice(0,2);
    v = Y0.slice(3,5);
    VPRINT("Conditions:\n\tx = %e m, y = %e m, z = %e m\n\tdx/dt = %e m/s, dy/dt = %e m/s, dz/dt = %e m/s\n",
	   r(0),r(1),r(2),v(0),v(1),v(2));
    
    a=Accel(Mjd0_TT,r,v,Area,mass,CR,CR,n,m,qSun,qMoon,qSRad,qDrag);
    fprintf(stdout,"Acceleration = %.17e, %.17e, %.17e\n",a(0),a(1),a(2));

    exit(0);
  }

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

    //THESE LINES DOES NOT HAVE ANY EFFECT
    /*Latitude*/xo=0.0;xe=90.0;qx=4;
    /*Longitude*/xo=0.0;xe=180.0;qx=5;
    /*Day of year*/xo=1.0;xe=120.0;qx=1;
    /*Second of the day*/xo=0.0;xe=86400;qx=2;

    //VVVV--- PUT HERE THE VARIABLE YOU WANT TO TEST
    /*Altitude*/xo=0.0;xe=1000.0;qx=3;

    fprintf(f,"%-25s %-25s %-25s\n",
	    "#1:x","2:rho(90)","3:T(90)");
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
