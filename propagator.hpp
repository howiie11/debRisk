////////////////////////////////////////////////////////////////////////
//CONFIGURATION
////////////////////////////////////////////////////////////////////////
#define VERBOSE 0

////////////////////////////////////////////////////////////////////////
//BASIC LIBRARIES
////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdexcept>
#include <iostream>

////////////////////////////////////////////////////////////////////////
//ADDITIONAL LIBRARIES
////////////////////////////////////////////////////////////////////////
#include <nrlmsise-00.h>
#include <SpiceUsr.h>
#include <eph_manager.h>
#include <SAT_Const.h>
#include <SAT_VecMat.h>
#include <SAT_RefSys.h>
#include <SAT_Force.h>

////////////////////////////////////////////////////////////////////////
//GSL
////////////////////////////////////////////////////////////////////////
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_legendre.h>

////////////////////////////////////////////////////////////////////////
//BASIC MACROS
////////////////////////////////////////////////////////////////////////
#define RAD 180/M_PI
#define DEG M_PI/180
#define EXTMET 1
#define TOLERANCE 1E-10
#define ATTEMPTS 12 /*SEE NUMBER_OF_STEPS*/
#define VPRINT if(VERBOSE) printf
#define PAUSE if(VERBOSE) getchar

////////////////////////////////////////////////////////////////////////
//PHYSICAL CONSTANTS
////////////////////////////////////////////////////////////////////////
#define GCONST GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT //m^3/kg s^2

//WGS84
#define REARTH 6.378137E6 //m
#define FEARTH 1/298.257223563
#define MUEARTH 3.986004418E14 //m^3/s^2
#define WEARTH 7.292115E-5 //rad/s

#define HOUR 3600.0
#define DAY 86400.0

////////////////////////////////////////////////////////////////////////
//GLOBAL VARIABLES
////////////////////////////////////////////////////////////////////////
static int NUMBER_OF_STEPS[]={2,4,6,8,12,16,24,32,48,64,96,128};
gsl_rng* RAND;

#define MAXN 4
#define MAXM 4
double C[][MAXM+1]={{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}};
double S[][MAXM+1]={{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}};
double J[][MAXM+1]={{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}};
double Q[][MAXM+1]={{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}};
double UL,UM,UT,UV;
double MEARTH;

////////////////////////////////////////////////////////////////////////
//INITIALIZE
////////////////////////////////////////////////////////////////////////
//void errorGSL(const char * reason,const char * file,int line,int gsl_errno){throw(1);}
int initPropagator(void)
{
  SpiceInt i,n;
  
  //CONFIGURATION
  #include<propagator.conf>

  //INITIALIZE GSL ERROR HANDLER
  //gsl_set_error_handler(&errorGSL); 

  //KERNELS
  furnsh_c("util/kernels/kernels.txt");

  //RANDOM NUMBERS
  RAND=gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(RAND,time(NULL));
  gsl_rng_set(RAND,3);

  //INIT ZONAL AND TESSERAL HARMINICS
  #include<harmonics-earth.hpp>
  //#include<harmonics-ideal.hpp>
  //#include<harmonics-spherical.hpp>
  //#include<harmonics-artificial.hpp>
  for(int n=MAXM+1;n-->0;){
    J[n][0]=C[n][0];
    for(int m=MAXN+1;m-->1;){
      J[n][m]=sqrt(C[n][m]*C[n][m]+S[n][m]*S[n][m]);
      if(m>0)
	Q[n][m]=(1/m)*atan2(S[n][m],C[n][m]);
    }
  }

  //DERIVED CONSTANTS
  MEARTH=MUEARTH/GCONST;

  //UNITS
  UT=sqrt(UL*UL*UL/(GCONST*UM));
  UV=UL/UT;
  
  return 0;
}

////////////////////////////////////////////////////////////////////////
//UTIL ROUTINES
////////////////////////////////////////////////////////////////////////
int copyVec(double tgt[],double src[],int n)
{
  memcpy(tgt,src,n*sizeof(double));
  return 0;
}

int copyVeci(int tgt[],int src[],int n)
{
  memcpy(tgt,src,n*sizeof(int));
  return 0;
}

//c=ca*a+cb*b
int sumVec(double c[],
	   double ca,double a[],
	   double cb,double b[],
	   int n=3)
{
  int i;
  for(i=n;i-->0;) c[i]=ca*a[i]+cb*b[i];
  return 0;
}

//Find the maximum in a
double maxAbsVec(double a[],int n)
{
  int i;
  double max=-1E100;
  for(i=n;i-->0;) if(fabs(a[i])>max) max=fabs(a[i]);
  return max;
}

double *newVector(int n)
{
  double *v;
  
  v=(double*)malloc(n*sizeof(double));
  
  return v;
}

double **newMatrix(int n,int m)
{
  double **M;
  
  M=(double**)malloc(n*sizeof(double*));
  for(int i=0;i<n;i++) M[i]=(double*)malloc(m*sizeof(double));
  
  return M;
}

char* dec2sex(double dec)
{
  double d,m,s;
  int id,im,sgn;
  char *str=(char*)calloc(sizeof(char),100); 
  d=fabs(dec);
  sgn=dec/d;
  id=floor(d);
  m=(d-id)*60;
  im=floor(m);
  s=(m-im)*60;
  sprintf(str,"%+d:%02d:%.3f",sgn*id,im,s);
  return str;
}

double sex2dec(double d,double m,double s)
{
  double s2d;
  s2d=d+m/60.0+s/3600.0;
  return s2d;
}

char* vec2str(double vec[],char frm[]="%.8e ")
{
  char format[100];
  char *str=(char*)calloc(sizeof(char),100); 
  sprintf(format,"%s %s %s",frm,frm,frm);
  sprintf(str,format,vec[0],vec[1],vec[2]);
  return str;
}

char* vec2strn(double vec[],int n,char frm[]="%.8e ")
{
  int i;
  char format[100];
  char *str=(char*)calloc(sizeof(char),100*n);
  sprintf(format,"%ss%s","%",frm);
  for(i=0;i<n;i++) sprintf(str,format,str,vec[i]);
  return str;
}

/*
  Save matriz into file
 */
int savetxt(char* fname,double** x,int n,int m,double* zcol=NULL)
{
    //SAVE OUTPUT
    FILE* f=fopen(fname,"w");
    for(int i=0;i<n;i++){
      if(zcol!=NULL)
	fprintf(f,"%.17e ",zcol[i]);
      fprintf(f,"%s\n",vec2strn(x[i],m,"%.17e "));
    }
    fclose(f);
    return 0;
}

/*
  Difference between files
*/
int difftxt(char* fname1,char *fname2,int n,int m)
{
  double v1,v2;
  char fdiff[1000];
  FILE *f1=fopen(fname1,"r");
  FILE *f2=fopen(fname2,"r");
  sprintf(fdiff,"%s_%s.diff",fname1,fname2);
  FILE *fd=fopen(fdiff,"w");
  for(int i=0;i<n;i++){
    fprintf(fd,"%d ",i);
    for(int j=0;j<m;j++){
      fscanf(f1,"%lf",&v1);
      fscanf(f2,"%lf",&v2);
      fprintf(fd,"%.17e ",fabs(v1-v2));
      //VPRINT("i=%d,j=%d: v1 = %e, v2=%e, v1-v2 = %e\n",i,j,v1,v2,fabs(v1-v2));
    }
    fprintf(fd,"\n");
  }
  fclose(f1);
  fclose(f2);
  fclose(fd);
  return 0;
}

/*
  Transform vector state from cartesian to spherical

  c: x,y,z,vx,vy,vz
  s: r,f(respect to plane xy),q(azimutal),vr,vq,vf

  In spice longitude = q, latitude = f

  See: http://www.astrosurf.com/jephem/library/li110spherCart_en.htm

  NOTES: Signs of the formulas must be corrected.
 */
int car2sph(double s[6],double c[6])
{
  //Convert cartesian coordinates to spherical
  reclat_c(c,&s[0],&s[2],&s[1]);

  //Convert velocity
  double rho2=c[0]*c[0]+c[1]*c[1];

  /*r.v*/double rv=vdot_c(c,c+3);
  /*dr*/s[3]=rv/s[0];
  /*df*/s[4]=-(c[2]*(rv-c[2]*c[5])-rho2*c[5])/(s[0]*s[0]*sqrt(rho2));
  /*dq*/s[5]=-(c[3]*c[1]-c[0]*c[4])/rho2;

  return 0;
}
int sph2car(double c[6],double s[6])
{
  latrec_c(s[0],s[2],s[1],c);
  double rho=sqrt(c[0]*c[0]+c[1]*c[1]);

  //Convert velocity
  /*Vx*/c[3]=c[0]/s[0]*s[3]-c[1]*s[5]-c[2]*c[0]/rho*s[4];
  /*Vy*/c[4]=c[1]/s[0]*s[3]+c[0]*s[5]-c[2]*c[1]/rho*s[4];
  /*Vz*/c[5]=c[2]/s[0]*s[3]+rho*s[4];

  return 0;
}
int vec_car2sph(double vs[3],double vc[3],double s[3])
{
  double q=M_PI/2-s[1];
  double f=s[2];
  double sinq=sin(q),cosq=cos(q),sinf=sin(f),cosf=cos(f);
  double T[][3]={{sinq*cosf,sinq*sinf,cosq},
		 {cosq*cosf,cosq*sinf,-sinq},
		 {-sinf,cosf,0}};
  mxv_c(T,vc,vs);
  double vr=vs[0],vq=vs[2],vf=-vs[1];
  vs[0]=vr;vs[1]=vf;vs[2]=vq;
  return 0;
}
int vec_sph2car(double vc[3],double vs[3],double s[3])
{
  double q=M_PI/2-s[1];
  double f=s[2];
  double sinq=sin(q),cosq=cos(q),sinf=sin(f),cosf=cos(f);
  double T[][3]={{sinq*cosf,cosq*cosf,-sinf},
		 {sinq*sinf,cosq*sinf,+cosf},
		 {cosq,-sinq,0}};
  mxv_c(T,vs,vc);
  return 0;
}

////////////////////////////////////////////////////////////////////////
//CORE ROUTINES
////////////////////////////////////////////////////////////////////////
/*
Adapted from: 
http://www.mymathlib.com/diffeq/bulirsch_stoer.html
 */
static int Rational_Extrapolation_to_Zero(double *fzero,double tableau[],
					  double x[],double f,int n) 
{
  double t, up, across, denominator, dum;
  int col;

  if (n==0) {  *fzero = f; tableau[0] = f; return 0; }
  if ( x[n] == 0.0 ) { *fzero = f; return -2; }
   
  across = 0.0;                                                        
  up = tableau[0];                                                    
  tableau[0] = f;                                               

  for (col = 1; col <= n; col++) {
    if(tableau[col-1]==0 && across==0){t=0;break;}
    denominator = tableau[col-1] - across;
    if (denominator == 0.0) return -1;
    dum = 1.0 - (tableau[col-1] - up) / denominator;
    denominator = (x[n - col] / x[n]) * dum - 1.0;
    if (denominator == 0.0) return -1;
    t = tableau[col-1] + ( tableau[col-1] - up ) / denominator;
    across = up;
    up = tableau[col];
    tableau[col] = t;
  }
  *fzero = t;
  return 0;
}

static int Polynomial_Extrapolation_to_Zero(double *fzero,double tableau[],
					    double x[], double f, int n )
{
  double back_two_columns;    //  T[row,col-2];
  double old_aux;             //  T[row-1,col];
  double new_value;           //  T[row,col];
  double vertical_diff;       //  T[row,col]-T[row-1,col]
  double backslant_diff;      //  T[row,col]-T[row,col-1]
  double forwardslant_diff;   //  T[row,col]-T[row-1,col-1];
  double denominator;        
  int i;

  if (n == 0) { tableau[0] = f; return 0; }
  if ( x[n] == 0.0 ) { *fzero = f; return -2; }

  back_two_columns = 0.0;
  old_aux = tableau[0];
  tableau[0] = f;
  for (i = 0; i < n; i++) {
    if(tableau[i]==0 && old_aux==0){tableau[n]=0.0;break;}
    vertical_diff = tableau[i] - old_aux;
    backslant_diff = tableau[i] - back_two_columns;
    forwardslant_diff = backslant_diff - vertical_diff;
    denominator = (x[n-i-1]/x[n]) * forwardslant_diff - backslant_diff;
    if (denominator == 0.0) return -1;
    back_two_columns = old_aux;
    old_aux = tableau[i+1];
    tableau[i+1] = tableau[i] + vertical_diff * backslant_diff / denominator;
  }
  *fzero = tableau[n];
  return 0;
}

static int Graggs_Method(int (*f)(double,double*,double*,void*),
			 double y0[],
			 double t0,double t,
			 int NUMBER_OF_STEPS,
			 void *params,
			 double yres[]) {
  
  double* pars=(double*)params;
  int order=(int)pars[0],i;
  double y1[order],dydt[order],y2[order],yaux[order];
  double h=(t-t0)/(double)NUMBER_OF_STEPS;
  double h2=h+h;

  copyVec(yaux,y0,order);
  (*f)(t0,yaux,dydt,params);
  sumVec(y1,1,yaux,h,dydt,order);

  while(--NUMBER_OF_STEPS) {
    t0+=h;
    (*f)(t0,y1,dydt,params);
    sumVec(y2,1,yaux,h2,dydt,order);
    copyVec(yaux,y1,order);
    copyVec(y1,y2,order);
  } 

  (*f)(t,y1,dydt,params);

  sumVec(yres,0.5,yaux,0.5,y1,order);
  sumVec(yres,1,yres,0.5*h,dydt,order);
  return 0;
}

int Gragg_Bulirsch_Stoer(int (*f)(double,double*,double*,void*), 
			 double y0[], double y1[],
			 double t, double h, double *h_new, 
			 double epsilon, double yscale,
			 void *params)
{
  double* pars=(double*)params;
  int order=(int)pars[0];
  double step_size2[ATTEMPTS];
  double tableau[order][ATTEMPTS+1];
  double dum;
  double est[order],dest[order],destmax;
  double old_est[order];
  
  int (*Extrapolate)(double*,double*,double*,double,int);
  int i,j;
  int err;

  if(yscale==0.0) return -3;
  if(EXTMET) Extrapolate=Rational_Extrapolation_to_Zero;
  else Extrapolate=Polynomial_Extrapolation_to_Zero;
 
  Graggs_Method(f,y0,t,t+h,NUMBER_OF_STEPS[0],params,est);
  step_size2[0]=(dum=h/(double)NUMBER_OF_STEPS[0],dum*dum);
  
  copyVec(y1,est,order);
  
  for(i=order;i-->0;){
    err=Extrapolate(&y1[i],tableau[i],step_size2,est[i],0);
    if(err<0) return err-1;
  }

  for(i = 1; i < ATTEMPTS; i++) {
    copyVec(old_est,y1,order);
    Graggs_Method(f,y0,t,t+h,NUMBER_OF_STEPS[i],params,est);
    step_size2[i]=(dum=h/(double)NUMBER_OF_STEPS[i],dum*dum);

    for(j=order;j-->0;){
      err=Extrapolate(&y1[j],tableau[j],step_size2,est[j],i);
      if(err<0) return err-1;
    }
    
    sumVec(dest,1.0/yscale,y1,-1.0/yscale,old_est,order);
    destmax=maxAbsVec(dest,order);

    if(destmax<epsilon){
      if(i>1) *h_new=8.0*h/(double)NUMBER_OF_STEPS[i-1];
      else *h_new=h;
      return 0;
    }
  }
  return -1;
}

////////////////////////////////////////////////////////////////////////
//INTEGRATION ROUTINE
////////////////////////////////////////////////////////////////////////
int NPS;
int integrateEoM(double tini,double X0[],double h,int npoints,double duration,
		 int nsys,int eom(double,double*,double*,void*),void *params,
		 double *ts,double** X)
{
  //GENERAL VARIABLES
  double deltat,h_next,h_adjust;
  int status;
  double x0[nsys],x[nsys];

  //DIRECTION OF INTEGRATION
  double direction=duration/fabs(duration);
  h*=direction;

  //PARAMETERS OF THE TIME INTERVAL
  double t_start=tini;
  double t_step=duration/(npoints-1);
  double tend=t_start+duration;
  double t_stop=tend;
  double t=t_start;
  double h_used=h;
  int p0=(int)*(double*)params;

  //INITIAL CONDITIONS
  copyVec(x0,X0,nsys);
  copyVec(x,x0,nsys);

  //CHECK NSYS
  if(nsys!=p0)
    fprintf(stderr,"System dimension (nsys = %d) and parameter vector (params[0] = %d) does not coincide. Please fix.\n",nsys,p0);

  //INTEGRATE
  NPS=0;
  for(int i=0;i<npoints;i++) {
    ts[i]=t;
    
    /*
    VPRINT("t=%e, h=%e, y=%e %e %e %e %e %e\n",
	    t,h_used,x0[0],x0[1],x0[2],x0[3],x0[4],x0[5]);
    getchar();
    //*/

    copyVec(X[i],x0,nsys);
    NPS++;

    deltat=t-tini;
    if(direction*((t_start+t_step)-tend)>0) t_step=(tend-t_start);
    t_stop=t_start+t_step;
    h_used=h;

    do {
      while(1){
	status=Gragg_Bulirsch_Stoer(eom,x0,x,t,h_used,&h_next,1.0,
				    TOLERANCE,params);

	/*
	VPRINT("t=%e, y=%e %e %e %e %e %e, hnext = %e, %d\n",
		t,x[0],x[1],x[2],x[3],x[4],x[5],h_next,direction);
	getchar();
	//*/
	
	if(status){
	  h_used/=4.0;
	  if(direction*h_used<0) h_used=h;
	}
	else break;
      }
      t+=h_used;
      copyVec(x0,x,nsys);

      if(direction*(t+h_next-t_stop)>0) h_used=h_next+(t_stop-(t+h_next));
      else h_used=h_next;

    }while(direction*(t-(t_stop-direction*fabs(t_step)*1.e-7))<0);

    if(direction*(t-t_stop)>0){
      h_adjust=(t_stop-t);
      status=Gragg_Bulirsch_Stoer(eom,x0,x,t,h_adjust,&h_next,1.0,
				  TOLERANCE,params);
      copyVec(x0,x,nsys);
      t=t_stop;
    }
    t_start = t;
    if(direction*(t_start-tend)>0) break;
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////
//ATMOSPHERIC MODEL NRLMSISE
////////////////////////////////////////////////////////////////////////
/*
  Sources:
    https://www.brodo.de/space/nrlmsise/
    https://git.linta.de/?p=~brodo/nrlmsise-00.git;a=tree

 *   Switches: to turn on and off particular variations use these switches.
 *   0 is off, 1 is on, and 2 is main effects off but cross terms on.
 *
 *   Standard values are 0 for switch 0 and 1 for switches 1 to 23. The 
 *   array "switches" needs to be set accordingly by the calling program. 
 *   The arrays sw and swc are set internally.
 *
 *   switches[i]:
 *    i - explanation
 *   -----------------
 *    0 - output in meters and kilograms instead of centimeters and grams
 *    1 - F10.7 effect on mean
 *    2 - time independent
 *    3 - symmetrical annual
 *    4 - symmetrical semiannual
 *    5 - asymmetrical annual
 *    6 - asymmetrical semiannual
 *    7 - diurnal
 *    8 - semidiurnal
 *    9 - daily ap [when this is set to -1 (!) the pointer
 *                  ap_a in struct nrlmsise_input must
 *                  point to a struct ap_array]
 *   10 - all UT/long effects
 *   11 - longitudinal
 *   12 - UT and mixed UT/long
 *   13 - mixed AP/UT/LONG
 *   14 - terdiurnal
 *   15 - departures from diffusive equilibrium
 *   16 - all TINF var
 *   17 - all TLB var
 *   18 - all TN1 var
 *   19 - all S var
 *   20 - all TN2 var
 *   21 - all NLB var
 *   22 - all TN3 var
 *   23 - turbo scale height var

 * Array containing the following magnetic values:
 *   0 : daily AP
 *   1 : 3 hr AP index for current time
 *   2 : 3 hr AP index for 3 hrs before current time
 *   3 : 3 hr AP index for 6 hrs before current time
 *   4 : 3 hr AP index for 9 hrs before current time
 *   5 : Average of eight 3 hr AP indicies from 12 to 33 hrs 
 *           prior to current time
 *   6 : Average of eight 3 hr AP indicies from 36 to 57 hrs 
 *           prior to current time 

    int year;      / year, currently ignored /
    int doy;       / day of year /
    double sec;    / seconds in day (UT) /
    double alt;    / altitude in kilometers /
    double g_lat;  / geodetic latitude /
    double g_long; / geodetic longitude /
    double lst;    / local apparent solar time (hours), see note below /
    double f107A;  / 81 day average of F10.7 flux (centered on doy) /
    double f107;   / daily F10.7 flux for previous day /
    double ap;     / magnetic index(daily) /
    struct ap_array *ap_a; / see above /

 *   NOTES ON INPUT VARIABLES: 
 *      UT, Local Time, and Longitude are used independently in the
 *      model and are not of equal importance for every situation.  
 *      For the most physically realistic calculation these three
 *      variables should be consistent (lst=sec/3600 + g_long/15).
 *      The Equation of Time departures from the above formula
 *      for apparent local time can be included if available but
 *      are of minor importance.
 *
 *      f107 and f107A values used to generate the model correspond
 *      to the 10.7 cm radio flux at the actual distance of the Earth
 *      from the Sun rather than the radio flux at 1 AU. The following
 *      site provides both classes of values:
 *      ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/
 *
 *      f107, f107A, and ap effects are neither large nor well
 *      established below 80 km and these parameters should be set to
 *      150., 150., and 4. respectively.


 *   OUTPUT VARIABLES:
 *      d[0] - HE NUMBER DENSITY(CM-3)
 *      d[1] - O NUMBER DENSITY(CM-3)
 *      d[2] - N2 NUMBER DENSITY(CM-3)
 *      d[3] - O2 NUMBER DENSITY(CM-3)
 *      d[4] - AR NUMBER DENSITY(CM-3)                       
 *      d[5] - TOTAL MASS DENSITY(GM/CM3) [includes d[8] in td7d]
 *      d[6] - H NUMBER DENSITY(CM-3)
 *      d[7] - N NUMBER DENSITY(CM-3)
 *      d[8] - Anomalous oxygen NUMBER DENSITY(CM-3)
 *      t[0] - EXOSPHERIC TEMPERATURE
 *      t[1] - TEMPERATURE AT ALT
 * 
 *
 *      O, H, and N are set to zero below 72.5 km
 *
 *      t[0], Exospheric temperature, is set to global average for
 *      altitudes below 120 km. The 120 km gradient is left at global
 *      average value for altitudes below 72 km.
 *
 *      d[5], TOTAL MASS DENSITY, is NOT the same for subroutines GTD7 
 *      and GTD7D
 *
 *        SUBROUTINE GTD7 -- d[5] is the sum of the mass densities of the
 *        species labeled by indices 0-4 and 6-7 in output variable d.
 *        This includes He, O, N2, O2, Ar, H, and N but does NOT include
 *        anomalous oxygen (species index 8).
 *
 *        SUBROUTINE GTD7D -- d[5] is the "effective total mass density
 *        for drag" and is the sum of the mass densities of all species
 *        in this model, INCLUDING anomalous oxygen.
  
  
 */
#define KAP_REF {150,22,22,22,22,22,22,22}
#define FLAGS_ALL {1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}
int NRLMSISE(int day,double sec,
	     double alt,double glat,double glon,
	     double KAP[],
	     double *rho,double *T)
{
  struct nrlmsise_output output;
  struct nrlmsise_input input;
  struct nrlmsise_flags flags;
  struct ap_array aph;
  int aflags[]=FLAGS_ALL;

  int k=0,j=2;
  aph.a[k++]=KAP[j++];aph.a[k++]=KAP[j++];aph.a[k++]=KAP[j++];
  aph.a[k++]=KAP[j++];aph.a[k++]=KAP[j++];aph.a[k++]=KAP[j++];
  aph.a[k++]=KAP[j++];
  copyVeci(flags.switches,aflags,24);

  input.doy=day;
  input.year=0; /* without effect */
  input.sec=sec;
  input.alt=alt;
  input.g_lat=glat;
  input.g_long=glon;
  input.lst=sec/3600.0+glon/15.0;
  input.f107A=KAP[0];
  input.f107=KAP[1];
  input.ap=40.0;

  gtd7(&input,&flags,&output);  

  *rho=output.d[5];
  *T=output.t[1];
  
  return 0;
}

////////////////////////////////////////////////////////////////////////
//PHYSICAL ROUTINES
////////////////////////////////////////////////////////////////////////
/*
 */
double earthRadius(double lon,double lat)
{
  double epos[3];
  georec_c(lon,lat,0.0,REARTH/UL,FEARTH,epos);
  double rs=vnorm_c(epos);
  return rs;
}

/*
  Geopotential

  Parameteres:
     1: mu: Gravitational parameter of the central planet (GM)
*/
double geoPotential(double et,double r,double q,double f,void *params)
{
  double *ps=(double*)params;
  double mu=ps[1];
  double phi,phic;

  //SPHERICAL COMPONENT
  phi=-mu/r;
  VPRINT("Spherical potential: phi = %.17e\n",phi);
  
  double cosm,Pnm,Ror,Rorn,muor,u;
  int inm;
  //COMPUTE HARMONIC CORRECTION
  Ror=1.0/r;
  muor=mu/r;
  u=sin(f);
  phic=0.0;  
  VPRINT("Coordinates:\n\tr = %e, f = %e, q = %e\n\tmu/r = %e, sin(f) = %e\n",r,f,q,muor,u);
  VPRINT("Components:\n");
  for(int n=MAXN+1;n-->0;){
    VPRINT("\tn = %d\n",n);
    for(int m=0;m<=n;m++){
      VPRINT("\t\tm = %d\n",m);
      if(J[n][m]==0) continue;
      cosm=cos(m*(q-Q[n][m]));
      Pnm=gsl_sf_legendre_Plm(n,m,u);
      Rorn=gsl_pow_int(Ror,n);
      VPRINT("\t\t\tJnm = %e, Pnm = %e, (R/r)^n = %e, Qnm = %e, cosm = %e\n",
	     J[n][m],Pnm,Rorn,Q[n][m],cosm);
      phic+=muor*J[n][m]*Rorn*Pnm*cosm;
    }
  }
  phi+=phic;

  return phi;
}

/*
  Gradient of geopotential in spherical coordinates

  Parameteres:
     1: mu: Gravitational parameter of the central planet (GM)
 */
int gradGeoPotential(double et,double r,double f,double q,double dphidy[],void *params)
{
  double *ps=(double*)params;
  double mu=ps[1];

  //SPHERICAL
  /*dphi/dr*/dphidy[0]=mu/(r*r);
  /*dphi/df*/dphidy[1]=0;
  /*dphi/dq*/dphidy[2]=0;
  VPRINT("Spherical force: %s\n",vec2strn(dphidy,3,"%.17e "));

  //HARMONIC EXPANSION
  double dphidyc[]={0.0,0.0,0.0};
  double cosm,sinm,Pnm,Pnmp,Ror,Rorn,muor,norm;
  double u=sin(f);
  int inm;
  
  //ASSOCIATED LEGENDRE POLYNOMIALS
  int nmax=gsl_sf_legendre_array_n(MAXN);
  double aPnm[nmax+1],aPnmp[nmax+1];
  gsl_sf_legendre_deriv_array((gsl_sf_legendre_t)1.0,MAXN,u,aPnm,aPnmp);

  //COMPUTE HARMONIC CORRECTION
  Ror=1.0/r;
  muor=mu/r;
  VPRINT("Coordinates:\n\tr = %e, f = %e, q = %e\n\tmu/r = %e, sin(f) = %e\n",r,f,q,muor,u);
  VPRINT("Components:\n");
  for(int n=MAXN+1;n-->0;){
    VPRINT("\tn = %d\n",n);
    for(int m=0;m<=n;m++){
      if(J[n][m]==0) continue;
      VPRINT("\t\tm = %d\n",m);

      inm=gsl_sf_legendre_array_index(n,m);
      cosm=cos(m*(q-Q[n][m]));
      sinm=sin(m*(q-Q[n][m]));
      Pnm=gsl_sf_legendre_Plm(n,m,u);
      norm=aPnm[inm]/Pnm;
      Pnmp=aPnmp[inm]/norm;
      Rorn=gsl_pow_int(Ror,n);

      VPRINT("\t\t\tJnm = %e, Pnm = %e, P'nm = %e, (R/r)^n = %e, Qnm = %e, cosm = %e, sinm = %e, \n",
	     J[n][m],Pnm,Pnmp,Rorn,Q[n][m],cosm,sinm);

      /*dphi/dr*/dphidyc[0]+=-muor/r*(n+1)*J[n][m]*Rorn*Pnm*cosm;
      /*dphi/df*/dphidyc[1]+=+muor*J[n][m]*Rorn*cos(f)*Pnmp*cosm;
      /*dphi/dq*/dphidyc[2]+=-muor*m*J[n][m]*Rorn*Pnm*sinm;
    }
  }
  VPRINT("Harmonic correction: %s\n",vec2strn(dphidyc,3,"%.17e "));
  PAUSE();
  
  dphidy[0]+=dphidyc[0];
  dphidy[1]+=dphidyc[1];
  dphidy[2]+=dphidyc[2];

  return 0;
}

double totalEnergy(double et,double x[],void *params)
{
  double r=x[0],f=x[1],q=x[2];
  double dr=x[3],df=x[4],dq=x[5];
  double K=0.5*(dr*dr+r*r*df*df+r*r*cos(f)*cos(f)*dq*dq);
  double V=geoPotential(et,r,q,f,params);
  double E=K+V;
  return E;
}

/*
  Tidal effects

  Parameteres:
     1: mu: Gravitational parameter of the central planet (GM)
 */
int tidalForce(double et,double r,double f,double q,double Ftid[],void *params)
{
  double *ps=(double*)params;
  
  return 0;
}

/*
  Atmospheric friction

  Parameteres:
     1: mu: Gravitational parameter of the central planet (GM)


  NOTE: Remember to include the rotation of Earth when calculating the
  relative velocity:

  v_e = r w_e cos b [ (sin f cos q - cos f cos i sin q) a_r - 
                     (cos f cos i cos q + sin f sin q) a_q +
		     (cos f sin i) a_f ]

  See: Wu (1965)

  Where:
  f: 
  
 */
int atmosDrag(double et,double r,double q,double f,double Fdrag[],void *params)
{
  double *ps=(double*)params;
  
  return 0;
}

/*
  Solar pressure

  Parameteres:
     1: mu: Gravitational parameter of the central planet (GM)
 */
int solarPressure(double r,double q,double f,double Ftid[],void *params)
{
  double *ps=(double*)params;
  
  return 0;
}

////////////////////////////////////////////////////////////////////////
//OTHER ROUTINES
////////////////////////////////////////////////////////////////////////
/*
  This is a prototype routine for integrate

  MAS:
    dx/dt = v
    dv/dt = -w^2 x
    
 */
int EoM_MAS(double t,double y[],double dydt[],void *params) 
{ 
  double *ps=(double*)params;
  int nsys=(int)ps[0];
  
  //PHYSICAL PARAMETERS
  double w=ps[1];

  //EQUATIONS
  dydt[0]=y[1];
  dydt[1]=-w*w*y[0];

  return 0;
}

/*
  This is a prototype routine for integrate

  2B:
    dr/dt = v
    dv/dt = -GM/r^3 r

  r,v vectors
    
 */
int EoM_2B(double t,double y[],double dydt[],void *params) 
{ 
  double *ps=(double*)params;
  int nsys=(int)ps[0];
  
  //PHYSICAL PARAMETERS
  double G=ps[1];
  double M=ps[2];

  //EQUATIONS
  dydt[0]=y[3];
  dydt[1]=y[4];
  dydt[2]=y[5];
  
  double r=sqrt(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]);
  dydt[3]=-G*M/(r*r*r)*y[0];
  dydt[4]=-G*M/(r*r*r)*y[1];
  dydt[5]=-G*M/(r*r*r)*y[2];

  /*
  VPRINT("%e %e %e %e %e %e\n%e %e %e %e %e %e\n",
	  y[0],y[1],y[2],y[3],y[4],y[5],
	  dydt[0],dydt[1],dydt[2],dydt[3],dydt[4],dydt[5]);
  getchar();
  */

  return 0;
}

/*
  Satellite equations in spherical coordinates
  
  Coordinates:
  y0: r
  y1: theta (q)
  y2: phi (f)
  
  Parameters:

     1: mu: Gravitational parameter of the central planet (GM)

  NOTES:

  - Equations come from Praire (2008).  However a slight error exists
    in his expressions in the equation for phi (the term proportional
    to sin phi cos phi have a factor dphi^2 that in fact is dtheta^2
    (see http://mathworld.wolfram.com/SphericalCoordinates.html)

 */
int EoM_Full(double t,double y[],double dydt[],void *params) 
{ 
  double *ps=(double*)params;
  int nsys=(int)ps[0];

  //PHYSICAL PARAMETERS
  double mu=ps[1];

  //VARIABLES
  double r=y[0],f=y[1],q=y[2],dr=y[3],df=y[4],dq=y[5];
  double dphidy[3];
  double cosf=cos(f),cosf2=cosf*cosf,sinf=sin(f),tanf=sinf/cosf,sinq=sin(q);

  //EQUATIONS
  /*dr/dt*/dydt[0]=dr;
  /*df/dt*/dydt[1]=df;
  /*dq/dt*/dydt[2]=dq;
  
  //GRADIENT OF POTENTIAL
  gradGeoPotential(t,r,f,q,dphidy,params);

  //ACCELERATION
  /*d(dr/dt)/dt*/dydt[3]=-dphidy[0]+r*dq*dq*cosf2+r*df*df;
  /*d(df/dt)/dt*/dydt[4]=-dphidy[1]/(r*r)-2*dr*df/r-dq*dq*sinf*cosf;
  /*d(dq/dt)/dt*/dydt[5]=-dphidy[2]/(r*r*cosf2)-2*dr*dq/r+2*dq*df*tanf;

  /*
  VPRINT("grad Potential = %e %e %e\n",
	  dphidy[0],dphidy[1],dphidy[2]);
  VPRINT("r = %e, f = %e, q = %e, dr/dt = %e, df/dt = %e, dq/dt = %e\n",y[0],y[1],y[2],y[3],y[4],y[5]);
  VPRINT("dr/dt = %e, df/dt = %e, dq/dt = %e, d2r/dt2 = %e, d2f/dt2 = %e, d2q/dt2 = %e\n",dydt[0],dydt[1],dydt[2],dydt[3],dydt[4],dydt[5]);
  getchar();
  //*/

  return 0;
}

int EoM_Fulla(double t,double y[],double dydt[],void *params) 
{ 
  double *ps=(double*)params;
  int nsys=(int)ps[0];

  ////////////////////////////////////////////////////
  //PHYSICAL PARAMETERS
  ////////////////////////////////////////////////////
  double mu=ps[1];

  ////////////////////////////////////////////////////
  //GET COORDINATES
  ////////////////////////////////////////////////////
  double r=y[0],f=y[1],q=y[2],dr=y[3],df=y[4],dq=y[5];
  VPRINT("Distance to corresponding earth's surface: r = %.17e, rs = %.17e\n",r,earthRadius(q,f));

  if(r<earthRadius(q,f)){
    fprintf(stderr,"Satellite has collided with the central body\n");
    throw(1);
  }
  
  double dphidy[3];
  double cosf=cos(f),cosf2=cosf*cosf,sinf=sin(f),tanf=sinf/cosf,sinq=sin(q);
  //VPRINT("r = %+e, f = %+e, q = %+e\n",r,f,q);

  ////////////////////////////////////////////////////
  //ACCELERATIONS
  ////////////////////////////////////////////////////
  double a[]={0,0,0};

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //GRAVITATIONAL POTENTIAL
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  gradGeoPotential(t,r,f,q,dphidy,params);
  double g[]={-dphidy[0],
	      -dphidy[1]/(r*r),
	      -dphidy[2]/(r*r*cosf2)};
  vadd_c(g,a,a);
  //VPRINT("ar = %+e, af = %+e, aq = %+e\n",g[0],g[1],g[2]);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //N-BODY
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  /*TEST
  //CONVERT TO CARTESIAN
  double rc[6];
  sph2car(rc,y);
  //VPRINT("x = %+e, y = %+e, z = %+e\n",rc[0],rc[1],rc[2]);
  double gc[]={
    -mu/(r*r*r)*rc[0],
    -mu/(r*r*r)*rc[1],
    -mu/(r*r*r)*rc[2]
  };
  //VPRINT("ax = %+e, ay = %+e, az = %+e\n",gc[0],gc[1],gc[2]);
  //VPRINT("|ac|=%e\n",vnorm_c(gc));
  double gs[3];
  //CONVERT TO SPHERICAL
  vec_car2sph(gs,gc,y);
  vadd_c(gs,a,a);
  //VPRINT("ar = %+e, af = %+e, aq = %+e\n",gs[0],gs[1],gs[2]);
  //VPRINT("|as|=%e\n",vnorm_c(gs));
  //getchar();
  //*/

  ////////////////////////////////////////////////////
  //EQUATIONS
  ////////////////////////////////////////////////////
  /*dr/dt*/dydt[0]=dr;
  /*df/dt*/dydt[1]=df;
  /*dq/dt*/dydt[2]=dq;
  /*d(dr/dt)/dt*/dydt[3]=a[0]+r*dq*dq*cosf2+r*df*df;
  /*d(df/dt)/dt*/dydt[4]=a[1]-2*dr*df/r-dq*dq*sinf*cosf;
  /*d(dq/dt)/dt*/dydt[5]=a[2]-2*dr*dq/r+2*dq*df*tanf;

  /*
  VPRINT("grad Potential = %e %e %e\n",
	  dphidy[0],dphidy[1],dphidy[2]);
  VPRINT("r = %e, f = %e, q = %e, dr/dt = %e, df/dt = %e, dq/dt = %e\n",y[0],y[1],y[2],y[3],y[4],y[5]);
  VPRINT("dr/dt = %e, df/dt = %e, dq/dt = %e, d2r/dt2 = %e, d2f/dt2 = %e, d2q/dt2 = %e\n",dydt[0],dydt[1],dydt[2],dydt[3],dydt[4],dydt[5]);
  getchar();
  //*/

  return 0;
}
