////////////////////////////////////////////////////////////////////////
//BASIC LIBRARIES
////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////
//ADDITIONAL LIBRARIES
////////////////////////////////////////////////////////////////////////
#include <nrlmsise-00.h>

////////////////////////////////////////////////////////////////////////
//BASIC MACROS
////////////////////////////////////////////////////////////////////////
#define EXTMET 1
#define TOLERANCE 1E-10
#define ATTEMPTS 12 /*SEE NUMBER_OF_STEPS*/
static int NUMBER_OF_STEPS[]={2,4,6,8,12,16,24,32,48,64,96,128};

////////////////////////////////////////////////////////////////////////
//PHYSICAL CONSTANTS
////////////////////////////////////////////////////////////////////////
#define GCONST 6.67E-11

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
  for(int i=0;i<npoints;i++) {
    ts[i]=t;
    
    /*
    fprintf(stdout,"t=%e, y=%e %e %e %e %e %e\n",
	    t,x0[0],x0[1],x0[2],x0[3],x0[4],x0[5]);
    getchar();
    //*/

    copyVec(X[i],x0,nsys);
    deltat=t-tini;
    if(direction*((t_start+t_step)-tend)>0) t_step=(tend-t_start);
    t_stop=t_start+t_step;
    h_used=h;

    do {
      while(1){
	status=Gragg_Bulirsch_Stoer(eom,x0,x,t,h_used,&h_next,1.0,
				    TOLERANCE,params);

	/*
	fprintf(stdout,"t=%e, y=%e %e %e %e %e %e, hnext = %e, %d\n",
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
//ATMOSPHERIC MODEL MSIS90E
////////////////////////////////////////////////////////////////////////
/*
  Source: ftp://ftp.irf.se/pub/perm/ESRAD/SPECTRUM/gtd6.f (JGR, 96, 1159-1172, 1991)
 */
typedef float real;
real *newVectorf(int n)
{
  real *v;
  v=(real*)malloc(n*sizeof(real));
  return v;
}
extern "C" void gtd6_(int* IYD,real* SEC,real* ALT,real* GLAT,real* GLONG,
		      real* STL,real* F107A,real* F107,real AP[],int* MASS,
		      real D[],real T[]);
/*
       IYD - YEAR AND DAY AS YYDDD or DDD (day of year from 1 to 365)
       SEC - UT(SEC)
       ALT - ALTITUDE(KM)
       GLAT - GEODETIC LATITUDE(DEG)
       GLONG - GEODETIC LONGITUDE(DEG)
       STL - LOCAL APPARENT SOLAR TIME(HRS)
       AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
          - ARRAY CONTAINING:
	    (1) F107A - 3 MONTH AVERAGE OF F10.7 FLUX
	    (2) F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
            (3) DAILY AP
            (4) 3 HR AP INDEX FOR CURRENT TIME
            (5) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
            (6) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
            (7) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
            (8) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
                TO CURRENT TIME
            (9) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 59 HRS PRIOR
                TO CURRENT TIME
       MASS - MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS
                CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL.

    Note:  Ut, Local Time, and Longitude are used independently in the
           model and are not of equal importance for every situation.  
           For the most physically realistic calculation these three
           variables should be consistent (STL=SEC/3600+GLONG/15).
           F107, F107A, and AP effects are not large below 80 km 
           and these can be set to 150., 150., and 4. respectively.

    OUTPUT:
       D(1) - HE NUMBER DENSITY(CM-3)
       D(2) - O NUMBER DENSITY(CM-3)
       D(3) - N2 NUMBER DENSITY(CM-3)
       D(4) - O2 NUMBER DENSITY(CM-3)
       D(5) - AR NUMBER DENSITY(CM-3)                       
       D(6) - TOTAL MASS DENSITY(GM/CM3)
       D(7) - H NUMBER DENSITY(CM-3)
       D(8) - N NUMBER DENSITY(CM-3)
       T(1) - EXOSPHERIC TEMPERATURE
       T(2) - TEMPERATURE AT ALT
 */
/*
int MSIS90E(int day,double sec,
	    double alt,double glat,double glon,
	    double KAP[],
	    double *rho,double *T)
{
  int IYD=day;
  real SEC=(real)sec;
  real ALT=(real)alt;
  real GLAT=(real)glat;
  real GLONG=(real)glon;
  real STL=sec/3600.0+GLONG/15.0;
  real F107A=KAP[0];
  real F107=KAP[1];
  int MASS=48;
  real* AP=newVectorf(7);
  int k=0,j=2;
  AP[k++]=KAP[j++];AP[k++]=KAP[j++];AP[k++]=KAP[j++];
  AP[k++]=KAP[j++];AP[k++]=KAP[j++];AP[k++]=KAP[j++];
  AP[k++]=KAP[j++];
  real* D=newVectorf(8);
  real* Ts=newVectorf(2);
  gtd6_(&IYD,&SEC,&ALT,&GLAT,&GLONG,&STL,&F107A,&F107,AP,&MASS,D,Ts);
  *rho=(double)D[5]*1e3;
  *T=(double)Ts[1];
  return 0;
}
*/

#define KAP_REF {150,22,22,22,22,22,22,22}
#define FLAGS_ALL {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}
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
  fprintf(stdout,"%e %e %e %e %e %e\n%e %e %e %e %e %e\n",
	  y[0],y[1],y[2],y[3],y[4],y[5],
	  dydt[0],dydt[1],dydt[2],dydt[3],dydt[4],dydt[5]);
  getchar();
  */

  return 0;
}

