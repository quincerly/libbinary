#include <stdio.h>
#include <math.h>

#include "binary.h"
#include "misc.h"

#define G 6.672*1e-11
#define Msun 1.9891*1e30
#define TOL 5.e-12
#define PI 3.1415926536

/* Variables needed by wrapper functions */
static double q_roche_xl1;
static double find_potl_xp,find_potl_yp,find_potl_zp;
static double find_potl_dx,find_potl_dy,find_potl_dz;
static double find_potl_potl,find_potl_q;

/* Protos for wrapper functions */
double roche_xl1_func(double x);
double find_potl_func(double l);

/* Calculate orbital separation */
double binary_sep(double q, double M1, double P)
{
  return pow(G*M1*Msun*(1.+q)*P*P/4./PI/PI,1./3.);
}

/* Calculate Eggleton radius of primary, RL1 */
double eggleton_radius(double q)
{
  return 0.49*pow(1./q,2./3.)/(0.6*pow(1./q,2./3.)+log(1+pow(1./q,1./3.)));
}

/* Find position of L1 point (xl1,0,0) */
double roche_xl1_func(double x)
{
  return -roche_potl(x,0.,0.,q_roche_xl1);
}
double roche_xl1(double q)
{
  
  double xa,xb,xc;
  double xp1=-q/(1.+q);
  double xp2=1./(1.+q);
  double potl1,xl1;
  
  xa=xp1+1.e-6;
  xc=xp2-1.e-6;
  xb=xa+0.38197*(xc-xa);
  q_roche_xl1=q;
  potl1=brentmin(xa,xb,xc,roche_xl1_func,TOL,&xl1);
  
  return xl1;
  
}

/* Calculate Roche potl at position [x,y,z] in J/kg / (G*M1*Msun*a) */
/* x,y and z in units of a */
double roche_potl(double x, double y, double z, double q)
{

  double r1,r2,rsqr,potl;
  double xp1=-q/(1.+q);
  double xp2=1./(1.+q);
  
  r1=sqrt((x-xp1)*(x-xp1)+y*y+z*z);
  r2=sqrt((x-xp2)*(x-xp2)+y*y+z*z);
  rsqr=(x*x+y*y);

  potl=(-1./r1-q/r2-rsqr*(1.+q)/2.);

  return potl;

}

// Find where roche_potl/(G*M1*Msun*a)=pot between P and Q
double find_potl_func(double l)
{
  double diff;

  //  printf("%g\n",l);

  diff=roche_potl(find_potl_xp+l*find_potl_dx,
		  find_potl_yp+l*find_potl_dy,
		  find_potl_zp+l*find_potl_dz,
		  find_potl_q)-find_potl_potl;

  return diff*diff;
}
double find_potl(double pot, double q, double xp, double yp, double zp, double xq, double yq, double zq)
{

  double la,lb,lc,r,dx,dy,dz;
  double dummy,lmin;

  r=sqrt((xp-xq)*(xp-xq)+(yp-yq)*(yp-yq)+(zp-zq)*(zp-zq));
  dx=(xq-xp)/r;
  dy=(yq-yp)/r;
  dz=(zq-zp)/r;
  la=0.;
  lc=r;
  lb=0.5*(la+lc);

  find_potl_xp=xp;
  find_potl_yp=yp;
  find_potl_zp=zp;
  find_potl_dx=dx;
  find_potl_dy=dy;
  find_potl_dz=dz;
  find_potl_q=q;
  find_potl_potl=pot;

  dummy=brentmin(la,lb,lc,find_potl_func,TOL,&lmin);
  
  return lmin/r;
  
}

/* Keplerian velocity vector vk[0,1,2] at position [x,y,z] */
void vkep(double *vk,double x,double y,double z,double M1)
{
  double rmag,r[3],omega_kep[3];
  
  r[0]=x;
  r[1]=y;
  r[2]=z;

  rmag=sqrt(x*x+y*y+z*z);

  omega_kep[0]=0.;
  omega_kep[1]=0.;
  omega_kep[2]=sqrt(G*M1*Msun/rmag/rmag/rmag);
  cross_vector(vk,omega_kep,r);
}

void omega_kep(double *omegak,double x,double y,double z,double M1)
{
  double rmag,r[3];
  
  r[0]=x;
  r[1]=y;
  r[2]=z;

  rmag=sqrt(x*x+y*y+z*z);

  omegak[0]=0.;
  omegak[1]=0.;
  omegak[2]=sqrt(G*M1*Msun/rmag/rmag/rmag);

}
