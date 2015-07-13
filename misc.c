#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "misc.h"

#define mu0 1.256637061*1e-6
#define PI 3.1415926536
#define ITMAX 1000
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

/* Return vector cross product of v and x to result */
void cross_vector(double *result, double *v, double *x)
{
  
  double cross[3];
  
  cross[0]=v[1]*x[2]-x[1]*v[2];
  cross[1]=x[0]*v[2]-v[0]*x[2];
  cross[2]=v[0]*x[1]-x[0]*v[1];

  result[0]=cross[0];
  result[1]=cross[1];
  result[2]=cross[2];
  
}

/* Return unit vector of x to result */
void unit_vector(double *result, double *x)
{
  
  double mag;
  
  mag=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);

  result[0]=x[0]/mag;
  result[1]=x[1]/mag;
  result[2]=x[2]/mag;
  
}

/* Return magnetic field B at r due to dipole [0,0,m] at [0,0,0]. */
/* If m=0 then just return direction of B */
void b_dipole(double *B, double *r, double m)
{

  double rxy,fac,theta;

  rxy=sqrt(r[0]*r[0]+r[1]*r[1]);

  theta=angle(r[2],rxy);
  B[0]=3.*sin(theta)*cos(theta)*r[0]/rxy;
  B[1]=3.*sin(theta)*cos(theta)*r[1]/rxy;
  B[2]=3.*cos(theta)*cos(theta)-1.;
  if (m!=0.) fac=mu0*m/4./PI/pow(r[0]*r[0]+r[1]*r[1]+r[2]*r[2],1.5);
  else fac=1./sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);

  B[0]*=fac;
  B[1]*=fac;
  B[2]*=fac;

}

/* For angle in triangle with opposite o and adjacent a */
/* this function returns the angle in the correct quadrant */
double angle(double a, double o)
{

  double theta;

  if (a!=0.) theta=atan(o/a);
  else {
    if (o>0.) theta=PI/2.;
    if (o<0.) theta=PI/-2.;
  }

if (a<0.) theta=PI+theta;

 if ((a==0.)&&(o==0.)) {
   printf("Angle undefined.\n");
   exit(1);
 }

 return theta;

}

/* Brent minimization from Numerical recipes */
double brentmin(double ax, double bx, double cx,
	     double (*f)(double), double tol, double *xmin)
{
	int iter;
	double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double e=0.0;

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=(*f)(x);
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*xmin=x;
			return fx;
		}
		if (fabs(e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=(*f)(u);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		}
	}
	printf("\nToo many iterations in brentmin\n");
	*xmin=x;
	return fx;
}

// Fold phase onto range 0-1
double fold(double phase)
{
  
  while (phase<0.) phase=phase+1.;
  
  return(phase-(double)((long int)phase));
}

// Return inverse and determinant of a 3*3 matrix
// Matrix input/output format is double precision
// 9 element matrix [a,b,c,d,e,f,g,h,i] representing
//
//     [ a b c ]
//     [ d e f ]
//     [ g h i ]
//
// Return value is determinant
double inverse_matrix(double *input, double *inverse)
{

  double a,b,c,d,e,f,g,h,i;
  double det,coa,cob,coc,cod,coe,cof,cog,coh,coi;

  a=input[0];
  b=input[1];
  c=input[2];
  d=input[3];
  e=input[4];
  f=input[5];
  g=input[6];
  h=input[7];
  i=input[8];

  coa=(e*i-h*f);
  cob=-(d*i-g*f);
  coc=(d*h-g*e);
  cod=-(b*i-h*c);
  coe=(a*i-g*c);
  cof=-(a*h-g*b);
  cog=(b*f-e*c);
  coh=-(a*f-d*c);
  coi=(a*e-d*b);

  det=a*coa+b*cob+c*coc;

  if (det!=0.) {
    inverse[0]=coa/det;
    inverse[1]=cod/det;
    inverse[2]=cog/det;
    inverse[3]=cob/det;
    inverse[4]=coe/det;
    inverse[5]=coh/det;
    inverse[6]=coc/det;
    inverse[7]=cof/det;
    inverse[8]=coi/det;
  }
  else {
      inverse[0]=0.;
      inverse[1]=0.;
      inverse[2]=0.;
      inverse[3]=0.;
      inverse[4]=0.;
      inverse[5]=0.;
      inverse[6]=0.;
      inverse[7]=0.;
      inverse[8]=0.;
  }

  return(det);
}

double min(double *data, long int ndata)
{

  long int i;
  double mindata;

  mindata=data[0];
  for (i=1;i<ndata;i++) if (data[i]<mindata) mindata=data[i];

  return(mindata);
}

double max(double *data, long int ndata)
{

  long int i;
  double maxdata;

  maxdata=data[0];
  for (i=1;i<ndata;i++) if (data[i]>maxdata) maxdata=data[i];

  return(maxdata);
}

double gaussprofile(double x, double fwhm)
  // Normalized gaussian profile
  // Peak value is 2 * sqrt(ln(2)) / fwhm / sqrt(pi)
{
 
  double ewidth,normfac;
  //ewidth= 0.5 * fwhm / sqrt(ln(2)) = 0.6005612043932249 fwhm
  
  ewidth=0.6005612043932249*fwhm;
  normfac=0.5641895835477563/ewidth;

  return(normfac*exp(-(x*x/(ewidth*ewidth))));
}

double gaussprofileamp(double x, double fwhm, double A)
  // Normalized gaussian profile
  // Peak value is A
{
 
  double ewidth;
  //ewidth= 0.5 * fwhm / sqrt(ln(2)) = 0.6005612043932249 fwhm
  
  ewidth=0.6005612043932249*fwhm;

  return(A*exp(-(x*x/(ewidth*ewidth))));
}
  
double lorentzprofile(double x, double fwhm)
  // Normalized Lorentzian profile
  // Peak value is 2 / fwhm
{

  return(2.*fwhm/PI/(fwhm*fwhm+4.*x*x));
}

double lorentzprofileamp(double x, double fwhm, double A)
     // Normalized Lorentzian profile
     // Peak value is A
{

  return(A*fwhm*fwhm/(fwhm*fwhm+4.*x*x));
}

double randomn(long *idum)
     // gasdev() from Numerical Recipes
{
	double randomu(long *idum);
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;

	if (*idum < 0) iset=0;
	if  (iset == 0) {
		do {
			v1=2.0*randomu(idum)-1.0;
			v2=2.0*randomu(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}

double randomu(long *idum)
     // ran1() from Numerical Recipes
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

long int randomseed()
     // Calculate seed from current time
{

  time_t currenttime;

  currenttime=time(NULL);

  return -(long int)currenttime;
}
