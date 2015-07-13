#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "stream.h"
#include "binary.h"
#include "misc.h"

#define G 6.672*1e-11
#define Msun 1.9891*1e30
#define c 2.99792e+08
#define PI 3.1415926536
#define NSUBDEF 500
#define Bphidisc 0.
#define Brdisc 0.
#define Bmagdisc sqrt(1.+Bphidisc*Bphidisc+Brdisc*Brdisc)

#ifdef __EXEC
#include <cpgplot.h>
/* Redundant main() function used when testing */
void main()
{

  double q=0.133333333;
  double M1=0.45;
  double P=82.*60.;
  double Pb=27.87;
  double k0=0.8e-5;
  double t=1.;
  long n=500;

  double x[n],y[n],z[n];
  double vix[n],viy[n],viz[n];
  double vkx[n],vky[n],vkz[n];
  double d[n];
  double a;
  long i;
  int plotdev;
  float xpl[n],ypl[n];

  a=binary_sep(q,M1,P);
  stream_calc(q,M1,P,Pb,k0,t,n,x,y,z,vix,viy,viz,vkx,vky,vkz,d);

  plotdev=cpgopen("/xw");
  if (plotdev<=0)
    {
      printf("** Couldn't open PGPLOT device.**\n\n");
      return;
    }

  cpgenv(-2.8,2.8,-2.,2.,0,1);
  cpglab("x/a","y/a","Stream trajectory");

  for (i=0; i<n; i++) 
    {
      xpl[i]=(float)(x[i]/a);
      ypl[i]=(float)(y[i]/a);
    }       

  cpgline(n,xpl,ypl);
  cpgend();

  exit(0);
}
#endif

/* Calculate acceleration of particle at r with velocity v */
void stream_accel(double *accel, double *vin, double *r, double *v, double dist,
		  double *k0, double kindex, double M1, double q, double rdisc, double Pb, double discbrfac,
	   double xp1, double xp2, double *omega, double *omegab,
	   double r0, double magdist, double *magdb, double *lmag)
{

  double s1,s2,s,s1xy;
  double vec1[3],vec2[3],vec[3];
  double vi[3],vxy[3],rxy[3],vdiff[3],r1[3];
  double coriolis[3],centrifugal[3],gravity[3],magnetic[3];
  double dirmod,direction[3];
  double k,vb[3],B[3],omegabdisc[3],vdiffproj;

  /* Distance of r from p1,p2 and c.o.m. */
  s1=sqrt((r[0]-xp1)*(r[0]-xp1)+r[1]*r[1]+r[2]*r[2]);
  s2=sqrt((r[0]-xp2)*(r[0]-xp2)+r[1]*r[1]+r[2]*r[2]);
  s=sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
  s1xy=sqrt((r[0]-xp1)*(r[0]-xp1)+r[1]*r[1]);

  /* Unit vectors in directions towards M1 and M2 and c.o.m. */
  vec1[0]=(xp1-r[0])/s1;
  vec1[1]=-r[1]/s1;
  vec1[2]=-r[2]/s1;
  vec2[0]=(xp2-r[0])/s2;
  vec2[1]=-r[1]/s2;
  vec2[2]=-r[2]/s2;
  vec[0]=r[0]/s;
  vec[1]=r[1]/s;
  vec[2]=r[2]/s;

  /* Velocity in instantaneous inertial frame */
  cross_vector(vi,omega,r);
  vi[0]+=v[0];
  vi[1]+=v[1];
  vi[2]+=v[2];
  vin[0]=vi[0];
  vin[1]=vi[1];
  vin[2]=vi[2];

  /* Calculate coriolis acceleration */
  vxy[0]=v[0];
  vxy[1]=v[1];
  vxy[2]=0.;
  cross_vector(coriolis,omega,vxy);
  coriolis[0]*=-2.;
  coriolis[1]*=-2.;
  coriolis[2]*=-2.;

  /* Calculate centrifugal acceleration */
  rxy[0]=r[0];
  rxy[1]=r[1];
  rxy[2]=0.;
  cross_vector(centrifugal,omega,rxy);
  cross_vector(centrifugal,omega,centrifugal);
  centrifugal[0]*=-1.;
  centrifugal[1]*=-1.;
  centrifugal[2]*=-1.;

  /* Calculate acceleration due to gravitation force of both stars */
  gravity[0]=(G*M1*Msun/s1/s1)*vec1[0]+(G*q*M1*Msun/s2/s2)*vec2[0];
  gravity[1]=(G*M1*Msun/s1/s1)*vec1[1]+(G*q*M1*Msun/s2/s2)*vec2[1];
  gravity[2]=(G*M1*Msun/s1/s1)*vec1[2]+(G*q*M1*Msun/s2/s2)*vec2[2];

  /* Add together all these contributions */
  accel[0]=gravity[0]+centrifugal[0]+coriolis[0];
  accel[1]=gravity[1]+centrifugal[1]+coriolis[1];
  accel[2]=gravity[2]+centrifugal[2]+coriolis[2];

  /* Include effect of magnetic drag. Magnetic field is */
  /* turned on smoothly from rdisc to 0.1*a along trajectory.*/
  /* See Wynn, King and Horne, MNRAS 286, 436-446. */
  *lmag=0.;
  if ((*k0<0.)&&(s1<rdisc)) {
    *k0*=-1.;
    *magdb=dist;
  }
  if (*k0>0.) {
      
    r1[0]=r[0]-xp1;
    r1[1]=r[1];
    r1[2]=r[2];
        
    /* k, Velocity and direction of magnetic field at r */
    if (Pb!=0.) {
      // For rotating dipole at primary
      cross_vector(vb,omegab,r1);
      b_dipole(B,r1,0.); // Get direction of B at r1
    }
    else {
      // Disc anchored field
      omega_kep(omegabdisc,r1[0]*discbrfac,r1[1]*discbrfac,r1[2]*discbrfac,M1);
      cross_vector(vb,omegabdisc,r1);
      B[0]=(-Bphidisc*r[1]/s1xy+Brdisc*(r[0]-xp1)/s1xy)/Bmagdisc;
      B[1]=(Bphidisc*(r[0]-xp1)/s1xy+Brdisc*r[1]/s1xy)/Bmagdisc;
      B[2]=1./Bmagdisc;
      //b_dipole(B,r1,0.); // Get direction of B at r1
    }	

    /* Switch off magnetic effects outside light cylinder */
    if (omegab[2]*s1>c) k=0.; else k=*k0*pow(s1/r0,kindex);

    /* Outside disc when using disc propeller */
    if ((Pb==0.)&&(s1>rdisc)) k=0.;
    
    /* Smooth transition between mag drag within disc and no drag at l1 */
    if (dist-*magdb<magdist) k=k*(dist-*magdb)*(dist-*magdb)/magdist/magdist;

    /* Difference between field velocity and velocity of particle in inertial frame */
    vdiff[0]=vi[0]-vb[0];
    vdiff[1]=vi[1]-vb[1];
    vdiff[2]=vi[2]-vb[2];

    /* Define magnetic drag as -k*(component of v-vb perpendicular to B) */
    cross_vector(direction,B,vdiff);
    cross_vector(direction,direction,B);
    dirmod=sqrt(direction[0]*direction[0]+direction[1]*direction[1]+direction[2]*direction[2]);
    if (dirmod!=0.) {
      direction[0]=direction[0]/dirmod;
      direction[1]=direction[1]/dirmod;
      direction[2]=direction[2]/dirmod;
    }
    vdiffproj=direction[0]*vdiff[0]+direction[1]*vdiff[1]+direction[2]*vdiff[2];
    /* vdiff resolved along direction */
    magnetic[0]=-k*direction[0]*vdiffproj; 
    magnetic[1]=-k*direction[1]*vdiffproj; 
    magnetic[2]=-k*direction[2]*vdiffproj; 
    
    /* Magnetic dissipation rate */
    *lmag=-k*vdiffproj*vdiffproj;

    /* Much quicker. 2-D only */
    /*magnetic[0]=-k*vdiff[0];
    magnetic[1]=-k*vdiff[1];
    magnetic[2]=-k*vdiff[2];*/

    /* Add magnetic drag to the acceleration */
    accel[0]+=magnetic[0];
    accel[1]+=magnetic[1];
    accel[2]+=magnetic[2];
    
  }
  
}

/* Calculate path of accretion stream in a close binary. */
/* calculation done in corotating frame. */
void stream_calc(double q, double M1, double P,
		 double Pb, double kk0, double kindex, double prdisc, double discbrfac,
		 double t, long n, long NSUB,
		 double *x, double *y, double *z,
		 double *vix, double *viy, double *viz,
		 double *vkx, double *vky, double *vkz,
		 double *d, double *lmag)
{

  double magdist; /* Distance over which magnetic field is turned on */
  double magdb; /* Distance travelled when B is turned on */
  double xp1,xp2; /* x positions of primary and secondary */
  double omegab[3]; /* Angular velocity of magnetic field */
  double omega[3]; /* Orbital angular velocity vector */
  double a; /* The orbital separation */
  double xl1; /* x position of L1 */
  double rl1; /* Distance from primary to L1 */
  double regg; /* The eggleton radius */
  double rdisc; /* Distance from primary at which B is turned on */
  double r0; /* Distance from primary corresponding to k0 */
  //double sx,sy,sz; /* Initial position of stream */
  //double svx,svy,svz; /* Initial velocity of stream */
  double dt; /* Time step */
  double k0,dlmag;
  double *vx,*vy,*vz; /* Velocity of stream in corotating frame */
  double acc[3],pos[3],vel[3],vin[3],vkep[3],dr[3],rxy[3],omegakep[3],dd;
  long i,j;

  vx=(double *)malloc(n*sizeof(double));
  vy=(double *)malloc(n*sizeof(double));
  vz=(double *)malloc(n*sizeof(double));

  /***************** INITIALIZE VARIABLES *************************/

  if (NSUB==0) NSUB=NSUBDEF;

  /* Calculate angular frequency of orbit */
  omega[0]=0.;
  omega[1]=0.;
  omega[2]=2.*PI/P;

  /* Calculate angular velocity of magnetic field */
  omegab[0]=0.;
  omegab[1]=0.;
  if ((kk0!=0.)&&(Pb>0.)) omegab[2]=2.*PI/Pb; else omegab[2]=0.;

  /* Calculate orbital separation */
  a=binary_sep(q,M1,P);

  /* Calculate positions of M1 and M2 */
  xp1=-q/(1.+q)*a;
  xp2=1./(1.+q)*a;

  /* Find position of L1 point (xl1,0,0) */
  xl1=roche_xl1(q)*a;
  rl1=xl1-xp1;

  /* Find the Eggleton radius */
  regg=eggleton_radius(q)*a;

  /* Set disc edge */
  if (prdisc==0.) rdisc=regg; else rdisc=prdisc*a;
  
  /* Set parameters for magnetic drag */
  k0=kk0;
  r0=1.e8;
  magdist=0.;
  //magdist=0.1*rdisc;
  //magdist=0.1*a;
  magdb=-1.;

  if (k0<0.) {
    printf("libbinary.so : k0 must be >= 0\n");
    k0=0.;
  }
  if ((k0!=0.)&&(rdisc<rl1)) k0*=-1.; // Start with mag drag off unless rdisc>=rl1
  if (t==0.) t=0.25;

  dt=t*P/(double)(n*NSUB);

  /* Set initial conditions */
  vx[0]=vix[0];
  vy[0]=viy[0];
  vz[0]=viz[0];
  d[0]=0.;
  
  /*************************************************************************/

  /* Advance using Euler method */
  for (i=0; i<n-1; i++) {
    pos[0]=x[i];
    pos[1]=y[i];
    pos[2]=z[i];
    vel[0]=vx[i];
    vel[1]=vy[i];
    vel[2]=vz[i];
    dd=d[i];
    lmag[i]=0.;

    /* Do NSUB substeps */
    for(j=0; j<NSUB; j++) {
      stream_accel(acc,vin,pos,vel,dd,&k0, kindex, M1,q,rdisc,Pb, discbrfac, xp1,xp2,omega,omegab,r0,magdist,&magdb,&dlmag);

      dr[0]=vel[0]*dt;
      dr[1]=vel[1]*dt;
      dr[2]=vel[2]*dt;

      dd+=sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

      vel[0]+=acc[0]*dt;
      vel[1]+=acc[1]*dt;
      vel[2]+=acc[2]*dt;

      pos[0]+=dr[0];
      pos[1]+=dr[1];
      pos[2]+=dr[2];

      lmag[i]+=dlmag/NSUB; // Add to average lmag in this step

    }    

    /* Update the arrays */
    x[i+1]=pos[0];
    y[i+1]=pos[1];
    z[i+1]=pos[2];
    vx[i+1]=vel[0];
    vy[i+1]=vel[1];
    vz[i+1]=vel[2];
    vix[i]=vin[0];
    viy[i]=vin[1];
    viz[i]=vin[2];
    d[i+1]=dd;
    lmag[i+1]=dlmag; // Use last calculated lmag for last point

    /* Calculate the Keplerian velocity of this point */
    rxy[0]=x[i]-xp1;
    rxy[1]=y[i];
    rxy[2]=0.;    
    omegakep[0]=0.;
    omegakep[1]=0.;
    omegakep[2]=sqrt(G*M1*Msun/pow(rxy[0]*rxy[0]+rxy[1]*rxy[1],3./2.));
    cross_vector(vkep,omegakep,rxy);
    vkx[i]=vkep[0];
    vky[i]=vkep[1]-q/(1.+q)*a*omega[2];
    vkz[i]=vkep[2];

  }

  pos[0]=x[n-1];
  pos[1]=y[n-1];
  pos[2]=z[n-1];
  cross_vector(vin,omega,pos);
  vix[n-1]=vx[n-1]+vin[0];
  viy[n-1]=vy[n-1]+vin[1];
  viz[n-1]=vz[n-1]+vin[2];

  rxy[0]=x[n-1]-xp1;
  rxy[1]=y[n-1];
  rxy[2]=0.;    
  omegakep[0]=0.;
  omegakep[1]=0.;
  omegakep[2]=sqrt(G*M1*Msun/pow(rxy[0]*rxy[0]+rxy[1]*rxy[1],3./2.));
  cross_vector(vkep,omegakep,rxy);
  vkx[n-1]=vkep[0];
  vky[n-1]=vkep[1]-q/(1.+q)*a*omega[2];
  vkz[n-1]=vkep[2];

  free(vx);
  free(vy);
  free(vz);

}
