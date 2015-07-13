#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "binary.h"
#include "misc.h"
#include "visibility.h"

#define TOL 5.e-6
#define PI 3.1415926536

/* Global variables for wrapper functions */
static double visibility_x,visibility_y,visibility_z;
static double visibility_dx,visibility_dy,visibility_dz;
static double visibility_q;

/* Protos for wrapper functions */
double visibility_func(double l);

/* Calculate R1 parameter required by visibility */
double vis_R1(double q)
{
  double xl1,rl1,rl2,x1,x2,potl1;

  xl1=roche_xl1(q);
  potl1=roche_potl(xl1,0.,0.,q);
  x1=-q/(1.+q);
  x2=1./(1.+q);
  rl1=xl1-x1;
  rl2=x2-xl1;

  return find_potl(potl1,q,x2,0.,0.01*rl2,x2,0.,rl2)*rl2;

}

/* Calculate R2 parameter required by visibility */
double vis_R2(double q)
{
  double xl1,x2;

  xl1=roche_xl1(q);
  x2=1./(1.+q);

  return x2-xl1;
}

// Calculate visibility grid at n phase values in orbphase
void calc_vis(int *vis_array, long DIM, double pix, double *orbphase, double *thetaphase, long n, double q, double incdeg, int display)
{

  double R1,R2; /* inner and outer circle radii contraining secondary Roche lobe */
  double xl1, potl1; /* x position of L1 and roche potential at L1 */
  double x1,x2;
  double xp,yp,xxp,yyp,x,y;
  double dx,dy,dz,IMAGESIZEA;
  double theta,costheta,sintheta;
  long p,i,j,ii,jj;
  int vis,v;

  x1=-q/(1.+q);
  x2=1./(1.+q);
  xl1=roche_xl1(q);
  potl1=roche_potl(xl1,0.,0.,q);
  R1=vis_R1(q);
  R2=vis_R2(q);
  IMAGESIZEA=pix*(double)DIM*0.5;

  if (display==1 ) printf("Calculating visibility    0%%");
  fflush(stdout);

  for (p=0;p<n;p++) {
    // (dx,dy,dx) is direction of observer measured in corotation frame at this phase 
    dxdydz(&dx,&dy,&dz,incdeg,orbphase[p]);

    // Go through each point in image calculating visibility. 
    // Inrease accuracy by calculating at 9 points in each pixel.
    // Frame of image is centred on M1, and rotated about M1 by angle
    // theta clockwise from line of centres. Image x-axis (xp) is at
    // angle theta to line of centres.
    theta=thetaphase[p]*2.*PI;
    costheta=cos(theta);
    sintheta=sin(theta);
    for(i=0;i<DIM;i++) {
      xp=(-(double)(DIM-1)/2.+(double)i)*pix;
      for(j=0;j<DIM;j++) {
	yp=(-(double)(DIM-1)/2.+(double)j)*pix;
	vis=0;
	if (xp*xp+yp*yp<=IMAGESIZEA*IMAGESIZEA) {
	  for(ii=-1;ii<2;ii++) {
	    xxp=pix/3.*ii;
	    for(jj=-1;jj<2;jj++) {	  
	      yyp=pix/3.*jj;
	      x=(xp+xxp)*costheta+(yp+yyp)*sintheta+x1;
	      y=-(xp+xxp)*sintheta+(yp+yyp)*costheta;
	      v=visibility(x,y,0.,dx,dy,dz,R1,R2,q,potl1);
/*	      printf("%f %f %f %f %f %f\n",yp+yyp,xp+xxp,sintheta,costheta,x,y);
	      printf("%f %f\n",x,y);*/
	      if (v==-1) exit(1);
	      vis=vis+v;
	    }
	  }
	}
	vis_array[i*DIM+j+DIM*DIM*p]=vis;
      }
    }
    if (display==1 ) printf("\b\b\b\b%3ld%%",100*(p+1)/n);
    fflush(stdout);
  }

}

// Calculate (dx,dy,dx) - direction of observer measured in corotating frame
void dxdydz(double *dx, double *dy, double *dz, double incdeg, double phase)
{
  double inc=incdeg/180.*PI;  
  double phi=phase*2.*PI;

  *dx=sin(inc)*cos(phi);
  *dy=-sin(inc)*sin(phi);
  *dz=cos(inc);

}

// Calculate full phase width of eclipse of centre of M1 for q,i
double eclipsewidth(double q, double incdeg)
{

  // Set up variables
  double xl1,potl1,R1,R2;
  double dx,dy,dz;
  double pa,pb,pc;
  double va,vb,vc;

  xl1=roche_xl1(q);
  potl1=roche_potl(xl1,0.,0.,q);
  R1=vis_R1(q);
  R2=vis_R2(q);

  pa=0.;
  pc=0.5;
  dxdydz(&dx,&dy,&dz,incdeg,pa);
  va=visibility(-q/(1.+q),0.,0.,dx,dy,dz,R1,R2,q,potl1);
  dxdydz(&dx,&dy,&dz,incdeg,pc);
  vc=visibility(-q/(1.+q),0.,0.,dx,dy,dz,R1,R2,q,potl1);
  if (va!=0) return 0.; // No eclipse of centre of M1  
  if (vc==0) return 1.; // Centre of M1 always eclipse!!?
  do {
    pb=0.5*(pa+pc);
    dxdydz(&dx,&dy,&dz,incdeg,pb);
    vb=visibility(-q/(1.+q),0.,0.,dx,dy,dz,R1,R2,q,potl1);
    if (vb==0) {
      va=vb;
      pa=pb;
      
    }
    else {
      vc=vb;
      pc=pb;
    }
  } while (pc-pa>TOL);

  return pa+pc;  

}

// Calculate full phase width of eclipse of centre of (x,y,z) for q,i
double eclipsecontact(double x, double y, double z, double q, double incdeg,double *pingress, double *pegress)
{

  // Set up variables
  double xl1,potl1,R1,R2;
  double dx,dy,dz;
  double pa,pb,pc;
  double va,vb,vc;
  double phase,pi1,pi2,pe1,pe2;
  double v;
  int i;

  xl1=roche_xl1(q);
  potl1=roche_potl(xl1,0.,0.,q);
  R1=vis_R1(q);
  R2=vis_R2(q);

  i=0;
  do {
    phase=(double)i/200-0.5;
    dxdydz(&dx,&dy,&dz,incdeg,(double)i/200-0.5);
    v=visibility(x,y,z,dx,dy,dz,R1,R2,q,potl1);
    i++;
  }
  while (v==1&&phase<0.5);
  i--;
  if (v==1) return 0; // No eclipse of (x,y,z) 
  if (i==0) {
    // phase -0.5 is eclipsed. Work back.
    do {
      phase=(double)i/200-0.5;
      dxdydz(&dx,&dy,&dz,incdeg,(double)i/200-0.5);
      v=visibility(x,y,z,dx,dy,dz,R1,R2,q,potl1);
      i--;
    }
    while (v==0&&phase>-1.);
    i++;
    if (v==0) return 1; // (x,y,z) always eclipsed!?
    i++;
  }
  pi1=(double)(i-1)/200-0.5;
  pi2=(double)i/200-0.5;
  // Eclipse ingress is between pi1 and pi2
  pe1=pi2;
  pe2=pi1+1.;
  // Eclipse egress is between pe1 and pe2

  // Bracket to find exact ingress phase
  pa=pi1;
  pc=pi2;
  dxdydz(&dx,&dy,&dz,incdeg,pa);
  va=visibility(x,y,z,dx,dy,dz,R1,R2,q,potl1);
  dxdydz(&dx,&dy,&dz,incdeg,pc);
  vc=visibility(x,y,z,dx,dy,dz,R1,R2,q,potl1);
  do {
    pb=0.5*(pa+pc);
    dxdydz(&dx,&dy,&dz,incdeg,pb);
    vb=visibility(x,y,z,dx,dy,dz,R1,R2,q,potl1);
    if (vb==1) {
      va=vb;
      pa=pb;
      
    }
    else {
      vc=vb;
      pc=pb;
    }
  } while (pc-pa>TOL*0.5);
  *pingress=0.5*(pa+pc);  

  // Bracket to find exact egress phase
  pa=pe1;
  pc=pe2;
  dxdydz(&dx,&dy,&dz,incdeg,pa);
  va=visibility(x,y,z,dx,dy,dz,R1,R2,q,potl1);
  dxdydz(&dx,&dy,&dz,incdeg,pc);
  vc=visibility(x,y,z,dx,dy,dz,R1,R2,q,potl1);
  do {
    pb=0.5*(pa+pc);
    dxdydz(&dx,&dy,&dz,incdeg,pb);
    vb=visibility(x,y,z,dx,dy,dz,R1,R2,q,potl1);
    if (vb==0) {
      va=vb;
      pa=pb;
      
    }
    else {
      vc=vb;
      pc=pb;
    }
  } while (pc-pa>TOL*0.5);
  *pegress=0.5*(pa+pc);  
  
  return *pegress-*pingress;

}

/* Calculate inclination of system for which */
/* mass ratio is q and full eclipse width of */
/* white dwarf is w (in phase units) */
double inclination(double q, double w)
{

  double max_width,incdeg;

  // Maximum with for q
  max_width=eclipsewidth(q,90.);
  if (w>max_width) {
    //    printf("Maximum width for mass ratio %g is %g when i = 90 degrees\n",q,max_width);
    return -max_width;
  }
  if (w<=0.) {
    printf("Don't be awkward!\n");
    return -1.;
  }

  // Find i to give width w for mass ratio q
  incdeg=90.;
  while (eclipsewidth(q,incdeg)>w) incdeg=incdeg-0.01;

  return incdeg;
}

/* Check to see if line from (x,y,z) in direction (dx,dy,dz) */
/* passes through secondary */                                 
double visibility_func(double l) /* Function used only by visibility */
{

  return roche_potl(visibility_x+l*visibility_dx,
		    visibility_y+l*visibility_dy,
		    visibility_z+l*visibility_dz,
		    visibility_q);
  
}
int visibility(double x, double y, double z,
		double dx, double dy, double dz,
		double R1, double R2,
		double q, double potl1)
{

  double lcl,rc,xc,yc,zc;
  double pot,potmin,lmin;
  double A,B,C,D,l1,l2;
  double x1,x2;

  x1=-q/(1.+q);
  x2=1./(1.+q);

  /* Find closest approach of line of sight with centre of secondary */
  lcl=((x2-x)*dx-y*dy-z*dz)/sqrt(dx*dx+dy*dy+dz*dz);
  xc=x+dx*lcl;
  yc=y+dy*lcl;
  zc=z+dz*lcl;
  rc=sqrt((xc-x2)*(xc-x2)+yc*yc+zc*zc); /* Closest approach distance */

  /* Easy checks */
  if (rc>R2) {
    return 1; /* Point fully visible */
  }
  if ((rc<R1)&&(lcl>0.)) {
    return 0; /* Point fully eclipsed */
  }
  if (lcl<0.) {
    /* Check to see if (x,y,z) is outside of and in front of secondary */
    if ((x-x2)*(x-x2)+y*y+z*z>R2*R2) {
      return 1; /* We're in front of the secondary */
    }
    pot=roche_potl(x,y,z,q);
    if (pot>potl1) {
      return 1; /* We're in front of the secondary */
    }
    if ((pot<potl1)&&(x>0.)) {
      return 0; /* We're inside the secondary */
    }
  }
  
  /* Find where line of sight intersects outer sphere */ 
  A=dx*dx+dy*dy+dz*dz;
  B=2.*(dx*(x-x2)+dy*y+dz*z);
  C=(x-x2)*(x-x2)+y*y+z*z-R2*R2;
  D=sqrt(B*B-4.*A*C);
  l1=(-B-D)/2./A;
  l2=(-B+D)/2./A;

  /* Find minimum potential between l1 and l2 */
  visibility_x=x;
  visibility_y=y;
  visibility_z=z;
  visibility_dx=dx;
  visibility_dy=dy;
  visibility_dz=dz;
  visibility_q=q;

  potmin=brentmin(l1,l1+0.5*(l2-l1),l2,visibility_func,TOL,&lmin);

  if (potmin<=potl1) return 0; else return 1;

}
