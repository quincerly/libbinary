/* Function declarations for visibility.c functions */

#ifndef __VIS
#define __VIS

#ifdef __cplusplus
extern "C" {
#endif
  
  double vis_R1(double q);

  double vis_R2(double q);

  int visibility(double x, double y, double z,
		 double dx, double dy, double dz,
		 double R1, double R2,
		 double q, double potl1);

  void calc_vis(int *vis_array, long DIM, double pix,
		double *orbphase, double *thetaphase, long n,
		double q, double incdeg, int display);

  void dxdydz(double *dx, double *dy, double *dz,
	      double incdeg, double phase);    

  double eclipsewidth(double q, double i);

  double eclipsecontact(double x, double y, double z,
		     double q, double incdeg,
		     double *pingress, double *pegress);

  double inclination(double q, double w);

#ifdef __cplusplus
}
#endif

#endif
