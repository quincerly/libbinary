/* Function declarations for binary.c functions */

#ifndef __BINARY
#define __BINARY

#ifdef __cplusplus
extern "C" {
#endif

double binary_sep(double q, double M1, double P);
double eggleton_radius(double q);
double roche_xl1(double q);
double roche_potl(double x, double y, double z, double q);
double find_potl(double pot, double q,
		 double xp, double yp, double zp,
		 double xq, double yq, double zq);
void vkep(double *vk,double x,double y,double z,double M1);
  void omega_kep(double *omegak,double x,double y,double z,double M1);

#ifdef __cplusplus
}
#endif

#endif

