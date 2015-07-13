/* Function declarations for miscellaneous functions */

#ifndef __MISC
#define __MISC

#ifdef __cplusplus
extern "C" {
#endif

  void cross_vector(double *result, double *v, double *x);
  void unit_vector(double *result, double *x);
  void b_dipole(double *B, double *r, double m);
  double angle(double a, double o);
  double brentmin(double ax, double bx, double cx,
		  double (*f)(double), double tol, double *xmin);
  double fold(double phase);
  double inverse_matrix(double *input, double *inverse);
  double min(double *data, long int ndata);
  double max(double *data, long int ndata);
  double gaussprofile(double x, double fwhm);
  double lorentzprofile(double x, double fwhm);
  double gaussprofileamp(double x, double fwhm, double A);
  double lorentzprofileamp(double x, double fwhm, double A);
  double randomn(long *idum);
  double randomu(long *idum);
  long int randomseed();

#ifdef __cplusplus
}
#endif

#endif
