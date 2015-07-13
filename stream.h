/* Function declarations for stream.c functions */

#ifndef __STREAM
#define __STREAM

#ifdef __cplusplus
extern "C" {
#endif
  
void stream_accel(double *accel, double *vin, double *r, double *v, double dist,
		  double *k0, double kindex, double M1, double q, double rdisc, double Pb, double discbrfac,
		  double xp1, double xp2, double *omega, double *omegab,
		  double r0, double magdist, double *magdb, double *lmag);
  
void stream_calc(double q, double M1, double P,
		 double Pb, double kk0, double kindex, double prdisc, double discbrfac,
		 double t, long n, long NSUB,
		 double *x, double *y, double *z,
		 double *vix, double *viy, double *viz,
		 double *vkx, double *vky, double *vkz,
		 double *d, double *lmag);

#ifdef __cplusplus
}
#endif

#endif
