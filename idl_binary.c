#include "stream.h"
#include "binary.h"
#include "visibility.h"
#include "misc.h"
#include "malloc.h"
#include "stdio.h"
#include "export.h"

int idl_stream(int argc, void *argv[])
{
  stream_calc(*((double *)argv[0]), *((double *)argv[1]), *((double *)argv[2]), /* q, M1, P    */
	      *((double *)argv[3]), *((double *)argv[4]), *((double *)argv[5]), /* Pb, k0, kindex */
	      *((double *)argv[6]), *((double *)argv[7]),                       /* rdisc, discbrfac*/
	      *((double *)argv[8]),                                             /* t              */
	      (long)(*(IDL_LONG *)argv[9]),(long)(*(IDL_LONG *)argv[10]),       /* n, NSUB        */
	      (double *)argv[11],(double *)argv[12],(double *)argv[13],         /* x,y,z       */
	      (double *)argv[14],(double *)argv[15],(double *)argv[16],         /* vix,viy,viz */
	      (double *)argv[17],(double *)argv[18],(double *)argv[19],         /* vkx,vky,vkz */
	      (double *)argv[20],(double *)argv[21]);                           /* d,lmag      */ 
  return 1;
}

double idl_eggleton_radius(int argc, void *argv[])
{ return eggleton_radius(*((double *)argv[0])); }

double idl_binary_sep(int argc, void *argv[])
{ return binary_sep(*((double *)argv[0]),*((double *)argv[1]),*((double *)argv[2])); }

double idl_roche_xl1(int argc, void *argv[])
{ return roche_xl1(*((double *)argv[0])); }

double idl_roche_potl(int argc, void *argv[])
{ return roche_potl(*((double *)argv[0]),*((double *)argv[1]),*((double *)argv[2]),*((double *)argv[3])); }

double idl_find_roche(int argc, void *argv[])
{
  return find_potl(*((double *)argv[0]),*((double *)argv[1]),
		   *((double *)argv[2]),*((double *)argv[3]),*((double *)argv[4]),            /* xp,yp,zp */
		   *((double *)argv[5]),*((double *)argv[6]),*((double *)argv[7]));           /* xq,yq,zq */
}

double idl_eclipsewidth(int argc, void *argv[])
{
  return eclipsewidth(*((double *)argv[0]),*((double *)argv[1]));
}

double idl_eclipsecontact(int argc, void *argv[])
{
  return eclipsecontact(*((double *)argv[0]),*((double *)argv[1]),*((double *)argv[2]),
		      *((double *)argv[3]),*((double *)argv[4]),
		      (double *)argv[5],(double *)argv[6]);
}

double idl_inclination(int argc, void *argv[])
{
  return inclination(*((double *)argv[0]),*((double *)argv[1]));
}
