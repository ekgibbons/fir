/* original version found at http://www.ews.uiuc.edu/%7Ejanovetz/DSP/remez/ */
/* added static modifiers before most functions */

/**************************************************************************
 * Parks-McClellan algorithm for FIR filter design (C version)
 *-------------------------------------------------
 * This code is adapted from the scipysignalmodule.  

 *************************************************************************/
#ifndef __REMEZ_SCIPY_H__
#define __REMEZ_SCIPY_H__

#define BANDPASS       1
#define DIFFERENTIATOR 2
#define HILBERT        3

#define NEGATIVE       0
#define POSITIVE       1

#define GRIDDENSITY    16
#define MAXITERATIONS  90
#define M_2PI          6.28318530717958623200

/* Function prototype for remez() - the only function that should need be
 * called from external code
 */
int remez_scipy(double *h2, int numtaps, int numbands, double *bands,
	      double *response, double *weight, int type, int maxiter,
	      int grid_density);


#endif /* __REMEZ_SCIPY_H__ */
