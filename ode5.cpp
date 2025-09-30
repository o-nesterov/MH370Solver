/************************************************************************************
* Disclaimer:                                                                       *
*                                                                                   *
* This file is a part of the modelling software developed to assist in finding      *
* the terminal location of the Malaysian MH370, which went missing on 8 March 2014. *
* This project is not intended for profit; it is not owned by, affiliated to,       *
* or sponsored by any governmental or private business entity. The author grants    * 
* permission to use, redistribute, and modify this software, subject to that it     *
* will serve its original purpose.                                                  *
*                                                                                   *
* The software is provided "as is", without warranty of any kind, explicit or       *
* implied, including but not limited to the warranties of merchantability,          *
* fitness for a particular purpose, including the intended one, and                 *
* non-infringement. In no event shall the author be liable for any claim,           *
* damages, or losses, including third-party liability, whether in an action         *
* of contract, tort or otherwise, arising from, out of or in connection with        *
* this software, or the use, or other dealings in the software.                     *
*                                                                                   *
* Author: O. Nesterov, PhD, Independent Researcher/Consultant, September 2025       *
************************************************************************************/

#include <math.h>
#include <stdio.h>
#include <malloc.h>
#include <limits.h>
#include <memory.h>
#include "ode5.h"


double sign (double a, double b)
{
  return (b > 0.0) ? fabs(a) : -fabs(a);

} /* sign */


double min_d (double a, double b)
{
  return (a < b)?a:b;

} /* min_d */


double max_d (double a, double b)
{
  return (a > b)?a:b;

} /* max_d */


////////////////////////////////////////////////////////////////////
// set initial step
////////////////////////////////////////////////////////////////////
double hinit(unsigned n, FcnEqDiff fcn, double x, double* y, void* param,
			 double posneg, double* f0, double* f1, double* yy1, int iord,
			 double hmax, double* atoler, double* rtoler)
{
  double   dnf, dny, atoli, rtoli, sk, h, h1, der2, der12, sqr;
  unsigned i;

  dnf = 0.0;
  dny = 0.0;
  atoli = atoler[0];
  rtoli = rtoler[0];

  for (i = 0; i < n; i++)
  {
      sk = atoler[i] + rtoler[i] * fabs(y[i]);
      sqr = f0[i] / sk;
      dnf += sqr*sqr;
      sqr = y[i] / sk;
      dny += sqr*sqr;
  }

  if ((dnf <= 1.0E-10) || (dny <= 1.0E-10))
    h = 1.0E-6;
  else
    h = sqrt (dny/dnf) * 0.01;

  h = min_d (h, hmax);
  h = sign (h, posneg);

  /* perform an explicit Euler step */
  for (i = 0; i < n; i++)
    yy1[i] = y[i] + h * f0[i];
  fcn (n, x+h, yy1, f1, param);

  /* estimate the second derivative of the solution */
  der2 = 0.0;
  for (i = 0; i < n; i++)
  {
      sk = atoler[i] + rtoler[i] * fabs(y[i]);
      sqr = (f1[i] - f0[i]) / sk;
      der2 += sqr*sqr;
  }
  der2 = sqrt (der2) / h;

  /* step size is computed such that h**iord * max_d(norm(f0),norm(der2)) = 0.01 */
  der12 = max_d (fabs(der2), sqrt(dnf));
  if (der12 <= 1.0E-15)
    h1 = max_d (1.0E-6, fabs(h)*1.0E-3);
  else
    h1 = pow (0.01/der12, 1.0/(double)iord);
  h = min_d (100.0 * h, min_d (h1, hmax));

  return sign (h, posneg);

} /* hinit */


/////////////////////////////////////////////////////////////////////////
// dense output function
// i is the index of desired component
// returns interpolated solution at x of the desired component i 
/////////////////////////////////////////////////////////////////////////

double contd5(unsigned long i, double x, 
			  double xold, double hout,
			  double* rcont1, double* rcont2, double* rcont3, double* rcont4, double* rcont5,
			  FILE* fileout)
{
  double   theta, theta1;

  if (i == UINT_MAX)
  {
	  if (fileout) fprintf (fileout, "No dense output available for %uth component\n", i);
	  return 0.0;
  }

  theta = (x - xold) / hout;
  theta1 = 1.0 - theta;

  return rcont1[i] + theta*(rcont2[i] + theta1*(rcont3[i] + theta*(rcont4[i] + theta1*rcont5[i])));

} // contd5


/////////////////////////////////////////////////////////////////////////
// function to prepare output
/////////////////////////////////////////////////////////////////////////

void solout(long nr, double xold, double hout, double x, double* y, unsigned long neq,
			double* pXout, double** pYout, unsigned long nOut, unsigned long& currentOutIdx,
			double* rcont1, double* rcont2, double* rcont3, double* rcont4, double* rcont5,
			FILE* fileout)
{
	if (nOut<=0 || pYout==NULL || (currentOutIdx>=nOut)) return;

	if (nr == 1)
	{
		currentOutIdx = 0;
		while (x > pXout[currentOutIdx] && currentOutIdx < nOut) currentOutIdx++; // skip outputs at points, which are requested before xstart
		if (x == pXout[currentOutIdx])
		{
			for (unsigned int i=0; i<neq; i++) 
			{
				double*p = pYout[i];
				if (p) p[currentOutIdx] = y[i];
			}
		}
	}
	else
	{

		while (x >= pXout[currentOutIdx])
		{
			for (unsigned int i=0; i<neq; i++) 
			{
				double*p = pYout[i];
				if (p) p[currentOutIdx] = contd5(i,pXout[currentOutIdx],
												 xold, hout,
												 rcont1, rcont2, rcont3, rcont4, rcont5,
												 fileout);
			}

			currentOutIdx++;

			if (currentOutIdx >= nOut) break;
		}
	}
} /* solout */



////////////////////////////////////////////////////////////////////
// core integrator
////////////////////////////////////////////////////////////////////

int dopcor(unsigned long n, FcnEqDiff fcn, 	FcnTerminate chkTerm, FcnGetMaxStep getMaxStep,
		   double* pXout, double** pYout, unsigned long nOut, void* param,
		   double* pWorkspace, double h, double hmax, double uround, unsigned long nmax,
		   double safe, double fac1, double fac2, double beta, unsigned long nstiff, 
		   unsigned long& nfcn, unsigned long& nstep, unsigned long& naccpt, unsigned long& nrejct, unsigned long& currentOutIdx, 
		   double& xout, double* yout, double& xoldout, double* yoldout, double& hout,
		   int& termcode, FILE* fileout)
{
  double   facold, expo1, fac, facc1, facc2, fac11, posneg, xph;
  double   hlamb, err, sk, hnew, yd0, ydiff, bspl;
  double   stnum, stden, sqr;
  int      iasti, iord, reject, last, nonsti;
  unsigned i;
  double   c2, c3, c4, c5, e1, e3, e4, e5, e6, e7, d1, d3, d4, d5, d6, d7;
  double   a21, a31, a32, a41, a42, a43, a51, a52, a53, a54;
  double   a61, a62, a63, a64, a65, a71, a73, a74, a75, a76;

  /* initialisations */
  c2=0.2, c3=0.3, c4=0.8, c5=8.0/9.0;
  a21=0.2, a31=3.0/40.0, a32=9.0/40.0;
  a41=44.0/45.0, a42=-56.0/15.0; a43=32.0/9.0;
  a51=19372.0/6561.0, a52=-25360.0/2187.0;
  a53=64448.0/6561.0, a54=-212.0/729.0;
  a61=9017.0/3168.0, a62=-355.0/33.0, a63=46732.0/5247.0;
  a64=49.0/176.0, a65=-5103.0/18656.0;
  a71=35.0/384.0, a73=500.0/1113.0, a74=125.0/192.0;
  a75=-2187.0/6784.0, a76=11.0/84.0;
  e1=71.0/57600.0, e3=-71.0/16695.0, e4=71.0/1920.0;
  e5=-17253.0/339200.0, e6=22.0/525.0, e7=-1.0/40.0;
  d1=-12715105075.0/11282082432.0, d3=87487479700.0/32700410799.0;
  d4=-10690763975.0/1880347072.0, d5=701980252875.0/199316789632.0;
  d6=-1453857185.0/822651844.0, d7=69997945.0/29380423.0;


  // initializations
  // assign arrays in the provided heap of memory (pWorkspace is either allocated by dopri5 or user-provided)
  // supposedly it is more efficient for small systems than allocating each array individually because of contiguous memory 
  // and only one call to allocate and deallocate the heap
  double* y =      &pWorkspace[n*0];
  double* rtoler = &pWorkspace[n*1];
  double* atoler = &pWorkspace[n*2];
  double* rcont1 = &pWorkspace[n*3];
  double* rcont2 = &pWorkspace[n*4];
  double* rcont3 = &pWorkspace[n*5];
  double* rcont4 = &pWorkspace[n*6];
  double* rcont5 = &pWorkspace[n*7];
  double* yy1 =    &pWorkspace[n*8];
  double* k1 =     &pWorkspace[n*9];
  double* k2 =     &pWorkspace[n*10];
  double* k3 =     &pWorkspace[n*11];
  double* k4 =     &pWorkspace[n*12];
  double* k5 =     &pWorkspace[n*13];
  double* k6 =     &pWorkspace[n*14];
  double* ysti =   &pWorkspace[n*15];

  double x = pXout[0];			// start x 
  double xend =  pXout[nOut-1];	// end x
  double hmax_ = hmax;

  nfcn = 0;
  nstep = 0;
  naccpt = 0;
  nrejct = 0;
  currentOutIdx = 0;
  nonsti = 0;

  facold = 1.0E-4;
  expo1 = 0.2 - beta * 0.75;
  facc1 = 1.0 / fac1;
  facc2 = 1.0 / fac2;
  posneg = sign (1.0, xend-x);


  /* initial preparations */
  last  = 0;
  hlamb = 0.0;
  iasti = 0;
  fcn(n, x, y, k1, param);
  hmax_ = fabs (hmax_);
  if (getMaxStep)
  {
	  // custom maximum step
	  double hmax_custom = getMaxStep(n, x, y, param);
	  if (hmax_>abs(hmax_custom)) hmax_ = abs(hmax_custom);
  };
  iord = 5;
  if (h == 0.0) h = hinit(n, fcn, x, y, param, posneg, k1, k2, k3, iord, hmax_, atoler, rtoler);
  nfcn += 2;
  reject = 0;
  double xold = x;

  hout = h;
  xout = x;
  if (yout) for (i = 0; i < n; i++) yout[i] = y[i];

  solout(naccpt+1, xold, hout, x, y, n,
		 pXout, pYout, nOut, currentOutIdx,
		 rcont1, rcont2, rcont3, rcont4, rcont5,
		 fileout);
 

  xoldout = x;
  if (yoldout) for (i = 0; i < n; i++) yoldout[i] = yout[i];


  // basic integration step
  while (1)
  {
    if (nstep > nmax)
    {
      if (fileout) fprintf (fileout, "Exit of dopri5 at x = %.16e, more than nmax = %li are needed\r\n", x, nmax);
      hout = h;
	  xout = x;
	  if (yout) for (i = 0; i < n; i++) yout[i] = y[i];
      return ODE_STATUS_NMAX;
    }

    if (0.1 * fabs(h) <= fabs(x) * uround)
    {
		if (fileout) fprintf (fileout, "Exit of dopri5 at x = %.16e, step size too small h = %.16e\r\n", x, h);
		hout = h;
		xout = x;
		if (yout) for (i = 0; i < n; i++) yout[i] = y[i];
		return ODE_STATUS_TOOSMALLSTEP;
    }

	if ((x + 1.01*h - xend) * posneg > 0.0)
	{
		h = xend - x;
		last = 1;
	}

	nstep++;

    /* the first 6 stages */
    for (i = 0; i < n; i++) yy1[i] = y[i] + h * a21 * k1[i];
    fcn(n, x+c2*h, yy1, k2, param);

    for (i = 0; i < n; i++) yy1[i] = y[i] + h * (a31*k1[i] + a32*k2[i]);
    fcn(n, x+c3*h, yy1, k3, param);

    for (i = 0; i < n; i++) yy1[i] = y[i] + h * (a41*k1[i] + a42*k2[i] + a43*k3[i]);
    fcn(n, x+c4*h, yy1, k4, param);

    for (i = 0; i < n; i++) yy1[i] = y[i] + h * (a51*k1[i] + a52*k2[i] + a53*k3[i] + a54*k4[i]);
    fcn(n, x+c5*h, yy1, k5, param);

    for (i = 0; i < n; i++) ysti[i] = y[i] + h * (a61*k1[i] + a62*k2[i] + a63*k3[i] + a64*k4[i] + a65*k5[i]);
    xph = x + h;
    fcn(n, xph, ysti, k6, param);

    for (i = 0; i < n; i++) yy1[i] = y[i] + h * (a71*k1[i] + a73*k3[i] + a74*k4[i] + a75*k5[i] + a76*k6[i]);
    fcn (n, xph, yy1, k2, param);

	for (i = 0; i < n; i++) rcont5[i] = h * (d1*k1[i] + d3*k3[i] + d4*k4[i] + d5*k5[i] + d6*k6[i] + d7*k2[i]);


	for (i = 0; i < n; i++) k4[i] = h * (e1*k1[i] + e3*k3[i] + e4*k4[i] + e5*k5[i] + e6*k6[i] + e7*k2[i]);
	nfcn += 6;

    /* error estimation */
	err = 0.0;
    for (i = 0; i < n; i++)
    {
		sk = atoler[i] + rtoler[i] * max_d(fabs(y[i]), fabs(yy1[i]));
		sqr = k4[i] / sk;
		err += sqr*sqr;
	}

	err = sqrt (err / (double)n);

	/* computation of hnew */
	fac11 = pow (err, expo1);
	/* Lund-stabilization */
	fac = fac11 / pow(facold,beta);
	/* we require fac1 <= hnew/h <= fac2 */
	fac = max_d (facc2, min_d (facc1, fac/safe));
	hnew = h / fac;

    if (err <= 1.0)
    {
		// step accepted
		facold = max_d (err, 1.0E-4);
		naccpt++;

		// stiffness detection
		if (!(naccpt % nstiff) || (iasti > 0))
		{
			stnum = 0.0;
			stden = 0.0;
			for (i = 0; i < n; i++)
			{
				sqr = k2[i] - k6[i];
				stnum += sqr*sqr;
				sqr = yy1[i] - ysti[i];
				stden += sqr*sqr;
			};

			if (stden > 0.0) hlamb = h * sqrt (stnum / stden);
			
			if (hlamb > 3.25)
			{
				nonsti = 0;
				iasti++;
				if (iasti == 15)
					if (fileout)
						fprintf (fileout, "The problem seems to become stiff at x = %.16e\r\n", x);
					else
					{
						hout = h;
						xout = x;
						if (yout) for (i = 0; i < n; i++) yout[i] = y[i];
						return ODE_STATUS_STIFF;
					}
			}
			else
			{
				nonsti++;
				if (nonsti == 6) iasti = 0;
			}
		}


		for (i = 0; i < n; i++)
		{
			yd0 = y[i];
			ydiff = yy1[i] - yd0;
			bspl = h * k1[i] - ydiff;
			rcont1[i] = y[i];
			rcont2[i] = ydiff;
			rcont3[i] = bspl;
			rcont4[i] = -h * k2[i] + ydiff - bspl;
		}

		memcpy (k1, k2, n * sizeof(double)); 
		memcpy (y, yy1, n * sizeof(double));
		xold = x;
		x = xph;


		hout = h;
		xoldout = xout;
		xout = x;
		if (yoldout) for (i = 0; i < n; i++) yoldout[i] = yout[i];
		if (yout) for (i = 0; i < n; i++) yout[i] = y[i];

		solout(naccpt+1, xold, hout, x, y, n,
			   pXout, pYout, nOut, currentOutIdx,
			   rcont1, rcont2, rcont3, rcont4, rcont5,
			   fileout);

		if (chkTerm)
		{
			// check termination condition
			termcode = chkTerm(n, x, y, param);
			if (termcode!=0)
			{
				hout = hnew;
				xout = x;
				if (yout) for (i = 0; i < n; i++) yout[i] = y[i];
				return ODE_STATUS_INTERRUPTED;
			}
		}

		// normal exit
		if (last)
		{
			hout = hnew;
			xout = x;
			if (yout) for (i = 0; i < n; i++) yout[i] = y[i];
			return ODE_STATUS_COMPLETED;
		}

		hmax_ = hmax;
		if (getMaxStep)
		{
			// custom maximum step
		  double hmax_custom = getMaxStep(n, x, y, param);
		  if (hmax_>abs(hmax_custom)) hmax_ = abs(hmax_custom);
		};

		if (fabs(hnew) > hmax_) hnew = posneg * hmax_;
		if (reject) hnew = posneg * min_d (fabs(hnew), fabs(h));

		reject = 0;
	}
	else
	{
		// step rejected
		hnew = h / min_d (facc1, fac11/safe);
		reject = 1;
		if (naccpt >= 1) nrejct++;
		last = 0;
	}
	h = hnew;
}

} /* dopcor */






/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Main Interfaces
// integration function; allocates and returns array of solution of size [neq x xlen] 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

double** dopri5(
				unsigned long neq,		/* number of equations */
				FcnEqDiff fnc,			/* ODE callback function to evaluate y'(x,y) */
				double* x,				/* vector of x where solution is required */
				unsigned long xlen,		/* length of vector x */
				double* ystart,			/* initial values for y - vector of neq length*/
				void* param,			/* pointer to parameter, which can be converted to an object, structure, array, etc.*/
				ODE45Options& options,	/* integration options*/
				ODE45Stats& stats		/* output statistics*/
			   )
{
	// prepare output statistics in case inputs are wrong
	stats.nfcn = 0;
	stats.nstep = 0;
	stats.naccpt = 0;
	stats.nrejct = 0;
	stats.xout = x[0];
	stats.hout = 0.0;
	stats.yout = NULL;
	stats.xoldout = 0.0;
	stats.yoldout = NULL;
	stats.status = ODE_STATUS_NOTSTARTED;	// not started status
	stats.termcode = 0; // caller-defined termination code in case termination is requested externally due to certain conditions met

	FILE* fileout = options.fileout;	// output message stream (can be NULL)

	// Allocate workspace if needed
	double* pWorkspace_ = NULL; // buffer to be 64-byte aligned
	double* pWorkspace = NULL;

	if (options.pWorkspace == NULL) 
	{
		pWorkspace_ = new double[16*neq+9]; // allocate heap workspace if it was not externally allocated to speed-up calculations in case of multiple calls, with extra 72 bytes
		if (!pWorkspace_)
		{
			if (fileout) fprintf (fileout, "Not enough free memory\n");
			return NULL;
		}
		int pWorkspace__ = (int) pWorkspace_;
		pWorkspace__ = pWorkspace__ >> 6;
		pWorkspace__ = pWorkspace__ << 6;
		pWorkspace__ = pWorkspace__ + 64;
		pWorkspace = (double*)pWorkspace__;
	}
	else
	{
		pWorkspace_ = NULL;
		pWorkspace = options.pWorkspace;
	}

	// assign rtoler and atoler arrays in the heap
	double* y = &pWorkspace[neq*0];
	double* rtoler = &pWorkspace[neq*1];
	double* atoler = &pWorkspace[neq*2];

	unsigned long n = neq;					// dimension of the system <= UINT_MAX-1
	unsigned long nmax = options.nmax;		// maximal number of allowed steps
	unsigned long nstiff = options.nstiff;	// test for stiffness
	double uround = options.uround;			// rounding unit
	double hmax = options.hmax;				// maximal step size
	double h = options.h;					// initial step size

	// copy arrays of the relative and absolute accuracies
	for (unsigned int i=0; i<neq; i++) rtoler[i] = options.rtoler[i];
	for (unsigned int i=0; i<neq; i++) atoler[i] = options.atoler[i];

	double safe = options.safe;		// safety factor
	double fac1 = options.fac1;		// parameters for step size selection
	double fac2 = options.fac2;
	double beta = options.beta;		// for stabilized step size control

	double* yout = options.yout;		// last y[neq], for which solution is computed (can be NULL)
	double* yoldout = options.yoldout;	// second last y[neq], for which solution is computed (can be NULL)
	// ------------------------------------------------------------------------
	// sanity checks

	bool arret = false; 

	// n, the dimension of the system
	if (n == UINT_MAX)
	{
		if (fileout) fprintf (fileout, "System too big, max. n = %u\r\n", UINT_MAX-1);
		arret = true;
	}

	// nmax, the maximal number of steps
	if (!nmax)
	{
		nmax = 100000;
	}
	else if (nmax <= 0)
	{
		if (fileout) fprintf (fileout, "Wrong input, nmax = %li\r\n", nmax);
		arret = true;
	}


	// nstiff, parameter for stiffness detection
	if (!nstiff)
	{
		nstiff = 1000;
	}
	else if (nstiff < 0)
	{
		nstiff = nmax + 10;
	}


	// uround, smallest number satisfying 1.0+uround > 1.0
	if (uround == 0.0)
	{
		uround = 2.3E-16;
	}
	else if ((uround <= 1.0E-35) || (uround >= 1.0))
	{
		if (fileout) fprintf (fileout, "Which machine do you have ? Your uround was : %.16e\r\n", uround);
		arret = true;
	}

	// safety factor
	if (safe == 0.0)
	{
			safe = 0.9;
	}
	else if ((safe >= 1.0) || (safe <= 1.0E-4))
	{
		if (fileout) fprintf (fileout, "Curious input for safety factor, safe = %.16e\r\n", safe);
		arret = true;
	}

	// fac1, fac2, parameters for step size selection
	if (fac1 == 0.0) fac1 = 0.2;
	if (fac2 == 0.0) fac2 = 10.0;

	// beta for step control stabilization */
	if (beta == 0.0) 
		beta = 0.04;
	else if (beta < 0.0)
		beta = 0.0;
	else if (beta > 0.2)
	{
		if (fileout) fprintf (fileout, "Curious input for beta : beta = %.16e\r\n", beta);
		arret = true;
	}

	// allocate output arrays
	double** pYout = new double*[neq+1];
	if (pYout)
	{
		for (unsigned long i=0; i<neq; i++) 
		{
			pYout[i] = new double [xlen+1]; // allocate output arrays of required length per component
			if (!pYout[i]) 
			{
				if (fileout) fprintf (fileout, "Not enough free memory to allocate output space\n");
				arret = true;
			}
		}
	}
	else
	{
		if (fileout) fprintf (fileout, "Not enough free memory to allocate output space\n");
		arret = true;
	}


	// if a failure has occured...
	if (arret)
	{
		stats.status =  ODE_STATUS_INCONSISTENT;
		if (pWorkspace_) delete pWorkspace_;
		for (unsigned long i=0; i<neq; i++) if (pYout[i]) delete pYout[i];
		delete pYout;
		return NULL;
	}

	// --------------------------------------

	// maximal step size
	if (hmax == 0.0) hmax = x[xlen-1] - x[0];

	// copy initial conditions
	for (unsigned int i=0; i<neq; i++) y[i] = ystart[i];

	// call core integration
	unsigned long nfcnRead;		// number of ODE function calls
	unsigned long nstepRead;	// number of used steps
	unsigned long naccptRead;	// number of accepted steps
	unsigned long nrejctRead;	// number of rejected steps
	unsigned long nout;
	double hout;		// predicted step size of the last accepted step (useful for a subsequent call to dopri5)
	double xout;		// x value for which the solution has been computed (x=xend after successful return)
	double xoldout;		// second last x value for which the solution has been computed - used to assess interruption conditions

	for (unsigned long i=0; i<neq; i++) y[i] = ystart[i];	// copy initial conditions (y modified by dopcor)

	stats.status = dopcor(neq, fnc, options.chkTerm, options.getMaxStep, x, pYout, xlen, param,
						  pWorkspace, h, hmax, uround, nmax,
						  safe, fac1, fac2, beta, nstiff, 
						  nfcnRead, nstepRead, naccptRead, nrejctRead, nout, 
						  xout, yout, xoldout, yoldout, hout,
						  stats.termcode, fileout);


	stats.nfcn = nfcnRead;
	stats.nstep = nstepRead;
	stats.naccpt = naccptRead;
	stats.nrejct = nrejctRead;
	stats.nout = nout;
	stats.hout = hout; 
	stats.xout = xout;
	stats.yout = yout;
	stats.xoldout = xoldout;
	stats.yoldout = yoldout;

	// clear workspace
	if (pWorkspace_) delete pWorkspace_;

	return pYout;
} // dopri5

