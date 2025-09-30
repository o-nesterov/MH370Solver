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

/*     
This code computes the numerical solution of a system of first order ordinary
differential equations y'=f(x,y). It uses an explicit Runge-Kutta method of
order (4)5 due to Dormand & Prince with step size control and dense output.

Authors : E. Hairer & G. Wanner
	  Universite de Geneve, dept. de Mathematiques
	  CH-1211 GENEVE 4, SWITZERLAND
	  E-mail : HAIRER@DIVSUN.UNIGE.CH, WANNER@DIVSUN.UNIGE.CH

The code is described in : E. Hairer, S.P. Norsett and G. Wanner, Solving
ordinary differential equations I, nonstiff problems, 2nd edition,
Springer Series in Computational Mathematics, Springer-Verlag (1993).



INPUT PARAMETERS
----------------

n        Dimension of the system (n < UINT_MAX).

fcn      A pointer the the function definig the differential equation, this
	 function must have the following prototype

	   void fcn (unsigned n, double x, double *y, double *f)

	 where the array f will be filled with the function result.

x        Initial x value.

*y       Initial y values (double y[n]).

xend     Final x value (xend-x may be positive or negative).

*rtoler  Relative and absolute error tolerances. They are both
*atoler  vectors of length n (in the scalar case pass the addresses of
	     variables where you have placed the tolerance values).
	     The code keeps the local error of y[i] below
		 rtoler[i]*abs(y[i])+atoler[i].

solout   Output function called during integration.
	 Continuous output : during the calls to solout, a continuous solution
	 for the interval (xold,x) is available through the function

	   contd5(i,s)

	 which provides an approximation to the i-th component of the solution
	 at the point s (s must lie in the interval (xold,x)).


fileout  A pointer to the stream used for messages, if you do not want any
	 message, just pass NULL.


Sophisticated setting of parameters
-----------------------------------

	 Several parameters have a default value (if set to 0) but, to better
	 adapt the code to your problem, you can specify particular initial
	 values.

uround   The rounding unit, default 2.3E-16 (this default value can be
	 replaced in the code by DBL_EPSILON providing float.h defines it
	 in your system).

safe     Safety factor in the step size prediction, default 0.9.

fac1     Parameters for step size selection; the new step size is chosen
fac2     subject to the restriction  fac1 <= hnew/hold <= fac2.
	 Default values are fac1=0.2 and fac2=10.0.

beta     The "beta" for stabilized step size control (see section IV.2 of our
	 book). Larger values for beta ( <= 0.1 ) make the step size control
	 more stable. dopri5 needs a larger beta than Higham & Hall. Negative
	 initial value provoke beta=0; default beta=0.04.

hmax     Maximal step size, default xend-x.

h        Initial step size, default is a guess computed by the function hinit.

nmax     Maximal number of allowed steps, default 100000.

nstiff   Test for stiffness is activated when the current step number is a
	 multiple of nstiff. A negative value means no test and the default
	 is 1000.


Memory requirements
-------------------

	 The function dopri5 allocates dynamically 8*n doubles for the method
	 stages, 5*n doubles for the interpolation


OUTPUT PARAMETERS
-----------------

y       numerical solution at x=xRead() (see below).

dopri5 returns the following values

	 1 : computation successful,
	 2 : computation successful interrupted by solout,
	-1 : input is not consistent,
	-2 : larger nmax is needed,
	-3 : step size becomes too small,
	-4 : the problem is probably stff (interrupted).


Several functions provide access to different values :

xRead   x value for which the solution has been computed (x=xend after
	successful return).

hRead   Predicted step size of the last accepted step (useful for a
	subsequent call to dopri5).

nstepRead   Number of used steps.
naccptRead  Number of accepted steps.
nrejctRead  Number of rejected steps.
nfcnRead    Number of function calls.


*/


#include <stdio.h>
#include <limits.h>

// not started
#define ODE_STATUS_NOTSTARTED 0
// completed
#define ODE_STATUS_COMPLETED 1
// interrupted
#define ODE_STATUS_INTERRUPTED 2
//input is not consistent
#define ODE_STATUS_INCONSISTENT -1
// larger nmax is needed
#define ODE_STATUS_NMAX -2
// step size becomes too small
#define ODE_STATUS_TOOSMALLSTEP -3
// the problem is probably stiff (interrupted)
#define ODE_STATUS_STIFF -4


// ODE callback function definition to evaluate y'(x,y)
typedef void (*FcnEqDiff)(unsigned long neq, double x, double *y, double *f, void* param);

// callback function definition to prematurely terminate integration if user-specified conditions are met (e.g., target altitude or waypoint)
// it should return code, which identifies reason for termination
typedef int (*FcnTerminate)(unsigned long neq, double x, double *y, void* param);

// callback function definition to dynamically impose user-specified maximum integration time step 
typedef double (*FcnGetMaxStep)(unsigned long neq, double x, double *y, void* param);


// Integration options
struct ODE45Options
{
	long nmax;					// maximal number of allowed steps
	long nstiff;				// test for stiffness
	double uround;				// rounding unit
	double hmax;				// maximal step size
	double h;					// initial step size
	double* rtoler;				// relative error tolerance for every component
	double* atoler;				// absolute error tolerance for every component
	double safe;				// safety factor
	double fac1;				// parameters for step size selection
	double fac2;
	double beta;				// for stabilized step size control
	FILE* fileout;				// messages stream - can be NULL
	double* pWorkspace;			// reusalble user-provided workspace of the minimum size: 16 x number of ODE components (if NULL, then allocated internally)
	double* yout;				// buffer to contain last values of y[neq], for which solution is computed (can be NULL)
	double* yoldout;			// buffer to contain secnd last values of y[neq], for which solution is computed (can be NULL)
	FcnTerminate chkTerm;		// custom function to check user-defined termination condition (can be NULL)
	FcnGetMaxStep getMaxStep;	// custom function to get dynamically updated maximum time step (can be NULL)
};

// -----------------
// Output statistics
struct ODE45Stats
{
	long nfcn;		// Number of ODE function calls.
	long nstep;		// Number of used steps.
	long naccpt;	// Number of accepted steps.
	long nrejct;	// Number of rejected steps.
	long nout;		// Number of successful outputs.
	double hout;	// Predicted step size of the last accepted step (useful for a subsequent call to dopri5).
	double xout;	// final x value for which the solution has been computed (x=xend after successful return).
	double* yout;	// final array of y[neq], for which the solution has been computed
	double xoldout;	// the second last x value for which the solution has been computed (used in case of interruption)
	double* yoldout;// the second last array of y[neq], for which the solution has been computed (used in case of interruption)
	int status;		// status flag returned by dopri5
	int termcode;	// termination code returned by external termination request (which conditions are met)
};



// integration function; allocates and returns array of solution vectors [neq x xlen] 
double** dopri5(
				unsigned long neq,		/* number of equations */
				FcnEqDiff fnc,			/* ODE callback function to evaluate y'(x,y) */
				double* x,				/* vector of x where solution is required */
				unsigned long xlen,		/* length of vector x */
				double* ystart,			/* initial values for y - vector of neq length*/
				void* param,			/* pointer to parameter to be passed to function (can be converted to an object, structure, array, etc.)*/
				ODE45Options& options,	/* integration options*/
				ODE45Stats& stats		/* output statistics*/
			   );

