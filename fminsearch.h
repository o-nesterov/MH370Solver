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

#include <stdlib.h>
#include <stdio.h>

// Integration options
struct FMINSearchOptions
{
	long MaxFunEvals;	// maximum number of allowed function evaluations.
	long MaxIter;		// maximum number of iterations.
	double* TolX;		// relative error tolerance for every component.
	double TolFun;		// absolute error tolerance for every component.
	double* pWorkspace;	// workspace of minimum size n*(9+n)+1. Can be null - then allocated dynamically.
	FILE* outstream;	// output message stream (can be NULL - then no output).
};


// Output statistics
struct FMINSearchStats
{
	long iterations;	// Number of ODE function calls.
	long funcCount;		// Number of used steps.
	double fVal;		// Current minimal value of function.
};



//////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Handle to function to minimize: 
// Input: nvars - number of variables, x[..] - input vector, params - pointer to parameters (can be any)
typedef double (*Fcn2Min)(unsigned long nvars, double* x, void* param);



//////////////////////////////////////////////////////////////////////////////////////////////////////////////


//   Minimizes objective function
//   Input:  double funfcn(unsigned n, double* x, void* params) - objective function;
//           nvars - length of x
//           X_in - array of input arguments
//           params - pointer on parameters to be passed to the objective function (can be any) 
//           options - optimization options
//   Output: fval - minimum value of the value of the objective function at Xout.
//           X_out - array of X that minimizes the objective function

bool fminsearch(Fcn2Min funfcn, unsigned long nvars, double* X_in, double* X_out, void* params, FMINSearchOptions& options, FMINSearchStats& searchStats);
