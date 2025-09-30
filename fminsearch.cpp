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
#include <math.h>
#include "fminsearch.h"


//////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool fminsearch(Fcn2Min funfcn, unsigned long nvars, double* X_in, double* X_out, void* params, FMINSearchOptions& options, FMINSearchStats& searchStats)
{
	//   Minimizes objective function
	//   Input:  double funfcn(unsigned n, double* x, void* params) - objective function;
	//           nvars - length of x
	//           X_in - array of input arguments
	//           params - pointer on parameters to be passed to the objective function (can be any) 
	//           options - optimization options
	//   Output: fval - minimum value of the value of the objective function at Xout.
	//           X_out - array of X that minimizes the objective function

	// Return true if FMINSEARCH converged to a solution X.
	//        false if the number of function evaluations or iterations reached.
	//
	//   FMINSEARCH uses the Nelder-Mead simplex (direct search) method.
	//
	//  Reference: Jeffrey C. Lagarias, James A. Reeds, Margaret H. Wright,
	//   Paul E. Wright, "Convergence Properties of the Nelder-Mead Simplex
	//   Method in Low Dimensions", SIAM Journal of Optimization, 9(1):
	//   p.112-147, 1998.

	// defaultopt = struct('Display','notify','MaxIter','200*numberOfVariables',...
	// 'MaxFunEvals','200*numberOfVariables','TolX',1e-4,'TolFun',1e-4, ...'FunValCheck','off','OutputFcn',[]);

	unsigned long n = nvars;
	double* tolx = options.TolX;
	double tolf = options.TolFun;
	unsigned long maxfun = options.MaxFunEvals;
	unsigned long maxiter = options.MaxIter;
	double* pWorkspace = options.pWorkspace;

	if (pWorkspace == NULL) pWorkspace = new double[(n+1)*(n+1)+7*n+2];


	// Initialize parameters
	const double rho = 1.0;
	const double chi = 2.0;
	const double psi = 0.5; 
	const double sigma = 0.5;

	// Set up a simplex near the initial guess.
	// xin = x(:); % Force xin to be a column vector
	double* v = &pWorkspace[0];             // array of size n*(n+1); linear memory will be treated here column-vise {v[:,1],v[:,2], ... }
	for (unsigned long i=0; i<(n*(n+1)); i++) v[i] = 0.0;

	double* fv = &pWorkspace[n*(1+n)];      // aray of size n+1
	for (unsigned long i=0; i<n+1; i++) fv[i] = 0.0;

	double* x =    &pWorkspace[(n+1)*(n+1)+1];      // aray of size n
	double* y =    &pWorkspace[(n+1)*(n+1)+n+1];    // aray of size n
	double* xbar = &pWorkspace[(n+1)*(n+1)+2*n+1];	// aray of size n
	double* xr =   &pWorkspace[(n+1)*(n+1)+3*n+1];	// aray of size n
	double* xe =   &pWorkspace[(n+1)*(n+1)+4*n+1];	// aray of size n
	double* xc =   &pWorkspace[(n+1)*(n+1)+5*n+1];	// aray of size n
	double* xcc =  &pWorkspace[(n+1)*(n+1)+6*n+1];	// aray of size n


	for (unsigned long i=0; i<n; i++) v[i] = X_in[i]; // Place input guess in the simplex! (credit L.Pfeffer at Stanford)
	fv[0] = funfcn(n,X_in,params);

	unsigned long func_evals = 1;
	unsigned long itercount = 0;

	// Following improvement suggested by L.Pfeffer at Stanford
	const double usual_delta = 0.05;			// 5 percent deltas for non-zero terms
	const double zero_term_delta = 0.00025;	// Even smaller delta for zero elements of x

	for (unsigned long j=0; j<n; j++)
	{
		for (unsigned long i=0; i<n; i++) y[i]=X_in[i];
		if (y[j] != 0)
		{
			y[j] = (1 + usual_delta)*y[j];
		}
		else
		{
			y[j] = zero_term_delta;
		}
		for (unsigned long i=0; i<n; i++) v[i+n*(j+1)] = y[i];
		for (unsigned long i=0; i<n; i++) x[i] = y[i]; 
		fv[j+1] = funfcn(n, x, params);
	}

	// sort so v(1,:) has the lowest function value
	// [fv,j] = sort(fv);
	// v = v(:,j);

	// use bubble method to sort
	for (unsigned long i=0; i<n; i++)
	{
		for (unsigned long j=0; j<n; j++)
		{
			double tmp = fv[j+1];
			if (tmp<fv[j])
			{
				// swap
				fv[j+1] = fv[j];
				fv[j] = tmp;

				for (unsigned long k=0; k<n; k++)
				{
					tmp = v[k+n*(j+1)];
					v[k+n*(j+1)] = v[k+n*j];
					v[k+n*j] = tmp;
				}

			}
		}
	}

	itercount = itercount + 1;
	func_evals = n+1;


	// Main algorithm
	// Iterate until the diameter of the simplex is less than tolx
	// AND the function values differ from the min by less than tolf,
	// or the max function evaluations are exceeded. (Cannot use OR instead of AND.)
	
	while (func_evals < maxfun && itercount < maxiter)
	{
		// check iteration break conditions 

		bool breakFlag = true;

		double maxf = 0.0;
		for (unsigned long i=0; i<n; i++)
		{
			double tmp = abs(fv[0]-fv[i+1]);
			if (maxf<tmp) maxf = tmp;
		}
		if (maxf > tolf) breakFlag = false;

		for (unsigned long i=0; i<n; i++)
		{
			double maxv = 0.0;
			for (unsigned long j=0; j<n; j++)
			{
				double tmp = abs(v[i+n*(j+1)]-v[i]);
				if (maxv<tmp) maxv = tmp;
			}
			if (maxv > tolx[i])  breakFlag = false;
		}

        if (breakFlag) break;


    
		// -----------------------------
		// Compute the reflection point
    
		// xbar = average of the n (NOT n+1) best points

		for (unsigned long i=0; i<n; i++)
		{
			double tmp = 0.0;
			for (unsigned long j=0; j<n; j++) tmp += v[i+n*j];
			xbar[i] = tmp/((double)n);
		}

		for (unsigned long i=0; i<n; i++) xr[i] = (1 + rho)*xbar[i] - rho*v[i+n*n];
		for (unsigned long i=0; i<n; i++) x[i] = xr[i];
		double fxr = funfcn(n, x, params);
		func_evals = func_evals+1;
    
		if (fxr < fv[0])
		{
			// Calculate the expansion point
			for (unsigned long i=0; i<n; i++) xe[i] = (1.0 + rho*chi)*xbar[i] - rho*chi*v[i+n*n];
			for (unsigned long i=0; i<n; i++)  x[i] = xe[i];
			double fxe = funfcn(n, x, params);
			func_evals++;
			if (fxe < fxr)
			{
				// expand
				for (unsigned long i=0; i<n; i++) v[i+n*n]=xe[i];
				fv[n] = fxe;
			}
			else
			{
				// reflect
				for (unsigned long i=0; i<n; i++) v[i+n*n]=xr[i];
				fv[n] = fxr;
			}
		}
		else // fv(:,1) <= fxr
		{
			if (fxr < fv[n-1])
			{
				// reflect
				for (unsigned long i=0; i<n; i++) v[i+n*n]=xr[i];
				fv[n] = fxr;
			}
			else // fxr >= fv(:,n)
			{
				//Perform contraction
				if (fxr < fv[n])
				{
					// Perform an outside contraction
					for (unsigned long i=0; i<n; i++) xc[i] = (1.0 + psi*rho)*xbar[i] - psi*rho*v[i+n*n];
					for (unsigned long i=0; i<n; i++) x[i] = xc[i]; 
					double fxc = funfcn(n, x, params);
					func_evals++;
                
					if (fxc <= fxr)
					{
						// contract outside'
	                   for (unsigned long i=0; i<n; i++) v[i+n*n] = xc[i];
					   fv[n] = fxc;
					}
					else
					{
						// perform a shrink
						for (unsigned long j=1; j<=n; j++)
						{
							for (unsigned long i=0; i<n; i++) v[i+j*n]=v[i]+sigma*(v[i+j*n] - v[i]);
							for (unsigned long i=0; i<n; i++) x[i] = v[i+j*n];
							fv[j] = funfcn(n, x, params);
						};
						func_evals += n;
					}
				}
				else
				{
					// Perform an inside contraction
					for (unsigned long i=0; i<n; i++) xcc[i] = (1-psi)*xbar[i] + psi*v[i+n*n];
					for (unsigned long i=0; i<n; i++) x[i] = xcc[i]; 
					double fxcc = funfcn(n, x, params);
					func_evals++;
                
					if (fxcc < fv[n])
					{
						// contract inside
						for (unsigned long i=0; i<n; i++) v[i+n*n]=xcc[i];
						fv[n] = fxcc;
					}
					else
					{
						// perform a shrink
						for (unsigned long j=1; j<=n; j++)
						{
							for (unsigned long i=0; i<n; i++) v[i+j*n]=v[i]+sigma*(v[i+j*n] - v[i]);
							for (unsigned long i=0; i<n; i++) x[i] = v[i+j*n];
							fv[j] = funfcn(n, x, params);
						};
						func_evals += n;
					}
				}
			}
		}

//		[fv,j] = sort(fv);
//		v = v(:,j);

		// use bubble method to sort
		for (unsigned long i=0; i<n; i++)
		{
			for (unsigned long j=0; j<n; j++)
			{
				double tmp = fv[j+1];
				if (tmp<fv[j])
				{
					// swap
					fv[j+1] = fv[j];
					fv[j] = tmp;

					for (unsigned long k=0; k<n; k++)
					{
						tmp = v[k+n*(j+1)];
						v[k+n*(j+1)] = v[k+n*j];
						v[k+n*j] = tmp;
					}

				}
			}
		}

		itercount++;

	}; //while

	for (unsigned long i=0; i<n; i++) X_out[i] = v[i];

	searchStats.iterations = itercount;
	searchStats.funcCount = func_evals;

	double minval = fv[0];
	for (unsigned long i=1; i<=n; i++) if(minval>fv[i]) minval = fv[i];
	searchStats.fVal = minval;

	// free workspace
	if (pWorkspace!=options.pWorkspace) delete pWorkspace;

	if (func_evals >= maxfun)
	{
		if (options.outstream)
		{
			fprintf(options.outstream, "Exiting: Maximum number of function evaluations has been exceeded.\n");
			fprintf(options.outstream, "Increase MaxFunEvals option.\n");
			fprintf(options.outstream, "Current function value: %f \n", minval);
		}
		return false;
	}

	if (itercount >= maxiter)
	{
		if (options.outstream)
		{
			fprintf(options.outstream, "Exiting: Maximum number of iterations has been exceeded.\n");
			fprintf(options.outstream, "Increase MaxIter option.\n");
			fprintf(options.outstream, "Current function value: %f \n", minval);
		}
		return false;
	}

	if (options.outstream)
	{
		fprintf(options.outstream, "Optimization successfully finished:\n");
		fprintf(options.outstream, "the current x satisfies the termination criteria OPTIONS.TolX\n");
		fprintf(options.outstream, "and F(X) satisfies the convergence criteria OPTIONS.TolFun\n");
		fprintf(options.outstream, "Function value: %f \n", minval);
	}
    return true;
}
