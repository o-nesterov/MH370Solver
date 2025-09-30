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


// Spline interpolator on the data inpX(:), inpY(:) at required X(:) 
// Based on Matlab spline(.,.,.) function (the same as interp1(...,'spline'))
// inpX and inpY must be of the same length n_in, which must be greater than 3
// pX is the array of X of the lenghth n_out, where interpolation is required;
// pY is the buffer of the lenght n_out to store result

void splineInterp(double* pY, double* pX, unsigned long n_out, double* pInpX, double* pInpY, unsigned long n_in)
{
	// Spline[i] = f[i] + b[i]*(x - x[i]) + c[i]*(x - x[i])^2 + d[i]*(x - x[i])^3

	// Find coefficients by solving linear system A*s=b, where
	// A[N, N] - Tridiagonal Matrix:
	// diagc(1)    diagu(1)     0            0           0   ...                                           | b(1)               
	// diagl(2)    diagc(2)    diagu(2)      0           0   ...                                           | b(2)              
	// 0           diagl(3)    diagc(3)    diagu(3)      0   ...                                           | b(3)              
	// 0           0           diagl(4)    diagl(5)      0   ...                                           | b(4)              
	// ...                                                                                                 | ...                                                                  
	// 0           0           0            0                    diagl(N-1)   diagc(N-1)   diagu(N-1)      | b(N-1) 
	// 0           0           0            0           0    ...              diagl(N)     diagc(N)        | b(N) 

	// Use Gauss algorithm to calculate x = (A^-1)b.
	// Then find coefficients f[i], b[i], c[i], d[i]

	// Sanity check:
	if (n_in <= 3 || n_out<1 || pY==NULL || pX==NULL || pInpX==NULL || pInpY==NULL) return;

	// allocate workspace
	double* pWorkspace = new double[5*n_in+1];

    double* b = &pWorkspace[0];				// RHS, must be of size n x nvars 
    double* diagu = &pWorkspace[n_in];		// upper diagonal
    double* diagc = &pWorkspace[2*n_in];	// central diagonal
    double* diagl = &pWorkspace[3*n_in];	// lower diagonal
    double* dx = &pWorkspace[4*n_in];		// dx

	for (unsigned long j=0; j+1<n_in; j++)
	{
		dx[j] = pInpX[j+1]-pInpX[j];		// note the actual length of dx is n-1
	}

    for (unsigned long j = 2; j<=n_in-1; j++)
	{
			b[j-1] = 3.0 * (dx[j-1] * ((pInpY[j-1] - pInpY[j-2]) / dx[j-2]) + dx[j-2] * ((pInpY[j] - pInpY[j-1]) / dx[j-1]));
	}
    
	b[0]  = ((dx[0] + 2.0 * (pInpX[2] - pInpX[0])) * dx[1] * ((pInpY[1] - pInpY[0]) / dx[0]) +
               (dx[0]*dx[0]) * ((pInpY[2] - pInpY[1]) / dx[1])) / (pInpX[2] - pInpX[0]);

	b[n_in-1] = ((dx[n_in-2]*dx[n_in-2]) * ((pInpY[n_in-2] - pInpY[n_in-3]) / dx[n_in-3]) + (2.0 * (pInpX[n_in-1] - pInpX[n_in-3])
               + dx[n_in-2]) * dx[n_in-3] * ((pInpY[n_in-1] - pInpY[n_in-2]) / dx[n_in-2])) / (pInpX[n_in-1] - pInpX[n_in-3]);

	for (unsigned long j = 2; j<=n_in-1; j++)
	{
        diagc[j-1] = 2.0 * (dx[j-1] + dx[j-2]);
        diagl[j-1] = dx[j-1];
        diagu[j-1] = dx[j-2];
	}

	diagc[0] = dx[1];
    diagu[0] = pInpX[2] - pInpX[0];
    diagc[n_in-1] = dx[n_in-3];
    diagl[n_in-1] = pInpX[n_in-1] - pInpX[n_in-3];   
    

	// Solve:
	// diagc(1)    diagu(1)     0            0           0   ...                                           | b(1)               
	// diagl(2)    diagc(2)     diagu(2)     0           0   ...                                           | b(2)              
	// 0           diagl(3)     diagc(3)     diagu(3)    0   ...                                           | b(3)              
	// 0           0            diagl(4)     diagl(5)    0   ...                                           | b(4)              
	// ...                                                                                                 | ... 
	// 0           0           0            0                      diagl(N-1)  diagc(N-1)   diagu(N-1)     | b(N-1) 
	// 0           0           0            0           0    ...               diagl(N)     diagc(N)       | b(N)  


	// Forward iterations
	for (unsigned long j = 1; j<=n_in-1; j++)
	{ 
        double c = 1.0 / diagc[j-1];
		b[j-1] *= c;
		diagc[j-1] = 1.0;
		diagu[j-1] = diagu[j-1] * c;

        c = diagl[j];
		b[j] -= (c*b[(j-1)]);
        diagc[j] -=  (c * diagu[j-1]);
	}
    
    b[n_in-1] /= diagc[n_in-1];
    diagc[n_in-1] = 1.0;
    
	
	// Backward iterations; main diagonals diagc is all ones here
    for (unsigned long j = n_in-1; j>=1; j--)
	{
		double c = diagu[j-1];
		b[j-1] -= (c * b[j]);
	}

	// ---------------------
	// interpolate

	// iterate for each pX[..]
	for (unsigned long n=0; n<n_out; n++)
	{
		double X = pX[n];

		// find index j that inpX[j-1] <= X < inpX[j]
		unsigned long  L = 0;

		if (X < pInpX[0])
		{
			pY[n] = pInpY[0];
			continue;
		}
		else
		{
			while (X >= pInpX[L])
			{
				L++;
				if (L>=n_in) break;
			}
		}

		if (L==n_in)
		{
			pY[n] = pInpY[n_in-1];
			continue;
		};

		double h = X - pInpX[L-1];
    
		// calculate spline coefficients
		double dzzdx = (pInpY[L] - pInpY[L-1]) / (dx[L-1]* dx[L-1]) - b[L-1] / dx[L-1];
		double dzdxdx = b[L] / dx[L-1] - (pInpY[L] - pInpY[L-1]) / (dx[L-1]*dx[L-1]);
		double coefs1 = (dzdxdx - dzzdx) / dx[L-1];
		double coefs2 = (2.0 * dzzdx - dzdxdx);
		double coefs3 = b[L-1];
		double coefs4 = pInpY[L-1];
        
		pY[n] = coefs4 + h * (coefs3 + h * (coefs2 + h * coefs1));
	};

	delete pWorkspace;
};


