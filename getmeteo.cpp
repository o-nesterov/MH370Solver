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
#include <string.h>
#include "environdata.h"

float cubicInterpolate (float p[4], double a);
float bicubicInterpolate (float p[4][4], double a, double b) ;
float tricubicInterpolate (float p[4][4][4], double a, double b, double c);
void splineInterpMeteo(double* pY, double* pInpX, double** pInpY, double X, unsigned long n, unsigned long nvars, double* pWorkspace);

////////////////////////////////////////////////////////////////////////////////////////
// Class to handle environmental and Inmarsat data
////////////////////////////////////////////////////////////////////////////////////////


// This function performs 4D interpolation and returns meteorological data at given time and position

bool CMH370EnvironData::GetMeteo(double& u, double& v, double& t, double& p, double& r,
								 double ts, double lon, double lat, double alt)
{
	if (meteo_dataset_ == METEO_DATASET_NULL || pWIND_U_ == NULL || pWIND_V_ == NULL || pLAYERH_ == NULL || pTEMPER_ == NULL || pRHUMID_ == NULL)
	{
		u = default_u_;
		v = default_v_;
		t = default_t_;
		p = default_p_;
		r = default_r_;
		return false;
	}

	int nx = meteo_nx_; // number of cells in longitudinal direction
	int ny = meteo_ny_; // number of cells in latitudinal direction
	int nz = meteo_nz_; // number of vertical layers
	int nt = meteo_nt_; // number of time steps available

	// ------------
	double dlon = (meteo_LON_MAX_-meteo_LON_MIN_)/(nx-1); // 0.25 deg for ERA5, 1 deg for GDAS1
	double dlat = (meteo_LAT_MAX_-meteo_LAT_MIN_)/(ny-1);

	int i0 = (int)floor((lon-meteo_LON_MIN_)/dlon); // nearest integer index
	int j0 = (int)floor((lat-meteo_LAT_MIN_)/dlat);

	double kx = (lon - meteo_LON_MIN_ - i0*dlon)/dlon; // grid fraction between 0 and 1
	double ky = (lat - meteo_LAT_MIN_ - j0*dlat)/dlat;

	// ------------
	// time:
	double ti = (ts - meteo_Start_Time_)/meteo_dt_;
	int ti_int = (int)floor(ti);
	double kt = ti-ti_int; // fraction of timestep

	// column array assignments (instead of memory allocation)
	double* h_column = &p_meteo_workspace_[0];
	double* u_column = &p_meteo_workspace_[nz];
	double* v_column = &p_meteo_workspace_[2*nz];
	double* t_column = &p_meteo_workspace_[3*nz];
	double* r_column = &p_meteo_workspace_[4*nz];

	// ------------
	// apply tricubic interpolation in horizontal-time space to each element of
	// a vertical column if possible; otherise bilinear (at the edges of the domain)
	// use custom interpolant to speed up calculation as the horizontal-time grid spacing is uniform

	if (i0>=1 && i0+2<=nx && j0>=1 && j0+2<=ny && ti_int>=1 && ti_int+2<=nt && (meteo_interp_method_==1 || meteo_interp_method_==3))
	{
		// tricubic interpolation
		//  4 (time) x 4 (lon) x 4 (lat) chunks of data
		float LAYERH_chunk[4][4][4];
		float WIND_U_chunk[4][4][4];
		float WIND_V_chunk[4][4][4];
		float TEMPER_chunk[4][4][4];
		float RHUMID_chunk[4][4][4];

		for (int l=0; l<nz; l++)
		{
			for (int m=0; m<4; m++)
			{
				for (int i=0; i<4; i++)
				{
					for (int j=0; j<4; j++)
					{
						LAYERH_chunk[m][i][j] = pLAYERH_[ti_int-1+m][l][(j0+j-1)*nx+(i0+i-1)];
						WIND_U_chunk[m][i][j] = pWIND_U_[ti_int-1+m][l][(j0+j-1)*nx+(i0+i-1)];
						WIND_V_chunk[m][i][j] = pWIND_V_[ti_int-1+m][l][(j0+j-1)*nx+(i0+i-1)];
						TEMPER_chunk[m][i][j] = pTEMPER_[ti_int-1+m][l][(j0+j-1)*nx+(i0+i-1)];
						RHUMID_chunk[m][i][j] = pRHUMID_[ti_int-1+m][l][(j0+j-1)*nx+(i0+i-1)];
					}
				}
			}

			h_column[l] = (double)tricubicInterpolate(LAYERH_chunk, kt, kx, ky);
			u_column[l] = (double)tricubicInterpolate(WIND_U_chunk, kt, kx, ky);
			v_column[l] = (double)tricubicInterpolate(WIND_V_chunk, kt, kx, ky);
			t_column[l] = (double)tricubicInterpolate(TEMPER_chunk, kt, kx, ky);
			r_column[l] = (double)tricubicInterpolate(RHUMID_chunk, kt, kx, ky);

			// trim humidity, which can exceed 100% or become negative in case of tricubic interpolation
			if (r_column[l]<0.0) r_column[l] = 0.0;
			if (r_column[l]>100.0) r_column[l] = 100.0;
		}
	}
	
	else
	{
		// use bilinear interpolation-extrapolation to avoid failure during minimization procedure

		int i1 = i0 + 1;
		int j1 = j0 + 1;
    
		if (i0 >= nx-1)
		{
			i0 = nx-1;
			i1 = nx-1;
			kx = 0.0;
		};

		if (i1 <= 0)
		{
			i0 = 0;
			i1 = 0;
			kx = 0.0;
		};

		if (j0 >= ny-1)
		{
			j0 = ny-1;
			j1 = ny-1;
			ky = 0.0;
		}

		if (j1 <= 0)
		{
			j0 = 0;
			j1 = 0;
			ky = 0.0;
		}

		int ti_int0 = ti_int;
		int ti_int1 = ti_int+1;

		if (ti_int<0)
		{
			ti_int0 = 0;
			ti_int1 = 0;
			kt = 0.0;
		}

        if (ti_int1 >= meteo_nt_)
		{
			ti_int0 = meteo_nt_-1;
			ti_int1 = meteo_nt_-1;
			kt = 0.0;
		}

		for (int l=0; l<nz; l++)
		{
			// horizontal & time interpolation for H
			float vals1 = pLAYERH_[ti_int0][l][j0*nx+i0];
			float vals2 = pLAYERH_[ti_int0][l][j0*nx+i1];
			float vals3 = pLAYERH_[ti_int0][l][j1*nx+i0];
			float vals4 = pLAYERH_[ti_int0][l][j1*nx+i1];
			double vals_ = (1.0-ky)*(1.0-kx)*vals1 + (1.0-ky)*kx*vals2 + ky*(1.0-kx)*vals3 + kx*ky*vals4;
			vals1 = pLAYERH_[ti_int1][l][j0*nx+i0];
			vals2 = pLAYERH_[ti_int1][l][j0*nx+i1];
			vals3 = pLAYERH_[ti_int1][l][j1*nx+i0];
			vals4 = pLAYERH_[ti_int1][l][j1*nx+i1];
			double vals__ = (1.0-ky)*(1.0-kx)*vals1 + (1.0-ky)*kx*vals2 + ky*(1.0-kx)*vals3 + kx*ky*vals4;
			h_column[l] = (double)((1.0-kt)*vals_ + kt*vals__);


			// horizontal & time interpolation for U
			vals1 = pWIND_U_[ti_int0][l][j0*nx+i0];
			vals2 = pWIND_U_[ti_int0][l][j0*nx+i1];
			vals3 = pWIND_U_[ti_int0][l][j1*nx+i0];
			vals4 = pWIND_U_[ti_int0][l][j1*nx+i1];
			vals_ = (1.0-ky)*(1.0-kx)*vals1 + (1.0-ky)*kx*vals2 + ky*(1.0-kx)*vals3 + kx*ky*vals4;
			vals1 = pWIND_U_[ti_int1][l][j0*nx+i0];
			vals2 = pWIND_U_[ti_int1][l][j0*nx+i1];
			vals3 = pWIND_U_[ti_int1][l][j1*nx+i0];
			vals4 = pWIND_U_[ti_int1][l][j1*nx+i1];
			vals__= (1.0-ky)*(1.0-kx)*vals1 + (1.0-ky)*kx*vals2 + ky*(1.0-kx)*vals3 + kx*ky*vals4;
			u_column[l] = (double)((1.0-kt)*vals_ + kt*vals__);


			// horizontal & time interpolation for V
			vals1 = pWIND_V_[ti_int0][l][j0*nx+i0];
			vals2 = pWIND_V_[ti_int0][l][j0*nx+i1];
			vals3 = pWIND_V_[ti_int0][l][j1*nx+i0];
			vals4 = pWIND_V_[ti_int0][l][j1*nx+i1];
			vals_ = (1.0-ky)*(1.0-kx)*vals1 + (1.0-ky)*kx*vals2 + ky*(1.0-kx)*vals3 + kx*ky*vals4;
			vals1 = pWIND_V_[ti_int1][l][j0*nx+i0];
			vals2 = pWIND_V_[ti_int1][l][j0*nx+i1];
			vals3 = pWIND_V_[ti_int1][l][j1*nx+i0];
			vals4 = pWIND_V_[ti_int1][l][j1*nx+i1];
			vals__= (1.0-ky)*(1.0-kx)*vals1 + (1.0-ky)*kx*vals2 + ky*(1.0-kx)*vals3 + kx*ky*vals4;
			v_column[l] = (double)((1.0-kt)*vals_ + kt*vals__);


			// horizontal & time interpolation for T
			vals1 = pTEMPER_[ti_int0][l][j0*nx+i0];
			vals2 = pTEMPER_[ti_int0][l][j0*nx+i1];
			vals3 = pTEMPER_[ti_int0][l][j1*nx+i0];
			vals4 = pTEMPER_[ti_int0][l][j1*nx+i1];
			vals_ = (1.0-ky)*(1.0-kx)*vals1 + (1.0-ky)*kx*vals2 + ky*(1.0-kx)*vals3 + kx*ky*vals4;
			vals1 = pTEMPER_[ti_int1][l][j0*nx+i0];
			vals2 = pTEMPER_[ti_int1][l][j0*nx+i1];
			vals3 = pTEMPER_[ti_int1][l][j1*nx+i0];
			vals4 = pTEMPER_[ti_int1][l][j1*nx+i1];
			vals__= (1.0-ky)*(1.0-kx)*vals1 + (1.0-ky)*kx*vals2 + ky*(1.0-kx)*vals3 + kx*ky*vals4;
			t_column[l] = (double)((1.0-kt)*vals_ + kt*vals__);


			// horizontal & time interpolation for RH
			vals1 = pRHUMID_[ti_int0][l][j0*nx+i0];
			vals2 = pRHUMID_[ti_int0][l][j0*nx+i1];
			vals3 = pRHUMID_[ti_int0][l][j1*nx+i0];
			vals4 = pRHUMID_[ti_int0][l][j1*nx+i1];
			vals_ = (1.0-ky)*(1.0-kx)*vals1 + (1.0-ky)*kx*vals2 + ky*(1.0-kx)*vals3 + kx*ky*vals4;
			vals1 = pRHUMID_[ti_int1][l][j0*nx+i0];
			vals2 = pRHUMID_[ti_int1][l][j0*nx+i1];
			vals3 = pRHUMID_[ti_int1][l][j1*nx+i0];
			vals4 = pRHUMID_[ti_int1][l][j1*nx+i1];
			vals__= (1.0-ky)*(1.0-kx)*vals1 + (1.0-ky)*kx*vals2 + ky*(1.0-kx)*vals3 + kx*ky*vals4;
			r_column[l] = (double)((1.0-kt)*vals_ + kt*vals__);

		};// l - vertical layers loop

	}; // cubic or linear interpolation


	// ------------------------------------------------------
	// altitude interpolation:
	double* p_column = pPressureLevels_;

	int l=0;
	while (l<nz)
	{
		if (alt<=h_column[l]) break;
		l++;
	}; // end while
                          
	if (l==0)
	{
		u = u_column[0];
		v = v_column[0];
		t = t_column[0];
		r = r_column[0];
		p = p_column[0];
		return false;
	}

	if (l==nz)
	{
		u = u_column[nz-1];
		v = v_column[nz-1];
		t = t_column[nz-1];
		r = r_column[nz-1];
		p = p_column[nz-1];
		return false;
	}

	double tau = double((alt-h_column[l-1])/(h_column[l]-h_column[l-1])); // coef. between 0.0 and 1.0

	if (meteo_interp_method_==3 || meteo_interp_method_==2)
	{
		// spline interpolation along vertical column
		// interpolate u,v,t and r in one round to speed up calculations; alt is a scalar
		// (air pressure is interpolated separately with exponential function)
		// this is equivalent of calling spline() separately for each of these variables
		double* uvtr_column[4];
		double arr[4];

		uvtr_column[0] = u_column;
		uvtr_column[1] = v_column;
		uvtr_column[2] = t_column;
		uvtr_column[3] = r_column;

		splineInterpMeteo(arr, h_column, uvtr_column, alt, nz, 4, &p_meteo_workspace_[5*nz]);

		u = arr[0];	// interpolated wind u - velocity component
		v = arr[1];	// interpolated wind v - velocity component
		t = arr[2];	// interpolated air temperature (C)
		r = arr[3];	// interpolated relative humidity (%)

		if (r<0.0) r = 0.0;
		if (r>100.0) r = 100.0;

	}
	else
	{
		// interpolate linearly 
		u = u_column[l-1] + (u_column[l]-u_column[l-1])*tau;
		v = v_column[l-1] + (v_column[l]-v_column[l-1])*tau;
		t = t_column[l-1] + (t_column[l]-t_column[l-1])*tau;
		r = r_column[l-1] + (r_column[l]-r_column[l-1])*tau;
	}

	// interpolate pressure exponentialy in either case:
	p = p_column[l-1] * exp(log(p_column[l]/p_column[l-1])*tau);
	return true;
};




// -------------------------------------------------------------------------
// A set of 1D, 2D and 3D cubic interpolants of a scalar on a regular grid; a,b,c are scalar coefficients in the range [0...1)
// Adapted from https://www.paulinternet.nl/?page=bicubic

float cubicInterpolate (float p[4], double a) 
{
	double res = p[1] + 0.5 * a*(p[2] - p[0] + a*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + a*(3.0*(p[1] - p[2]) + p[3] - p[0])));
	return (float) res;
}

float bicubicInterpolate (float p[4][4], double a, double b) 
{
	float arr[4];
	arr[0] = cubicInterpolate(p[0], b);
	arr[1] = cubicInterpolate(p[1], b);
	arr[2] = cubicInterpolate(p[2], b);
	arr[3] = cubicInterpolate(p[3], b);
	return cubicInterpolate(arr, a);
}

float tricubicInterpolate (float p[4][4][4], double a, double b, double c) 
{
	float arr[4];
	arr[0] = bicubicInterpolate(p[0], b, c);
	arr[1] = bicubicInterpolate(p[1], b, c);
	arr[2] = bicubicInterpolate(p[2], b, c);
	arr[3] = bicubicInterpolate(p[3], b, c);
	return cubicInterpolate(arr, a);
}



// -------------------------------------------------------------------------

// Customized spline interpolator on the data inpX(:), inpY(:,nvars) of nvars quantities at required X (scalar) 
// This is a customized version based on Matlab spline(.,.,.) function (the same as interp1(...,'spline'))
// inpX and inpY must be of the same length, which must be greater than 3
// pWorkspace is a user-supplied chunk of memory sufficient to accomodate (nvars+3) x n doubles (this is to avoid dynamic allocation/dealloction upon every call)

void splineInterpMeteo(double* pY, double* pInpX, double** pInpY, double X, unsigned long n, unsigned long nvars, double* pWorkspace)
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

	// nvars is the number of variables (u,v,temperature,...)

    double* b = &pWorkspace[0]; // RHS, must be of size n x nvars 
    double* diagu = &pWorkspace[n*nvars]; // upper diagonal
    double* diagc = &pWorkspace[n*nvars+n]; // central diagonal
    double* diagl = &pWorkspace[n*nvars+2*n]; // lower diagonal
    double* dx = &pWorkspace[n*nvars+3*n]; // dx

	for (unsigned long j=0; j+1<n; j++)
	{
		dx[j] = pInpX[j+1]-pInpX[j]; // note the actual length of dx is n-1
	}

    for (unsigned long j = 2; j<=n-1; j++)
	{
		for (unsigned long i = 0; i<nvars; i++)
		{
			b[(j-1)*nvars+i] = 3.0 * (dx[j-1] * ((pInpY[i][j-1] - pInpY[i][j-2]) / dx[j-2]) + dx[j-2] * ((pInpY[i][j] - pInpY[i][j-1]) / dx[j-1]));
		}
	}
    
	for (unsigned long i = 0; i<nvars; i++)
	{
//		b[i]  = ((dx[0] + 2.0 * (pInpX[2] - pInpX[0]))) * dx[1] * ((pInpY[i][1] - pInpY[i][0]) / dx[0]) + 
//			    (dx[0]*dx[0]) * ((pInpY[i][2] - pInpY[i][1]) / dx[1]) / (pInpX[2] - pInpX[0]);
		b[i]  = ((dx[0] + 2.0 * (pInpX[2] - pInpX[0])) * dx[1] * ((pInpY[i][1] - pInpY[i][0]) / dx[0]) +
                (dx[0]*dx[0]) * ((pInpY[i][2] - pInpY[i][1]) / dx[1])) / (pInpX[2] - pInpX[0]);

		b[(n-1)*nvars+i] = ((dx[n-2]*dx[n-2]) * ((pInpY[i][n-2] - pInpY[i][n-3]) / dx[n-3]) + (2.0 * (pInpX[n-1] - pInpX[n-3])
                + dx[n-2]) * dx[n-3] * ((pInpY[i][n-1] - pInpY[i][n-2]) / dx[n-2])) / (pInpX[n-1] - pInpX[n-3]);
	}


	for (unsigned long j = 2; j<=n-1; j++)
	{
        diagc[j-1] = 2.0 * (dx[j-1] + dx[j-2]);
        diagl[j-1] = dx[j-1];
        diagu[j-1] = dx[j-2];
	}

	diagc[0] = dx[1];
    diagu[0] = pInpX[2] - pInpX[0];
    diagc[n-1] = dx[n-3];
    diagl[n-1] = pInpX[n-1] - pInpX[n-3];   
    

	// Solve:
	// diagc(1)    diagu(1)     0            0           0   ...                                           | b(1)               
	// diagl(2)    diagc(2)     diagu(2)     0           0   ...                                           | b(2)              
	// 0           diagl(3)     diagc(3)     diagu(3)    0   ...                                           | b(3)              
	// 0           0            diagl(4)     diagl(5)    0   ...                                           | b(4)              
	// ...                                                                                                 | ... 
	// 0           0           0            0                      diagl(N-1)  diagc(N-1)   diagu(N-1)     | b(N-1) 
	// 0           0           0            0           0    ...               diagl(N)     diagc(N)       | b(N)  


	// Forward iterations
	for (unsigned long j = 1; j<=n-1; j++)
	{ 
        double c = 1.0 / diagc[j-1];
		for (unsigned long i=0; i<nvars; i++) b[(j-1)*nvars+i] *= c;
		diagc[j-1] = 1.0;
		diagu[j-1] = diagu[j-1] * c;

        c = diagl[j];
		for (unsigned long i=0; i<nvars; i++) b[j*nvars+i] -= (c*b[(j-1)*nvars+i]);
        diagc[j] -=  (c * diagu[j-1]);
	}
    
    for (unsigned long i=0; i<nvars; i++) b[(n-1)*nvars+i] /= diagc[n-1];
    diagc[n-1] = 1.0;
    
	
	// Backward iterations; main diagonals diagc is all ones here
    for (unsigned long j = n-1; j>=1; j--)
	{
		double c = diagu[j-1];
		for (unsigned long i=0; i<nvars; i++) b[(j-1)*nvars+i] -= (c * b[j*nvars+i]);
	}

	// ---------------------
	// interpolate
	// find index j that inpX[j-1] <= X < inpX[j]
	unsigned long  L = 0;

	if (X < pInpX[0])
	{
        for (unsigned long i=0; i<nvars; i++) pY[i] = pInpY[0][i];
        return;
	}
	else
	{
		while (X >= pInpX[L])
		{
			L++;
			if (L>=n) break;
		}
	}

	if (L==n)
	{
		for (unsigned long i=0; i<nvars; i++) pY[(n-1)*nvars+i] = pInpY[n-1][i];
        return;
	};

	double h = X - pInpX[L-1];
    
	// calculate spline coefficients
	for (unsigned long i=0; i<nvars; i++)
	{
		double dzzdx = (pInpY[i][L] - pInpY[i][L-1]) / (dx[L-1]* dx[L-1]) - b[(L-1)*nvars+i] / dx[L-1];
		double dzdxdx = b[L*nvars+i] / dx[L-1] - (pInpY[i][L] - pInpY[i][L-1]) / (dx[L-1]*dx[L-1]);
		double coefs1 = (dzdxdx - dzzdx) / dx[L-1];
		double coefs2 = (2.0 * dzzdx - dzdxdx);
		double coefs3 = b[(L-1)*nvars+i];
		double coefs4 = pInpY[i][L-1];
        
		pY[i] = coefs4 + h * (coefs3 + h * (coefs2 + h * coefs1));
	}
};


