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


////////////////////////////////////////////////////////////////////////////////////////
// Class to handle environmental:
// 1). Meteorological 
// 2). Magnetic declination
////////////////////////////////////////////////////////////////////////////////////////

// Constructor: various initializations to handle environmental data

CMH370EnvironData::CMH370EnvironData()
{
	// Meteorology:
	meteo_dataset_ = METEO_DATASET_NULL;	// (not initialized; GDAS1, or ERA5)

	meteo_Start_Time_ = 0.0;	// dataset start time since 2014-03-07 00:00:00 UTC
	meteo_dt_ = 0.0;			// time step (seconds)

	// dimensions:
	meteo_nt_ = 0;				// number of time steps 
	meteo_nx_ = 0;				// number of grid nodes in lon-direction
	meteo_ny_ = 0;				// number of grid nodes in lat-direction
	meteo_nz_ = 0;				// number of layers

	// bounding box
	meteo_LON_MIN_ = 0.0;		// western boundary longitude
	meteo_LON_MAX_ = 0.0;		// eastern boundary longitude
	meteo_LAT_MIN_ = 0.0;		// southern boundary latitude
	meteo_LAT_MAX_ = 0.0;		// northern boundary latitude

	// arrays of pointers on data [NumTs,nz,nx,ny]
	pWIND_U_ = NULL;			// Wind U-component (W->E, m/s)
	pWIND_V_ = NULL;			// Wind V-component (S->N, m/s)
	pLAYERH_ = NULL;			// Geopotential height (m)
	pTEMPER_ = NULL;			// Air temperature (deg C)
	pRHUMID_ = NULL;			// Relative humidity (%)
	pPressureLevels_ = NULL;	// Pressure levels

	meteo_interp_method_ = 3;	// default tricubic (horizontal x time) and spline (vertical)

	p_meteo_workspace_ = NULL;	// memory workspace for interpolations (to avoid dynamic allocations)

	// default values to return if data cannot be loaded:
	default_u_ = 0.0;			// m/s
	default_v_ = 0.0;			// m/s
	default_t_ = 15.0;			// deg C 
	default_p_ = 101325.0;		// Pa
	default_r_ = 0.0;			// %

	// ------------------------------------------------------
	// Magnetic declination (deg)
	magn_nx_ = 0;				// number of grid nodes in lon-direction
	magn_ny_ = 0;				// number of grid nodes in lat-direction
	pMagneticDecl_ =NULL;		// magnetic declination (deg)
	decl_interp_method_ = 1;	// available methods: 0 - bilinear; 1 - bicubic (default)


	// ------------
	// Verbose diganostic stream
	fileout_ = NULL;			 // message stream for diagnostic output (can be NULL) 
}


////////////////////////////////////////////////////////////////////////////////////////
// Destructor: deallocate all arrays
////////////////////////////////////////////////////////////////////////////////////////
CMH370EnvironData::~CMH370EnvironData()
{
	// deallocate U-arrays
	if (pWIND_U_)
	{
		for (int n=0; n<meteo_nt_; n++)
		{
			if (pWIND_U_[n])
			{
				for (int k=0; k<meteo_nz_; k++)
				{
					if (pWIND_U_[n][k]) delete pWIND_U_[n][k];
				}
				delete pWIND_U_[n];
			}
		}
	}

	// deallocate V-arrays
	if (pWIND_V_)
	{
		for (int n=0; n<meteo_nt_; n++)
		{
			if (pWIND_V_[n])
			{
				for (int k=0; k<meteo_nz_; k++)
				{
					if (pWIND_V_[n][k]) delete pWIND_V_[n][k];
				}
				delete pWIND_V_[n];
			}
		}
	}


	// deallocate geopotential height arrays
	if (pLAYERH_)
	{
		for (int n=0; n<meteo_nt_; n++)
		{
			if (pLAYERH_[n])
			{
				for (int k=0; k<meteo_nz_; k++)
				{
					if (pLAYERH_[n][k]) delete pLAYERH_[n][k];
				}
				delete pLAYERH_[n];
			}
		}
	}

	// deallocate temperature arrays
	if (pTEMPER_)
	{
		for (int n=0; n<meteo_nt_; n++)
		{
			if (pTEMPER_[n])
			{
				for (int k=0; k<meteo_nz_; k++)
				{
					if (pTEMPER_[n][k]) delete pTEMPER_[n][k];
				}
				delete pTEMPER_[n];
			}
		}
	}

	// deallocate humidity arrays
	if (pRHUMID_)
	{
		for (int n=0; n<meteo_nt_; n++)
		{
			if (pRHUMID_[n])
			{
				for (int k=0; k<meteo_nz_; k++)
				{
					if (pRHUMID_[n][k]) delete pRHUMID_[n][k];
				}
				delete pRHUMID_[n];
			}
		}
	}

	// deallocate pressure levels array
	if (pPressureLevels_) delete pPressureLevels_;

	// deallocate memory workspace used for interpolations
	if (p_meteo_workspace_) delete p_meteo_workspace_;

	// deallocate magnetic declination array
	if (pMagneticDecl_) delete pMagneticDecl_;

}


////////////////////////////////////////////////////////////////////////////////////////
// Overwrite default constant values to be returned by GetMeteo if data loading fails.
////////////////////////////////////////////////////////////////////////////////////////
void CMH370EnvironData::SetDefault(double u, double v, double t, double p, double r)
{
	// Input
	// u - wind x component (m/s)
	// v - wind y component (m/s)
	// t - air temperature (deg C)
	// p - air pressure (Pa)
	// r - relative humidity (%)

	default_u_ = (float)u;
	default_v_ = (float)v;
	default_t_ = (float)t;
	default_p_ = (float)p;
	default_r_ = (float)r;

	// sanity check
	if (r>100.0) default_r_ = 100.0;
	if (r<0.0) default_r_  =0.0;

}


////////////////////////////////////////////////////////////////////////////////////////
// Set 4-D interpolation method for meteorology.
////////////////////////////////////////////////////////////////////////////////////////

// Currently supported methods:
// 0 - linear in time, bilinear horizontally, linear vertically
// 1 - tricubic in time x horizontally, linear vertically
// 2 - linear in time, bilinear horizontally, cubic spline vertically
// 3 - tricubic in time x horizontally, cubic spline vertically

void CMH370EnvironData::SetMeteoInterpMethod(int method)
{
	if (!(method==0 || method==1 || method==2 || method==3))
	{
		if (fileout_) fprintf(fileout_,"Invalid interpolation method for meteorology.\n");
		return;
	}

	meteo_interp_method_ = method;
}

////////////////////////////////////////////////////////////////////////////////////////
// Set interpolation method for magnetic declination.
////////////////////////////////////////////////////////////////////////////////////////

// Currently supported methods:
// 0 - bilinear
// 1 - bicubic

void CMH370EnvironData::SetDeclInterpMethod(int method)
{
	if (!(method==0 || method==1))
	{
		if (fileout_) fprintf(fileout_,"Invalid interpolation method for magnetic declination.\n");
		return;
	}

	decl_interp_method_ = method;
}


////////////////////////////////////////////////////////////////////////////////////////
//  Verbose for diagnostic (e.g., if data loading fails)
////////////////////////////////////////////////////////////////////////////////////////

void CMH370EnvironData::SetVerbose(FILE* outputstream, int verboselevel)
{
	// set verbose level (0 - silent, 1 - make diagnostic output during data load)
	fileout_ = outputstream;
	verboseLevel_ = verboselevel;
}
