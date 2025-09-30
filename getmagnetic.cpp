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
//
// Class to handle environmental and Inmarsat data - get magnetic declination
//
////////////////////////////////////////////////////////////////////////////////////////

extern float bicubicInterpolate (float p[4][4], double a, double b) ;

// Returns magnetic declination (deg) at given location longitude lon (deg E), latitutde (deg N)
// This is the direction where compass norh points to with respect to the true north

double CMH370EnvironData::GetMagneticDeclination(double lon, double lat)
{
	if (pMagneticDecl_ == NULL) return 0.0; // magnetic declination data was not loaded

	int nx = magn_nx_; // number of cells in longitudinal direction
	int ny = magn_ny_; // number of cells in latitudinal direction

	double dlon = (magn_LON_MAX_-magn_LON_MIN_)/(nx-1); // 0.25 deg
	double dlat = (magn_LAT_MAX_-magn_LAT_MIN_)/(ny-1);

	int i0 = (int)floor((lon-magn_LON_MIN_)/dlon); // nearest integer index
	int j0 = (int)floor((lat-magn_LAT_MIN_)/dlat);

	double kx = (lon - magn_LON_MIN_ - i0*dlon)/dlon; // grid fraction between 0 and 1
	double ky = (lat - magn_LAT_MIN_ - j0*dlat)/dlat;

	// 2D space interpolation
	if (i0>=1 && i0+2<=nx && j0>=1 && j0+2<=ny && decl_interp_method_ == 1)
	{
		// bicubic interpolation
		//  4 (lon) x 4 (lat) chunks of data
		float MAGN_DECL_chunk[4][4];

		for (int i=0; i<4; i++)
		{
			for (int j=0; j<4; j++)
			{
				MAGN_DECL_chunk[i][j] = pMagneticDecl_[(j0+j-1)*nx+(i0+i-1)];
			}
		}

		return (double)bicubicInterpolate(MAGN_DECL_chunk, kx, ky);
	}
	else
	{
		// fall back to bilinear interpolation-extrapolation

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

		double vals1 = pMagneticDecl_[j0*nx+i0];
		double vals2 = pMagneticDecl_[j0*nx+i1];
		double vals3 = pMagneticDecl_[j1*nx+i0];
		double vals4 = pMagneticDecl_[j1*nx+i1];

		return ((1.0-ky)*(1.0-kx)*vals1 + (1.0-ky)*kx*vals2 + ky*(1.0-kx)*vals3 + kx*ky*vals4);
	}

}
