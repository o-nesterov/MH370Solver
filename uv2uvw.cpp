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
#include <math.h>

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Convert horizontal velocity components {u,v} in the local W->E, S->N system into ECEF 
// velocity components {U,V,W} in WGS'84 ellipsoid system at the surface level
// Note that u^2+v^2 == U^2+V^2+W^2
////////////////////////////////////////////////////////////////////////////////////////////////////////////

void uv2uvw(double& XYZ_U, double& XYZ_V, double& XYZ_W, double u, double v, double lon, double lat)
{
	// Input: u,v,lon, lat - u,v, longitude, and latitude
	// Output: XYZ_U,XYZ_V,XYZ_W in (X,Y,Z) ECEF system

	const double Earth_Rp = 6356752.314245; // polar radius
	const double Earth_Re = 6378137.0; // equatorial radius
	const double Rad = 3.1415926535897932384626433832795/180.0; // pi/180

	// --------------------------------------------
	double coslat = cos(Rad*lat);
	double sinlat = sin(Rad*lat);
	double coslon = cos(Rad*lon);
	double sinlon = sin(Rad*lon);

	double R = Earth_Rp*Earth_Re/sqrt((Earth_Rp*coslat)*(Earth_Rp*coslat) + (Earth_Re*sinlat)*(Earth_Re*sinlat));
	double mu = R*R*(Earth_Re*Earth_Re - Earth_Rp*Earth_Rp)/(Earth_Rp*Earth_Rp)/(Earth_Re*Earth_Re)*sinlat*coslat;

	double coef = 1.0 / sqrt(1.0+mu*mu);

	XYZ_U = -u*sinlon - v*coef*(sinlat+mu*coslat)*coslon;
	XYZ_V =  u*coslon - v*coef*(sinlat+mu*coslat)*sinlon;
	XYZ_W =             v*coef*(coslat-mu*sinlat);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Convert horizontal velocity components {u,v} in the local W->E, S->N system into spherical ECEF 
// velocity components {U,V,W}
// Note that u^2+v^2 == U^2+V^2+W^2
////////////////////////////////////////////////////////////////////////////////////////////////////////////

void uv2uvw_sph(double& XYZ_U, double& XYZ_V, double& XYZ_W, double u, double v, double lon, double lat)
{
	// Input: u,v,lon, lat - u,v, longitude, and latitude
	// Output: XYZ_U,XYZ_V,XYZ_W in (X,Y,Z) ECEF sperical system

	const double Rad = 3.1415926535897932384626433832795/180.0; // pi/180

	// --------------------------------------------
	double coslat = cos(Rad*lat);
	double sinlat = sin(Rad*lat);
	double coslon = cos(Rad*lon);
	double sinlon = sin(Rad*lon);

	XYZ_U = -u*sinlon - v*sinlat*coslon;
	XYZ_V =  u*coslon - v*sinlat*sinlon;
	XYZ_W =             v*coslat;
}