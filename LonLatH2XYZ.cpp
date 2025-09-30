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
// convert geocentric longitude lon, latitdue lat and geodetic height (perpendicular to ellipsoid's surface)
// in WGS'84 ellipsoid system into ECEF Cartesian coordinates (X,Y,Z)
////////////////////////////////////////////////////////////////////////////////////////////////////////////
void LonLatH2XYZ(double &X, double &Y, double &Z, double lon, double lat, double hgt)
{
	// Input: triplets of [longitude, latitude, height]
	// Output: triplets of [X, Y, Z]

	// Earth radii, WGS'84 ellipsoid, m:
	const double Earth_Rp = 6356752.314245; // polar radius
	const double Earth_Re = 6378137.0; // equatorial radius
	const double Rad = 3.1415926535897932384626433832795/180.0; // pi/180

	double coslon = cos(Rad*lon);
	double sinlon = sin(Rad*lon);
	double coslat = cos(Rad*lat);
	double sinlat = sin(Rad*lat);

	// Vector normal to the surface {nx,ny,nz}
	double R = Earth_Rp*Earth_Re/sqrt((Earth_Rp*coslat)*(Earth_Rp*coslat) + (Earth_Re*sinlat)*(Earth_Re*sinlat));
	double mu = R*R*(Earth_Re*Earth_Re - Earth_Rp*Earth_Rp)/(Earth_Rp*Earth_Rp)/(Earth_Re*Earth_Re)*coslat*sinlat;
	double coef = 1.0/sqrt(1.0+mu*mu);
	double nx = coef*(coslat-mu*sinlat)*coslon;
	double ny = coef*(coslat-mu*sinlat)*sinlon;
	double nz = coef*(sinlat+mu*coslat);

	X = R*coslat*coslon + nx*hgt;
	Y = R*coslat*sinlon + ny*hgt;
	Z = R*sinlat        + nz*hgt;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////
// convert geocentric longitude lon, latitdue lat and height (perpendicular to sphere's surface)
// in spherical system into ECEF Cartesian coordinates (X,Y,Z)
////////////////////////////////////////////////////////////////////////////////////////////////////////////
void LonLatH2XYZ_sph(double &X, double &Y, double &Z, double lon, double lat, double hgt)
{
	// Input: triplets of [longitude, latitude, height]
	// Output: triplets of [X, Y, Z]

	// Earth radius (m):
	const double Earth_R = 6378137.0; //6378137.0; // equatorial radius
	const double Rad = 3.1415926535897932384626433832795/180.0; // pi/180

	double coslon = cos(Rad*lon);
	double sinlon = sin(Rad*lon);
	double coslat = cos(Rad*lat);
	double sinlat = sin(Rad*lat);

	// Vector normal to the surface {nx,ny,nz}
	double nx = coslat*coslon;
	double ny = coslat*sinlon;
	double nz = sinlat;

	X = (Earth_R + hgt) * nx;
	Y = (Earth_R + hgt) * ny;
	Z = (Earth_R + hgt) * nz;
}
