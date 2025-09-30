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
// Convert ECEF Cartesian coordinates (X,Y,Z) to geocentric longitude lon, latitdue lat and 
// geodetic height (perpendicular to ellipsoid's surface) in WGS'84 ellipsoid system
////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Declaration of the auxilary function needed to find latitude by Newton-Raphson iterations 
void findlatfunc(double& funcval, double& derval, double lat, double D2, double Z);


void XYZ2LonLatH(double& lon_out, double& lat_out, double& hgt_out, double X, double Y, double Z)
{
	// Input: triplets of [longitude, latitude, height]
	// Output: triplets of [X, Y, Z]

	// Earth radii, WGS'84 ellipsoid, m:
	const double Earth_Rp = 6356752.314245; // polar radius
	const double Earth_Re = 6378137.0; // equatorial radius
	const double Rad = 3.1415926535897932384626433832795/180.0; // pi/180
	const double RRad = 180.0/3.1415926535897932384626433832795; // 180/pi


	double D2 = X*X + Y*Y + Z*Z;

	// initial guess (approximation) for latitude (degrees)
	double lat = RRad * asin(Z*Earth_Rp/sqrt(Earth_Rp*Earth_Rp*Earth_Re*Earth_Re+Z*Z*Earth_Rp*Earth_Rp-Z*Z*Earth_Re*Earth_Re));

	// params = [Earth_Rp + Z*0.0, Earth_Re + Z*0.0, D2, Z];

    // Newton-Raphson iterations to solve non-linear equation for latitude;
	// According to conducted tests 2-3 iterations are sufficient in the study domain
	for (int iter=0; iter<10; iter++)
	{
		double f, df;
		findlatfunc(f,df,lat,D2,Z);
		lat -= f/df;
    
		if (abs(f)<0.00001) break; // function value is less than 0.00001 m (error in geodetic altitude using found latitude)
	}


	// Calculate altitude in ellipsoidal system:
	double coslat = cos(Rad*lat);
	double sinlat = sin(Rad*lat);

	double R = Earth_Rp*Earth_Re/sqrt((Earth_Rp*coslat)*(Earth_Rp*coslat) + (Earth_Re*sinlat)*(Earth_Re*sinlat));
	double mu = R*R*(Earth_Re*Earth_Re - Earth_Rp*Earth_Rp)/(Earth_Rp*Earth_Rp)/(Earth_Re*Earth_Re)*coslat*sinlat;
	double coef = 1.0 / sqrt(1.0+mu*mu);
	// find altitude from the relation (h+R*coef)^2 + R^2*mu^2/(1+mu^2) = X^2+Y^2+Z^2
	double hgt = sqrt(D2-R*R*mu*mu/(1.0+mu*mu))-R*coef;


	// find longitude from equation:
	// X = R*cosd(lat)*cosd(lon) + nx*alt == R*cosd(lat)*cosd(lon) + alt*coef*(cosd(lat)-mu*sind(lat))*cosd(lon) ==
	// {R*cosd(lat)+alt*coef*(cosd(lat)-mu*sind(lat))}*cosd(lon)

	double lon = RRad*acos(X/(R*coslat+hgt*coef*(coslat-mu*sinlat)));
	if (Y < 0.0) lon = (360.0-lon); // change to 360-... if location is at the other side of Earth

	lon_out = lon;
	lat_out = lat;
	hgt_out = hgt;

}


// ------------------------------------------------------------------------
// Function and its derivative to find latitude by Newton-Raphson iterations

void findlatfunc(double& funcval, double& derval, double lat, double D2, double Z)
{
	// Input: lat - column of latitudes
	// Output: function to be zero and its derivative [funval, derivative]

	// Earth radii, WGS'84 ellipsoid, m:
	const double Earth_Rp = 6356752.314245; // polar radius
	const double Earth_Re = 6378137.0; // equatorial radius
	const double Rad = 3.1415926535897932384626433832795/180.0; // pi/180

	double coslat = cos(Rad*lat);
	double sinlat = sin(Rad*lat);

    double R = Earth_Rp*Earth_Re/sqrt((Earth_Rp*coslat)*(Earth_Rp*coslat) + (Earth_Re*sinlat)*(Earth_Re*sinlat));
	double mu = R*R*(Earth_Re*Earth_Re - Earth_Rp*Earth_Rp)/(Earth_Rp*Earth_Rp)/(Earth_Re*Earth_Re)*coslat*sinlat;
	double coef2 = 1.0/(1.0+mu*mu);
	double nz = sqrt(coef2)*(sinlat+mu*coslat); 
	double D = sqrt(D2 - R*R*mu*mu*coef2);

	// calculate value of the function ->0 
	funcval = Z + R*mu*coef2*(coslat-mu*sinlat) - nz*D;

	// calculate derivative:
	double dR = ((Earth_Rp*coslat)*(Earth_Rp*coslat) + (Earth_Re*sinlat)*(Earth_Re*sinlat));
	dR = -Earth_Rp*Earth_Re/sqrt(dR*dR*dR)*(Earth_Re*Earth_Re-Earth_Rp*Earth_Rp)*sinlat*coslat;
	double dmu = (Earth_Re*Earth_Re - Earth_Rp*Earth_Rp)/(Earth_Rp*Earth_Rp)/(Earth_Re*Earth_Re) *R*(2*dR*coslat*sinlat + R*(coslat*coslat-sinlat*sinlat));
	double dnz = sqrt(coef2)/(1+mu*mu)*(1+mu*mu+dmu)*(coslat-mu*sinlat);
	double dcoef2 = -2*mu*dmu/((1+mu*mu)*(1+mu*mu));
	double dD = 0.5/D*(-2*R*dR*mu*mu*coef2 - 2*R*R*mu*dmu*coef2 - R*R*mu*mu*dcoef2);

	derval = (dR*mu*coef2*(coslat-mu*sinlat) + 
			  R*dmu*coef2*(coslat-mu*sinlat) +
			  R*mu*dcoef2*(coslat-mu*sinlat) +
			  R*mu*coef2*(-sinlat-dmu*sinlat-mu*coslat) -
			  dnz*D - nz*dD) * Rad;
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Convert ECEF Cartesian coordinates (X,Y,Z) to longitude lon, latitdue lat and 
// height in spherical system
////////////////////////////////////////////////////////////////////////////////////////////////////////////

void XYZ2LonLatH_sph(double& lon_out, double& lat_out, double& hgt_out, double X, double Y, double Z)
{
	// Input: triplets of [longitude, latitude, height]
	// Output: triplets of [X, Y, Z]

	// Earth radius (m):
	const double Earth_R = 6378137.0;
	const double Rad = 3.1415926535897932384626433832795/180.0; // pi/180
	const double RRad = 180.0/3.1415926535897932384626433832795; // 180/pi


	double D = sqrt(X*X + Y*Y + Z*Z);

	// latitude (degrees)
	double lat = RRad * asin(Z/Earth_R);

	// Calculate altitude in ellipsoidal system:
	double coslat = cos(Rad*lat);
	double sinlat = sin(Rad*lat);

	// find altitude = sqrt(X^2+Y^2+Z^2) - radius
	double hgt = D - Earth_R;

	// find longitude from equation:
	// X = R*cosd(lat)*cosd(lon) + nx*alt == R*cosd(lat)*cosd(lon) + alt*coef*(cosd(lat)-mu*sind(lat))*cosd(lon) ==
	// {R*cosd(lat)+alt*coef*(cosd(lat)-mu*sind(lat))}*cosd(lon)

	double lon = RRad*acos(X/D);
	if (Y < 0.0) lon = (360.0-lon); // change to 360-... if location is at the other side of Earth

	lon_out = lon;
	lat_out = lat;
	hgt_out = hgt;

}

