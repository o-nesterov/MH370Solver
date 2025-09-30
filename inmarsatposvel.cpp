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
#include "inmarsatdata.h"


////////////////////////////////////////////////////////////////////////////////////////
// Class to handle Inmarsat data
////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////
//
// Return Inmarsat's F3 satellite position in ECEF system
//
////////////////////////////////////////////////////////////////////////////////////////
void CMH370InmarsatData::GetInmarsatPos(double& X, double& Y, double& Z, double t_inp)
{
	// Input t_inp is the time (s) since since UTC 2014-03-07 00:00:00
	// Output: {X,Y,Z} - position of Inmarsat F3 satellite in meters in ECEF system
	// This is based on fitting of Table 4: Satellite Location and Velocity (ECEF) in
	// "mr052_MH370_Definition_of_Sea_Floor_Wide_Area_Search.pdf" by ATSB, 2014

	const double pi = 3.1415926535897932384626433832795;

	const double start_t = 59400.0;		// 59400 is 16:30:00 offset (s) - first value provided by ATSB

	double t = 2.0*pi*(t_inp - start_t)/86400.0;	

	const double a0x = 18119.47816;
	const double b0x = 0.001044687659;
	const double a1x = 28.86238622;
	const double b1x = -6.2831853;
	const double a2x = -9.18690232;
	const double b2x = 5.838509104;
	const double a3x = 1.162423398;
	const double b3x = 5.804735151;

	X = (a0x + b0x*(t_inp-start_t) + a1x*sin(t+b1x) + a2x*sin(2*t+b2x) + a3x*sin(3*t+b3x))*1000.0; // X coordinate (m)

	const double a0y = 38077.61366;
	const double b0y = -0.001173120311;
	const double a1y = 5.77388945;
	const double b1y = -6.056387291;
	const double a2y = -0.8155913874;
	const double b2y = 6.2831853;
	const double a3y = -1.445978767;
	const double b3y = 5.431131899;

	Y = (a0y + b0y*(t_inp-start_t) + a1y*sin(t+b1y) + a2y*sin(2*t+b2y) + a3y*sin(3*t+b3y))*1000.0; // Y coordinate (m)

	const double a1z = 1209.666621;
	const double b1z = 0.7587541017;
	const double a2z = 0.4082203063;
	const double b2z = 5.952583427;
	const double c0z = -3.634127047;


	Z = (c0z + a1z*sin(t+b1z) + a2z*sin(2*t+b2z))*1000.0; // Z coordinate (m)

};



////////////////////////////////////////////////////////////////////////////////////////
//
// Return Inmarsat's F3 satellite velocity in ECEF system
//
////////////////////////////////////////////////////////////////////////////////////////
void CMH370InmarsatData::GetInmarsatVel(double& U, double& V, double& W, double t_inp)
{
	// Input t_inp is the time (s) since since UTC 2014-03-07 00:00:00
	// Output is {U,V,W} - velocity of Inmarsat satellite in m/s in ECEF system
	// Based on fitting Table 4: Satellite Location and Velocity (ECEF) from
	// "mr052_MH370_Definition_of_Sea_Floor_Wide_Area_Search.pdf" by ATSB

	const double pi = 3.1415926535897932384626433832795;

	double t = 2.0*pi*(t_inp - 59400.0)/86400.0;	// 59400 is 16:30:00 offset (s) - first value provided by ATSB

	const double a0u = 0.1505177265;
	const double a1u = 2.876429598;
	const double b1u = 1.020516292;
	const double a2u = 1.147506149;
	const double b2u = -2.727698203;
	const double a3u = 0.02234532061;
	const double b3u = -4.268483808;

	U = a0u + a1u*sin(t+b1u) + a2u*sin(2*t+b2u) + a3u*sin(3*t+b3u);

	const double a0v = 0.2330510348;
	const double a1v = 2.396458238;
	const double b1v = 3.173795694;
	const double a2v = 1.38063347;
	const double b2v = -1.141452099;
	const double a3v = 0.02964700134;
	const double b3v = 1.375217103;

	V = a0v + a1v*sin(t+b1v) + a2v*sin(2*t+b2v) + a3v*sin(3*t+b3v);

	const double a1w = 87.97752146;
	const double b1w = 2.329502275;
	const double a2w = 0.05147551934;
	const double b2w = 1.285845504;

	W = a1w*sin(t+b1w) + a2w*sin(2*t+b2w);
}
