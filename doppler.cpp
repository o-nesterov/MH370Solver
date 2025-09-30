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
#include "inmarsatdata.h" // some constants

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This function calculates Doppler shift (Hz)
////////////////////////////////////////////////////////////////////////////////////////////////////////////
double doppler(double XYZs[3], double XYZr[3], double UVWs[3], double UVWr[3], double Frq, int Method)
{
	// Input: 
	// XYZs: triplets of [Xs, Ys, Zs] source coordinates
	// XYZr: triplets of [Xr, Yr, Zr] receiver coordinates
	// UVWs: triplets of [Us, Vs, Ws] velocity components of the source
	// UVWr: triplets of [Ur, Vr, Wr] velocity components of the receiver
	// Frq: base frequency
	// Method - Doppler computation methods:
	// - DOPPLER_METHOD_SIMPLE:Simple Doppler shift
	// - DOPPLER_METHOD_REL: Doppler shift with relativistic effect
	// - DOPPLER_METHOD_ADV: Advanced: regular + relativistic + gravitational

	// Output: Doppler shift (Hz)

	// Speed of Light, m/s
	const long double Speed_of_light = 299792500.0;

	// --------------------------------------------

	// relative velocity of the receiver point w.r.t. source
	double UVW[3];
	UVW[0] = UVWr[0] - UVWs[0];
	UVW[1] = UVWr[1] - UVWs[1];
	UVW[2] = UVWr[2] - UVWs[2];
	long double r = sqrt(UVW[0]*UVW[0] + UVW[1]*UVW[1] + UVW[2]*UVW[2]) / Speed_of_light;

	// compute normilized orientation vector
	double XYZ[3];
	XYZ[0] = XYZr[0] - XYZs[0];
	XYZ[1] = XYZr[1] - XYZs[1];
	XYZ[2] = XYZr[2] - XYZs[2];

	double D = sqrt(XYZ[0]*XYZ[0] + XYZ[1]*XYZ[1] + XYZ[2]*XYZ[2]);

	XYZ[0] = XYZ[0] / D;
	XYZ[1] = XYZ[1] / D;
	XYZ[2] = XYZ[2] / D;

	// compute ratio of the projection of the relative velocity on the orientation vector to speed of light
	// note that computing projection is equivalent to multiplying by cosine of the angle 
	long double rp = (XYZ[0]*UVW[0] + XYZ[1]*UVW[1] + XYZ[2]*UVW[2])/Speed_of_light; // rp > 0 if source and receiver are moving away from each other; 
	                                                                                 // rp < 0 if they move toward each other  

	// Compute Doppler shift accosiated with relative velocities, e.g.
	// https://www.feynmanlectures.caltech.edu/I_34.html, 
	// https://en.wikipedia.org/wiki/Relativistic_Doppler_effect#Motion_in_an_arbitrary_direction
	// https://www.fourmilab.ch/etexts/einstein/specrel/www/, "ON THE ELECTRODYNAMICS OF MOVING BODIES By A. Einstein June 30, 1905"
	// http://www.mrelativity.net/CompGravRelDopEffect/Complete%20Gravitational%20Relativistic%20Doppler%20Effect.htm

	// long double c1 = sqrt(1-r*r)/(1.0+rp); // including relativistic transformation effect
	// note: Intel and GCC support long double and use FPU; Microsoft equates long double to double and uses SSE
	// accuracy may deteriorate if r and rp are small, so we use approximate sqrt(1-r*r) ~= 1 - r^2/2 - r^4/8 - r^6/16 where needed

	// 1/(1+rp) == 1 - rp/(1+rp) - exact
	long double c1rpm1 = -rp/(1.0+rp);
	long double c1rm1, c1r, c1f, R0, R1, F0, F1;

	// here c1 = (1+c1rm1)*(1+c1rpm1) = 1+c1rm1+c1rpm1+c1rm1*c1rpm1

	switch (Method)
	{
	case DOPPLER_METHOD_SIMPLE:
		return (double) (Frq*c1rpm1);
		break;
		
	case DOPPLER_METHOD_REL:

		// sqrt(1-r*r) ~= 1 - r^2/2 - r^4/8 - r^6/16 - approximate; if r<0.001, then error of this approximation would be smaller than 1.0E-18
		c1rm1 = -r*r*(0.5 + 0.125*r*r*(1.0 + 0.5*r*r));
		c1r = c1rm1+c1rpm1+c1rm1*c1rpm1;
		return (double) (Frq*c1r);
		break;

	case DOPPLER_METHOD_ADV:
		// Gravitational constant
		const long double Gravit = 6.67384E-11; // m^3 kg^-1 s^-2

		// Earth mass
		const long double Earth_mass = 5.97219E+24; // Earth mass, kg

		c1rm1 = -r*r*(0.5 + 0.125*r*r*(1.0 + 0.5*r*r));
		c1r = c1rm1+c1rpm1+c1rm1*c1rpm1;

		R0 = sqrt(XYZs[0]*XYZs[0] + XYZs[1]*XYZs[1] + XYZs[2]*XYZs[2]);	// distance from the source to the Earth center
		R1 = sqrt(XYZr[0]*XYZr[0] + XYZr[1]*XYZr[1] + XYZr[2]*XYZr[2]);	// distance from the receiver to the Earth center
		F0 = Gravit*Earth_mass / (Speed_of_light*Speed_of_light*R0);		// source's gravitational potential
		F1 = Gravit*Earth_mass / (Speed_of_light*Speed_of_light*R1);		// receiver's gravitational potential

		// c1 = c1*(1.0-F0)*(1.0+F1);
		// here c1 ~= (1+c1r) => (1+c1r)*(1-F0)*(1+F1) = (1+c1r)*(1-F0+F1-F0*F1) = 1 + c1r + (-F0+F1-F0*F1) + c1r*(-F0+F1-F0*F1)

		c1f = -F0+F1-F0*F1;

		//return (double) (Frq*(c1 - 1.0));
		return (double) (Frq*(c1r + c1f + c1r*c1f));
		break;
	};

	return 0.0;
}

