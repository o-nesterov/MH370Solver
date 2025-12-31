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
#include "trajectory.h"

extern void uv2uvw(double& XYZ_U, double& XYZ_V, double& XYZ_W, double u, double v, double lon, double lat);
extern void w2uvw(double& XYZ_U, double& XYZ_V, double& XYZ_W, double w, double lon, double lat);

////////////////////////////////////////////////////////////////////////////////////////
//
// Auxilary function, which calculates TAS based on input IAS
//
////////////////////////////////////////////////////////////////////////////////////////

inline double get_tas_from_ias(double ias, double air_t, double air_pressure)
{
	// Input: 
	// ias - IAS speed (note m/s, not knots)
	// air_t - air temperature (deg C)
	// air_pressure - air pressure (Pa)

	// Output:
	// TAS (m/s)

	double air_dens = air_pressure/(287.058*(air_t+273.16)); // air density
	double sound_sp = 331.3*sqrt(1.0+air_t/273.16); // sound speed 
	// (it can also be corrected for humidity, but effect is likely to be minimal) 

	// Use q/P = (1+0.2*M^2)^3.5 - 1, where q is impact pressure, M is Mach
	// Reference:
	// Nguyen X. Vinh (1993). Flight Mechanics of High-Performance Aircraft
	// Volume 4 of Cambridge Aerospace Series, Cambridge University Press, 1993
	// ISBN	0521478529, 9780521478526
	// 382 pages
	// Equation (3.42)

	// By definition IAS = CAS = TAS = EAS at the sea level standard atmospere
	// q/P0 = (1+0.2*(IAS/a0)^2)^3.5 - 1, where IAS - indicated airspeed, 
	// a0 = 340.27 m/s - sound speed at the sea level, standard atmospere +15.0C (288.15K)

	double dp = 1.0+0.2*(ias*ias/(340.27*340.27));
	dp = 101325.0*(dp*dp*dp*sqrt(dp) - 1.0); // Pa

	double mach = sqrt(5.0*(pow((1.0+dp/air_pressure),(1.0/3.5))-1.0));

	return mach*sound_sp;

};



////////////////////////////////////////////////////////////////////////////////////////
//
// Auxilary function, which calculates TAS based on input Mach
//
////////////////////////////////////////////////////////////////////////////////////////

inline double get_tas_from_mach(double mach, double air_t)
{
	// Input: 
	// ias - IAS speed (note m/s, not knots)
	// air_t - air temperature (deg C)
	// air_pressure - air pressure (Pa)

	// Output:
	// TAS (m/s)

	double sound_sp = 331.3*sqrt(1.0+air_t/273.16); // sound speed
	return mach*sound_sp;

};




////////////////////////////////////////////////////////////////////////////////////////
//
// Auxilary function, which calculates fuel flow rate
//
////////////////////////////////////////////////////////////////////////////////////////

inline double get_ff(double TAS, double m, double air_t, double rho_a, bool isSingleEngine = false)
{
	// Input: 
	// TAS - true airspeed (note m/s, not knots)
	// m - weigth (kg)
	// air_t - air temperature (deg C)
	// rho_a - air density

	// Output:
	// fuel flow (kg/s)

	// Fitting constants are based on fitting tables PI.21.3 (both engines) / PI.23.8 (one engine) and PI.21.4 / PI.23.9
	// 777-200/-200LR/-300/-300ER/-200F
	// Flight Crew Operations Manual, the Boeing Company, 1994, Rev. June 16, 2008,
	// Doc. no. D632W001-TBC

	double f0, f1, f2, Ct, Cd10, Cd11, Cd12;

	if (!isSingleEngine)
	{
		// Both engines operational
		f0 = 4.876818835801012e-006;
		f1 = 4.532449681392403e-007;
		f2 = 1.705274362111052e-008;
		Ct = 0.00105994131967;
		Cd10 = 8.17913242456180;
		Cd11 = 6.95174565307222;
		Cd12 = 2.066849788581953e+003;
	}
	else
	{
		// Single engine operational
		f0 = 3.101827008210732e-006;
		f1 = 6.802506429445600e-007;
		f2 = 2.081182444241502e-008;
		Cd10 = 6.42843142211523;
		Cd11 = 1.88747236450717;
		Cd12 = 5.077231209547657e+003;
		Ct = 3.014011529708891e-004;
	};

	// Calculate angle of attack, rads
	// References: Richardson T.S., Beaverstock C., Isikveren A., Meheri A., Badcock K, Ronch A.D.: Analysis of the Boeing 747-100 using CEASIOM. 
	// Progress in Aerospace Sciences, V.47(8), 2011, 660--673. https://doi.org/10.1016/j.paerosci.2011.08.009;
	// Also: https://aerospaceweb.org/question/aerodynamics/q0252.shtml#:~:text=We%20can%20again%20use%20the,provided%20in%20the%20original%20data.

	double CLa = 5.5; // B747 CL linear model: CL = CL0 + CLa*alpha
	double CL0 = 0.29;
	double S_ref = 427.8; // wings area, m^2
	double alpha = (2.0*m*9.81/(rho_a*TAS*TAS*S_ref) - CL0)/CLa; //angle of attack, rads

	// calculate drag coefficient
	double Cd1 = Cd10 + Cd11*alpha + Cd12*alpha*alpha;

	// calculate fuel flow (kg/s)
	double ff = (f0+f1*sqrt(273.16+air_t)+f2*TAS)*0.5*Cd1*rho_a*TAS*TAS + Ct*(15.0-air_t); // fuel flow (kg/s)

	return ff;
};




inline double min_dbl(double a, double b)
{
  return (a < b)?a:b;

}



////////////////////////////////////////////////////////////////////////////////////////
//
// Calculate u (W->E), v (S->N), w (vertical) velocity components in the local system, and fuel flow rate q
// based on maneuver type
//
////////////////////////////////////////////////////////////////////////////////////////
void CManeuver::Calc_uvwq(double& u, double& v, double& w, double& q, double& dphi,
						  double t, double lon, double lat, double alt, double m, double phi)
{
	const double Rad = 3.1415926535897932384626433832795/180.0; // pi/180

	// ---------------------------------------------
	// Get atmospheric conditions:
	double Wind_X, Wind_Y, air_t, air_p, rel_hum;

	if (pEnvironData_)
	{
		pEnvironData_->GetMeteo(Wind_X, Wind_Y, air_t, air_p, rel_hum, t+t0_, lon, lat, alt);
	}
	else
	{
		Wind_X = 0.0;
		Wind_Y = 0.0;
		air_t = 25.0 - alt*6.5/1000.0;
		air_p = 101325.0 * exp(-9.81*0.02896968*alt/288.15/8.314462618); // e.g., https://en.wikipedia.org/wiki/Atmospheric_pressure
		rel_hum = 0.0;
	}

	double air_dens = air_p/(287.058*(air_t+273.16)); // air density

	// --------------------------------------------
	// Calculate true air speed (TAS) from indicated air speed (IAS) or Mach
	double Max_TAS_VFC = get_tas_from_ias(AIRCRAFT_VFC*0.514444, air_t, air_p); // max TAS for stability based on IAS
	double Max_TAS_MFC = get_tas_from_mach(AIRCRAFT_MFC, air_t); // max TAS for stability based on Mach
	double MAX_TAS = min_dbl(Max_TAS_VFC, Max_TAS_MFC); // max TAS based on B777 envelope protection
	double TAS = 0.0;

	if (SPEED_MODE_ == SPEED_MODE_IAS)
	{
		double airspeed = min_dbl(airspeed_, AIRCRAFT_VFC*0.514444);
		TAS = get_tas_from_ias(airspeed, air_t, air_p);
		TAS = min_dbl(TAS, MAX_TAS);
	}
	else
	{
		double airspeed = min_dbl(airspeed_, AIRCRAFT_MFC);
		TAS = get_tas_from_mach(airspeed, air_t);
		TAS = min_dbl(TAS, MAX_TAS);
	};

	// --------------------------------------------
	// Magnetic declination (deg) if heading or track mode is MAGN
	double decl = 0.0;

	if (HEADING_MODE_ == HEADING_MODE_MAGNHDG || HEADING_MODE_ == HEADING_MODE_MAGNTRK)
	{
		if (pEnvironData_)
		{
			 decl = pEnvironData_->GetMagneticDeclination(lon,lat);
		}
		else
		{
			 decl = 0.0;
		}
	}

	// --------------------------------------------
	// Gyroscopic heading (HEADING_MODE_GYROHDG mode)
	double n_lon_gyro = 0.0;
	double n_lat_gyro = 0.0;

	if (HEADING_MODE_ == HEADING_MODE_GYROHDG)
	{
		// Define gyroscopic North vector in the non-spinning reference frame (note ECEF system is spinning). 
		// Consider the two parameters being longitude and latitude in the absolute reference frame. ECEF coincides with this frame on 2014-03-07 00:00:00.
		// Note: this vector remains constant over the integration period
		// Presumably it may match the time when SAARU's north was last set, and is may not necessarily be when the SAARU was started up in the KLIA

		double p0x, p0y, p0z;
		//uv2uvw(p0x, p0y, p0z, 0.0, 1.0, gyro_heading_param1_, gyro_heading_param2_);
		uv2uvw(p0x, p0y, p0z, 0.0, 1.0, gyro_heading_param1_-360.0/86400.0*(t+t0_), gyro_heading_param2_); // unit vector, which defines gyroscopic "north" in the ECEF frame

		// Now calculate {u,v} direction vectors of the local tangential plane in the non-rotating reference frame
		// use uv2uvw to compute W->E and S->N direction vectors, and then projections of the "gyroscopic north" vector in ECEF (assuming that SAARU thinks it is the true north)
		// on those to determine the deviation of the "gyroscopic reference heading" from the true north
		double n_we_x, n_we_y, n_we_z;
		//uv2uvw(n_we_x, n_we_y, n_we_z, 1.0, 0.0, lon-360.0/86400.0*(t+t0_), lat);	// unit vector west to east in the non-spinning reference frame
		uv2uvw(n_we_x, n_we_y, n_we_z, 1.0, 0.0, lon, lat);	// unit vector west to east in the ECEF reference frame

		double n_sn_x, n_sn_y, n_sn_z;
		//uv2uvw(n_sn_x, n_sn_y, n_sn_z, 0.0, 1.0, lon-360.0/86400.0*(t+t0_), lat);	// unit vector south to north in the non-spinning reference frame
		uv2uvw(n_sn_x, n_sn_y, n_sn_z, 0.0, 1.0, lon, lat);	// unit vector south to north in the ECEF reference frame

		// Now we need to present gyroscopic heading vector as the sum of vectors of the local reference frame, which moves with the airplane
		// To do that calculate projections of the gyroscopic N vector on the axes of the tangential plane (vertical component is ignored)
		n_lon_gyro = p0x*n_we_x + p0y*n_we_y + p0z*n_we_z;
		n_lat_gyro = p0x*n_sn_x + p0y*n_sn_y + p0z*n_sn_z;

		// When the plane moves along a meridian, this vector becomes not parallel to the tangential plane, so its length becomes < 1.0. 
		// Scale it to have a unit vector in the local system.
		double s = sqrt(n_lon_gyro*n_lon_gyro + n_lat_gyro*n_lat_gyro);
		n_lon_gyro = n_lon_gyro/s;
		n_lat_gyro = n_lat_gyro/s;
	}

	// --------------------------------------------
	// Horizontal velocity components w.r.t. air
	double heading = hdg0_ + phi; // heading ('straight' and turn modes)

	if (HEADING_MODE_==HEADING_MODE_TRUEHDG || HEADING_MODE_==HEADING_MODE_MAGNHDG)
	{
		// True or magnetic heading mode:
		u = TAS * cos((90.0-(heading+decl))*Rad) + Wind_X;
		v = TAS * sin((90.0-(heading+decl))*Rad) + Wind_Y;
	}

	if (HEADING_MODE_==HEADING_MODE_TRUETRK|| HEADING_MODE_==HEADING_MODE_MAGNTRK)
	{
		// True or magnetic track (course) mode. Note:
		// ground_spd * cosd(90.0-(heading+decl)) == TAS*cosd(alpha) + Wind_X;
		// ground_spd * sind(90.0-(heading+decl)) == TAS*sind(alpha) + Wind_Y; 
		// ground_spd and alha are unknown
    
		// wind projection on required local track of aircraft:
		double nx = cos((90.0-(heading+decl))*Rad);
		double ny = sin((90.0-(heading+decl))*Rad);
    
		// calculate projection of wind component on the direction perpendicular to track (i.e., in direction {ny, -nx})
		double wind_p = Wind_X*ny - Wind_Y*nx;
		double sin_delta = wind_p/TAS; // sine of the angle, which the aircraft needs to deviate to compensate for 
    
		// ground speed: sum of the projection of aircraft velocity and wind velocity on {nx, ny} 
		double ground_spd = TAS*sqrt(1.0-sin_delta*sin_delta) + Wind_X * nx + Wind_Y * ny;
    
		// velocity components
		u = ground_spd * nx;
		v = ground_spd * ny;
	}

	if (HEADING_MODE_ == HEADING_MODE_GYROHDG)
	{
		// Turn "gyroscopic north" unit-vector in the local system CLOCKWISE by specified heading degrees
		double n_lon =  n_lon_gyro*cos(heading*Rad) + n_lat_gyro*sin(heading*Rad);
		double n_lat = -n_lon_gyro*sin(heading*Rad) + n_lat_gyro*cos(heading*Rad);

		// compute local velocity components
		u = TAS * n_lon + Wind_X;
		v = TAS * n_lat + Wind_Y;
	}


    //  --------------------------------------------
	// vertical velocity component
	double ws = 0.0;

	if (pTrajectory_) ws = pTrajectory_->GetVerticalSpd_(t+t0_,alt,TAS);

	w = ws;

	// ---------------------------------------------
	// rate of angular velocity
	dphi = 0.0; // deg/s

    //  --------------------------------------------
	// fuel flow rate (kg/s) based on parameterization of tabular values
	bool isSingleEngine = pTrajectory_ -> IsEngineModeSingle();
	q = get_ff(TAS, m, air_t, air_dens, isSingleEngine);

};






////////////////////////////////////////////////////////////////////////////////////////
//
// Calculate u (W->E), v (S->N), w (vertical) velocity components in the local system, and fuel flow rate q
// Turn maneuver
//
////////////////////////////////////////////////////////////////////////////////////////
void CManeuverTurn::Calc_uvwq(double& u, double& v, double& w, double& q, double& dphi,
							  double t, double lon, double lat, double alt, double m, double phi)
{
	const double Rad = 3.1415926535897932384626433832795/180.0; // pi/180

	// ---------------------------------------------
	// Get atmospheric conditions:
	double Wind_X, Wind_Y, air_t, air_p, rel_hum;

	if (pEnvironData_)
	{
		pEnvironData_->GetMeteo(Wind_X, Wind_Y, air_t, air_p, rel_hum, t+t0_, lon, lat, alt);
	}
	else
	{
		Wind_X = 0.0;
		Wind_Y = 0.0;
		air_t = 25.0 - alt*6.5/1000.0;
		air_p = 101325.0 * exp(-9.81*0.02896968*alt/288.15/8.314462618); // e.g., https://en.wikipedia.org/wiki/Atmospheric_pressure
		rel_hum = 0.0;
	}

	double air_dens = air_p/(287.058*(air_t+273.16)); // air density

	// --------------------------------------------
	// Calculate true air speed (TAS) from indicated air speed (IAS) or Mach
	double Max_TAS_VFC = get_tas_from_ias(AIRCRAFT_VFC*0.514444, air_t, air_p); // max TAS for stability based on IAS
	double Max_TAS_MFC = get_tas_from_mach(AIRCRAFT_MFC, air_t); // max TAS for stability based on Mach
	double MAX_TAS = min_dbl(Max_TAS_VFC, Max_TAS_MFC); // max TAS based on B777 envelope protection
	double TAS = 0.0;

	if (SPEED_MODE_ == SPEED_MODE_IAS)
	{
		double airspeed = min_dbl(airspeed_, AIRCRAFT_VFC*0.514444);
		TAS = get_tas_from_ias(airspeed, air_t, air_p);
		TAS = min_dbl(TAS, MAX_TAS);
	}
	else
	{
		double airspeed = min_dbl(airspeed_, AIRCRAFT_MFC);
		TAS = get_tas_from_mach(airspeed, air_t);
		TAS = min_dbl(TAS, MAX_TAS);
	};

	// --------------------------------------------
	// Magnetic declination (deg) if heading or track mode is MAGN
	double decl = 0.0;

	if (HEADING_MODE_ == HEADING_MODE_MAGNHDG || HEADING_MODE_ == HEADING_MODE_MAGNTRK)
	{
		if (pEnvironData_)
		{
			 decl = pEnvironData_->GetMagneticDeclination(lon,lat);
		}
		else
		{
			 decl = 0.0;
		}
	}

	// --------------------------------------------
	// Gyroscopic heading (HEADING_MODE_GYROHDG mode)
	double n_lon_gyro = 0.0;
	double n_lat_gyro = 0.0;

	if (HEADING_MODE_ == HEADING_MODE_GYROHDG)
	{
		// Define gyroscopic North vector in the non-spinning reference frame (note ECEF system is spinning). 
		// Consider the two parameters being longitude and latitude in the absolute reference frame. ECEF coincides with this frame on 2014-03-07 00:00:00.
		// Note: this vector remains constant over the integration period
		// Presumably it may match the time when SAARU's north was last set, and is may not necessarily be when the SAARU was started up in the KLIA

		double p0x, p0y, p0z;
		// uv2uvw(p0x, p0y, p0z, 0.0, 1.0, gyro_heading_param1_, gyro_heading_param2_);
		uv2uvw(p0x, p0y, p0z, 0.0, 1.0, gyro_heading_param1_-360.0/86400.0*(t+t0_), gyro_heading_param2_);

		// Now calculate {u,v} direction vectors of the local tangential plane in the non-rotating reference frame
		// use uv2uvw to compute W->E and S->N direction vectors, and then projections of the "gyroscopic north" vector in ECEF (assuming that SAARU thinks it is the true north)
		// on those to determine the deviation of the "gyroscopic reference heading" from the true north
		double n_we_x, n_we_y, n_we_z;
		//uv2uvw(n_we_x, n_we_y, n_we_z, 1.0, 0.0, lon-360.0/86400.0*(t+t0_), lat);	// unit vector west to east in the non-spinning reference frame
		uv2uvw(n_we_x, n_we_y, n_we_z, 1.0, 0.0, lon, lat);	// unit vector west to east in the ECEF reference frame

		double n_sn_x, n_sn_y, n_sn_z;
		//uv2uvw(n_sn_x, n_sn_y, n_sn_z, 0.0, 1.0, lon-360.0/86400.0*(t+t0_), lat);	// unit vector south to north in the non-spinning reference frame
		uv2uvw(n_sn_x, n_sn_y, n_sn_z, 0.0, 1.0, lon, lat);	// unit vector south to north in the ECEF reference frame

		// Now we need to present gyroscopic heading vector as the sum of vectors of the local reference frame, which moves with the airplane
		// To do that calculate projections of the gyroscopic N vector on the axes of the tangential plane (vertical component is ignored)
		n_lon_gyro = p0x*n_we_x + p0y*n_we_y + p0z*n_we_z;
		n_lat_gyro = p0x*n_sn_x + p0y*n_sn_y + p0z*n_sn_z;

		// When the plane moves along a meridian, this vector becomes not parallel to the tangential plane, so its length becomes < 1.0. 
		// Scale it to have a unit vector in the local system.
		double s = sqrt(n_lon_gyro*n_lon_gyro + n_lat_gyro*n_lat_gyro);
		n_lon_gyro = n_lon_gyro/s;
		n_lat_gyro = n_lat_gyro/s;
	}

	// --------------------------------------------
	// Horizontal velocity components w.r.t. air
	double heading = hdg0_ + phi; // heading ('straight' and turn modes)

	if (HEADING_MODE_==HEADING_MODE_TRUEHDG || HEADING_MODE_==HEADING_MODE_MAGNHDG)
	{
		// True or magnetic heading mode:
		u = TAS * cos((90.0-(heading+decl))*Rad) + Wind_X;
		v = TAS * sin((90.0-(heading+decl))*Rad) + Wind_Y;
	}

	if (HEADING_MODE_==HEADING_MODE_TRUETRK|| HEADING_MODE_==HEADING_MODE_MAGNTRK)
	{
		// True or magnetic track (course) mode. Note:
		// ground_spd * cosd(90.0-(heading+decl)) == TAS*cosd(alpha) + Wind_X;
		// ground_spd * sind(90.0-(heading+decl)) == TAS*sind(alpha) + Wind_Y; 
		// ground_spd and alha are unknown
    
		// wind projection on required local track of aircraft:
		double nx = cos((90.0-(heading+decl))*Rad);
		double ny = sin((90.0-(heading+decl))*Rad);
    
		// calculate projection of wind component on the direction perpendicular to track (i.e., in direction {ny, -nx})
		double wind_p = Wind_X*ny - Wind_Y*nx;
		double sin_delta = wind_p/TAS; // sine of the angle, which the aircraft needs to deviate to compensate for 
    
		// ground speed: sum of the projection of aircraft velocity and wind velocity on {nx, ny} 
		double ground_spd = TAS*sqrt(1.0-sin_delta*sin_delta) + Wind_X * nx + Wind_Y * ny;
    
		// velocity components
		u = ground_spd * nx;
		v = ground_spd * ny;
	}

	if (HEADING_MODE_ == HEADING_MODE_GYROHDG)
	{
		// Turn "gyroscopic north" unit-vector in the local system CLOCKWISE by specified heading degrees
		double n_lon =  n_lon_gyro*cos(heading*Rad) + n_lat_gyro*sin(heading*Rad);
		double n_lat = -n_lon_gyro*sin(heading*Rad) + n_lat_gyro*cos(heading*Rad);

		// compute local velocity components
		u = TAS * n_lon + Wind_X;
		v = TAS * n_lat + Wind_Y;
	}


    //  --------------------------------------------
	// vertical velocity component
	double ws = 0.0;

	if (pTrajectory_) ws = pTrajectory_->GetVerticalSpd_(t+t0_,alt,TAS);

	w = ws;

	// ---------------------------------------------
	// rate of angular velocity
	// TAS^2/R = g x tan(bank_angle_*pi/180.0) => R = TAS^2/[g x tan(bank_angle_*pi/180.0)] - turning radius
	// => T = 2*pi*R/TAS = 2*pi*TAS/[g * tan(bank_angle_*pi/180.0)] - time required to compete one full rotation of 360 deg
	const double pi = 3.1415926535897932384626433832795;
	dphi = turn_sign_ * 180.0/(pi*TAS) * 9.81 * tan(bank_angle_*pi/180.0); // deg/s


    //  --------------------------------------------
	// fuel flow rate (kg/s) based on parameterization of tabular values
	bool isSingleEngine = pTrajectory_ -> IsEngineModeSingle();
	q = get_ff(TAS, m, air_t, air_dens, isSingleEngine);

};

