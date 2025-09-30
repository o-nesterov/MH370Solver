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
#include "maneuver.h"
#include "trajectory.h"
#include "ode5.h"



extern void XYZ2LonLatH(double& lon, double& lat, double& hgt, double X, double Y, double Z);
extern void LonLatH2XYZ(double &X, double &Y, double &Z, double lon, double lat, double hgt);
extern void uv2uvw(double& XYZ_U, double& XYZ_V, double& XYZ_W, double u, double v, double lon, double lat);
extern void w2uvw(double& XYZ_U, double& XYZ_V, double& XYZ_W, double w, double lon, double lat);
extern double get_tas_from_ias(double ias, double air_t, double air_pressure);
extern double get_tas_from_mach(double mach, double air_t);
extern double min_dbl(double a, double b);

// Integration termination codes:
// VNAV - change in vertical maneuver type
// LNAV - reached waypoint
#define TERMINATION_CODE_NONE 0
#define TERMINATION_CODE_VNAV 1
#define TERMINATION_CODE_LNAV 2

////////////////////////////////////////////////////////////////////////////////////////
// Maneuver class
// Integrates a primitive element of trajectrory
// The whole trajectory is comprised of the sequence of these maneuvers (legs, turns, etc.)
////////////////////////////////////////////////////////////////////////////////////////

CManeuver::CManeuver(double t0, double tmax, double timestep, double* pSTimes, unsigned long nSTimes,
					 double lon0, double lat0, double alt0, double M0, 
					 int speed_mode, double airspeed, int heading_mode, double hdg0,
					 CMH370EnvironData* pEnvironData, CTrajectory* pTrajectory)
{
	// Airspeed mode
	SPEED_MODE_ = speed_mode;		// SPEED_MODE_IAS or SPEED_MODE_MACH
	airspeed_ = airspeed;			// air speed parameter (IAS or MACH)

	// Heading mode
	HEADING_MODE_ = heading_mode;	// HEADING_MODE_TRUEHDG - true heading, HEADING_MODE_MAGNHDG - magnetic heading, HEADING_MODE_TRUETRK - true track, HEADING_MODE_MAGNTRK - magnetic track
	hdg0_ = hdg0;					// heading parameter (true/magnetic heading or true/magnetic track); it may changes through integration

	// store initial/reference parameters 
	pEnvironData_ = pEnvironData;		// object providing environmental data 
	pTrajectory_ = pTrajectory;			// pointer to whole trajectory providing vertical flight profile
	t0_ = t0;							// initial time (s) since 2014-03-07 00:00:00 UTC of the maneuver
	lon0_ = lon0;						// initial longitude
	lat0_ = lat0;						// initial latitude
	alt0_ = alt0;						// initial altitude
	M0_ = M0;							// initial weigth

	// integration parameters
	max_time_step_ = timestep;			// maximum time step
	duration_ = tmax - t0;				// initialize duration (it can be redefined later - shorter then this for turns)
	last_hdg_ = hdg0;					// initialize last heading (it can be redefined later, e.g. for turns)

	// gyroscopic parameters (used in HEADING_MODE_GYROHDG mode)
	gyro_heading_param1_ = 0.0;
	gyro_heading_param2_ = 0.0;

	// solution output
	nPoints_ = 0;			// number of output trajectory points, which depend on integration time step and duration
	pTime_ = NULL;			// array of time stamps (s) since 2014-03-07 00:00:00 UTC
	pLon_ = NULL;			// array or longitudes
	pLat_ = NULL;			// array of latitudes
	pAlt_ = NULL;			// array of altitudes
	pM_ = NULL;				// array of weights
	pu_ = NULL;				// u-velocity component in local system
	pv_ = NULL;				// v-velocity component in local system
	pw_ = NULL;				// w-velocity component in local system
	pq_ = NULL;				// fuel flow in local system

	// specific output times
	if (pSTimes!=NULL && nSTimes>0)
	{
		nSTimes_ = nSTimes;
		pSTimes_ = new double[nSTimes_+1]; // allocate array

		for (unsigned long i=0; i<nSTimes; i++) pSTimes_[i] = pSTimes[i] - t0; // set specific output times relative to the integration time of this maneuver
	}
	else
	{
		nSTimes_ = 0;
		pSTimes_ = NULL;
	}

	// message stream for diagnostic output (can be NULL) 
	fileout_ = NULL; 
};


////////////////////////////////////////////////////////////////////////////////////////
// 'Straight' leg maneuver class of specified duration
////////////////////////////////////////////////////////////////////////////////////////

CManeuverLeg::CManeuverLeg(double t0, double duration, double timestep,  double* pSTimes, unsigned long nSTimes,
						   double lon0, double lat0, double alt0, double M0, 
						   int speed_mode, double airspeed, int heading_mode, double hdg0,
						   CMH370EnvironData* pEnvironData, CTrajectory* pTrajectory):
CManeuver(t0, t0+duration, timestep, pSTimes, nSTimes,
		  lon0, lat0, alt0, M0, 
		  speed_mode, airspeed, heading_mode, hdg0,
		  pEnvironData, pTrajectory)
{
	duration_ = duration;		// duration must be specified for straight leg
	last_hdg_ = hdg0;			// last heading
}


////////////////////////////////////////////////////////////////////////////////////////
// Gyroscopic 'straight' leg maneuver class of specified duration
////////////////////////////////////////////////////////////////////////////////////////

CManeuverLeg::CManeuverLeg(double t0, double duration, double timestep,  double* pSTimes, unsigned long nSTimes,
						   double lon0, double lat0, double alt0, double M0, 
						   int speed_mode, double airspeed, double hdg0, double gyro_hdg_param1, double gyro_hdg_param2,
						   CMH370EnvironData* pEnvironData, CTrajectory* pTrajectory):
CManeuver(t0, t0+duration, timestep, pSTimes, nSTimes,
		  lon0, lat0, alt0, M0, 
		  speed_mode, airspeed, HEADING_MODE_GYROHDG, hdg0,
		  pEnvironData, pTrajectory)
{
	duration_ = duration;		// duration must be specified for straight leg
	last_hdg_ = hdg0;			// last heading

	gyro_heading_param1_ = gyro_hdg_param1; 
	gyro_heading_param2_ = gyro_hdg_param2;
}


////////////////////////////////////////////////////////////////////////////////////////
// Turn maneuver class of specified duration
////////////////////////////////////////////////////////////////////////////////////////


CManeuverTurn::CManeuverTurn(double t0, double max_duration, double timestep,  double* pSTimes, unsigned long nSTimes,
							 double lon0, double lat0, double alt0, double M0, 
							 int speed_mode, double airspeed, int heading_mode, double hdg0,
							 double turn_angle, double target_hdg, double bank_angle,
							 CMH370EnvironData* pEnvironData, CTrajectory* pTrajectory):
CManeuver(t0, t0+max_duration, timestep, pSTimes, nSTimes,
		  lon0, lat0, alt0, M0, 
		  speed_mode, airspeed, heading_mode, hdg0,
		  pEnvironData, pTrajectory)
{
	ref_hdg0_ = 0.0;  // initial reference heading to N (==initial magnetic declination or gyro)
	turn_type_ = TURN_TYPE_BY_ANGLE;
	turn_sign_ = 1.0;
	isCompleted_ = false;
	bank_angle_ = bank_angle;
	max_duration_ = max_duration;

	// change in heading
	pPhi_ = NULL;

	// derive turn rate and duration based on TAS at the beginning of maneuver
	init_(hdg0, turn_angle, target_hdg);
}


////////////////////////////////////////////////////////////////////////////////////////
// Turn maneuver class of specified duration
////////////////////////////////////////////////////////////////////////////////////////


CManeuverTurn::CManeuverTurn(double t0, double max_duration, double timestep,  double* pSTimes, unsigned long nSTimes,
							 double lon0, double lat0, double alt0, double M0, 
							 int speed_mode, double airspeed, double hdg0, double gyro_hdg_param1, double gyro_hdg_param2,
							 double turn_angle, double target_hdg, double bank_angle,
							 CMH370EnvironData* pEnvironData, CTrajectory* pTrajectory):
CManeuver(t0, t0+max_duration, timestep, pSTimes, nSTimes,
		  lon0, lat0, alt0, M0, 
		  speed_mode, airspeed, HEADING_MODE_GYROHDG, hdg0,
		  pEnvironData, pTrajectory)
{
	ref_hdg0_ = 0.0;  // initial reference heading to N (==initial magnetic declination or gyro)
	turn_type_ = TURN_TYPE_BY_ANGLE;
	turn_sign_ = 1.0;
	isCompleted_ = false;
	bank_angle_ = bank_angle;
	gyro_heading_param1_ = gyro_hdg_param1;
	gyro_heading_param2_ = gyro_hdg_param2;
	max_duration_ = max_duration;

	// change in heading
	pPhi_ = NULL;

	// derive turn rate and duration based on TAS at the beginning of maneuver
	init_(hdg0, turn_angle, target_hdg);
}


// common initializations for turn maneuver
void CManeuverTurn::init_(double current_hdg, double turn_angle, double target_hdg)
{
	// Make sure hdg (current heading) is in [0,360) range
	hdg0_ = current_hdg;
	while (hdg0_ <0.0) hdg0_ += 360.0;
	while (hdg0_>=360.0) hdg0_ -= 360.0;

	if (abs(target_hdg)>1.0E+10)
	{
		// this is turn by specified angle, which can be more than 360 deg
		target_turn_angle_ = turn_angle;
		turn_type_ =  TURN_TYPE_BY_ANGLE;
		if (target_turn_angle_>0.0)
		{
			turn_sign_ = 1.0;
		}
		else
		{
			turn_sign_ = -1.0;
		}
		return;
	}

	// ---------

	turn_type_ =  TURN_TYPE_TO_HDG;

	// Make sure target_hdg_ is in [0,360) range
	target_hdg_ = target_hdg;
	while (target_hdg_ <0.0) target_hdg_ += 360.0;
	while (target_hdg_>=360.0) target_hdg_ -= 360.0;

	double hdg_diff = target_hdg_ - hdg0_;	// difference

	if (abs(hdg_diff)<=0.001)
	{
		// assume hdg and hdg_new identical, so do nothing
		last_hdg_ = target_hdg_;
		isCompleted_ = true;
		nPoints_ = 0;
	}
	else
	{
		// new heading is different from the previous one, so execute a "Turn" maneuver
		if (hdg_diff>0.0)
		{
			// target_hdg_ > hdg0_
			if (hdg_diff<180.0)
			{
				// this is R-turn (CW) by angle within (0.001,180)
				turn_sign_ = 1.0;					// positive
			}
			else
			{
				// this is L-turn (CCW) by negative angle in the range [-180,-0.001)
				turn_sign_ = -1.0;
				hdg0_ += 360.0;	// add 360 to the current heading, so it can turn CCW till it reaches target_hdg < hdg0_+360;
			}
		}
		else
		{
			// target_hdg_ < hdg0_
			if (hdg_diff>=-180.0)
			{
				// this is L-turn (CWW) by negative angle in the range [-180,-0.001)
				turn_sign_ = -1.0;
			}
			else
			{
				// this is R-turn (CCW) by 360.0+hdg_diff deg
				turn_sign_ = 1.0;
				hdg0_ -= 360.0; // deduct 360 from the current heading, so it can turn CW till it reaches target_hdg > hdg0_-360;
			}
		}
	}
}




////////////////////////////////////////////////////////////////////////////////////////
// Destructors: deallocate all arrays
////////////////////////////////////////////////////////////////////////////////////////
CManeuver::~CManeuver()
{
	// deallocate arrays (free memory)
	if (pTime_) delete pTime_;
	if (pLon_) delete pLon_;
	if (pLat_) delete pLat_;
	if (pAlt_) delete pAlt_;
	if (pM_) delete pM_;
	if (pu_) delete pu_;
	if (pv_) delete pv_;
	if (pw_) delete pw_;
	if (pq_) delete pq_; 

	if (pSTimes_) delete  pSTimes_;
}


CManeuverTurn::~CManeuverTurn()
{
	if (pPhi_) delete pPhi_;
}



////////////////////////////////////////////////////////////////////////////////////////
//  Verbose for diagnostic
////////////////////////////////////////////////////////////////////////////////////////

void CManeuver::SetVerbose(FILE* outputstream)
{
	// set verbose stream (can be NULL - no output)
	fileout_ = outputstream;
}


////////////////////////////////////////////////////////////////////////////////////////
//  Interface functions
////////////////////////////////////////////////////////////////////////////////////////

// return duration of the maneuver
double CManeuver::GetDuration()
{
	return duration_;
}

// return last heading
double CManeuver::GetLastHdg()
{
	return last_hdg_;
}


// return last time
double CManeuver::GetLastTime()
{
	if (pTime_) 
	{
		return pTime_[nPoints_-1];
	}
	else
	{
		return t0_;
	}
}


// return last longitude
double CManeuver::GetLastLon()
{
	if (pLon_) 
	{
		return pLon_[nPoints_-1];
	}
	else
	{
		return lon0_;
	}
}

// return last latitude
double CManeuver::GetLastLat()
{
	if (pLat_) 
	{
		return pLat_[nPoints_-1];
	}
	else
	{
		return lat0_;
	}
}

// return last altitude
double CManeuver::GetLastAlt()
{
	if (pAlt_) 
	{
		return pAlt_[nPoints_-1];
	}
	else
	{
		return alt0_;
	}
}

// return last weight
double CManeuver::GetLastM()
{
	if (pM_) 
	{
		return pM_[nPoints_-1];
	}
	else
	{
		return M0_;
	}
}


// return number of output points (nPoints_)
unsigned long CManeuver::GetNOut()
{
	return nPoints_;
}

// return pointer to array of times (of length nPoints_)
double* CManeuver::GetTimeOut()
{
	return pTime_;
}

// return pointer to array of longitudes (of length nPoints_)
double* CManeuver::GetLonOut()
{
	return pLon_;
}

// return pointer to array of latitudes (of length nPoints_)
double* CManeuver::GetLatOut()
{
	return pLat_;
}

// return pointer to array of altitudes (of length nPoints_)
double* CManeuver::GetAltOut()
{
	return pAlt_;
}

// return pointer to array of weights (of length nPoints_)
double* CManeuver::GetMOut()
{
	return pM_;
}


// return pointer to array of u-velocity (W->E, m/s, of length nPoints_)
double* CManeuver::GetuOut()
{
	return pu_;
}

// return pointer to array of v-velocity (S->N, m/s, of length nPoints_)
double* CManeuver::GetvOut()
{
	return pv_;
}


// return pointer to array of w-velocity (vertical, m/s, of length nPoints_)
double* CManeuver::GetwOut()
{
	return pw_;
}

// return pointer to array of fuel flow rate (kg/s, of length nPoints_)
double* CManeuver::GetqOut()
{
	return pq_;
}


// check conditions and update VNAV maneuver: returns true if update takes place or false if current maneuver stays 
bool CManeuver::UpdateVNavManeuver(double t, double *XYZMPhi)
{
	return pTrajectory_->UpdateVNavManeuver_(t+t0_, XYZMPhi);
}


// update LNAV maneuver - default return false - current mode stays
bool CManeuver::UpdateLNavManeuver(double t, double *XYZMPhi)
{
	if (t > duration_) return true;
	return false;
}



double CManeuver::GetMaxIntegrationTimeStep(double t, double *XYZM)
{
	// get maximum integration time step based on LNAV profile
	double max_vnav_step = 	pTrajectory_-> GetMaxVNAVStep_(t+t0_, XYZM);

	// return the smallest limit
	return (max_vnav_step < max_time_step_)?max_vnav_step:max_time_step_;
}



// Check if to terminate current turn maneuver (return true), or continue (false) 
bool CManeuverTurn::UpdateLNavManeuver(double t, double *XYZMPhi)
{
	double phi = XYZMPhi[4];
	double ref_hdg = 0.0; // reference heading

	if (t > max_duration_) return true;

	switch(turn_type_)
	{
	case TURN_TYPE_BY_ANGLE:
		if (turn_sign_>0.0 && phi>target_turn_angle_) return true;
		if (turn_sign_<0.0 && phi<target_turn_angle_) return true;
		break;

	case TURN_TYPE_TO_HDG:
		if (turn_sign_>0.0 && phi+hdg0_+ref_hdg0_>target_hdg_+ref_hdg) return true;
		if (turn_sign_<0.0 && phi+hdg0_+ref_hdg0_<target_hdg_+ref_hdg) return true;
		break;
	}
	return false;
}




///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Estimate end time of the maneuver (starting from 0), which is equal to the actual duration
// This is directly related to the termination conditions in UpdateLNavManeuver
//
///////////////////////////////////////////////////////////////////////////////////////////////////

double CManeuverTurn::EstimateEndTurnTime_(double t0, double* XYZMPhi0, 
										   double t1, double* XYZMPhi1)
{
	// Input: t0 - the second last time, when current maneuver is still valid
	//        pXYZM0 - ECEF position, weight, and rotation angle at t0
	//        t1 - the fist time, when current maneuver becomes invalid
	//        pXYZM1 - ECEF position, weight and rotation angle at t1
	//        
	// Output: end time of the maneuver (starting from 0), which is equal to the actual duration

	double phi0 = XYZMPhi0[4];
	double phi1 = XYZMPhi1[4];

	double ref_hdg = 0.0; // reference heading

	if (t1 > max_duration_) return max_duration_;

	switch(turn_type_)
	{
	case TURN_TYPE_BY_ANGLE:
		return t0+(t1-t0)*(target_turn_angle_-phi0)/(phi1-phi0);
		break;

	case TURN_TYPE_TO_HDG:
		return t0+(t1-t0)*(target_hdg_+ref_hdg - (phi0+hdg0_+ref_hdg0_))/(phi1-phi0);
		break;
	}
	return false;
}




///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Integrate trajectory - define derivatives in ECEF system
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// ODE derivatives
void dxyzm(unsigned long n, double t, double *XYZM, double *dXYZM, void* pParam)
{
	CManeuver* pManeuver = (CManeuver*)pParam;	// pointer to maneuver object, which contains all needed parameters

	double X = XYZM[0];		// X in ECEF system
	double Y = XYZM[1];		// Y in ECEF system
	double Z = XYZM[2];		// Z in ECEF system
	double m = XYZM[3];		// weight

	double lon, lat, alt;
	double u, v, w, q, dphi;

	// ---------------
	// convert X, Y, Z to lon, lat, alt
	XYZ2LonLatH(lon, lat, alt, X, Y, Z);

	// ---------------
	// calculate {u,v,w} velocity components (m/s) and angular rotation rate dphi (deg/s) in local coordinate system, and fuel flow rate q (kg/s) 
	pManeuver -> Calc_uvwq(u, v, w, q, dphi, t, lon, lat, alt, m, 0.0);

	// ----------------
	// convert horizontal velocity components to  (X, Y, Z) ECEF system
	double ECEF_HU, ECEF_HV, ECEF_HW;
	uv2uvw(ECEF_HU, ECEF_HV, ECEF_HW, u, v, lon, lat);

	// ----------------
	// convert vertical velocity component to (X, Y, Z) ECEF system
    double ECEF_VU, ECEF_VV, ECEF_VW;
	w2uvw(ECEF_VU, ECEF_VV, ECEF_VW, w, lon, lat);	

	// ---------------
	// failsafe check and final output
	if (alt>=0 && alt<=20000.0)
	{
		dXYZM[0] = ECEF_HU + ECEF_VU; // sum of the horizontal and vertical vector components in ECEF system
		dXYZM[1] = ECEF_HV + ECEF_VV;
		dXYZM[2] = ECEF_HW + ECEF_VW;
		dXYZM[3] = -q;
	}
	else
	{
		dXYZM[0] = 0.0;
		dXYZM[1] = 0.0;
		dXYZM[2] = 0.0;
		dXYZM[3] = 0.0;
	}
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Integrate trajectory - define derivatives in ECEF system for turn maneuvers,which include heading
// Explicit integration allows for changing TAS during turns (e.g., spiral dives), and termination based on reached bearing
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// ODE derivatives
void dxyzmphi(unsigned long n, double t, double *XYZMPhi, double *dXYZMPhi, void* pParam)
{
	CManeuver* pManeuver = (CManeuver*)pParam;	// pointer to maneuver object, which contains all needed parameters

	double X = XYZMPhi[0];		// X in ECEF system (m)
	double Y = XYZMPhi[1];		// Y in ECEF system (m)
	double Z = XYZMPhi[2];		// Z in ECEF system (m)
	double m = XYZMPhi[3];		// weight (kg)
	double phi = XYZMPhi[4];	// change in heading (deg) starting from initial 0.0; positive clockwise, negative counterclockwise

	double lon, lat, alt;
	double u, v, w, q, dphi;

	// ---------------
	// convert X, Y, Z to lon, lat, alt
	XYZ2LonLatH(lon, lat, alt, X, Y, Z);

	// ---------------
	// calculate {u,v,w} velocity components (m/s) and angular rotation rate dphi (deg/s) in local coordinate system, and fuel flow rate q (kg/s) 
	pManeuver -> Calc_uvwq(u, v, w, q, dphi, t, lon, lat, alt, m, phi);

	// ----------------
	// convert horizontal velocity components to  (X, Y, Z) ECEF system
	double ECEF_HU, ECEF_HV, ECEF_HW;
	uv2uvw(ECEF_HU, ECEF_HV, ECEF_HW, u, v, lon, lat);

	// ----------------
	// convert vertical velocity component to (X, Y, Z) ECEF system
    double ECEF_VU, ECEF_VV, ECEF_VW;
	w2uvw(ECEF_VU, ECEF_VV, ECEF_VW, w, lon, lat);	

	// ---------------
	// failsafe check and final output
	if (alt>=0 && alt<=20000.0)
	{
		dXYZMPhi[0] = ECEF_HU + ECEF_VU; // sum of the horizontal and vertical vector components in ECEF system
		dXYZMPhi[1] = ECEF_HV + ECEF_VV;
		dXYZMPhi[2] = ECEF_HW + ECEF_VW;
		dXYZMPhi[3] = -q;
		dXYZMPhi[4] = dphi;
	}
	else
	{
		dXYZMPhi[0] = 0.0;
		dXYZMPhi[1] = 0.0;
		dXYZMPhi[2] = 0.0;
		dXYZMPhi[3] = 0.0;
		dXYZMPhi[4] = 0.0;
	}
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Callback function to terminate integration
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int terminateIntegration(unsigned long n, double t, double *XYZMPhi, void* pParam)
{
	CManeuver* pManeuver = (CManeuver*)pParam;			// pointer to maneuver object

	// check if vertical maneuver changes (e.g., time limit or target level reached)
	// returns true if change took place, or false if current maneuver stays
	bool vmaneuver_change_flag = pManeuver -> UpdateVNavManeuver(t, XYZMPhi);

	if (vmaneuver_change_flag) return TERMINATION_CODE_VNAV;	// terminated integration with code VNAV

	bool hmaneuver_change_flag = pManeuver -> UpdateLNavManeuver(t, XYZMPhi);

	if (hmaneuver_change_flag) return TERMINATION_CODE_LNAV;	// terminated integration with code LNAV

	return TERMINATION_CODE_NONE;	// continue integration
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Callback function to dynamically update maximum time step
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double getMaxTimeStep(unsigned long n, double t, double *XYZM, void* pParam)
{
	CManeuver* pManeuver = (CManeuver*)pParam;	// pointer to maneuver object
	return pManeuver->GetMaxIntegrationTimeStep(t, XYZM);
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Integrate trajectory
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void CManeuver::Integrate()
{
	// this function integrates an element of trajectory based on provided maneuver's (leg, turn, etc.) parameters

	ODE45Options options;
	ODE45Stats stats;

	unsigned neq = 4; // number of components in the ODE system

	double rtoler[] = {1.0E-9, 1.0E-9, 1.0E-9, 1.0E-9};	// relative error tolerance (neq components)
	double atoler[] = {1.0, 1.0, 1.0, 1.0};				// absolute error tolerance (neq components)

	double final_sol[4];								// buffer to contain final computed solution
	double second_last_sol[4];							// buffer to contain second last computed solution

	// ---------------------------
	// Set solver options

	options.nmax = 100000;				// maximal number of allowed steps
	options.nstiff = 0;					// test for stiffness
	options.uround = 0.0;				// rounding unit - set automatically
	options.hmax = max_time_step_;		// maximum step size
	options.h = 0.0;					// initial integration step
	options.safe = 0.0;					// safety factor
	options.fac1 = 0.0;					// parameters for step size selection
	options.fac2 = 0.0;
	options.beta = 0.0;					// for stabilized step size control
	options.rtoler = rtoler;			// relative error tolerance (must be array of length neq)
	options.atoler = atoler;			// absolute error tolerance (must be array of length neq)
	options.fileout = stdout;			// output message stream (can be NULL)
	options.pWorkspace = NULL;			// internally allocate and deallocate workspace 
	options.chkTerm = NULL;				// callback function to check termination conditions (can be NULL)
	options.getMaxStep = NULL;			// callback function to get maximum time step (can be NULL)

	options.yout = final_sol;			// buffer to store final computed solution
	options.yoldout = second_last_sol;	// buffer to store second last computed solution

	options.chkTerm = terminateIntegration;	// callback function to interrupt integration
	options.getMaxStep = getMaxTimeStep;	// callback function to dynamically impose external limit on integration time step

	// ---------------------------
	// deallocate old arrays if needed
	if (pTime_) {delete pTime_; pTime_=NULL;};
	if (pLon_) {delete pLon_; pLon_=NULL;};
	if (pLat_) {delete pLat_; pLat_=NULL;};
	if (pAlt_) {delete pAlt_; pAlt_=NULL;};
	if (pM_) {delete pM_; pM_=NULL;};
	if (pu_) {delete pu_; pu_=NULL;};
	if (pv_) {delete pv_; pv_=NULL;};
	if (pw_) {delete pw_; pw_=NULL;};
	if (pq_) {delete pq_; pq_=NULL;};

	// ---------------------------
	// Set initial conditions (convert lon, lat, alt to ECEF X, Y, Z)
	double solstart[4];					 // initial conditions (neq components)
	LonLatH2XYZ(solstart[0], solstart[1], solstart[2], lon0_, lat0_, alt0_);
	solstart[3] = M0_;

	// ---------------------------
	// Create array of times for output, when solution is desired
	int tlen_regular = (unsigned int) floor(duration_/max_time_step_)+2; // desired time steps for the maneuver
	pTime_ = new double[tlen_regular+nSTimes_+1];	// allocate array to accomodate 'regular' times and specific caller-defined ones

	// set time from zero to duration_, whilst embedding specific ones
	int idx0 = 0;	// regular time index
	int idx1 = 0;	// specific time index
	int tlen = 0;	// total length, tlen_regular <= tlen <= tlen_regular+nSTimes_

	if (nSTimes_>0 && pSTimes_!=NULL)
	{
		// find index of the first specific time >= epsilon_time_, which may need to be added
		// all specific times < 0.0 are before this maneuver, and therefore they are to be ignored
		// 0.0 is added as a regular time (beginning of this maneuver)
		while ((pSTimes_[idx1] < EPS_TIME) && (idx1 < (int)nSTimes_)) idx1++;
	}

	// first time
	pTime_[0] = 0.0;
	idx0++;
	tlen++;

	while (idx0 < tlen_regular)
	{
		double reg_time_next = duration_*((double)idx0)/((double)(tlen_regular-1));		// next regular output time
		double reg_time_prev = duration_*((double)(idx0-1))/((double)(tlen_regular-1));	// previous regular output time

		if (idx1 < (int)nSTimes_)
		{
			double sp_time = pSTimes_[idx1];	// specific time

			if (abs(sp_time-reg_time_prev) < EPS_TIME)
			{
				idx1++;	// just skip; no need to add it because previous regular time output was already close to caller-specified time output
				continue;
			}

			if (abs(sp_time-reg_time_next) < EPS_TIME)
			{
				pTime_[tlen] = reg_time_next;	// do not differentiate between regular and specific times; in a particular case they can be equal
				idx0++;							// increase both the counters
				idx1++;
			}
			else
			{
				if (reg_time_next < sp_time)		// regular output occurs before caller-specified one
				{
					pTime_[tlen] = reg_time_next;	// choose regular time, and increase regular time counter
					idx0++;
				}
				else
				{
					pTime_[tlen] = sp_time;			// choose caller-defined time, and increase specific time counter
					idx1++;
				}
			}
		}
		else
		{
			// no  times left in caller-specidied list, so add regular time and increase regular counter
			pTime_[tlen] = reg_time_next;
			idx0++;
		}

		tlen++;
	}


	// ---------------------------
	double** ppSolution; // array of pointers
	double** ppLon = NULL;
	double** ppLat = NULL;
	double** ppAlt = NULL;
	double** ppM = NULL;
	double** ppu = NULL;
	double** ppv = NULL;
	double** ppw = NULL;
	double** ppq = NULL;

	int* pIdxS = NULL;	// start index
	int* pIdxL = NULL;	// length of each section
	int IDXs = 0;

//	double lastX, lastY, lastZ, lastM;				// last computed values of ECEF X,Y,Z and weight for specified times pTime[]

	double* pTime = new double[tlen+1];
	for (int i=0; i<tlen; i++) pTime[i] = pTime_[i];	// duplicate because this new array could be modified

	int nSection = 0;

	// Iteratively call ODE solver
	while (1)
	{
		// ----------------------------
		// Reallocate all arrays
		double** ppLon_new = new double*[nSection+1];
		double** ppLat_new = new double*[nSection+1];
		double** ppAlt_new = new double*[nSection+1];
		double** ppM_new = new double*[nSection+1];
		double** ppu_new = new double*[nSection+1];
		double** ppv_new = new double*[nSection+1];
		double** ppw_new = new double*[nSection+1];
		double** ppq_new = new double*[nSection+1];
		int* pIdxS_new = new int[nSection+1];
		int* pIdxL_new = new int[nSection+1];

		for (int i=0; i<nSection; i++)
		{
			ppLon_new[i] = ppLon[i];
			ppLat_new[i] = ppLat[i];
			ppAlt_new[i] = ppAlt[i];
			ppM_new[i] = ppM[i];
			ppu_new[i] = ppu[i];
			ppv_new[i] = ppv[i];
			ppw_new[i] = ppw[i];
			ppq_new[i] = ppq[i];
			pIdxS_new[i] = pIdxS[i];
			pIdxL_new[i] = pIdxL[i];
		}

		if (ppLon) delete ppLon;
		if (ppLat) delete ppLat;
		if (ppAlt) delete ppAlt;
		if (ppM) delete ppM;
		if (ppu) delete ppu;
		if (ppv) delete ppv;
		if (ppw) delete ppw;
		if (ppq) delete ppq;
		if (pIdxS) delete pIdxS;
		if (pIdxL) delete pIdxL;

		ppLon = ppLon_new;
		ppLat = ppLat_new;
		ppAlt = ppAlt_new;
		ppM = ppM_new;
		ppu = ppu_new;
		ppv = ppv_new;
		ppw = ppw_new;
		ppq = ppq_new;
		pIdxS = pIdxS_new;
		pIdxL = pIdxL_new;


		// ----------------------------
		// Do integration
		pIdxS[nSection] = IDXs; // first index in this section

		int tlen_s = tlen - IDXs;	// attemp to integrate for the whole length

		ppSolution =  dopri5(neq, dxyzm, &pTime[IDXs], tlen_s, solstart, (void*)this, options, stats);

		// check status
		if (stats.status ==  ODE_STATUS_COMPLETED || ODE_STATUS_INTERRUPTED)
		{
			// End of integration or end of VNAV maneuver - 
			// Finalize output (convert (X,Y,Z) into (lon,lat,alt) and compute velocity)
			// Utilize the same arrays to store (X,Y,Z) and (lon,lat, alt) so as to avoid re-allocation: 
			// this is done by replacing triplets {X,Y,Z} with {lon,lat,alt}
			ppLon[nSection]  = ppSolution[0];
			ppLat[nSection]  = ppSolution[1];
			ppAlt[nSection]  = ppSolution[2];
			ppM[nSection]    = ppSolution[3];

			tlen_s = (int)stats.nout;
			pIdxL[nSection] = tlen_s;		// store length of this section; note last element may need to be discarded if integration was interrupted rather than completed

			for (int i=0; i<tlen_s; i++)
			{
				double lon, lat, alt;
				XYZ2LonLatH(lon,lat,alt, ppSolution[0][i], ppSolution[1][i], ppSolution[2][i]);
				ppLon[nSection][i] = lon;
				ppLat[nSection][i] = lat;
				ppAlt[nSection][i] = alt;
				ppM[nSection][i]= ppSolution[3][i];
			}

			delete ppSolution;	// this array is no longer needed. Note: ppSolution[0..3] arrays stay.

			// Compute velocity components and fuel flow (the former are needed for BFO calculation)
			ppu[nSection] = new double[tlen_s+1];
			ppv[nSection] = new double[tlen_s+1];
			ppw[nSection] = new double[tlen_s+1];
			ppq[nSection] = new double[tlen_s+1]; 

			double dphi; // change of heading rate is not required for output
			for (int i=0; i<tlen_s; i++) Calc_uvwq(ppu[nSection][i], ppv[nSection][i], ppw[nSection][i], ppq[nSection][i], dphi,
												   pTime[IDXs+i], ppLon[nSection][i], ppLat[nSection][i], ppAlt[nSection][i], ppM[nSection][i], 0.0);
		}


		if (stats.status == ODE_STATUS_COMPLETED)
		{
			pIdxL[nSection] = tlen_s;
			nSection++;
			break;	// break while
		}


		if (stats.status == ODE_STATUS_INTERRUPTED)
		{
			if (stats.termcode == TERMINATION_CODE_VNAV) // termination occured becasue of VNAV mode change
			{
				// --------------------------------------
				// Define new time and initial conditions
				double start_t;
				pTrajectory_ -> MoveToNextVNAVManeuver_(start_t, &solstart[0], 
												    stats.xoldout + t0_, &stats.yoldout[0], 
												    stats.xout + t0_, &stats.yout[0]);

				start_t -= t0_;	// adjust time to current LNAV maneuver's time

				IDXs += (tlen_s-1);	// (IDXs + tlen_s -1) is the last output element, after which termination occured

				// Exclude all outputs after start_t time (the latter corresponds to termination time, estimated upon reaching termination condition)
				while (pTime[IDXs]>start_t && tlen_s>0)
				{
					IDXs--;
					tlen_s--;
				}

				pIdxL[nSection] = tlen_s;
				pTime[IDXs] = start_t;	// replace old time with new time


				nSection++;
				continue;
			}
			else
			{
				// termiantion was invoked by LNAV condition, e.g. reached a waypoint
				pIdxL[nSection] = tlen_s;
				nSection++;
				break;	// break while
			}
		}

		// Handle integration failure
		if (stats.status == ODE_STATUS_INCONSISTENT || stats.status == ODE_STATUS_NMAX || stats.status == ODE_STATUS_NOTSTARTED ||
			stats.status == ODE_STATUS_TOOSMALLSTEP || stats.status == ODE_STATUS_STIFF)
		{
			if (fileout_) fprintf(fileout_, "Integration failed before the requested end.\n");
			break;	// break while
		}
	}



	// --------------------------------------------------------------------------------------
	// Finalize output

	if (pTime) delete pTime;	// delete modified array of times, which is no longer needed

	// Add initial time, so that output time is in s since 2014-03-07 00:00:00 UTC
	for (int i=0; i<tlen; i++) pTime_[i] += t0_; 

	// Compute how many output points were computed in the total, which can be less than tlen
	// Note: the first output elements from all the sections except the very first one are to be excluded
	// because these elements were inserted as initial conditions
	nPoints_ = 0;
	for (int n=0; n<nSection; n++)
	{
		int idx_s = pIdxS[n];
		int idx_e = pIdxS[n]+pIdxL[n];
		if (n>0) idx_s = idx_s+1;
		nPoints_ += (idx_e-idx_s);
	}

	pLon_ = new double [nPoints_+1];
	pLat_ = new double [nPoints_+1];
	pAlt_ = new double [nPoints_+1];
	pM_ = new double[nPoints_+1];
	pu_ = new double[nPoints_+1];
	pv_ = new double[nPoints_+1];
	pw_ = new double[nPoints_+1];
	pq_ = new double[nPoints_+1];

	// Copy sections into contiguous arrays and free temporal arrays
	int idx = 0;
	for (int n=0; n<nSection; n++)
	{
		int idx_s = pIdxS[n];
		int idx_e = pIdxS[n]+pIdxL[n];
		int idxs = idx_s;
		if (n>0) idxs = idx_s+1;

		for (int i=idxs; i<idx_e; i++)
		{
			pLon_[idx] = ppLon[n][i-idx_s];
			pLat_[idx] = ppLat[n][i-idx_s];
			pAlt_[idx] = ppAlt[n][i-idx_s];
			pM_[idx] = ppM[n][i-idx_s];
			pu_[idx] = ppu[n][i-idx_s];
			pv_[idx] = ppv[n][i-idx_s];
			pw_[idx] = ppw[n][i-idx_s];
			pq_[idx] = ppq[n][i-idx_s];
			idx++;
		}

		delete ppLon[n];
		delete ppLat[n];
		delete ppAlt[n];
		delete ppM[n];
		delete ppu[n];
		delete ppv[n];
		delete ppw[n];
		delete ppq[n];
	}

	delete ppLon;
	delete ppLat;
	delete ppAlt;
	delete ppM;
	delete ppu;
	delete ppv;
	delete ppw;
	delete ppq;
	delete pIdxS;
	delete pIdxL;

	// last heading does not change
	last_hdg_ = hdg0_;

	// adjust last heading to make it in [0, 360.0)
	while (last_hdg_<0.0) last_hdg_ += 360.0;
	while (last_hdg_>=360.0) last_hdg_ -= 360.0;

};






/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Integrate trajectory - Turn maneuvers,
// which have an additional integration variable - change in heading
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void CManeuverTurn::Integrate()
{
	if (isCompleted_) {nPoints_=0; return;};

	// this function integrates an element of trajectory based on provided maneuver's (leg, turn, etc.) parameters

	ODE45Options options;
	ODE45Stats stats;

	unsigned neq = 5; // number of components in the ODE system

	double rtoler[] = {1.0E-9, 1.0E-9, 1.0E-9, 1.0E-9, 1.0E-9};	// relative error tolerance (neq components)
	double atoler[] = {1.0, 1.0, 1.0, 1.0, 1.0E-3};				// absolute error tolerance (neq components)

	double final_sol[5];										// buffer to contain final computed solution
	double second_last_sol[5];									// buffer to contain second last computed solution

	// ---------------------------
	// Set solver options

	options.nmax = 100000;				// maximal number of allowed steps
	options.nstiff = 0;					// test for stiffness
	options.uround = 0.0;				// rounding unit - set automatically
	options.hmax = max_time_step_;		// maximum step size
	options.h = 0.0;					// initial integration step
	options.safe = 0.0;					// safety factor
	options.fac1 = 0.0;					// parameters for step size selection
	options.fac2 = 0.0;
	options.beta = 0.0;					// for stabilized step size control
	options.rtoler = rtoler;			// relative error tolerance (must be array of length neq)
	options.atoler = atoler;			// absolute error tolerance (must be array of length neq)
	options.fileout = stdout;			// output message stream (can be NULL)
	options.pWorkspace = NULL;			// internally allocate and deallocate workspace 
	options.chkTerm = NULL;				// callback function to check termination conditions (can be NULL)
	options.getMaxStep = NULL;			// callback function to get maximum time step (can be NULL)

	options.yout = final_sol;			// buffer to store final computed solution
	options.yoldout = second_last_sol;	// buffer to store second last computed solution

	options.chkTerm = terminateIntegration;	// callback function to interrupt integration
	options.getMaxStep = getMaxTimeStep;	// callback function to dynamically impose external limit on integration time step

	// ---------------------------
	// deallocate old arrays if needed
	if (pTime_) {delete pTime_; pTime_=NULL;};
	if (pLon_) {delete pLon_; pLon_=NULL;};
	if (pLat_) {delete pLat_; pLat_=NULL;};
	if (pAlt_) {delete pAlt_; pAlt_=NULL;};
	if (pM_) {delete pM_; pM_=NULL;};
	if (pu_) {delete pu_; pu_=NULL;};
	if (pv_) {delete pv_; pv_=NULL;};
	if (pw_) {delete pw_; pw_=NULL;};
	if (pq_) {delete pq_; pq_=NULL;};
	if (pPhi_) {delete pPhi_; pPhi_=NULL;};

	// ---------------------------
	// Set initial conditions (convert lon, lat, alt to ECEF X, Y, Z)
	double solstart[5];					 // initial conditions (neq components)
	LonLatH2XYZ(solstart[0], solstart[1], solstart[2], lon0_, lat0_, alt0_);
	solstart[3] = M0_;
	solstart[4] = 0.0;	// change in heading (total heading = hdg0 + Phi, with heading interpreation depending on mode, e.g. true or magnetic)

	// ---------------------------
	// Create array of times for output, when solution is desired
	int tlen_regular = (unsigned int) floor(duration_/max_time_step_)+2; // desired time steps for the maneuver
	pTime_ = new double[tlen_regular+nSTimes_+1];	// allocate array to accomodate 'regular' times and specific caller-defined ones

	// set time from zero to duration_, whilst embedding specific ones
	int idx0 = 0;	// regular time index
	int idx1 = 0;	// specific time index
	int tlen = 0;	// total length, tlen_regular <= tlen <= tlen_regular+nSTimes_

	if (nSTimes_>0 && pSTimes_!=NULL)
	{
		// find index of the first specific time >= epsilon_time_, which may need to be added
		// all specific times < 0.0 are before this maneuver, and therefore they are to be ignored
		// 0.0 is added as a regular time (beginning of this maneuver)
		while ((pSTimes_[idx1] < EPS_TIME) && (idx1 < (int)nSTimes_)) idx1++;
	}

	// first time
	pTime_[0] = 0.0;
	idx0++;
	tlen++;

	while (idx0 < tlen_regular)
	{
		double reg_time_next = duration_*((double)idx0)/((double)(tlen_regular-1));		// next regular output time
		double reg_time_prev = duration_*((double)(idx0-1))/((double)(tlen_regular-1));	// previous regular output time

		if (idx1 < (int)nSTimes_)
		{
			double sp_time = pSTimes_[idx1];	// specific time

			if (abs(sp_time-reg_time_prev) < EPS_TIME)
			{
				idx1++;	// just skip; no need to add it because previous regular time output was already close to caller-specified time output
				continue;
			}

			if (abs(sp_time-reg_time_next) < EPS_TIME)
			{
				pTime_[tlen] = reg_time_next;	// do not differentiate between regular and specific times; in a particular case they can be equal
				idx0++;							// increase both the counters
				idx1++;
			}
			else
			{
				if (reg_time_next < sp_time)		// regular output occurs before caller-specified one
				{
					pTime_[tlen] = reg_time_next;	// choose regular time, and increase regular time counter
					idx0++;
				}
				else
				{
					pTime_[tlen] = sp_time;			// choose caller-defined time, and increase specific time counter
					idx1++;
				}
			}
		}
		else
		{
			// no  times left in caller-specidied list, so add regular time and increase regular counter
			pTime_[tlen] = reg_time_next;
			idx0++;
		}

		tlen++;
	}


	// ---------------------------
	double** ppSolution; // array of pointers
	double** ppLon = NULL;
	double** ppLat = NULL;
	double** ppAlt = NULL;
	double** ppM = NULL;
	double** ppu = NULL;
	double** ppv = NULL;
	double** ppw = NULL;
	double** ppq = NULL;
	double** ppPhi = NULL;

	int* pIdxS = NULL;	// start index
	int* pIdxL = NULL;	// length of each section
	int IDXs = 0;


	double* pTime = new double[tlen+1];
	for (int i=0; i<tlen; i++) pTime[i] = pTime_[i];	// duplicate because this new array could be modified

	int nSection = 0;

	// Iteratively call ODE solver
	while (1)
	{
		// ----------------------------
		// Reallocate all arrays
		double** ppLon_new = new double*[nSection+1];
		double** ppLat_new = new double*[nSection+1];
		double** ppAlt_new = new double*[nSection+1];
		double** ppM_new = new double*[nSection+1];
		double** ppu_new = new double*[nSection+1];
		double** ppv_new = new double*[nSection+1];
		double** ppw_new = new double*[nSection+1];
		double** ppq_new = new double*[nSection+1];
		double** ppPhi_new = new double*[nSection+1];
		int* pIdxS_new = new int[nSection+1];
		int* pIdxL_new = new int[nSection+1];

		for (int i=0; i<nSection; i++)
		{
			ppLon_new[i] = ppLon[i];
			ppLat_new[i] = ppLat[i];
			ppAlt_new[i] = ppAlt[i];
			ppM_new[i] = ppM[i];
			ppu_new[i] = ppu[i];
			ppv_new[i] = ppv[i];
			ppw_new[i] = ppw[i];
			ppq_new[i] = ppq[i];
			ppPhi_new[i] = ppPhi[i];
			pIdxS_new[i] = pIdxS[i];
			pIdxL_new[i] = pIdxL[i];
		}

		if (ppLon) delete ppLon;
		if (ppLat) delete ppLat;
		if (ppAlt) delete ppAlt;
		if (ppM) delete ppM;
		if (ppu) delete ppu;
		if (ppv) delete ppv;
		if (ppw) delete ppw;
		if (ppq) delete ppq;
		if (ppPhi) delete ppPhi;
		if (pIdxS) delete pIdxS;
		if (pIdxL) delete pIdxL;

		ppLon = ppLon_new;
		ppLat = ppLat_new;
		ppAlt = ppAlt_new;
		ppM = ppM_new;
		ppu = ppu_new;
		ppv = ppv_new;
		ppw = ppw_new;
		ppq = ppq_new;
		ppPhi = ppPhi_new;
		pIdxS = pIdxS_new;
		pIdxL = pIdxL_new;


		// ----------------------------
		// Do integration
		pIdxS[nSection] = IDXs; // first index in this section

		int tlen_s = tlen - IDXs;	// attemp to integrate for the whole length

		ppSolution =  dopri5(neq, dxyzmphi, &pTime[IDXs], tlen_s, solstart, (void*)this, options, stats);

		// check status
		if (stats.status ==  ODE_STATUS_COMPLETED || ODE_STATUS_INTERRUPTED)
		{
			// End of integration or end of VNAV maneuver - 
			// Finalize output (convert (X,Y,Z) into (lon,lat,alt) and compute velocity)
			// Utilize the same arrays to store (X,Y,Z) and (lon,lat, alt) so as to avoid re-allocation: 
			// this is done by replacing triplets {X,Y,Z} with {lon,lat,alt}
			ppLon[nSection]  = ppSolution[0];
			ppLat[nSection]  = ppSolution[1];
			ppAlt[nSection]  = ppSolution[2];
			ppM[nSection]    = ppSolution[3];
			ppPhi[nSection]  = ppSolution[4];

			tlen_s = (int)stats.nout;
			pIdxL[nSection] = tlen_s;		// store length of this section; note last element may need to be discarded if integration was interrupted rather than completed

			for (int i=0; i<tlen_s; i++)
			{
				double lon, lat, alt;
				XYZ2LonLatH(lon,lat,alt, ppSolution[0][i], ppSolution[1][i], ppSolution[2][i]);
				ppLon[nSection][i] = lon;
				ppLat[nSection][i] = lat;
				ppAlt[nSection][i] = alt;
				ppM[nSection][i]   = ppSolution[3][i];
				ppPhi[nSection][i] = ppSolution[4][i];
			}

			delete ppSolution;	// this array is no longer needed. Note: ppSolution[0..3] arrays stay.

			// Compute velocity components and fuel flow (the former are needed for BFO calculation)
			ppu[nSection] = new double[tlen_s+1];
			ppv[nSection] = new double[tlen_s+1];
			ppw[nSection] = new double[tlen_s+1];
			ppq[nSection] = new double[tlen_s+1]; 

			double dphi; // it is auxilary parameter (rate of heading change, which is not required to be saved)

			for (int i=0; i<tlen_s; i++) Calc_uvwq(ppu[nSection][i], ppv[nSection][i], ppw[nSection][i], ppq[nSection][i], dphi, 
												   pTime[IDXs+i], ppLon[nSection][i], ppLat[nSection][i], ppAlt[nSection][i], ppM[nSection][i], ppPhi[nSection][i]);
		}


		if (stats.status == ODE_STATUS_COMPLETED)
		{
			pIdxL[nSection] = tlen_s;
			nSection++;
			break;	// break while
		}


		if (stats.status == ODE_STATUS_INTERRUPTED)
		{
			if (stats.termcode == TERMINATION_CODE_VNAV) // termination occured becasue of VNAV mode change
			{
				// --------------------------------------
				// Define new time and initial conditions
				double start_t;
				pTrajectory_ -> MoveToNextVNAVManeuver_(start_t, &solstart[0], 
													    stats.xoldout + t0_, &stats.yoldout[0], 
													    stats.xout + t0_, &stats.yout[0]);

				start_t -= t0_;	// adjust time to current LNAV maneuver's time

				IDXs += (tlen_s-1);	// (IDXs + tlen_s -1) is the last output element, after which termination occured

				// Exclude all outputs after start_t time (the latter corresponds to termination time, estimated upon reaching termination condition)
				while (pTime[IDXs]>start_t && tlen_s>0)
				{
					IDXs--;
					tlen_s--;
				}

				pIdxL[nSection] = tlen_s;
				pTime[IDXs] = start_t;	// replace old time with new time


				nSection++;
				continue;
			}
			else
			{
				// termiantion was invoked by LNAV condition, e.g. reached a waypoint
//				pIdxL[nSection] = tlen_s;
//				nSection++;

				// --------------------------------------
				// Define new time and initial conditions by linear interpolation at the time, when the termination condition is expected to be met
				// (otherwise non-linear equation needs to be solved)

				double t_end = EstimateEndTurnTime_(stats.xoldout, &stats.yoldout[0], 
													stats.xout, &stats.yout[0]);

				double coef = (t_end - stats.xoldout)/(stats.xout - stats.xoldout);

				double lon_end0, lat_end0, alt_end0, lon_end1, lat_end1, alt_end1;
				XYZ2LonLatH(lon_end0,lat_end0,alt_end0, stats.yoldout[0], stats.yoldout[1], stats.yoldout[2]);
				XYZ2LonLatH(lon_end1,lat_end1,alt_end1, stats.yout[0], stats.yout[1], stats.yout[2]);

				IDXs += (tlen_s-1);	// (IDXs + tlen_s -1) is the last output element, after which termination occured

				// Exclude all outputs after start_t time (the latter corresponds to termination time, estimated upon reaching termination condition)
				while (pTime[IDXs]>t_end && tlen_s>0)
				{
					IDXs--;
					tlen_s--;
				}

				pIdxL[nSection] = tlen_s;
				pTime[IDXs] = t_end;	// replace old time with new time

				ppLon[nSection][IDXs] = lon_end0 + (lon_end1-lon_end0)*coef;
				ppLat[nSection][IDXs] = lat_end0 + (lat_end1-lat_end0)*coef;
				ppAlt[nSection][IDXs] = alt_end0 + (alt_end1-alt_end0)*coef;
				ppM[nSection][IDXs]   = stats.yoldout[3] + (stats.yout[3] - stats.yoldout[3])*coef;
				ppPhi[nSection][IDXs] = stats.yoldout[4] + (stats.yout[4] - stats.yoldout[4])*coef;

				double dphi; // angular rotation rate - not required to save;

				Calc_uvwq(ppu[nSection][IDXs], ppv[nSection][IDXs], ppw[nSection][IDXs], ppq[nSection][IDXs], dphi, 
						  pTime[IDXs], ppLon[nSection][IDXs], ppLat[nSection][IDXs], ppAlt[nSection][IDXs], ppM[nSection][IDXs], ppPhi[nSection][IDXs]);

				nSection++;

				break;	// break while
			}
		}

		// Handle integration failure
		if (stats.status == ODE_STATUS_INCONSISTENT || stats.status == ODE_STATUS_NMAX || stats.status == ODE_STATUS_NOTSTARTED ||
			stats.status == ODE_STATUS_TOOSMALLSTEP || stats.status == ODE_STATUS_STIFF)
		{
			if (fileout_) fprintf(fileout_, "Integration failed before the requested end.\n");
			break;	// break while
		}
	}



	// --------------------------------------------------------------------------------------
	// Finalize output

	if (pTime) delete pTime;	// delete modified array of times, which is no longer needed

	// Add initial time, so that output time is in s since 2014-03-07 00:00:00 UTC
	for (int i=0; i<tlen; i++) pTime_[i] += t0_; 

	// Compute how many output points were computed in the total, which can be less than tlen
	// Note: the first output elements from all the sections except the very first one are to be excluded
	// because these elements were inserted as initial conditions
	nPoints_ = 0;
	for (int n=0; n<nSection; n++)
	{
		int idx_s = pIdxS[n];
		int idx_e = pIdxS[n]+pIdxL[n];
		if (n>0) idx_s = idx_s+1;
		nPoints_ += (idx_e-idx_s);
	}

	pLon_ = new double [nPoints_+1];
	pLat_ = new double [nPoints_+1];
	pAlt_ = new double [nPoints_+1];
	pM_ = new double[nPoints_+1];
	pu_ = new double[nPoints_+1];
	pv_ = new double[nPoints_+1];
	pw_ = new double[nPoints_+1];
	pq_ = new double[nPoints_+1];
	pPhi_ = new double[nPoints_+1];

	// Copy sections into contiguous arrays and free temporal arrays
	int idx = 0;
	for (int n=0; n<nSection; n++)
	{
		int idx_s = pIdxS[n];
		int idx_e = pIdxS[n]+pIdxL[n];
		int idxs = idx_s;
		if (n>0) idxs = idx_s+1;

		for (int i=idxs; i<idx_e; i++)
		{
			pLon_[idx] = ppLon[n][i-idx_s];
			pLat_[idx] = ppLat[n][i-idx_s];
			pAlt_[idx] = ppAlt[n][i-idx_s];
			pM_[idx] = ppM[n][i-idx_s];
			pu_[idx] = ppu[n][i-idx_s];
			pv_[idx] = ppv[n][i-idx_s];
			pw_[idx] = ppw[n][i-idx_s];
			pq_[idx] = ppq[n][i-idx_s];
			pPhi_[idx] = ppPhi[n][i-idx_s];
			idx++;
		}

		delete ppLon[n];
		delete ppLat[n];
		delete ppAlt[n];
		delete ppM[n];
		delete ppu[n];
		delete ppv[n];
		delete ppw[n];
		delete ppq[n];
		delete ppPhi[n];
	}

	delete ppLon;
	delete ppLat;
	delete ppAlt;
	delete ppM;
	delete ppu;
	delete ppv;
	delete ppw;
	delete ppq;
	delete ppPhi;
	delete pIdxS;
	delete pIdxL;

	// Save last heading
	last_hdg_ = hdg0_ + pPhi_[idx-1];

	// Adjust last heading to make it in [0, 360.0) range
	while (last_hdg_<0.0) last_hdg_ += 360.0;
	while (last_hdg_>=360.0) last_hdg_ -= 360.0;

};