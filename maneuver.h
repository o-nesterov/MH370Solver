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
#include "environdata.h"

#ifndef CMANEUVER_CLASS

#define CMANEUVER_CLASS

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Class to handle horizontal and vertical maneuvers: 'straight' (e.g., true/magnetic heading/tracks), turns, ...
// Each maneuver is characterized by constant flight settings, e.g., constant bank angle during turns, constant descent rate
// Each maneuver forms an element of trajectory
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define SPEED_MODE_IAS 0
#define SPEED_MODE_MACH 1

#define HEADING_MODE_TRUEHDG 0
#define HEADING_MODE_MAGNHDG 1
#define HEADING_MODE_TRUETRK 2
#define HEADING_MODE_MAGNTRK 3
#define HEADING_MODE_GYROHDG 4
#define HEADING_MODE_BEARING 5

#define VNAV_MODE_LEVEL 0
#define VNAV_MODE_VS 1
#define VNAV_MODE_FPA 2
#define VNAV_MODE_FLCH 3

#define ENGINE_MODE_SINGLE 0
#define ENGINE_MODE_DUAL 1

// -------------------------------
// The whole trajectory consists of maneuvers
// Currently there are 3 maneuvers supported: 
// leg/LEG/Leg - keep previous heading/track, true/magnetic 
// turn/Turn/TURN - turn by specified angle at specified bank angle
// hleg/HLEG/Hleg - keep specified heading/track, true/magnetic; insert turn maneuver if needed
// More can be added

#define MANEUVER_TYPE_LEG 0
#define MANEUVER_TYPE_TURN 1
#define MANEUVER_TYPE_HLEG 2

// Each maneuver is characterized by its own set of parameters:
// Leg: heading/track (deg from N, clockwise) and duration (seconds)
// Turn: turn angle (clockwise from the present heading/track), and bank angle

// Maximum number of parameters (predefined)
#define MAX_MANEUVER_IPARAMS 4
#define MAX_MANEUVER_FPARAMS 8

// Turn types
// TURN_TYPE_BY_ANGLE - turn by angle
// TURN_TYPE_TO_ANGLE - turn to specified heading
#define TURN_TYPE_BY_ANGLE 0
#define TURN_TYPE_TO_HDG 1



struct LNAVProfile
{
	int maneuverType;								// type
	int iparams[MAX_MANEUVER_IPARAMS];				// integer parameters to pass to maneuver path integrator
	double fparams[MAX_MANEUVER_FPARAMS];			// floating point parameters to pass to maneuver path integrator
	bool ofparams[MAX_MANEUVER_FPARAMS];			// optimization of floating point parameters: YES/NO (used only for optimization)
	bool isOrig;									// flag indicating if the maneuver is original (true if specified in lnav file) or automatically inserted (false) - needed to save optimized file
};


struct VNAVProfile
{
	int maneuverType;								// type: VNAV_MODE_FL, VNAV_MODE_VS, VNAV_MODE_FPA, VNAV_MODE_FLCH
	double startAltitude;							// starting altitude
	double targetAltitude;							// target altitude
	double param;									// parameter: VNAV_MODE_FL - duration (s); VNAV_MODE_VS - vertical speed (m/s); VNAV_MODE_FPA - flight path angle (deg)
	double startTime;								// start vmaneuver time (s since 2014-03-07 00:00:00) - start and end are dynamically updated because in FPA mode duration is unknown
	double endTime;									// end vmaneuver time (s since 2014-03-07 00:00:00) 
};


// -------------------------------
// Aircraft limits
// maximum operating speed (knots)
#define AIRCRAFT_VMO 330.0

// maximum speed for stability characteristics B777F - "Boeing 777F First Flight Workshop in Vienna, Austria 2009"  (knots)
#define AIRCRAFT_VFC 357.0

// maximum operating mach
#define AIRCRAFT_MMO 0.87

//maximum mach for stability characteristics B777F - "Boeing 777F First Flight Workshop in Vienna, Austria 2009"
#define AIRCRAFT_MFC 0.90


// difference to distinguish between regular and specific time outputs (s)
#define EPS_TIME 1.0/(1048576.0)



// Declaration of CTrajectory class (called from Maneuver)
class CTrajectory;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// CManeuver class
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class CManeuver
{
public:
	// Constructor
	// t0 - start time (since 2014-03-07 00:00:00 UTC)
	// tmax - maximum time (since 2014-03-07 00:00:00 UTC); maneuver time can be shorter than this
	// timestep - maximum integration time step;
	// pSTimes - array of specific times (s since 2014-03-07 00:00:00 UTC) when output is required (can be NULL)
	// nSTime - number of specific outputs required (can be 0)
	// {lon0, lat0, alt0, M0} - initial longitude (deg E), latitude (deg N); altitude (m) and weight (kg);
	// speed mode: either SPEED_MODE_IAS or SPEED_MODE_MACH
	// airspeed - air speed parameter (IAS or MACH)
	// heading_mode: HEADING_MODE_TRUEHDG - true heading, HEADING_MODE_MAGNHDG - magnetic heading, HEADING_MODE_TRUETRK - true track, HEADING_MODE_MAGNTRK - magnetic track
	// hdg0 - initial heading (degrees), with interpretation depending on mode
	// pEnvironData - pointer to object providing environmental data (meteorology and magnetic declination)
	// pTrajectory - pointer to the whole trajectory object, which is needed to provide VNAV 

	CManeuver(double t0, double tmax, double timestep, double* pSTimes, unsigned long nStimes,
			  double lon0, double lat0, double alt0, double M0, 
			  int speed_mode, double airspeed, int heading_mode, double hdg0,
			  CMH370EnvironData* pEnvironData, CTrajectory* pTrajectory);

	~CManeuver(); // destructor


	// Integrate path and weigth - core function
	virtual void Integrate(); 

	// ------------------------------------

	// Interface functions:

	// return duration of the maneuver
	double GetDuration();

	// return last time
	double GetLastTime();

	// return last heading
	double GetLastHdg();

	// return last longitude
	double GetLastLon();

	// return last latitude
	double GetLastLat();

	// return last altitude
	double GetLastAlt();

	// return last weight
	double GetLastM();

	// return number of output points (nPoints_)
	unsigned long GetNOut();

	// return pointer to array of times (s since 2014-03-07 00:00:00, of length nPoints_)
	double* GetTimeOut();

	// return pointer to array of longitudes (deg E, of length nPoints_)
	double* GetLonOut();

	// return pointer to array of latitudes (deg N, of length nPoints_)
	double* GetLatOut();

	// return pointer to array of altitudes (m, of length nPoints_)
	double* GetAltOut();

	// return pointer to array of weights (of length nPoints_, kg)
	double* GetMOut();

	// return pointer to array of u (W->E velocity component of length nPoints_, m/s)
	double* GetuOut();

	// return pointer to array of v (S->N velocity component of length nPoints_, m/s)
	double* GetvOut();

	// return pointer to array of w (vertical velocity component of length nPoints_, m/s)
	double* GetwOut();

	// return pointer to array of q (flow rate of length nPoints_, kg/s)
	double* GetqOut();

	// ------------------------------------
	// Get maximum temporal integration time step, which may improve transition between flight modes
	double GetMaxIntegrationTimeStep(double t, double *XYZM);

	// ------------------------------------
	// update VNAV maneuver: returns true if update took place or false if current mode stays 
	bool UpdateVNavManeuver(double t, double *XYZMPhi);

	// ------------------------------------
	// update LNAV maneuver: returns true if update took place or false if current mode stays 
	virtual bool UpdateLNavManeuver(double t, double *XYZMPhi);

	// ------------------------------------

	// calculate u (W->E), v (S->N), w (vertical) velocity components in the local system, change in heading (deg), and fuel flow rate q
	// depending on maneuver type (it can be default, or maneuver-specific overloaded function)
	virtual void Calc_uvwq(double& u, double& v, double& w, double& q, double& dphi,
						   double t, double lon, double lat, double alt, double m, double phi);

	// ------------------------------------
	// set diagnostic stream (can be NULL)
	void SetVerbose(FILE* outputstream);


protected:
	// internal parameters
	double duration_;				// maneuver duration, either specified or inferred 
	double last_hdg_;				// last heading

	// gyroscopic parameters (used in HEADING_MODE_GYROHDG mode)
	double gyro_heading_param1_;
	double gyro_heading_param2_;

protected:	
	// Airspeed mode
	int SPEED_MODE_;				// 0 - IAS, 1 - MACH
	double airspeed_;				// air speed parameter (IAS or MACH)

	// Heading mode
	int HEADING_MODE_;				// 0 - true heading, 1 - magnetic heading, 2 - true track, 3 - magnetic track

	// store initial/reference parameters 
	double t0_;						// initial time (s) since 2014-03-07 00:00:00 UTC 
	double tmax_;					// tmax - maximum time (since 2014-03-07 00:00:00 UTC); maneuver time can be shorter than this
	double max_time_step_;			// user-specified maximum integration time step (s)

	double lon0_;					// initial longitude
	double lat0_;					// initial latitude
	double alt0_;					// initial altitude
	double M0_;						// initial weight
	double hdg0_;					// initial heading parameter (true/magnetic heading or true/magnetic track)
	CMH370EnvironData* pEnvironData_; // pointer on object providing environmental data 
	CTrajectory* pTrajectory_;		// pointer to trajectory object, which provides VNAV profile independent on LNAV

	// output
	unsigned long nPoints_;			// number of output trajectory points, which depend on integration time step and duration
	double* pTime_;					// array of time stamps (s) since 2014-03-07 00:00:00 UTC
	double* pLon_;					// array or longitudes
	double* pLat_;					// array of latitudes
	double* pAlt_;					// array of altitudes
	double* pM_;					// array of weights

	double* pu_;					// u-velocity component in local system (m/s, W->E, output)
	double* pv_;					// v-velocity component in local system (m/s, S->N, output)
	double* pw_;					// w-velocity component in local system (m/s, vertical, output)
	double* pq_;					// fuel flow rate (kg/s, output)

	unsigned long nSTimes_;			// number of specific caller-defined time outputs
	double* pSTimes_;				// list of specific caller-defined output times

	// Verbose diganostic stream
	FILE* fileout_; // message stream for diagnostic output (can be NULL) 
};



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Straight leg maneuver - its properties can be handled by the parent class
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class CManeuverLeg : public CManeuver
{
public:
	// Constructor
	// t0 - start time (since 2014-03-07 00:00:00 UTC)
	// duration
	// timestep - maximum integration time step;
	// pSTimes - array of specific times (s since 2014-03-07 00:00:00 UTC) when output is required (can be NULL)
	// nSTime - number of specific outputs required (can be 0)
	// {lon0, lat0, alt0, M0} - initial longitude (deg E), latitude (deg N); altitude (m) and weight (kg);
	// speed mode: either SPEED_MODE_IAS or SPEED_MODE_MACH
	// airspeed - air speed parameter (IAS or MACH)
	// heading_mode: HEADING_MODE_TRUEHDG - true heading, HEADING_MODE_MAGNHDG - magnetic heading, HEADING_MODE_TRUETRK - true track, HEADING_MODE_MAGNTRK - magnetic track
	// hdg0 - heading (degrees), with interpretation depending on mode
	// pEnvironData - pointer to object providing environmental data (meteorology and magnetic declination)
	// pTrajectory - pointer to the whole trajectory object, which is needed to provide VNAV 

	CManeuverLeg(double t0, double duration, double timestep, double* pSTimes, unsigned long nStimes,
				 double lon0, double lat0, double alt0, double M0, 
				 int speed_mode, double airspeed, int heading_mode, double hdg0,
				 CMH370EnvironData* pEnvironData, CTrajectory* pTrajectory);

	// same as above, but in gyroscopic mode
	CManeuverLeg(double t0, double duration, double timestep, double* pSTimes, unsigned long nStimes,
				 double lon0, double lat0, double alt0, double M0, 
				 int speed_mode, double airspeed, double hdg0, double gyro_hdg_param1, double gyro_hdg_param2,
				 CMH370EnvironData* pEnvironData, CTrajectory* pTrajectory);
};



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Turn maneuver
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class CManeuverTurn : public CManeuver
{
public:
	// Constructor
	// t0 - start time (since 2014-03-07 00:00:00 UTC)
	// timestep - maximum integration time step;
	// pSTimes - array of specific times (s since 2014-03-07 00:00:00 UTC) when output is required (can be NULL)
	// nSTime - number of specific outputs required (can be 0)
	// {lon0, lat0, alt0, M0} - initial longitude (deg E), latitude (deg N); altitude (m) and weight (kg);
	// speed mode: either SPEED_MODE_IAS or SPEED_MODE_MACH
	// airspeed - air speed parameter (IAS or MACH)
	// heading_mode: HEADING_MODE_TRUEHDG - true heading, HEADING_MODE_MAGNHDG - magnetic heading, HEADING_MODE_TRUETRK - true track, HEADING_MODE_MAGNTRK - magnetic track
	// hdg0 - heading (degrees), with interpretation depending on mode
	// turn_angle - turn angle (degrees)
	// bank_angle - maximum bank angle (degrees), typically 5, 10, 15, 20 or 25 deg on selector
	// turn_type - turn by specified angle or to specified heading, or bearing
	// pEnvironData - pointer to object providing environmental data (meteorology and magnetic declination)
	// pTrajectory - pointer to the whole trajectory object, which is needed to provide VNAV 

	CManeuverTurn(double t0, double tmax, double timestep, double* pSTimes, unsigned long nStimes,
				  double lon0, double lat0, double alt0, double M0, 
				  int speed_mode, double airspeed, int heading_mode, double hdg0,
				  double turn_angle, double target_hdg, double bank_angle,
				  CMH370EnvironData* pEnvironData, CTrajectory* pTrajectory);

	// same as above, but in gyroscopic mode
	CManeuverTurn(double t0, double tmax, double timestep, double* pSTimes, unsigned long nStimes,
				  double lon0, double lat0, double alt0, double M0, 
				  int speed_mode, double airspeed, double hdg0, double gyro_hdg_param1, double gyro_hdg_param2,
				  double turn_angle, double target_hdg, double bank_angle,
				  CMH370EnvironData* pEnvironData, CTrajectory* pTrajectory);

	// destructor (to deallocate additional arrays)
	~CManeuverTurn();

	void Calc_uvwq(double& u, double& v, double& w, double& q, double& dphi,
				  double t, double lon, double lat, double alt, double m, double phi);

	// update Turn maneuver: returns true if end is reached or false if current mode stays 
	bool UpdateLNavManeuver(double t, double *XYZMPhi);

	void Integrate(); 

private:

	void init_(double current_hdg, double turn_angle, double target_hdg);

	// estimate end time (local maneuver's time starting from 0) when termination conditions are met, == actual maneuver's duration
	double EstimateEndTurnTime_(double t0, double* pXYZM0, 
								double t1, double* pXYZM1);

	// additional parameters
	double target_turn_angle_;		// turn by specified angle (deg)
	double target_hdg_;				// turn to specified heading (deg)
	double ref_hdg0_;				// initial reference heading to N (==initial magnetic declination or gyro)

	int turn_type_;					// turn type: TURN_TYPE_BY_ANGLE; TURN_TYPE_TO_HDG - turn to specified heading

	double bank_angle_;				// bank angle limit
	double turn_sign_;				// +1.0 turn clockwise, -1.0 turn counterclockwise

	double max_duration_;			// max duration

	bool isCompleted_;				// true if turn is completed

	double* pPhi_;					// heading change (deg, output) - array of size nPoints_; total heading = hdg0+Phi_[..], with interpretation depending on mode

};





#endif
