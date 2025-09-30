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
#include "maneuver.h"
#include "environdata.h"
#include "params_parser.h"


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Class to build trajectory consisting of the series of maneuvers
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef CTRAJECTORY_CLASS

#define CTRAJECTORY_CLASS


class CTrajectory
{
	friend class CManeuver;
	friend class CManeuverTurn;

public:

	CTrajectory();		// constructor: working folder and main setup file
	~CTrajectory(); // destructor

	// create/read lateral navigation profile from specified file
	bool CreateLNAVProfile(CInpParams* pInpParams);

	// create/read vertical navigation profile from specified file
	bool CreateVNAVProfile(CInpParams* pInpParams);

	// set specific times (s since 2014-03-07 00:00:00), when output is required, e.g. when pings are available
	void SetSpecificOutputTimes(double* pTimes, unsigned long n, double AES_Delay);

	// get output at specific time, which must be previously specified using SetSpecificOutputTimes
	void GetSpecificOutput(double& lon, double& lat, double& alt, double& M, double& u, double& v, double& w, double& q, unsigned long n, bool isDelayed);

	// overwrite LNAV optimization parameters (used for optimization only) that were marked with * in LNAV profile
	void ResetLNAVOptimizationParameters(double* params, int nparams);

	// get initial LNAV parameters for optimization (used for optimization only); also returns number of parameters nparams
	void GetLNAVOptimizationParameters(double* params, int& nparams, int maxparams);

	// integrate path - core function
	void BuildTrajectory(double t0, double endtime,
						 double lon0, double lat0, double alt0, double M0,
						 int speed_mode, double airspeed, int heading_mode, double hdg0, int hdg0_type,
						 CMH370EnvironData* pEnvironData,
						 double gyro_hdg_param1 = 71.0, double gyro_hdg_param2 = 0.0);							

	// ------------------------------------
	// Interface functions
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

	// set/reset dafault bank angle to be used to connect legs if turn is not explicitly specified
	// if force_default==true, override all bank angles explicitly specified in lnav file, which are not subjected to individual optimization (marked with *) 
	void SetDefaultBankAngle(double bank_angle, bool force_default);

	// set engine mode: single or dual (both) for fuel flow calculation
	void SetEngineModeSingle();
	void SetEngineModeDual();
	bool IsEngineModeSingle();	// returns true if single engine, or false if dual

	// ------------------------------------
	// set diagnostic stream (can be NULL)
	void SetVerbose(FILE* outputstream);

	// ------------------------------------
	// save trajectory; if filename = NULL, then print it into outputstream
	bool SaveTrajectory(char* filename);

	// ------------------------------------
	// save lateral navigation file (to be called after optimization with lnav flag)
	bool SaveLNAVProfile(char* filename);


private:

	// get vertical speed, which depends on current altitude and may depend on TAS in case of FPA mode 
	double GetVerticalSpd_(double t, double alt, double TAS);

	// update VNAV maneuver: returns true if update took place or false if current mode stays 
	bool UpdateVNavManeuver_(double t, double *XYZM);

	// return maximum time step based on VNAV profile 
	double GetMaxVNAVStep_(double t, double *XYZM);

	// move to next VNAV maneuver and return starting position and weight
	void MoveToNextVNAVManeuver_(double& start_t, double* pXYZMNew, 
								 double t0, double* pXYZM0, 
								 double t1, double* pXYZM1);


	unsigned long nManeuvers_;		// number of maneuvers
	CManeuver** ppManeuverList_;	// list of pointers to each maneuver

	// Airspeed mode
	int SPEED_MODE_;				// 0 - IAS, 1 - MACH
	double airspeed_;				// air speed parameter (IAS or MACH)

	// Heading mode
	int HEADING_MODE_;				// 0 - true heading, 1 - magnetic heading, 2 - true track, 3 - magnetic track
	double heading_;				// heading parameter (true/magnetic heading or true/magnetic track)

	// Default bank angle:
	double default_bank_angle_;		// set to 25 deg

	// Single or double engine flag
	bool isSingleEngine_;

	CMH370EnvironData* pEnvironData_; // pointer on object providing environmental data 
	CTrajectory* pTrajectory_;		// pointer to trajectory object, which provides VNAV profile independent on LNAV

	LNAVProfile* pLNAVProfile_;		// list of lateral maneuvers (of length nManeuvers_)

	// vertical navigation
	int nVManeuvers_;				// number of vertical maneuvers
	VNAVProfile* pVNAVProfile_;		// list of vertical maneuvers
	int cVmaneuver_;				// current vertical maneuver (0<=cVmaneuver_<nVManeuvers_)

	// add maneuver and associated with it parameters to the list
	void addManeuver2List_(int maneuvertype, unsigned long num_iparams, int iparams[MAX_MANEUVER_IPARAMS], 
											 unsigned long num_fparams, double fparams[MAX_MANEUVER_FPARAMS], bool ofparams[MAX_MANEUVER_FPARAMS],
											 bool isOriginal);

	// output
	unsigned long nPoints_;			// number of output trajectory points, which depend on integration time step and duration
	double* pTime_;					// array of times (s) since 2014-03-07 00:00:00 UTC
	double* pLon_;					// array or longitudes (deg E)
	double* pLat_;					// array of latitudes (deg N)
	double* pAlt_;					// array of altitudes (m)
	double* pM_;					// array of weight (kg)
	double* pu_;					// array of u- ground velocity component in local system (W->E)
	double* pv_;					// array of v- ground velocity component in local system (S->E)
	double* pw_;					// array of w velocity component in local system (vertical)
	double* pq_;					// array of fuel flow rates (kg/s)

	// output at specific times set by SetSpecificOutputTimes, e.g. when pings are available
	unsigned long nSPoints_;		// number of specific caller-requested output trajectory points
	double* pSTime_;				// array of user-specified times (s since 2014-03-07 00:00:00 UTC)
	double* pSLon_;					// array or longitudes (deg E)
	double* pSLat_;					// array of latitudes (deg N)
	double* pSAlt_;					// array of altitudes (m)
	double* pSM_;					// array of weight (kg)
	double* pSu_;					// array of u- ground velocity component in local system (W->E)
	double* pSv_;					// array of v- ground velocity component in local system (S->E)
	double* pSw_;					// array of w velocity component in local system (vertical)
	double* pSq_;					// array of fuel flow rates (kg/s)

	unsigned long nSPointsUnique_;	// number of sorted requested output trajectory points with duplicated values removed
	double* pSTimeUnique_;			// array of specified times (s since 2014-03-07 00:00:00 UTC; sorted in ascending order with duplicated and close values removed

	// Verbose diganostic stream
	FILE* fileout_; // message stream for diagnostic output (can be NULL - no output) 

};

#endif
