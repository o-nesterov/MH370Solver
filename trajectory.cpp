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
#include "trajectory.h"
#include "inmarsatdata.h"
#include "params_parser.h"


extern void splineInterp(double* pY, double* pX, unsigned long n_out, double* pInpX, double* pInpY, unsigned long n_in);
extern bool isPingIndexInList(CInpParams* pInpParams, int idx);
extern void print_time_formatted(char* buf, double t);
extern bool read_time_formatted(double& t, char* buf);


extern double True2MagneticHDG(double hdg, double lon, double lat, CMH370EnvironData* pEnvironData);
extern double Magnetic2TrueHDG(double hdg, double lon, double lat, CMH370EnvironData* pEnvironData);
extern double True2GyroHDG(double hdg, double lon, double lat, double t, double  gyro_heading_param1, double gyro_heading_param2);
extern double Gyro2TrueHDG(double hdg, double lon, double lat, double t, double gyro_heading_param1, double gyro_heading_param2);
extern double Magnetic2GyroHDG(double hdg, double lon, double lat, double t, double gyro_heading_param1, double gyro_heading_param2, CMH370EnvironData* pEnvironData);
extern double Gyro2MagneticHDG(double hdg, double lon, double lat, double t, double gyro_heading_param1, double gyro_heading_param2, CMH370EnvironData* pEnvironData);


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Build trajectory
//
void BuildTrajectory(CInpParams* pInpParams, CTrajectory* pTrajectory, CMH370EnvironData* pMH370EnvironData, CMH370InmarsatData* pInmarsatData)
{
	double t0 = pInpParams->ts;				// initial time (s since 2014-03-07 00:00:00)
	double t1 = pInpParams->te;				// end time (s since 2014-03-07 00:00:00)
	int spd_mode = pInpParams->SPD_MODE;	// speed mode
	int hdg_mode = pInpParams->HDG_MODE;	// heading mode

	double lon0 = pInpParams->lon0;			// initial longitude (deg E)
	double lat0 = pInpParams->lat0;			// initial latitude (deg N)
	double alt0 = pInpParams->alt0;			// altitude (m)
	double M0 = pInpParams->wgt0;			// weight (kg)
	double hdg0 = pInpParams->hdg0;			// deg
	double airspd = pInpParams->spd0;		// IAS m/s
	int hdg0_type = pInpParams->hdg0_type;	// initial heading value type: true, magnetic or gyroscopic

	double gyro_hdg_param1 = pInpParams->gyro_hdg_param1;	// gyroscopic heading parameters
	double gyro_hdg_param2 = pInpParams->gyro_hdg_param2;

	if (pInpParams->nLNAVOptParams>0) pTrajectory->ResetLNAVOptimizationParameters(pInpParams->LNAVOptParams, pInpParams->nLNAVOptParams);	// reset first leg duration if it was previously optimized
	pTrajectory->SetDefaultBankAngle(pInpParams->bank_angle, pInpParams->force_default_bank_angle);
	if (pInpParams->ENG_MODE == ENGINE_MODE_SINGLE) pTrajectory->SetEngineModeSingle();
	if (pInpParams->ENG_MODE == ENGINE_MODE_DUAL) pTrajectory->SetEngineModeDual();

	// Find ping in the interval [t0, t1];
	unsigned long nPings =  pInmarsatData->GetNumberOfPings();
	unsigned long nfirst = nPings;	// set an unrealistic number
	unsigned long nlast = nPings;

	for (unsigned long i=0; i<nPings; i++)
	{
		double ts;
		bool res = pInmarsatData->GetPingData(ts, i);
		if (res)
		{
			if (ts>=t0 && ts<=t1 && nfirst == nPings) {nfirst = i; break;} // break loop
		}
	};

	for (unsigned long i=0; i<nPings; i++)
	{
		double ts;
		bool res = pInmarsatData->GetPingData(ts, i);
		if (res)
		{
			if (ts>=t0 && ts<=t1) nlast = i;
		}
	};

	// ------------------------
	// set special output times
	double* pSpecialOutputTimes = new double[nlast-nfirst+2];
	for (unsigned long i=nfirst; i<=nlast; i++)
	{
		double ts;
		pInmarsatData->GetPingData(ts, i);
		pSpecialOutputTimes[i-nfirst] = ts;
	};
	pTrajectory->SetSpecificOutputTimes(pSpecialOutputTimes, nlast-nfirst+1, pInpParams->AES_Delay);

	// ------------------------
	// Build trajectory
	pTrajectory->BuildTrajectory(t0, t1, lon0, lat0, alt0, M0, spd_mode, airspd, hdg_mode, hdg0, hdg0_type, pMH370EnvironData, gyro_hdg_param1, gyro_hdg_param2);

	// ------------------------
	// Trajectory output
	if (pInpParams->pOutTrajFilename) pTrajectory->SaveTrajectory(pInpParams->pOutTrajFilename);

	// ------------------------
	//  Statistics output
	FILE* fid = NULL;
	char buf[256];

	if (pInpParams->pOutStatFilename!=NULL)
	{
		fid = fopen(pInpParams->pOutStatFilename,"wt");
	}

	if (fid)
	{
		// Output BTO comparison
		fprintf(fid,"***************************************************************************************************\n");
		fprintf(fid,"Time (UTC), Measured BTO (microsecond), Modelled BTO (microsecond), Difference BTO (microsecond)\n");

		for (unsigned long i=0; i<nlast-nfirst+1; i++)
		{
			// check if ping is in the custom list;  skip if not
			if (!isPingIndexInList(pInpParams, i+nfirst)) continue;

			double lon, lat, alt, M, u, v, w, q; // actual position and velocity data

			pTrajectory->GetSpecificOutput(lon, lat, alt, M, u, v, w, q, i, false);

			double ts, dist, bto, bfo;
			int channelType;
			pInmarsatData->GetPingData(ts, dist, bto, bfo, channelType, i+nfirst);

			double bto_computed = pInmarsatData->CalcBTO(pSpecialOutputTimes[i], lon, lat, alt, channelType);

			print_time_formatted(buf, pSpecialOutputTimes[i]);

			if (channelType!=CH_TYPE_C_RX1200) fprintf(fid, "%s %9.0f %10.1f %9.2f\n", buf, bto, bto_computed, bto_computed-bto);
		}

		// -------------------------------------------------
		// Output satellite-to-aircraft distance comparison
		fprintf(fid,"\n");
		fprintf(fid,"***************************************************************************************************\n");
		fprintf(fid,"Time (UTC), Inferred distance (km), Modelled distance (km), Difference distance (km)\n");
		for (unsigned long i=0; i<nlast-nfirst+1; i++)
		{
			// check if ping is in the custom list; skip if not
			if (!isPingIndexInList(pInpParams, i+nfirst)) continue;

			double lon, lat, alt, M, u, v, w, q;	// actual position and velocity data

			pTrajectory->GetSpecificOutput(lon, lat, alt, M, u, v, w, q, i, false);

			double ts, dist, bto, bfo;
			int channelType;
			pInmarsatData->GetPingData(ts, dist, bto, bfo, channelType, i+nfirst);

			double dist_computed = pInmarsatData->CalcDistToSat(pSpecialOutputTimes[i], lon, lat, alt);

			print_time_formatted(buf, pSpecialOutputTimes[i]);
			// fprintf(fid, "%11.3f %f %f %f %f\n", pSpecialOutputTimes[i], bto, bto_computed, bfo, bfo_computed);
			if (channelType!=CH_TYPE_C_RX1200) fprintf(fid, "%s %11.3f %11.3f %8.3f\n", buf, dist*0.001, dist_computed*0.001, (dist_computed-dist)*0.001);
		}

		// -----------------------------------------
		// Output BFO comparison
		fprintf(fid,"\n");
		fprintf(fid,"***************************************************************************************************\n");
		fprintf(fid,"Time (UTC), Measured BFO (Hz), Modelled BFO (Hz), Difference BFO (Hz)\n");
		for (unsigned long i=0; i<nlast-nfirst+1; i++)
		{
			// check if ping is in the custom list; skip if not
			if (!isPingIndexInList(pInpParams, i+nfirst)) continue;

			double lon, lat, alt, M, u, v, w, q;	// actual position and velocity data
			double aes_lon, aes_lat, aes_alt, aes_M, aes_u, aes_v, aes_w, aes_q;	// delayed position and velocity data used for AES compsation term calculation


			pTrajectory->GetSpecificOutput(lon, lat, alt, M, u, v, w, q, i, false);
			pTrajectory->GetSpecificOutput(aes_lon, aes_lat, aes_alt, aes_M, aes_u, aes_v, aes_w, aes_q, i, true);	// get delayed position and velocity data

			double ts, dist, bto, bfo;
			int channelType;
			pInmarsatData->GetPingData(ts, dist, bto, bfo, channelType, i+nfirst);

			double bfo_computed = pInmarsatData->CalcBFO(pSpecialOutputTimes[i], lon, lat, alt, u, v, w, aes_lon, aes_lat, aes_alt, aes_u, aes_v, channelType);

			print_time_formatted(buf, pSpecialOutputTimes[i]);

			fprintf(fid, "%s %6.0f %8.1f %8.2f\n", buf, bfo, bfo_computed, bfo_computed-bfo);
		}

		fclose(fid);
	}

	delete pSpecialOutputTimes;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Wrapper interface function called from main

void BuildTrajectory(char* inpfilename, bool isVerbose)
{
	if (isVerbose) printf("Integrating...\n");
	// Read-in setup file
	CInpParams* pInpParams = new CInpParams();	// input parameters
	if (isVerbose) pInpParams -> SetVerbose(stdout);
	if (!pInpParams -> parseParams(inpfilename))
	{
		delete pInpParams;
		return;
	}

	// Inmarsat data
	CMH370InmarsatData* pInmarsatData = new CMH370InmarsatData();
	if (isVerbose) pInmarsatData -> SetVerbose(stdout);
	// set parameters, which may affect Doppler calculation and AES Doppler compensation term
	if (pInpParams->DopplerMethod != DOPPLER_METHOD_UDF) pInmarsatData -> SetDopplerMethod(pInpParams->DopplerMethod);
	if (pInpParams->AES_EarthModel == AES_EARTH_MODEL_SPH || pInpParams->AES_EarthModel == AES_EARTH_MODEL_WGS)  pInmarsatData -> SetAESEarthModel(pInpParams->AES_EarthModel);
	if (pInpParams->AES_DopplerMethod != DOPPLER_METHOD_UDF)  pInmarsatData -> SetAESDopplerMethod(pInpParams->AES_DopplerMethod);
	if (pInpParams->AES_Earth_Radius > 0.0) pInmarsatData -> SetAESEarthRadius(pInpParams->AES_Earth_Radius);
	if (pInpParams->AES_GeoSat_Alt > 0.0) pInmarsatData -> SetAESGeoSatAlt(pInpParams->AES_GeoSat_Alt);


	// Trajectory
	CTrajectory* pTrajectory = new CTrajectory();
	if (isVerbose) pTrajectory->SetVerbose(stdout);

	if(!pTrajectory->CreateLNAVProfile(pInpParams))
	{
		if (isVerbose) printf("Integration failed. Missing or wrong LNAV file %s.\n", pInpParams->pLNAVfilename);
		delete pTrajectory;
		delete pInmarsatData;
		delete pInpParams;
		return;
	}

	if (pInpParams->bank_angle >0.0) pTrajectory->SetDefaultBankAngle(pInpParams->bank_angle, pInpParams->force_default_bank_angle);	// set default bank angle to connect path legs if turns are not explicitly specified	
	// if (pInpParams->first_leg_duration>0.0) pTrajectory->SetFirstLegDuration(pInpParams->first_leg_duration); // first leg duration should be specified in LNAV file

	// set engine mode 
	if (pInpParams->ENG_MODE == ENGINE_MODE_SINGLE) pTrajectory->SetEngineModeSingle();
	if (pInpParams->ENG_MODE == ENGINE_MODE_DUAL) pTrajectory->SetEngineModeDual();


	// Load environmental data:
	CMH370EnvironData* pMH370EnvironData = new CMH370EnvironData();
	if (isVerbose) pMH370EnvironData -> SetVerbose(stdout,1);
	pMH370EnvironData -> LoadMeteo(pInpParams->METEO_DATASET, pInpParams->pMeteoDir);
	pMH370EnvironData -> SetMeteoInterpMethod(pInpParams->MeteoInterpMethod);
	pMH370EnvironData -> LoadMagnetic(pInpParams->MAGN_DECL_DATASET, pInpParams->pMagnDeclDir);
	pMH370EnvironData -> SetDeclInterpMethod(pInpParams->MagnDeclMethod);

	BuildTrajectory(pInpParams, pTrajectory, pMH370EnvironData, pInmarsatData);

	delete pTrajectory;
	delete pMH370EnvironData;
	delete pInmarsatData;
	delete pInpParams;
	if (isVerbose) printf("Integration completed.\n");
}





////////////////////////////////////////////////////////////////////////////////////////
// Trajectory class
// Comprises a sequence of maneuvers
////////////////////////////////////////////////////////////////////////////////////////
CTrajectory::CTrajectory()
{
	nManeuvers_ = 0;		// number of maneuvers
	ppManeuverList_ = NULL;	// array of pointers on maneuvers

	pLNAVProfile_ = NULL;	// list of lateral maneuvers

	// vertical maneuvers
	nVManeuvers_ = 0;			// number of vertical maneuvers
	pVNAVProfile_ = NULL;		// list of vertical maneuvers
	cVmaneuver_ = 0;			// current vertical maneuver

	// output
	nPoints_ = 0;			// number of output trajectory points, which depend on integration time step and duration
	pTime_ = NULL;			// array of times (s) since 2014-03-07 00:00:00 UTC
	pLon_ = NULL;			// array or longitudes (deg E)
	pLat_ = NULL;			// array of latitudes (deg N)
	pAlt_ = NULL;			// array of altitudes (m)
	pM_ = NULL;				// array of weight (kg)
	pu_ = NULL;				// array of u- ground velocity component in local system (W->E)
	pv_ = NULL;				// array of v- ground velocity component in local system (S->E)
	pw_ = NULL;				// array of w velocity component in local system (vertical)
	pq_ = NULL;				// array of fuel flow rates (kg/s)

	// output at specific times
	nSPoints_ = 0;			// number of specific output trajectory points
	nSPointsUnique_ = 0;	// number of 'unique' specific output trajectory points (duplicated or close values are treated as single)
	pSTime_ = NULL;			// array of times (s) since 2014-03-07 00:00:00 UTC
	pSTimeUnique_ = NULL;	// array of sorted times with duplicated and close values removed
	pSLon_ = NULL;			// array or longitudes (deg E)
	pSLat_ = NULL;			// array of latitudes (deg N)
	pSAlt_ = NULL;			// array of altitudes (m)
	pSM_ = NULL;			// array of weight (kg)
	pSu_ = NULL;			// array of u- ground velocity component in local system (W->E)
	pSv_ = NULL;			// array of v- ground velocity component in local system (S->E)
	pSw_ = NULL;			// array of w velocity component in local system (vertical)
	pSq_ = NULL;			// array of fuel flow rates (kg/s)

	default_bank_angle_ = 25.0; // default bank angle

	isSingleEngine_ = false;	// default dual engine 

	pEnvironData_ = NULL;	// pointer to environmental data provider object

	// Verbose diganostic stream
	fileout_ = NULL; // message stream for diagnostic output (can be NULL) 
}


////////////////////////////////////////////////////////////////////////////////////////
// Destructor: deallocate all arrays
////////////////////////////////////////////////////////////////////////////////////////
CTrajectory::~CTrajectory()
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

	if (pSTime_) delete pSTime_;
	if (pSTimeUnique_) delete pSTimeUnique_;
	if (pSLon_) delete pSLon_;
	if (pSLat_) delete pSLat_;
	if (pSAlt_) delete pSAlt_;
	if (pSM_) delete pSM_;
	if (pSu_) delete pSu_;
	if (pSv_) delete pSv_;
	if (pSw_) delete pSw_;
	if (pSq_) delete pSq_; 

	if (ppManeuverList_)
	{
		for (unsigned long i=0; i<nManeuvers_; i++)	if (ppManeuverList_[i]) delete ppManeuverList_[i];
		delete ppManeuverList_;
	}

	if (pLNAVProfile_) delete pLNAVProfile_;
	if (pVNAVProfile_) delete pVNAVProfile_;
}


////////////////////////////////////////////////////////////////////////////////////////
//  Output interface functions
////////////////////////////////////////////////////////////////////////////////////////

// return number of output points (nPoints_)
unsigned long CTrajectory::GetNOut()
{
	return nPoints_;
}

// return pointer to array of times (of length nPoints_)
double* CTrajectory::GetTimeOut()
{
	return pTime_;
}

// return pointer to array of longitudes (of length nPoints_)
double* CTrajectory::GetLonOut()
{
	return pLon_;
}

// return pointer to array of latitudes (of length nPoints_)
double* CTrajectory::GetLatOut()
{
	return pLat_;
}

// return pointer to array of altitudes (of length nPoints_)
double* CTrajectory::GetAltOut()
{
	return pAlt_;
}

// return pointer to array of weights (of length nPoints_)
double* CTrajectory::GetMOut()
{
	return pM_;
}


// return pointer to array of u-velocity (W->E, m/s, of length nPoints_)
double* CTrajectory::GetuOut()
{
	return pu_;
}

// return pointer to array of v-velocity (S->N, m/s, of length nPoints_)
double* CTrajectory::GetvOut()
{
	return pv_;
}


// return pointer to array of w-velocity (vertical, m/s, of length nPoints_)
double* CTrajectory::GetwOut()
{
	return pw_;
}

// return pointer to array of fuel flow rate (kg/s, of length nPoints_)
double* CTrajectory::GetqOut()
{
	return pq_;
}

// set default bank angle to be used when legs are not explicitly connected by turns
void CTrajectory::SetDefaultBankAngle(double bank_angle, bool force_default)
{
	// bank_angle - new default bank angle
	// force_default -  override bank angle in all explicitly listed turns, which are not subjected to individual optimization (* in lnav file) 
	//                  it is used for optimization of all bank angles simultaneously

	if (bank_angle <= 0.0) return; // do nothing

	default_bank_angle_ = bank_angle;

	if (nManeuvers_ == 0) return;

	for (int i=0; i<nManeuvers_; i++)
	{
		if (pLNAVProfile_[i].maneuverType == MANEUVER_TYPE_TURN)
		{
			// override explicitly specified bank angle (pLNAVProfile_[i].fparams[1]) with a new value if the follwoing conditions are met:
			// 1). force_default flag is set
			// 2). positive bank angle was explicitly specified in lnav file (negative values imply the bank angle is not specified and it should be treated as default)
			// 3). the bank angle is not a subject to individual optimization (marked with * in lnav file) 
			if (force_default && (pLNAVProfile_[i].fparams[1]>0.0) && (!pLNAVProfile_[i].ofparams[1])) pLNAVProfile_[i].fparams[1] = default_bank_angle_;
		}
	}
}

void CTrajectory::SetEngineModeSingle()
{
	isSingleEngine_ = true;
}

void CTrajectory::SetEngineModeDual()
{
	isSingleEngine_ = false;
}

bool CTrajectory::IsEngineModeSingle()
{
	return isSingleEngine_;
}

////////////////////////////////////////////////////////////////////////////////////////
//  Verbose for diagnostic
////////////////////////////////////////////////////////////////////////////////////////

void CTrajectory::SetVerbose(FILE* outputstream)
{
	// set verbose stream (can be NULL - no output)
	fileout_ = outputstream;
}



////////////////////////////////////////////////////////////////////////////////////////
//  Set specific output times,  when output is required;
//  The list is duplicated to acount for delay in AES position/velocity data update
//  (ARINC transmits data at either 20, 40 or 50 millisecond intervals, and AES may wait till next update potentially doubling this;
//   GPS data could probably be updated at lower frequency (1s?), which may notably affect AES compensation term during turns)
////////////////////////////////////////////////////////////////////////////////////////
void CTrajectory::SetSpecificOutputTimes(double* pTimes, unsigned long n, double AES_Delay)
{
	// Input pTimes - array of times (s since 2014-03-07 00:00:00), when output is required

	if (pSTime_) delete pSTime_;
	if (pSTimeUnique_) delete pSTimeUnique_;
	if (pSLon_) delete pSLon_;
	if (pSLat_) delete pSLat_;
	if (pSAlt_) delete pSAlt_;
	if (pSM_) delete pSM_;
	if (pSu_) delete pSu_;
	if (pSv_) delete pSv_;
	if (pSw_) delete pSw_;
	if (pSq_) delete pSq_; 

	nSPoints_ = 0;			// number of specific output trajectory points
	nSPointsUnique_ = 0;	// number of 'unique' (undistinguishable) output trajectory points
	pSTime_ = NULL;			// array of times (s) since 2014-03-07 00:00:00 UTC as provided
	pSTimeUnique_ = NULL;	// array of times (s) since 2014-03-07 00:00:00 UTC sorted, with duplicted values eliminated
	pSLon_ = NULL;			// array or longitudes (deg E)
	pSLat_ = NULL;			// array of latitudes (deg N)
	pSAlt_ = NULL;			// array of altitudes (m)
	pSM_ = NULL;			// array of weight (kg)
	pSu_ = NULL;			// array of u- ground velocity component in local system (W->E)
	pSv_ = NULL;			// array of v- ground velocity component in local system (S->E)
	pSw_ = NULL;			// array of w velocity component in local system (vertical)
	pSq_ = NULL;			// array of fuel flow rates (kg/s)

	if (pTimes == NULL || n==0) return;

	nSPoints_ = 2*n;
	nSPointsUnique_ = 2*n;	// assume all input times unique here 

	// allocate arrays
	pSTime_ = new double[nSPoints_+1];
	pSTimeUnique_ = new double[nSPoints_+1];
	pSLon_ = new double[nSPoints_+1];
	pSLat_ = new double[nSPoints_+1];
	pSAlt_ = new double[nSPoints_+1];
	pSM_ = new double[nSPoints_+1];
	pSu_ = new double[nSPoints_+1];
	pSv_ = new double[nSPoints_+1];
	pSw_ = new double[nSPoints_+1];
	pSq_ = new double[nSPoints_+1];

	// here indexes 2*i match delayed times and 2*i+1 match actual times

	for (int i=0; i<(int)n; i++)
	{
		pSTime_[2*i] = pTimes[i] - AES_Delay; // as provided minus delay
		pSTime_[2*i+1] = pTimes[i]; // as provided

		pSLon_[2*i] = -1.0E+10;	// undefined value
		pSLon_[2*i+1] = -1.0E+10;

		pSLat_[2*i] = -1.0E+10;
		pSLat_[2*i+1] = -1.0E+10;

		pSAlt_[2*i] = -1.0E+10;
		pSAlt_[2*i+1] = -1.0E+10;

		pSM_[2*i] = -1.0E+10;
		pSM_[2*i+1] = -1.0E+10;

		pSu_[2*i] = -1.0E+10;
		pSu_[2*i+1] = -1.0E+10;

		pSv_[2*i] = -1.0E+10;
		pSv_[2*i+1] = -1.0E+10;

		pSw_[2*i] = -1.0E+10;
		pSw_[2*i+1] = -1.0E+10;

		pSq_[2*i] = -1.0E+10;
		pSq_[2*i+1] = -1.0E+10;
	}

	// -------------------------------------------------------------------------------------------
	// The list of required output times may not be in order, or it may contain duplicated records, or very close times 
	// so that sort them in ascending order, while simultaneously removing duplicated or close values


	// set the first value, which is minimum
	double t_min = pTimes[0];
	for (int i=1; i<(int)n; i++) if (t_min > pTimes[i]) t_min = pTimes[i];

	// find maximum
	double t_max = t_min;
	for (int i=0; i<(int)n; i++) if (t_max < pTimes[i]) t_max = pTimes[i];

	// if minimum and maximum are close to each other, then take average point and return
	if (t_max - t_min < 2.0*EPS_TIME)
	{
		if (AES_Delay>2.0*EPS_TIME)
		{
			pSTimeUnique_[0] = 0.5*(t_max + t_min) - AES_Delay;
			pSTimeUnique_[1] = 0.5*(t_max + t_min);
			nSPointsUnique_ = 2;
		}
		else
		{
			pSTimeUnique_[0] = 0.5*(t_max + t_min);
			nSPointsUnique_ = 1;
		}
		return;
	}

	pSTimeUnique_[0] = t_min - AES_Delay;


	for (int i=1; i<(int)nSPoints_; i++)
	{
		double next_min = pSTimeUnique_[i-1] + 2.0*EPS_TIME;	// next value should be greater than previous by at least EPS_TIME
		double t_next = t_max + 2.0*EPS_TIME;					// assume a value greater than maximum

		// ------
		// scan and pick the next value from the list, which should be not less than next_min

		// actual time list
		for (int j=0; j<(int)n; j++)
		{
			if ((pTimes[j]>=next_min) && (pTimes[j]<=t_next))
			{
				t_next = pTimes[j];
			}
		}

		// delayed time list
		for (int j=0; j<(int)n; j++)
		{
			if (((pTimes[j]-AES_Delay)>=next_min) && ((pTimes[j]-AES_Delay)<=t_next))
			{
				t_next = pTimes[j]-AES_Delay;
			}
		}

		// ----------
		// store found value or break the loop
		if (t_next > t_max+EPS_TIME)
		{
			nSPointsUnique_ = i;
			break;			// break i-loop because no new value can be found
		}
		else
		{
			pSTimeUnique_[i] = t_next;	// store value and search next one
		}
	}

}





////////////////////////////////////////////////////////////////////////////////////////
//
// Get output at specific time, which must be previously specified using SetSpecificOutputTimes
// if isDelayed = false, then return actual position and velocity data;
// if isDelayed = true,  then return data matching delayed position and velocity for AES compensation term calculation 
////////////////////////////////////////////////////////////////////////////////////////

void CTrajectory::GetSpecificOutput(double& lon, double& lat, double& alt, double& M, double& u, double& v, double& w, double& q, unsigned long n, bool isDelayed)
{
	if (nSPoints_==0 || 2*n>=nSPoints_)
	{
		lon = -1.0E+10;
		lat = -1.0E+10;
		alt = -1.0E+10;
		M = -1.0E+10;
		u = -1.0E+10;
		v = -1.0E+10;
		w = -1.0E+10;
		q = -1.0E+10;
		return;
	}

	if (isDelayed)
	{
		lon = pSLon_[2*n];
		lat = pSLat_[2*n];
		alt = pSAlt_[2*n];
		M = pSM_[2*n];
		u = pSu_[2*n];
		v = pSv_[2*n];
		w = pSw_[2*n];
		q = pSq_[2*n];
	}
	else
	{
		lon = pSLon_[2*n+1];
		lat = pSLat_[2*n+1];
		alt = pSAlt_[2*n+1];
		M = pSM_[2*n+1];
		u = pSu_[2*n+1];
		v = pSv_[2*n+1];
		w = pSw_[2*n+1];
		q = pSq_[2*n+1];
	}
}




////////////////////////////////////////////////////////////////////////////////////////
//
// Integrate entire trajectory
//
////////////////////////////////////////////////////////////////////////////////////////

void CTrajectory::BuildTrajectory(double t0, double endtime,
								  double lon0, double lat0, double alt0, double M0,
								  int speed_mode, double airspeed, int heading_mode, double hdg0, int hdg0_type,
								  CMH370EnvironData* pEnvironData,
								  double gyro_hdg_param1, double gyro_hdg_param2)
{
	// Input:
	// t0 - start time (s)
	// t1 - end time (s)
	// RefTime - reference time == offset of the start (s) since 2014-03-07 00:00:00 
	// lon0 - initial longitude (deg E)
	// lat0 - initial latitude (deg N)
	// alt0 - initial altitude (m)
	// M0 - initial weight (kg)
	// hdg0 - inidial heading (deg)
	// hdg0_type - initial heading value type: true, magnetic, or gyroscopic 
	// spd_mode - speed mode (SPEED_MODE_IAS or SPEED_MODE_MACH)
	// hdg_mode - heading mode (HEADING_MODE_TRUEHDG, HEADING_MODE_MAGNHDG, HEADING_MODE_TRUETRK, or HEADING_MODE_MAGNTRK)
	// bank_angle - maximum bank angle (deg) in turns


	// -------------------------------
	// clear maneuvers array if needed:
	if (ppManeuverList_)
	{
		for (unsigned long i=0; i<nManeuvers_; i++) if (ppManeuverList_[i]) delete ppManeuverList_[i];
		delete ppManeuverList_;
		ppManeuverList_ = NULL;
	}

	if (nManeuvers_==0)	// sanity check
	{
		if (fileout_) fprintf(fileout_,"No valid maneuvers to integrate. Check lateral navigation file.\n");
		return;
	}

	// -------------------------------
	// create new array of maneuvers
	ppManeuverList_ = new CManeuver*[nManeuvers_+1];

	// beginning of the section
	int current_maneuver = 0;		 // maneuver number - start from the straight path

	double current_time = t0;		// current time
	double current_lon = lon0;		// current longitude
	double current_lat = lat0;		// current latitude
	double current_alt = alt0;		// current altitude
	double current_m = M0;			// current mass

	double current_hdg = hdg0;		// current heading in the specified reference system: true, magnetic or gyroscopic
	// convert specified heading to heading in current reference system 
	switch (hdg0_type)
	{
		case HEADING_MODE_TRUEHDG:
		if (heading_mode == HEADING_MODE_MAGNHDG || heading_mode == HEADING_MODE_MAGNTRK) current_hdg = True2MagneticHDG(hdg0, current_lon, current_lat, pEnvironData_);
		if (heading_mode == HEADING_MODE_GYROHDG) current_hdg = True2GyroHDG(hdg0, current_lon, current_lat, current_time, gyro_hdg_param1, gyro_hdg_param2);
		break;

		case HEADING_MODE_MAGNHDG:
		if (heading_mode == HEADING_MODE_TRUEHDG || heading_mode == HEADING_MODE_TRUETRK) current_hdg = Magnetic2TrueHDG(hdg0, current_lon, current_lat, pEnvironData_);
		if (heading_mode == HEADING_MODE_GYROHDG) current_hdg = Magnetic2GyroHDG(hdg0, current_lon, current_lat, current_time, gyro_hdg_param1, gyro_hdg_param2, pEnvironData_);
		break;

		case HEADING_MODE_GYROHDG:
		if (heading_mode == HEADING_MODE_TRUEHDG || heading_mode == HEADING_MODE_TRUETRK) current_hdg = Gyro2TrueHDG(hdg0, current_lon, current_lat, current_time, gyro_hdg_param1, gyro_hdg_param2);
		if (heading_mode == HEADING_MODE_MAGNHDG || heading_mode == HEADING_MODE_MAGNTRK) current_hdg = Gyro2MagneticHDG(hdg0, current_lon, current_lat, current_time, gyro_hdg_param1, gyro_hdg_param2, pEnvironData_);
		break;
	}


	// -------------------------------------------------------------------------
	// maneuver iterations, including turns and constant heading/track sections

	for (unsigned long current_maneuver=0; current_maneuver<nManeuvers_; current_maneuver++)
	{
		// Currently implemented types of maneuvers:
		// CManeuverLeg - 'straight' leg (constant true/magnetic heading/track) of specified duration
		// CManeuverTurn - 'turn' with specified bank angle

		// create maneuver
		CManeuver* pManeuver;
		double max_duration = endtime - current_time;
		double timestep, duration, turn_by_angle, turn_to_hdg, bank_angle;

		// failsafe if the simulation time has reached the specified end, but there are still unfinished maneuvers in queue 
		if (max_duration<0.001)
		{
			ppManeuverList_[current_maneuver] = NULL;
			continue;
		}

		switch (pLNAVProfile_[current_maneuver].maneuverType)
		{

		// simple leg with two parameters  - duration and max integration time step and 
		// leg with specified heading - heading, duration and max integration time step; headings were adjusted in lnav.cpp with the help of turn inserted if needed 

		case MANEUVER_TYPE_LEG:
		case MANEUVER_TYPE_HLEG:

			if (pLNAVProfile_[current_maneuver].maneuverType==MANEUVER_TYPE_LEG)
			{
				duration = pLNAVProfile_[current_maneuver].fparams[0];	// duration
				timestep = pLNAVProfile_[current_maneuver].fparams[1];	// maximum integration time step (s) 
			}
			else
			{
				duration = pLNAVProfile_[current_maneuver].fparams[1];	// duration
				timestep = pLNAVProfile_[current_maneuver].fparams[2];	// maximum integration time step (s) 
			}

			if (duration>max_duration) duration = max_duration;

			if (heading_mode != HEADING_MODE_GYROHDG)
			{
				pManeuver = new CManeuverLeg(current_time, duration, timestep, pSTimeUnique_, nSPointsUnique_,
											 current_lon, current_lat, current_alt, current_m, 
											 speed_mode, airspeed, heading_mode, current_hdg,
											 pEnvironData, this);
			}
			else
			{
				pManeuver = new CManeuverLeg(current_time, duration, timestep, pSTimeUnique_, nSPointsUnique_,
											 current_lon, current_lat, current_alt, current_m, 
											 speed_mode, airspeed, current_hdg, gyro_hdg_param1, gyro_hdg_param2,
											 pEnvironData, this);
			}
			break;



		case MANEUVER_TYPE_TURN:

			turn_by_angle = -1.0E+20;
			turn_to_hdg = -1.0E+20;
			if (pLNAVProfile_[current_maneuver].iparams[0] == TURN_TYPE_BY_ANGLE)
			{
				turn_by_angle = pLNAVProfile_[current_maneuver].fparams[0];			// turn angle (deg)
			}
			else
			{
				turn_to_hdg = pLNAVProfile_[current_maneuver].fparams[0];

				// convert heading to current heading system 
				switch (pLNAVProfile_[current_maneuver].iparams[1])
				{
				case HEADING_MODE_TRUEHDG:
					if (heading_mode == HEADING_MODE_MAGNHDG || heading_mode == HEADING_MODE_MAGNTRK) turn_to_hdg = True2MagneticHDG(turn_to_hdg, current_lon, current_lat, pEnvironData_);
					if (heading_mode == HEADING_MODE_GYROHDG) turn_to_hdg = True2GyroHDG(turn_to_hdg, current_lon, current_lat, current_time, gyro_hdg_param1, gyro_hdg_param2);
					break;

				case HEADING_MODE_MAGNHDG:
					if (heading_mode == HEADING_MODE_TRUEHDG || heading_mode == HEADING_MODE_TRUETRK) turn_to_hdg = Magnetic2TrueHDG(turn_to_hdg, current_lon, current_lat, pEnvironData_);
					if (heading_mode == HEADING_MODE_GYROHDG) turn_to_hdg = Magnetic2GyroHDG(turn_to_hdg, current_lon, current_lat, current_time, gyro_hdg_param1, gyro_hdg_param2, pEnvironData_);
					break;

				case HEADING_MODE_GYROHDG:
					if (heading_mode == HEADING_MODE_TRUEHDG || heading_mode == HEADING_MODE_TRUETRK) turn_to_hdg = Gyro2TrueHDG(turn_to_hdg, current_lon, current_lat, current_time, gyro_hdg_param1, gyro_hdg_param2);
					if (heading_mode == HEADING_MODE_MAGNHDG || heading_mode == HEADING_MODE_MAGNTRK) turn_to_hdg = Gyro2MagneticHDG(turn_to_hdg, current_lon, current_lat, current_time, gyro_hdg_param1, gyro_hdg_param2, pEnvironData_);
					break;
				}

			}

			bank_angle = pLNAVProfile_[current_maneuver].fparams[1];	// bank angle (deg)
			if (bank_angle < 0.0) bank_angle = default_bank_angle_;		// bank angle was not specified, so assume default bank angle, which may vary during optimization

			timestep = pLNAVProfile_[current_maneuver].fparams[2];		// maximum integration time step (s) 

			if (heading_mode != HEADING_MODE_GYROHDG)
			{
				pManeuver = new CManeuverTurn(current_time, max_duration, timestep, pSTimeUnique_, nSPointsUnique_,
											  current_lon, current_lat, current_alt, current_m, 
											  speed_mode, airspeed, heading_mode, current_hdg,
											  turn_by_angle, turn_to_hdg, bank_angle, 
											  pEnvironData, this);
			}
			else
			{
				pManeuver = new CManeuverTurn(current_time, max_duration, timestep, pSTimeUnique_, nSPointsUnique_,
											  current_lon, current_lat, current_alt, current_m, 
											  speed_mode, airspeed, current_hdg, gyro_hdg_param1, gyro_hdg_param2,
											  turn_by_angle, turn_to_hdg, bank_angle, 
											  pEnvironData, this);
			}

		};

		// Integrate:
		pManeuver->Integrate();


		// Add maneuver to the list of maneuvers
		ppManeuverList_[current_maneuver] = pManeuver;

		current_time = pManeuver -> GetLastTime();
		current_lon  = pManeuver -> GetLastLon();
		current_lat  = pManeuver -> GetLastLat();
		current_alt  = pManeuver -> GetLastAlt();
		current_m    = pManeuver -> GetLastM();
		current_hdg  = pManeuver -> GetLastHdg();

	};

	// -------------------------------------------------------------------------
	//  Assemble the whole output trajectory (time, lon, lat, alt, u, v, w, M)

	unsigned long n_out = 0; // total number of output points
	unsigned long n_actual_maneuvers = 0;

	// recursively add numbers of outputs in each maneuver
	for (unsigned long current_maneuver=0; current_maneuver<nManeuvers_; current_maneuver++)
	{
		if (ppManeuverList_[current_maneuver])
		{
			n_out += ppManeuverList_[current_maneuver] -> GetNOut();
			n_actual_maneuvers++;
		}
		else
		{
			break; // maneuver is NULL because time reached the specified end time
		}
	}
	n_out -= (n_actual_maneuvers-1);	// beginning of each next maneuver is the same as the end of previous one

	// deallocate existing arrays, if any
	if (pTime_) delete pTime_;
	if (pLon_) delete pLon_;
	if (pLat_) delete pLat_;
	if (pAlt_) delete pAlt_;
	if (pM_) delete pM_;
	if (pu_) delete pu_;
	if (pv_) delete pv_;
	if (pw_) delete pw_;
	if (pq_) delete pq_; 

	// allocate new arrays
	pTime_ = new double[n_out+1];
	pLon_ = new double[n_out+1];
	pLat_ = new double[n_out+1];
	pAlt_ = new double[n_out+1];
	pM_ = new double[n_out+1];
	pu_ = new double[n_out+1];
	pv_ = new double[n_out+1];
	pw_ = new double[n_out+1];
	pq_ = new double[n_out+1];

	// recursively add pieces of trajectory
	unsigned long idx = 0;
	for (unsigned long current_maneuver=0; current_maneuver<n_actual_maneuvers; current_maneuver++)
	{
		unsigned long n = ppManeuverList_[current_maneuver] -> GetNOut(); // number of output in a maneuver

		double* pTime = ppManeuverList_[current_maneuver] -> GetTimeOut();
		double* pLon =  ppManeuverList_[current_maneuver] -> GetLonOut();
		double* pLat =  ppManeuverList_[current_maneuver] -> GetLatOut();
		double* pAlt =  ppManeuverList_[current_maneuver] -> GetAltOut();
		double* pM =    ppManeuverList_[current_maneuver] -> GetMOut();
		double* pu =    ppManeuverList_[current_maneuver] -> GetuOut();
		double* pv =    ppManeuverList_[current_maneuver] -> GetvOut();
		double* pw =    ppManeuverList_[current_maneuver] -> GetwOut();
		double* pq =    ppManeuverList_[current_maneuver] -> GetqOut();

		// copy arrays
		for (int i=0; i<(int)n; i++)
		{
			pTime_[idx] = pTime[i];
			pLon_[idx] =  pLon[i];
			pLat_[idx] =  pLat[i];
			pAlt_[idx] =  pAlt[i];
			pM_[idx] =  pM[i];
			pu_[idx] =  pu[i];
			pv_[idx] =  pv[i];
			pw_[idx] =  pw[i];
			pq_[idx] =  pq[i];
			idx++;
		}
		idx--;	// beginning of the next maneuver coincides with the end of the previous one
	}

	nPoints_ = n_out;

	// -----------------------------------------------------------------------------------------
	// Extract trajectory at specified times
	// Approximate times were already 'embedded' to address case when integration step is inadequately large 
	// (for example, if max integration step during turns is 10s, interpolation may not be sufficiently accurate)
	// EPS_TIME in maneuver.h controls the resolution of 'embedded' output times
	// However, input list may differ in that it may contain times in arbitrary order, or it may contain duplicated times
	// So, the purpose of the following is to "fill in" trajectory variables exactly at requested times exactly in requested order

	splineInterp(pSLon_, pSTime_, nSPoints_, pTime_, pLon_, nPoints_);
	splineInterp(pSLat_, pSTime_, nSPoints_, pTime_, pLat_, nPoints_);
	splineInterp(pSAlt_, pSTime_, nSPoints_, pTime_, pAlt_, nPoints_);
	splineInterp(pSM_, pSTime_, nSPoints_, pTime_, pM_, nPoints_);
	splineInterp(pSu_, pSTime_, nSPoints_, pTime_, pu_, nPoints_);
	splineInterp(pSv_, pSTime_, nSPoints_, pTime_, pv_, nPoints_);
	splineInterp(pSw_, pSTime_, nSPoints_, pTime_, pw_, nPoints_);
	splineInterp(pSq_, pSTime_, nSPoints_, pTime_, pq_, nPoints_);

	// set value to undefined if requested time is outside of the modelled interval
	for (int i=0; i<(int)nSPoints_; i++)
	{
		if (pSTime_[i]<t0 || pSTime_[i]>endtime)
		{
			pSLon_[i] = -1.0E+10;
			pSLat_[i] = -1.0E+10;
			pSAlt_[i] = -1.0E+10;
			pSM_[i] = -1.0E+10;
			pSu_[i] = -1.0E+10;
			pSv_[i] = -1.0E+10;
			pSw_[i] = -1.0E+10;
			pSq_[i] = -1.0E+10;
		}
	}

	// set value to be the same as actual one if AES value (delayed) is outside of the modelled interval
	for (int i=0; i<(int)nSPoints_/2; i++)
	{
		if (pSTime_[2*i]<t0 || pSTime_[2*i]>endtime)
		{
			pSLon_[2*i] = pSLon_[2*i+1];
			pSLat_[2*i] = pSLat_[2*i+1];
			pSAlt_[2*i] = pSAlt_[2*i+1];
			pSM_[2*i] = pSM_[2*i+1];
			pSu_[2*i] = pSu_[2*i+1];
			pSv_[2*i] = pSv_[2*i+1];
			pSw_[2*i] = pSw_[2*i+1];
			pSq_[2*i] = pSq_[2*i+1];
		}
	}


}

