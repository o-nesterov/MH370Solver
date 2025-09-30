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

extern void XYZ2LonLatH(double& lon, double& lat, double& hgt, double X, double Y, double Z);

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Create/read vertical navigation profile
//
///////////////////////////////////////////////////////////////////////////////////////////////////

bool CTrajectory::CreateVNAVProfile(CInpParams* pInpParams)
{
	// vertical velocity component
	// V/S: -8000 to +6000 fpm, 100 fpm increments - 777 Flight Manual, Continental Rev. 11/01/02 #9
	// FPA: -9.9 to +9.9 degrees, 0.1 deg increments - 777 Flight Manual, Continental Rev. 11/01/02 #9
	// FLCH, aimed to complete flight level change within 125 seconds, CAL_Training_Manual_B777-200_WB371.pdf - 
	//       Boeing 777 Training manual, Coninental Airlines Inc., Airplane General, WB371

	char* filename = pInpParams->pVNAVfilename;

	// If no filename specified, then assume constant flight level
	if (filename == NULL)
	{
		if (pVNAVProfile_) delete pVNAVProfile_;
		pVNAVProfile_ = NULL;
		nVManeuvers_ = 0;

		return true;	// input file was not specified, so constant altitude specified in pInpParams is sufficient to assume constant level
	}

	// Same check as above to cover the case of filename==""
	if (strlen(filename) == 0)
	{
		if (pVNAVProfile_) delete pVNAVProfile_;
		pVNAVProfile_ = NULL;
		nVManeuvers_ = 0;

		return true;	// input file was not specified, so constant altitude specified in pInpParams is sufficient to assume constant level
	}

	// ----------------------
	// read sequence of maneuvers from the specified file
	FILE* fid = fopen(filename,"rt");

	if (fid == NULL)
	{
		if (fileout_) fprintf(fileout_,"Error: cannot open LNAV file: %s\n", filename);
		return false; // file was specified, but it cannot be opened, so return false
	}

	// reset profile
	if (pVNAVProfile_) delete pVNAVProfile_;
	pVNAVProfile_ = NULL;
	nVManeuvers_ = 0;

	// parse vertical navigation file
	char buf[1024];
	int vmaneuver;
	double param, target_alt;

	while(!feof(fid))
	{
		if (fgets(buf, 1000, fid)!=NULL)
		{
			char* p=&buf[0];
			int n = (int)strlen(p);
			while((*p==' ') || (*p=='	') && p<(&buf[n])) p++;	// skip spaces
			if (p==&buf[n] || n<=3) continue;	// skip empty line
			if (*p=='#' || *p=='!' || *p=='%%') continue; // this is a comment line - skip it

			// initialize
			vmaneuver = -1;
			param = 0.0;
			target_alt = 0.0;

			// ----------------------
			// Keep level maneuver
			char* p1 = strstr(p,"FL");
			char* p2 = strstr(p,"LEVEL");
			char* p3 = strstr(p,"FLCH");
			if ((p1!=NULL || p2!=NULL) && p3==NULL)
			{
				if (p1)
				{
					p = p1+strlen("FL");
				}
				else
				{
					p = p2+strlen("LEVEL");
				}
				while((*p==' ') || (*p=='	') && p<(&buf[n])) p++;	// skip spaces etc.
				char* pp;
				param = strtod(p,&pp);	// duration
				if (p==pp) 
				{
					if (fileout_) fprintf(fileout_, "Duration of flight level not specified. Skeeping VNAV line.\n");
					continue;
				}
				vmaneuver =  VNAV_MODE_LEVEL;
			}


			// ----------------------
			// VS mode
			p1 = strstr(p,"VS");
			if (p1 && (vmaneuver == -1))
			{
				p+=strlen("VS");
				while((*p==' ') || (*p=='	') && p<(&buf[n])) p++;	// skip spaces etc.
				char* pp;
				param = strtod(p,&pp);	// vertical speed (m/s)
				if (p==pp) 
				{
					if (fileout_) fprintf(fileout_, "Erroneus V/S. Skipping VNAV line.\n");
					continue;
				}
				p = pp;
				while((*p==' ') || (*p=='	') && p<(&buf[n])) p++;	// skip spaces etc.
				target_alt = strtod(p,&pp);	// target altitude (m)
				if (p==pp) 
				{
					if (fileout_) fprintf(fileout_, "Erroneus target altitude. Skipping VNAV line.\n");
					continue;
				}
				vmaneuver =  VNAV_MODE_VS;
			}

			// ----------------------
			// FPA mode
			p1 = strstr(p,"FPA");
			if (p1 && (vmaneuver == -1))
			{
				p+=strlen("FPA");
				while((*p==' ') || (*p=='	') && p<(&buf[n])) p++;	// skip spaces etc.
				char* pp;
				param  = strtod(p,&pp);	// flight path angle (deg)
				if (p==pp) 
				{
					if (fileout_) fprintf(fileout_, "Erroneus FPA. Skipping VNAV line.\n");
					continue;
				}
				p = pp;
				while((*p==' ') || (*p=='	') && p<(&buf[n])) p++;	// skip spaces etc.
				target_alt = strtod(p,&pp);	// target altitude (m)
				if (p==pp) 
				{
					if (fileout_) fprintf(fileout_, "Erroneus target altitude. Skipping VNAV line.\n");
					continue;
				}
				vmaneuver =  VNAV_MODE_FPA;
			}

			// ----------------------
			// FLCH mode
			p1 = strstr(p,"FLCH");
			if (p1 && (vmaneuver == -1))
			{
				p+=strlen("FLCH");
				while((*p==' ') || (*p=='	') && p<(&buf[n])) p++;	// skip spaces etc.
				char* pp;
				target_alt = strtod(p,&pp);	// target altitude (m)
				if (p==pp) 
				{
					if (fileout_) fprintf(fileout_, "Erroneus FLCH altitude. Skipping VNAV line.\n");
					continue;
				}
				param = 125.0;	// seconds
				vmaneuver =  VNAV_MODE_FLCH;
			}

			if (vmaneuver == -1)
			{
				if (fileout_)
				{
					fprintf(fileout_,"Erroneus VNAV command. Must be: FL/LEVEL, VS, FPA or FLCH.\n");
					fprintf(fileout_,"Skipping VNAV line.\n");
					continue;
				}
			}

			// --------------------------------------
			// Sanity checks:
			if (vmaneuver == VNAV_MODE_LEVEL)
			{
				if (param<=0)
				{
					if (fileout_) fprintf(fileout_,"Zero or negative duration of level flight. Skipping VNAV line.\n");
					continue;
				}
			}


			if (vmaneuver == VNAV_MODE_VS || vmaneuver == VNAV_MODE_FPA || vmaneuver == VNAV_MODE_FLCH)
			{
				if (target_alt<0)
				{
					if (fileout_) fprintf(fileout_,"Warning: input target altitude < 0. Reset to 0.\n");
					target_alt = 0.0;
				}

				if (target_alt>15000)
				{
					if (fileout_) fprintf(fileout_,"Warning: too high altitude in VNAV file. Resetting to 15000 m.\n");
					target_alt = 15000;
				}
			}

			if (vmaneuver == VNAV_MODE_VS)
			{
				if (param < (-8000.0*0.3048/60.0))
				{
					if (fileout_) fprintf(fileout_,"Warning: too fast rate of descent. Resetting to -8000 ft/minute.\n");
					param = -8000.0*0.3048/60.0;
				}

				if (param > (6000.0*0.3048/60.0))
				{
					if (fileout_) fprintf(fileout_,"Warning: too fast rate of climb. Resetting to 6000 ft/minute.\n");
					param = 6000.0*0.3048/60.0;
				}
			}

			if (vmaneuver == VNAV_MODE_FPA)
			{
				if (param < -9.9)
				{
					if (fileout_) fprintf(fileout_,"Warning: too fast FPA descent. Resetting to -9.9 deg.\n");
					param = -9.9;
				}

				if (param > 9.9)
				{
					if (fileout_) fprintf(fileout_,"Warning: too fast FPA climb. Resetting to 9.9 deg.\n");
					param = 9.9;
				}
			}

			// -----------------------------------------
			// here everything is ok; reallocate new vnav profile
			VNAVProfile* pVNAVProfile_new = new VNAVProfile[nVManeuvers_+1];
			for (int i = 0; i<nVManeuvers_; i++)
			{
				pVNAVProfile_new[i].maneuverType = pVNAVProfile_[i].maneuverType;
				pVNAVProfile_new[i].targetAltitude = pVNAVProfile_[i].targetAltitude;
				pVNAVProfile_new[i].param = pVNAVProfile_[i].param;
				pVNAVProfile_new[i].startTime = 0.0;
				pVNAVProfile_new[i].endTime = 1.0E+10;
				pVNAVProfile_new[i].startAltitude = 0.0;
			}
			pVNAVProfile_new[nVManeuvers_].maneuverType = vmaneuver;
			pVNAVProfile_new[nVManeuvers_].targetAltitude = target_alt;
			pVNAVProfile_new[nVManeuvers_].param = param;
			pVNAVProfile_new[nVManeuvers_].startTime = 0.0;
			pVNAVProfile_new[nVManeuvers_].endTime = 1.0E+10;
			pVNAVProfile_new[nVManeuvers_].startAltitude = 0.0;

			if (pVNAVProfile_) delete pVNAVProfile_;
			pVNAVProfile_ = pVNAVProfile_new;
			nVManeuvers_++;

		} // fgets
	} //feof - while loop

	if (fid) fclose(fid);

	// ------------------------------------------
	// The last maneuver must be keeping constant level; if not, add it to the list with unlimited time
	if (nVManeuvers_>0)
	{
		if (pVNAVProfile_[nVManeuvers_-1].maneuverType != VNAV_MODE_LEVEL)
		{
			VNAVProfile* pVNAVProfile_new = new VNAVProfile[nVManeuvers_+1];
			for (int i = 0; i<nVManeuvers_; i++)
			{
				pVNAVProfile_new[i].maneuverType = pVNAVProfile_[i].maneuverType;
				pVNAVProfile_new[i].targetAltitude = pVNAVProfile_[i].targetAltitude;
				pVNAVProfile_new[i].param = pVNAVProfile_[i].param;
				pVNAVProfile_new[i].startTime = 0.0;
				pVNAVProfile_new[i].endTime = 1.0E+10;
				pVNAVProfile_new[i].startAltitude = 0.0;
			}
			pVNAVProfile_new[nVManeuvers_].maneuverType = VNAV_MODE_LEVEL;
			pVNAVProfile_new[nVManeuvers_].targetAltitude = pVNAVProfile_new[nVManeuvers_-1].targetAltitude;
			pVNAVProfile_new[nVManeuvers_].param = 1.0E+10;	// duration
			pVNAVProfile_new[nVManeuvers_].startTime = 0.0;
			pVNAVProfile_new[nVManeuvers_].endTime = 1.0E+10;
			pVNAVProfile_new[nVManeuvers_].startAltitude = pVNAVProfile_new[nVManeuvers_-1].targetAltitude;

			if (pVNAVProfile_) delete pVNAVProfile_;
			pVNAVProfile_ = pVNAVProfile_new;
			nVManeuvers_++;
		}
	}

	// Assign start time and start altitude
	if (nVManeuvers_>0)
	{
		pVNAVProfile_[0].startTime = pInpParams->ts;
		pVNAVProfile_[0].startAltitude = pInpParams->alt0;
	}



	return true;
}



///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Get vertical speed, which depends on current altitude and may depend on TAS in case of FPA mode 
// this function is called from CManeuver's calc_uvwq
// specified vertical profile of a flight path is independent on horizontal one.
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////

double CTrajectory::GetVerticalSpd_(double t, double alt, double TAS)
{
	// Input:
	// t  - time (s since 2017-03-07 00:00:00)
	// alt - current altitude (m)
	// TAS - current TAS (m/s)
	// dt_max is the maximum integration time step, which could be helphul to suppress overshots and undershots
	//        of the target altitude, which occurs in case of 

    //  --------------------------------------------
	// vertical velocity component
	// V/S: -8000 to +6000 fpm, 100 fpm increments - 777 Flight Manual, Continental Rev. 11/01/02 #9
	// FPA: -9.9 to +9.9 degrees, 0.1 deg increments - 777 Flight Manual, Continental Rev. 11/01/02 #9
	// FLCH, aimed to complete flight level change within 125 seconds, CAL_Training_Manual_B777-200_WB371.pdf - 
	//       Boeing 777 Training manual, Coninental Airlines Inc., Airplane General, WB371


	// keep flight level if there no vertical maneuvers
	if (nVManeuvers_ == 0) return 0.0;
	if (cVmaneuver_ >= nVManeuvers_) return 0.0; // this case is not supposed to happen - just for failsafe

	double ws = 0.0;

	int vnav_mode = pVNAVProfile_[cVmaneuver_].maneuverType;
	double param = pVNAVProfile_[cVmaneuver_].param; 
	double target_altitude = pVNAVProfile_[cVmaneuver_].targetAltitude;
	double start_altitude = pVNAVProfile_[cVmaneuver_].startAltitude;

	double diff_alt = target_altitude - start_altitude;	// altitude difference

	switch(vnav_mode)
	{

	case VNAV_MODE_LEVEL:
		// Keep flight level
		ws = 0.0;
		break;

	case VNAV_MODE_VS:
		// Ascent or descent at specified RoC or RoD till target altitude is reached 
		if ((param>0.0 && alt<target_altitude) || (param<0.0 && alt>target_altitude)) 
		{
			ws = param;
		}
		else
		{
			ws = 0.0;
		}

		break;


	case VNAV_MODE_FPA:
		// ascent or descent with specifdied flight path angle till target altitude is reached 
		if ((param>0.0 && alt<target_altitude) || (param<0.0 && alt>target_altitude)) 
		{
			ws = TAS*sin(param*3.1415926536/180.0);
		}
		else
		{
			ws = 0.0;
		}
		break;


	case VNAV_MODE_FLCH:
		// change in altitude occurs within 125 s per Continental manual
		if ((diff_alt>0.0 && alt<target_altitude) || (diff_alt<0.0 && alt>target_altitude)) 
		{
			ws = diff_alt/125.0;
		}
		else
		{
			ws = 0.0;
		}
		break;
	}

	return ws;
}



///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Check conditions to move to next vertical maneuver. If conditions are met, move to next maneuver
// and return true
//
///////////////////////////////////////////////////////////////////////////////////////////////////

bool CTrajectory::UpdateVNavManeuver_(double t, double *XYZM)
{
	// keep last maneuver if there are no other maneuvers
	if (nVManeuvers_ == 0) return false;
	if (cVmaneuver_ >= nVManeuvers_) return false;

	int vnav_mode = pVNAVProfile_[cVmaneuver_].maneuverType;
	double param = pVNAVProfile_[cVmaneuver_].param; 
	double target_altitude = pVNAVProfile_[cVmaneuver_].targetAltitude;
	double start_altitude = pVNAVProfile_[cVmaneuver_].startAltitude;
	double diff_alt = target_altitude - start_altitude;	// altitude difference

	const double alt_eps = 0.01;		// altitude margin

	// ---------------
	// convert X, Y, Z to lon, lat, alt
	double lon, lat, alt;
	XYZ2LonLatH(lon, lat, alt, XYZM[0], XYZM[1], XYZM[2]);


	switch(vnav_mode)
	{
	case VNAV_MODE_LEVEL:
		// Flight level
		if (t<pVNAVProfile_[cVmaneuver_].startTime + param)
		{
			return false;	// no changes required
		}
		else
		{
			return true;
		}
		break;


	case VNAV_MODE_VS:
		// Ascent or descent at specified RoC or RoD till target altitude is reached 
		if ((param>0.0 && alt<target_altitude-alt_eps) || (param<0.0 && alt>target_altitude+alt_eps)) 
		{
			return false;	// no changes required
		}
		else
		{
			return true;
		}

		break;



	case VNAV_MODE_FPA:
		// Ascent or descent with specifdied flight path angle till target altitude is reached 
		if ((param>0.0 && alt<target_altitude-alt_eps) || (param<0.0 && alt>target_altitude+alt_eps)) 
		{
			return false;	// no changes required
		}
		else
		{
			return true;
		}
		break;


	case VNAV_MODE_FLCH:
		// change in altitude occurs within 125 s per Continental manual
		if ((diff_alt>0.0 && alt<target_altitude) || (diff_alt<0.0 && alt>target_altitude)) 
		{
			return false;	// no changes required
		}
		else
		{
			return true;
		}
		break;
	}

	return false;
}



///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Return suggested maximum time step based on vertical flight profile to minimize error due to mode switching 
//
///////////////////////////////////////////////////////////////////////////////////////////////////

double CTrajectory::GetMaxVNAVStep_(double t, double *XYZM)
{
	// keep last maneuver if there are no other maneuvers
	if (nVManeuvers_ == 0) return 1.0E+10; // no limit
	if (cVmaneuver_ >= nVManeuvers_) return 1.0E+10; // no limit

	int vnav_mode = pVNAVProfile_[cVmaneuver_].maneuverType;
	double param = pVNAVProfile_[cVmaneuver_].param; 
	double target_altitude = pVNAVProfile_[cVmaneuver_].targetAltitude;
	double start_altitude = pVNAVProfile_[cVmaneuver_].startAltitude;
	double diff_alt = target_altitude - start_altitude;	// altitude difference

	const double alt_eps = 5.0;		// altitude margin (allowable overshots, m)

	// ---------------
	// convert X, Y, Z to lon, lat, alt
	double lon, lat, alt;
	XYZ2LonLatH(lon, lat, alt, XYZM[0], XYZM[1], XYZM[2]);

	double max_dt = 1.0E+10;

	switch(vnav_mode)
	{
	case VNAV_MODE_LEVEL:
		// Flight level - allow for 1s longer than time needed to reach the end
		if (t<pVNAVProfile_[cVmaneuver_].startTime + param) max_dt = pVNAVProfile_[cVmaneuver_].startTime + param - t + 0.1;
		break;

	case VNAV_MODE_VS:
		// Ascent or descent at specified RoC or RoD till target altitude is reached 
		if ((param>0.0 && alt<target_altitude) || (param<0.0 && alt>target_altitude)) max_dt = (abs(target_altitude-alt)+alt_eps)/(abs(param)+1.0E-6);
		break;


	case VNAV_MODE_FPA:
		// Ascent or descent with specifdied flight path angle till target altitude is reached; assume TAS = 250 m/s - close to maximum; otherwise need to pass TAS or calculate it here
		if ((param>0.0 && alt<target_altitude) || (param<0.0 && alt>target_altitude)) max_dt = (abs(target_altitude-alt)+alt_eps)/(250.0*abs(sin(param*3.1415926536/180.0))+1.0E-6);
		break;


	case VNAV_MODE_FLCH:
		// change in altitude occurs within 125 s per Continental manual
		if ((param>0.0 && alt<target_altitude) || (param<0.0 && alt>target_altitude)) max_dt = (abs(target_altitude-alt)+alt_eps)/(abs(diff_alt)/125.0+1.0E-6);
		break;
	}

	return max_dt;
}




///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Move to next VNAV maneuver
//
///////////////////////////////////////////////////////////////////////////////////////////////////

void CTrajectory::MoveToNextVNAVManeuver_(double& start_t, double* pXYZMNew, 
										  double t0, double* pXYZM0, 
										  double t1, double* pXYZM1)
{
	// Input: t0 - the second last time, when current VNAV maneuver is still valid
	//        pXYZM0 - ECEF position and weight at t0
	//        t1 - the fist time, when current VNAV maneuver becomes invalid
	//        pXYZM1 - ECEF position and weight at t1
	//        
	// Output: start_t - start time of new VNAV maneuver (since 2014-03-07 00:00:00)
	//         pXYZMNew - position in ECEF and weight at the beginning of next maneuver

	if (nVManeuvers_ == 0) return;
	if (cVmaneuver_ >= nVManeuvers_-1) return;

	int vnav_mode = pVNAVProfile_[cVmaneuver_].maneuverType;
	double param = pVNAVProfile_[cVmaneuver_].param; 
	double target_altitude = pVNAVProfile_[cVmaneuver_].targetAltitude;
	double start_altitude = pVNAVProfile_[cVmaneuver_].startAltitude;

	// ---------------
	// convert X, Y, Z to lon, lat, alt
	double lon0, lat0, alt0;
	XYZ2LonLatH(lon0, lat0, alt0, pXYZM0[0], pXYZM0[1], pXYZM0[2]);

	double lon1, lat1, alt1;
	XYZ2LonLatH(lon1, lat1, alt1, pXYZM1[0], pXYZM1[1], pXYZM1[2]);

	double coef;

	switch(vnav_mode)
	{
	case VNAV_MODE_LEVEL:

		// Flight level
		start_t = pVNAVProfile_[cVmaneuver_].startTime + param;
		coef = (start_t-t0)/(t1-t0);

		// update time and altitude to start next vertical maneuver
		pVNAVProfile_[cVmaneuver_].endTime = pVNAVProfile_[cVmaneuver_].startTime + param;
		pVNAVProfile_[cVmaneuver_].targetAltitude = pVNAVProfile_[cVmaneuver_].startAltitude;

		if (cVmaneuver_+1<nVManeuvers_)	
		{
			pVNAVProfile_[cVmaneuver_+1].startTime = pVNAVProfile_[cVmaneuver_].endTime;
			pVNAVProfile_[cVmaneuver_+1].startAltitude = pVNAVProfile_[cVmaneuver_].targetAltitude;
		}

		break;

	case VNAV_MODE_VS:
	case VNAV_MODE_FPA:
	case VNAV_MODE_FLCH:

		// Ascent or descent 
		coef = (target_altitude-alt0)/(alt1-alt0);
		start_t = (1.0-coef)*t0 + coef*t1;

		pVNAVProfile_[cVmaneuver_].endTime = start_t;
		if (cVmaneuver_+1<nVManeuvers_)	
		{
			pVNAVProfile_[cVmaneuver_+1].startTime = start_t;
			pVNAVProfile_[cVmaneuver_+1].startAltitude = pVNAVProfile_[cVmaneuver_].targetAltitude;
		}

		break;
	}

	// linearly interpolate initial conditions for next VNAV maneuver to avoid overshots and undershots
	pXYZMNew[0] = (1.0-coef)*pXYZM0[0] + coef*pXYZM1[0];
	pXYZMNew[1] = (1.0-coef)*pXYZM0[1] + coef*pXYZM1[1];
	pXYZMNew[2] = (1.0-coef)*pXYZM0[2] + coef*pXYZM1[2];
	pXYZMNew[3] = (1.0-coef)*pXYZM0[3] + coef*pXYZM1[3];


	cVmaneuver_++;
}
