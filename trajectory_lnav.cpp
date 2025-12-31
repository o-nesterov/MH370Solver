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



///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Create / read LNAV profile
//
///////////////////////////////////////////////////////////////////////////////////////////////////

// read num_iparams of integer parameters and num_fparams of fp paramaters, and accompanying optimization flag from string p
bool readManeuverParameters(char* p_in, int num_iparams, int iparams[MAX_MANEUVER_IPARAMS], int num_fparams, double fparams[MAX_MANEUVER_FPARAMS], bool ofparams[MAX_MANEUVER_FPARAMS])
{
	char* p = p_in;
	int n = (int) strlen(p_in);

	for (int i=0; i<MAX_MANEUVER_IPARAMS; i++) iparams[i] = 0;	// initialize
	for (int i=0; i<MAX_MANEUVER_FPARAMS; i++) fparams[i] = 0.0;
	for (int i=0; i<MAX_MANEUVER_FPARAMS; i++) ofparams[i] = false;

	// read integer parameters
	for (int i=0; i<num_iparams; i++)
	{
		while((*p==' ') || (*p==',') || (*p==';') && p<(&p_in[n])) p++;	// skip spaces, etc.
		if (p==&p_in[n]) return false;

		char* pp;
		iparams[i] = strtol(p,&pp,10);
		if (p==pp) return false;	// nothing was scanned
		p = pp;
	}

	// read fp parameters
	for (int i=0; i<num_fparams; i++)
	{
		while((*p==' ') || (*p==',') || (*p==';') && p<(&p_in[n])) p++;	// skip spaces, etc.
		if (p==&p_in[n]) return false;

		char* pp;

		if (*p=='*') {ofparams[i] = true; p++;};

		fparams[i] = strtod(p,&pp);
		if (p==pp) return false;	// nothing was scanned
		p = pp;

		if (*p=='*') {ofparams[i] = true; p++;};
	}

	return true;
}



// read num_iparams of integer parameters and num_fparams of fp paramaters with character attribute from string p
bool readManeuverParametersA(char* p_in, int num_iparams, int iparams[MAX_MANEUVER_IPARAMS], int num_fparams, double fparams[MAX_MANEUVER_FPARAMS], char cfparams[MAX_MANEUVER_FPARAMS], bool ofparams[MAX_MANEUVER_FPARAMS])
{
	char* p = p_in;
	int n = (int) strlen(p_in);

	for (int i=0; i<MAX_MANEUVER_IPARAMS; i++) iparams[i] = 0;	// initialize
	for (int i=0; i<MAX_MANEUVER_FPARAMS; i++) fparams[i] = 0.0;
	for (int i=0; i<MAX_MANEUVER_FPARAMS; i++) ofparams[i] = false;

	// read integer parameters
	for (int i=0; i<num_iparams; i++)
	{
		while((*p==' ') || (*p==',') || (*p==';') && p<(&p_in[n])) p++;	// skip spaces, etc.
		if (p==&p_in[n]) return false;

		char* pp;
		iparams[i] = strtol(p,&pp,10);
		if (p==pp) return false;	// nothing was scanned
		p = pp;
	}

	// read fp parameters
	for (int i=0; i<num_fparams; i++)
	{
		while((*p==' ') || (*p==',') || (*p==';') && p<(&p_in[n])) p++;	// skip spaces, etc.
		if (p==&p_in[n]) return false;

		if (*p=='*') {ofparams[i] = true; p++;}; // check optimization flag

		char* pp;
		fparams[i] = strtod(p,&pp);
		if (p==pp)
		{
			// try to read prefix character attribute first
			if (((*p)>='A' && (*p)<='Z') || ((*p)>='a' && (*p)<='z'))
			{
				cfparams[i] = *p;
				p++;
				fparams[i] = strtod(p,&pp);
				if (p==pp) return false; // incorrect format
				p = pp;	// move to the end of variable
			}
			else
			{
				return false; // nothing to read - incorrect format
			}
		}
		else
		{
			if (((*pp)>='A' && (*pp)<='Z') || ((*pp)>='a' && (*pp)<='z')) // try to read suffix attribute
			{
				cfparams[i] = *pp;
				p = pp+1;
			}
			else
			{
				cfparams[i] = 0; // fp variable was read, but there is neither prefix nor suffix
				p = pp;
			}
		}

		if (*p=='*') {ofparams[i] = true; p++;}; // check for optimization flag
	}

	return true;
}




////////////////////////////////////////////////////////////////////////////////////////
//
// Load lateral navigation profile from ASCII file or create a single leg
//
////////////////////////////////////////////////////////////////////////////////////////


bool CTrajectory::CreateLNAVProfile(CInpParams* pInpParams)
{
	// Input: filename - name of the file containing the list of maneuvers
	char* filename = pInpParams->pLNAVfilename;

	char buf[1024];
	int iparams[MAX_MANEUVER_IPARAMS];		// integer maneuver parameters
	char cfparams[MAX_MANEUVER_IPARAMS];	// auxilary character at fp parameter
	double fparams[MAX_MANEUVER_FPARAMS];	// fp maneuver parameters
	bool ofparams[MAX_MANEUVER_FPARAMS];	// fp maneuver parameters optimization flags
	double fparams1[MAX_MANEUVER_FPARAMS];	// fp maneuver parameters - second set for turns
	bool ofparams1[MAX_MANEUVER_FPARAMS];	// fp maneuver parameters - second set for turns optimization flags

	// if no filename specified, then create a simple leg
	if (filename == NULL)
	{
		for (int i=0; i<MAX_MANEUVER_IPARAMS; i++) iparams[i] = 0;	// initialize
		for (int i=0; i<MAX_MANEUVER_FPARAMS; i++) fparams[i] = 0.0;
		for (int i=0; i<MAX_MANEUVER_FPARAMS; i++) ofparams[i] = false;
		fparams[0] = pInpParams->hdg0;
		fparams[1] = 1.0E+10;	// fictitious duration indicating integration of the leg till the end
		fparams[2] = 30.0;		// maximum integration time step and regular output (s)

		addManeuver2List_(MANEUVER_TYPE_HLEG, 0, iparams, 3, fparams, ofparams, true); // add maneuver with 'infinite' duration
		return true;	// input file was not specified, but constants specified in pInpParams are sufficient
	}



	// ----------------------
	// read sequence of maneuvers from the specified file
	FILE* fid = fopen(filename,"rt");

	if (fid == NULL)
	{
		if (fileout_) fprintf(fileout_,"Cannot open LNAV file: %s\n", filename);
		return false; // file was specified, but it cannot be opened, so return false
	}

	// reset profile
	if (pLNAVProfile_) delete pLNAVProfile_;
	pLNAVProfile_ = NULL;

	double hdg = pInpParams->hdg0; // set initial heading or track; it must be in [0, 360) deg range

	// parse lateral navigation file
	while(!feof(fid))
	{
		if (fgets(buf, 1000, fid)!=NULL)
		{
			char* p=&buf[0];
			int n = (int)strlen(p);
			while((*p==' ') || (*p=='	') && p<(&buf[n])) p++;	// skip spaces
			if (p==&buf[n] || n<=3) continue;	// skip empty line
			if (*p=='#' || *p=='!' || *p=='%%') continue; // this is a comment line - skip it

			// ----------------------
			// Simple leg maneuver - keep previous heading
			char* p1 = strstr(p,"LEG");
			char* p2 = strstr(p,"leg");
			char* p3 = strstr(p,"Leg");
			if (p1!=NULL || p2!=NULL || p3!=NULL)
			{
				// maneuver name must come first
				if (!(p==p1 || p==p2 || p==p3))
				{
					if (fileout_) fprintf(fileout_,"Incorrect line format in LNAV file. Skipping line.\n");
					continue;
				}
				else
				{
					p+=3;
					// read maneuver parameters (paramsi - integer; paramsf - floating point) from the string p
					if(!readManeuverParameters(p, 0, iparams, 2, fparams, ofparams))
					{
						if (fileout_) fprintf(fileout_,"Incorrect line format in LNAV file. Skipping line.\n");
						continue;
					}
					else
					{
						ofparams[1] = false; // only duration can be optimized
						addManeuver2List_(MANEUVER_TYPE_LEG, 0, iparams, 3, fparams, ofparams, true);
					}
				}
			}




			// ----------------------
			// Leg with specified heading; if heading does not match the previous one, then insert turn maneuver
			p1 = strstr(p,"HLG");
			p2 = strstr(p,"hlg");
			p3 = strstr(p,"Hlg");
			if (p1!=NULL || p2!=NULL || p3!=NULL)
			{
				// maneuver name must come first
				if (!(p==p1 || p==p2 || p==p3))
				{
					if (fileout_) fprintf(fileout_,"Incorrect line format in LNAV file. Skipping line.\n");
					continue;
				}
				else
				{
					p+=3;
					// read maneuver parameters (paramsi - integer; paramsf - floating point, cfparams - character attribute at fp variable) from the string p
					if(!readManeuverParametersA(p, 0, iparams, 3, fparams, cfparams, ofparams))
					{
						if (fileout_) fprintf(fileout_,"Incorrect line format in LNAV file. Skipping line.\n");
						continue;
					}
					else
					{
						double hdg_new = fparams[0];	// heading of this leg
						if (hdg_new>=360.0 || hdg_new<0.0)
						{
							if (fileout_) fprintf(fileout_, "Error: heading or track must be in the range [0,360) deg.\n");
							if (fileout_) fprintf(fileout_, "Use Turn maneuver if a L or R turn exceeding 180 deg implied.\n");
							break;
						}

						// now both hdg and hdg_leg are in [0,360) range
						fparams1[0] = hdg_new;				// target heading
						ofparams1[0] = false;				// set optimize flag to false, but update target heading as needed in ResetLNAVOptimizationParameters

						fparams1[1] = -1.0;					// use default bank angle default_bank_angle_, which could be reset by specifying it in a setup file, and calling SetDefaultBankAngle 
						ofparams1[1] = false;				// do not optimize this bank angle individually (optimize the default bank angle still applies)

						fparams1[2] = 2.0;					// time step output (s)
						ofparams1[2] = false;				// do not optimize

						iparams[0] = TURN_TYPE_TO_HDG;		// indicate that parameter means turn to specified HDG
						
						// additional parameter - specified heading: true, magnetic or gyroscopic
						// they cannot be converted here right away because conversion varies in time and space
						// if mode is not specified, then assume that heading is specified as true
						iparams[1] = HEADING_MODE_TRUEHDG;
						switch (cfparams[0])
						{
						case 'T':
							iparams[1] = HEADING_MODE_TRUEHDG;
							break;

						case 'M':
							iparams[1] = HEADING_MODE_MAGNHDG;
							break;

						case 'G':
							iparams[1] = HEADING_MODE_GYROHDG;
							break;
						};

						ofparams1[2] = false;	// only heading and duration can be optimized
						ofparams[2] = false;

						addManeuver2List_(MANEUVER_TYPE_TURN, 2, iparams, 3, fparams1, ofparams1, false);
						addManeuver2List_(MANEUVER_TYPE_HLEG, 2, iparams, 3, fparams, ofparams, true);

					}
				}
			}





			// ----------------------
			// Turn maneuver
			p1 = strstr(p,"TRN");
			p2 = strstr(p,"turn");
			p3 = strstr(p,"TURN");
			char* p4 = strstr(p,"Turn");
			if (p1!=NULL || p2!=NULL || p3!=NULL || p4!=NULL)
			{
				// maneuver name must come first
				if (!(p==p1 || p==p2 || p==p3 || p==p4))
				{
					if (fileout_) fprintf(fileout_,"Incorrect line format in LNAV file. Skipping line.\n");
					continue;
				}
				else
				{
					if (p1)
					{
						p+=3;
					}
					else
					{
						p+=4;
					}

					// read maneuver parameters (paramsi - integer; paramsf - floating point) from string p
					if(!readManeuverParameters(p, 0, iparams, 3, fparams, ofparams))	// read 3 FP parameters (turn angle (deg), and bank angle (deg), time step(s)) from the string p
					{
						if (fileout_) fprintf(fileout_,"Incorrect line format in LNAV file. Skipping line.\n");
						continue;
					}
					else
					{
						// add maneuver to the list
						iparams[0] = TURN_TYPE_BY_ANGLE;		// indicate that parameter means turn by specified angle

						ofparams[2] = false;					// only heading and bank angle can be optimized

						addManeuver2List_(MANEUVER_TYPE_TURN, 1, iparams, 3, fparams, ofparams, true);

						hdg += fparams[0]; // add change in heading

						// use simple method to reduce heading to the range [0,360) deg.
						while (hdg<0.0) hdg += 360.0; // make sure that 360.0>hdg0>=0.0
						while (hdg>=360.0) hdg -=360.0;
					}
				}
			}
		} // fgets
		else
		{
			if (ferror(fid))
			{
				if (fileout_) fprintf(fileout_,"Failure to read LNAV file.\n");
				fclose(fid);
				return false;
			}
			else
			{
				if (nManeuvers_==0)
				{
					if (fileout_) fprintf(fileout_,"No valid maneuver was read. Wrong formatting?\n");
					fclose(fid);
					return false;
				}
				else
				{
					fclose(fid);
					return true; // end of file
				}
			}
		}
	}

	fclose(fid);

	if (nManeuvers_==0)
	{
		if (fileout_) fprintf(fileout_,"No valid maneuver was read. Wrong formatting?\n");
		return false;
	}

	return true; // everything is ok
}





// ----------------------------
// Add a maneuver to the list
void CTrajectory::addManeuver2List_(int maneuvertype, unsigned long num_iparams, int iparams[MAX_MANEUVER_IPARAMS], 
													  unsigned long num_fparams, double fparams[MAX_MANEUVER_FPARAMS],
													  bool ofparams[MAX_MANEUVER_FPARAMS],
													  bool isOriginal)
{
	// create a new pLNAVProfile_
	LNAVProfile* pLNAVProfile = new LNAVProfile[nManeuvers_+1];

	if (pLNAVProfile_)
	{
		for (int i=0; i<(int)nManeuvers_; i++)
		{
			pLNAVProfile[i].maneuverType = pLNAVProfile_[i].maneuverType;
			for (int j=0; j<MAX_MANEUVER_IPARAMS; j++) pLNAVProfile[i].iparams[j] = pLNAVProfile_[i].iparams[j];
			for (int j=0; j<MAX_MANEUVER_FPARAMS; j++) pLNAVProfile[i].fparams[j] = pLNAVProfile_[i].fparams[j];
			for (int j=0; j<MAX_MANEUVER_FPARAMS; j++) pLNAVProfile[i].ofparams[j] = pLNAVProfile_[i].ofparams[j];
			pLNAVProfile[i].isOrig = pLNAVProfile_[i].isOrig;
		}

		delete pLNAVProfile_; // delete old list
	}

	// add maneuver
	pLNAVProfile[nManeuvers_].maneuverType = maneuvertype;
	for (int j=0; j<MAX_MANEUVER_IPARAMS; j++) pLNAVProfile[nManeuvers_].iparams[j] = 0;		// initialize to 0
	for (int j=0; j<MAX_MANEUVER_FPARAMS; j++) pLNAVProfile[nManeuvers_].fparams[j] = 0.0;		// initialize to 0.0
	for (int j=0; j<(int)num_iparams; j++) pLNAVProfile[nManeuvers_].iparams[j] = iparams[j];	// copy provided integer parameters
	for (int j=0; j<(int)num_fparams; j++) pLNAVProfile[nManeuvers_].fparams[j] = fparams[j];	// copy provided fp parameters
	for (int j=0; j<(int)num_fparams; j++) pLNAVProfile[nManeuvers_].ofparams[j] = ofparams[j];	// copy provided optimization flags for fp parameters (yes/no)
	pLNAVProfile[nManeuvers_].isOrig = isOriginal;												// flag indicating if the maneuver was originally specified (true) or automatically inserted  - for output

	// replace old list with new one and increase counter of maneuvers
	pLNAVProfile_ = pLNAVProfile;
	nManeuvers_++;
}


// ----------------------------
// Overwrite LNAV optimization parameters (used for optimization only) that were marked with * in LNAV profile
void CTrajectory::ResetLNAVOptimizationParameters(double* params, int nparams)
{
	if (nManeuvers_ == 0) return;

	int i = 0;

	for (int n=0; n<nManeuvers_; n++)
	{
		if (i==nparams) break;

		switch (pLNAVProfile_[n].maneuverType)
		{
		case MANEUVER_TYPE_LEG:
			if (pLNAVProfile_[n].ofparams[0]) {if (params[i]>0.1) {pLNAVProfile_[n].fparams[0] = params[i];}; i++;}; // reset leg duration
			break;

		case MANEUVER_TYPE_HLEG:
			// reset heading
			if (pLNAVProfile_[n].ofparams[0]) 
			{
				pLNAVProfile_[n].fparams[0] = params[i]; 
				if (n>=1)
				{
					// check if preceeding turn needs to be updated to
					if (pLNAVProfile_[n-1].maneuverType == MANEUVER_TYPE_TURN)
					{
						if (pLNAVProfile_[n-1].iparams[0] = TURN_TYPE_TO_HDG)
						{
							pLNAVProfile_[n-1].fparams[0] = params[i];
						}
					}
				}
				i++;
			}; 
		
			// reset leg duration
			if (pLNAVProfile_[n].ofparams[1]) {if (params[i]>0.1) {pLNAVProfile_[n].fparams[1] = params[i];}; i++;};
			break;

		case MANEUVER_TYPE_TURN:
			if (pLNAVProfile_[n].ofparams[0]) {pLNAVProfile_[n].fparams[0] = params[i]; i++;}; // reset turn angle
			if (pLNAVProfile_[n].ofparams[1]) {if (params[i]>0.1) {pLNAVProfile_[n].fparams[1] = params[i];}; i++;}; // reset bank angle
			break;
		};
	};
};




// ----------------------------
// get initial LNAV parameters for optimization (used for optimization only); also returns number of parameters nparams
void CTrajectory::GetLNAVOptimizationParameters(double* params, int& nparams, int maxparams)
{
	nparams = 0;

	if (nManeuvers_ == 0) return;


	for (int n=0; n<nManeuvers_; n++)
	{
		if (nparams == maxparams) break;

		switch (pLNAVProfile_[n].maneuverType)
		{
		case MANEUVER_TYPE_LEG:
			if (pLNAVProfile_[n].ofparams[0]) {params[nparams] = pLNAVProfile_[n].fparams[0]; nparams++;}; // leg duration
			break;

		case MANEUVER_TYPE_HLEG:
			if (pLNAVProfile_[n].ofparams[0]) {params[nparams] = pLNAVProfile_[n].fparams[0]; nparams++;}; // heading
			if (pLNAVProfile_[n].ofparams[1]) {params[nparams] = pLNAVProfile_[n].fparams[1]; nparams++;}; // leg duration
			break;

		case MANEUVER_TYPE_TURN:
			if (pLNAVProfile_[n].ofparams[0]) {params[nparams] = pLNAVProfile_[n].fparams[0]; nparams++;}; // turn angle
			if (pLNAVProfile_[n].ofparams[1]) {params[nparams] = pLNAVProfile_[n].fparams[1]; nparams++;}; // bank angle
			break;
		};
	};
}




// ----------------------------
// save lateral navigation file (to be called after optimization with lnav flag)
bool CTrajectory::SaveLNAVProfile(char* filename)
{
	if (nManeuvers_ == 0) return true;		// nothing to save
	if (pLNAVProfile_ == NULL) return true;	//nothing to save

	FILE* fid = fopen(filename, "wt");
	if (fid == NULL) 
	{
		if (fileout_) fprintf(fileout_, "Unable to save file.\n");
		return false;
	}

	// header
	fprintf(fid,"#####################################################################################################\n");
	fprintf(fid,"# Lateral navigation file - list of maneuvers. Supported maneuvers:\n");
	fprintf(fid,"# 1). LEG (also LEG or LEG) - simple leg: duration (s),  and output time step\n");
	fprintf(fid,"# 2). HLG (also hlg or Hlg) - leg sith heading: heading/track (depending on mode),  duration (s),  and output time step\n");
	fprintf(fid,"# 3). TRN (also TURN, Turn or turn): turn by (deg, CW), bank angle (deg), time step (s)\n");
	fprintf(fid,"#\n");
	fprintf(fid,"# If LEG's heading in inconsistent with the previous one using HLG, then a turn is automatically inserted\n");
	fprintf(fid,"#####################################################################################################\n");
	fprintf(fid," \n");

	// ---------

	for (int i=0; i<nManeuvers_; i++)
	{
		if (!(pLNAVProfile_[i].isOrig)) continue;

		switch (pLNAVProfile_[i].maneuverType)
		{
		case MANEUVER_TYPE_LEG:
			fprintf(fid, "LEG %f", pLNAVProfile_[i].fparams[0]);	// print duration
			if (pLNAVProfile_[i].ofparams[0]) fprintf(fid,"*");		// print optimization mark if needed
			fprintf(fid, ", %f\n", pLNAVProfile_[i].fparams[1]);	// print time step
			break;
		
		case MANEUVER_TYPE_HLEG:
			fprintf(fid, "HLG ");
			if (pLNAVProfile_[i].iparams[1] == HEADING_MODE_TRUEHDG) fprintf(fid, "T"); // print heading reference (true, magnetic or gyroscopic)
			if (pLNAVProfile_[i].iparams[1] == HEADING_MODE_MAGNHDG) fprintf(fid, "M");
			if (pLNAVProfile_[i].iparams[1] == HEADING_MODE_GYROHDG) fprintf(fid, "G");
			fprintf(fid, "%f", pLNAVProfile_[i].fparams[0]);		// print heading value
			if (pLNAVProfile_[i].ofparams[0]) fprintf(fid,"*");		// print optimization mark if needed
			fprintf(fid, ", %f", pLNAVProfile_[i].fparams[1]);		// print duration
			if (pLNAVProfile_[i].ofparams[1]) fprintf(fid,"*");		// print optimization mark if needed
			fprintf(fid, ", %f\n", pLNAVProfile_[i].fparams[2]);	// print time step
			break;

		case MANEUVER_TYPE_TURN:
			fprintf(fid, "TRN ");
			fprintf(fid, "%f", pLNAVProfile_[i].fparams[0]);		// print turn angle value value
			if (pLNAVProfile_[i].ofparams[0]) fprintf(fid,"*");		// print optimization mark if needed
			fprintf(fid, ", %f", pLNAVProfile_[i].fparams[1]);		// print bank angle
			if (pLNAVProfile_[i].ofparams[1]) fprintf(fid,"*");		// print optimization mark if needed
			fprintf(fid, ", %f\n", pLNAVProfile_[i].fparams[2]);	// print time step
			break;
		}

	}

	fclose(fid);
	return true;
}