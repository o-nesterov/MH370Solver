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
#include <sys/stat.h>
#include <math.h>
#include <string.h>

#include "params_parser.h"
#include "maneuver.h"
#include "inmarsatdata.h"

extern void uv2uvw(double& XYZ_U, double& XYZ_V, double& XYZ_W, double u, double v, double lon, double lat);


extern void print_time_formatted(char* buf, double t);
extern bool read_time_formatted(double& t, char* buf);
// read value from a string
extern bool parseLine(void* var, char* varname, char* pline, char* vartype);
// scan FP value of varname in line along with attribute character
extern bool parseLineA(double& var, char& attr, char* varname, char* pline);

////////////////////////////////////////////////////////////////////////////////////////
// Contructor: set default values
////////////////////////////////////////////////////////////////////////////////////////
CInpParams::CInpParams()
{
	// Parser version (input parameters list version)
	parserVersion_ = 0;				// default (no file parser was called)


	ts = 70862.906;					// start time (s since 2014-03-07 00:00:00 UTC)
	te = 87569.416;					// end time (s since 2014-03-07 00:00:00 UTC)
	lon0 = 93.323106598760;			// initial longitude (deg E)
	lat0 = 5.997782469407;			// initial latitude (deg N)
	alt0 = 8300.0;					// initial altitude (m)
	wgt0 = 209000.0;				// initial weight (kg)
	hdg0 = 171.775059224434;		// initial heading (deg from N, clockwise) 
	spd0 = 151.478716126326;		// initial speed (IAS, m/s)
	bank_angle = -999.0;			// bank angle - optional parameter used to connect legs if turns are not explicitly specified in lnav file.
	force_default_bank_angle = false; // override turn bank angles specified in lnav file with the default value, unless those are subject to individual optimization (marked with *)
	hdg0_type = -1;					// assume initial heading type the same as heading type

	// gyroscopic heading; assume it parallel to west-to-east direction at KLIA, which is perpendicular to south-to-east direction, and normal 
	double U_gyro, V_gyro, W_gyro;
	uv2uvw(U_gyro, V_gyro, W_gyro, 1.0, 0.0, 101.7058, 2.7521);	// 101.7058E, 2.7521N is the KLIA coordinates; see "inmarsatdata.cpp" for details
	// now solve:
	// U_gyro = sin(gyro_hdg_param2)*cos(gyro_hdg_param1-360.0/86400.0*t);
	// V_gyro = sin(gyro_hdg_param2)*sin(gyro_hdg_param1-360.0/86400.0*t);
	// W_gyro = cos(gyro_hdg_param2); W_gyro == 0 here
	const double Rad = 180.0/3.1415926535897932384626433832795;
	gyro_hdg_param2 = acos(W_gyro)*Rad;
	gyro_hdg_param1 = acos(U_gyro/sqrt(1.0-W_gyro*W_gyro))*Rad;
	if (U_gyro<0.0) gyro_hdg_param1 = 360.0-gyro_hdg_param1;
	gyro_hdg_param1 += 360.0*16.0/24.0; // assume this alighnment was at 16:00:00 utc
	if (gyro_hdg_param1>360.0) gyro_hdg_param1-=360.0;

	bto_bias_rx1200 = -491144.851;	// bto bias for RX1200 channel (Hz)
	bto_bias_rx600 = -495724.681;	// bto bias for RX600 channel (Hz)
	bto_bias_tx = -500700.612;		// bto bias for TX channel (Hz)

	bfo_bias_rx = 153.18;			// bfo bias for RX channel
	bfo_bias_tx = 151.04;			// bfo bias for TX channel
	bfo_bias = 153.73;				// bfo bias for other channels (CX)

	SPD_MODE = HEADING_MODE_TRUEHDG;// speed mode
	HDG_MODE = SPEED_MODE_IAS;		// heading mode
	ENG_MODE = ENGINE_MODE_DUAL;	// engine mode

	METEO_DATASET = METEO_DATASET_ERA5;			// currently supported METEO_DATASET_GDAS1 or METEO_DATASET_ERA5
	MAGN_DECL_DATASET = MAGN_DECL_DATASET_NOAA;	// currently supported only MAGN_DECL_DATASET_NOAA

	MeteoInterpMethod = 3;			// interpolation method for meteorological forcing (0,1,2, or 3)
	MagnDeclMethod = 1;				// interpolation method for magnetic declination (0 or 1)

	pMeteoDir = NULL;				// directory path, where meteorological data are stored
	pMagnDeclDir = NULL;			// directory path, where magnetic declination data are stored

	pLNAVfilename = NULL;			// lateral navigation profile file
	pVNAVfilename = NULL;			// vertical navigation profile file

	pOutTrajFilename = NULL;		// output trajectory filename
	pOutStatFilename = NULL;		// output statistics filename

	pOutputPingIDX = NULL;			// list of indices of pings to be included in statistic file when possible
	nOutputPingIDX = 0;				// number of indices of pings in the list to be included in statistic file

	pOptPingsIDX = NULL;			// list of indices of pings to be used for path optimization 
	nOptPingsIDX = 0;				// number of pings used for path optimization


	DopplerMethod = DOPPLER_METHOD_UDF;		// Doppler computation method (DOPPLER_METHOD_UDF, DOPPLER_METHOD_SIMPLE,  DOPPLER_METHOD_REL,  DOPPLER_METHOD_ADV)
	AES_DopplerMethod = DOPPLER_METHOD_UDF;	// Doppler method used by AES (default - undefined)
	AES_EarthModel = AES_EARTH_MODEL_SPH;	// AES Earth model: SPH (spherical) or WGS (WGS'84 ellipsoid)
	AES_Earth_Radius = -999.0;				// AES Earth radius - optional parameter (affects AES doppler compensation term)
	AES_GeoSat_Alt = -999.0;				// AES geosynch. sat. orbit altitude - optional (affects AES doppler compensation term)
	AES_Delay = 0.1;						// Default delay in position/velocity data used to calculate AES doppler compensation term

	outputstream_ = NULL;			// output diagnostic stream (file); if NULL - no output messages

	pOptPingsIDX = new int[12];
	nOptPingsIDX = 11;
	pOptPingsIDX[0] = 499;
	pOptPingsIDX[1] = 501;
	pOptPingsIDX[2] = 502;
	pOptPingsIDX[3] = 503;
	pOptPingsIDX[4] = 504;
	pOptPingsIDX[5] = 505;
	pOptPingsIDX[6] = 511;
	pOptPingsIDX[7] = 512;
	pOptPingsIDX[8] = 553;
	pOptPingsIDX[9] = 562;
	pOptPingsIDX[10] = 563;


	nLNAVOptParams= 0;		// number of optimization parameters in LNAV profile 
};


////////////////////////////////////////////////////////////////////////////////////////
// Destructor: deallocate all arrays
////////////////////////////////////////////////////////////////////////////////////////
CInpParams::~CInpParams()
{
	if (pMeteoDir) delete pMeteoDir;
	if (pMagnDeclDir) delete pMagnDeclDir;
	if (pLNAVfilename) delete pLNAVfilename;
	if (pVNAVfilename) delete pVNAVfilename;
	if (pOutTrajFilename) delete pOutTrajFilename;
	if (pOutStatFilename) delete pOutStatFilename;
	if (pOutputPingIDX) delete pOutputPingIDX;
	if (pOptPingsIDX) delete pOptPingsIDX;
}


////////////////////////////////////////////////////////////////////////////////////////
//  Verbose (for diagnostic)
////////////////////////////////////////////////////////////////////////////////////////

void CInpParams::SetVerbose(FILE* outputstream)
{
	// set verbose stream (can be NULL - no output)
	outputstream_ = outputstream;
}






/////////////////////////////////////////////////////////////////////////////////
//
// Check parser version. Reserved for future use if list of parameters changes
//
/////////////////////////////////////////////////////////////////////////////////

int CInpParams::checkVersion(char* filename)
{
	// read in the whole file, skipping comment lines
	FILE* fid = fopen(filename,"rt");

	if (fid == NULL) 
	{
		if (outputstream_) fprintf(outputstream_, "Unable to open file: %s.\n", filename);
		return -1;
	}

	char buf[1024];

	// ----------------
	// load the whole file
	while(!feof(fid))
	{
		if (fgets(buf, 1000, fid)!=NULL)
		{
			char* p=&buf[0];
			int n = (int)strlen(p);
			while((*p==' ') || (*p=='	') && p<(&buf[n])) p++;	// skip spaces
			if (p==&buf[n] || n<=3) continue;	// skip empty line
			if (*p=='#' || *p=='!' || *p=='%%') continue; // this is a comment line - skip it

			int version = -1;
			if (parseLine((void*)(&version), "VERSION", buf,"int"))
			{
				parserVersion_ = version;
				fclose(fid);
				return version;
			}
		}
	}
	fclose(fid);
	if (outputstream_) fprintf(outputstream_, "Parser version is not specified.\n");
	return -1;
}



/////////////////////////////////////////////////////////////////////////////////
//
// Parse parameters
//
/////////////////////////////////////////////////////////////////////////////////


bool CInpParams::parseParams(char* filename)
{
	if (checkVersion(filename) == -1) 
	{
		if (outputstream_) fprintf(outputstream_, "Unable to check version or file does not exist.\n");
		return false; // check parser version
	}

	if (pMeteoDir) delete pMeteoDir;
	pMeteoDir = NULL;

	if (pMagnDeclDir) delete pMagnDeclDir;
	pMagnDeclDir = NULL;

	if (pLNAVfilename) delete pLNAVfilename;
	pLNAVfilename = NULL;

	if (pVNAVfilename) delete pVNAVfilename;
	pVNAVfilename = NULL;

	if (pOutTrajFilename) delete pOutTrajFilename;
	pOutTrajFilename = NULL;

	if (pOutStatFilename) delete pOutStatFilename;
	pOutStatFilename = NULL;

	if (pOutputPingIDX) delete pOutputPingIDX;
	pOutputPingIDX = NULL;
	nOutputPingIDX = 0;

	if (pOptPingsIDX) delete pOptPingsIDX;
	pOptPingsIDX = NULL;
	nOptPingsIDX = 0;

	// read in the whole file, skipping comment lines
	FILE* fid = fopen(filename,"rt");

	if (fid == NULL) return false;

	char buf[1024];

	int nmaxLines = 256;
	int nLines = 0;
	char** ppLines = new char*[nmaxLines];

	// ----------------
	// load the whole file
	while(!feof(fid))
	{
		if (fgets(buf, 1000, fid)!=NULL)
		{
			char* p=&buf[0];
			int n = (int)strlen(p);
			while((*p==' ') || (*p=='	') && p<(&buf[n])) p++;	// skip spaces
			if (p==&buf[n] || n<=3) continue;	// skip empty line
			if (*p=='#' || *p=='!' || *p=='%%') continue; // this is a comment line - skip it

			if (nLines+1 >= nmaxLines)
			{
				// increase array of string pointers
				nmaxLines += 256;
				char** ppLines_new = new char*[nmaxLines];
				for (int j=0; j<nLines; j++) ppLines_new[j] = ppLines[j];	// copy to new array
				delete ppLines;
				ppLines = ppLines_new;
			}

			ppLines[nLines] = new char[n+2];
			strcpy(ppLines[nLines], buf);
			nLines++;
		}
	}
	fclose(fid);

	// ------------------
	// now parse lines
	// FP variable names to read
	char* fvarnames[] = {"T_START", "T_END", "LON", "LAT", "ALT", "WGT", "HDG", "TRK", "SPD", "MACH", "BTO_BIAS_RX1200", "BTO_BIAS_RX600", "BTO_BIAS_TX", "BFO_BIAS_RX", "BFO_BIAS_TX", "BFO_BIAS", "GYRO_HDG_PARAM1", "GYRO_HDG_PARAM2", "BANK_ANGLE", "AES_EARTH_RADIUS", "AES_GEOSAT_ALT", "AES_DELAY"};
	double fvars[sizeof(fvarnames)/sizeof(char*)+1];
	bool isfvars[sizeof(fvarnames)/sizeof(char*)+1];	// flag indicating if variable was read
	bool isFPAttribute[sizeof(fvarnames)/sizeof(char*)+1] = {false, false, false, false, false, false, true, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false};	// read with or without attribute letter
	char cfvars[sizeof(fvarnames)/sizeof(char*)+1];		// optional attribute letters (used for heading type; can be used for units)

	// integer variable names to read
	char* ivarnames[] = {"METEO_INTERP_METHOD", "MAGNETIC_DECL_INTERP_METHOD"};
	int ivars[sizeof(ivarnames)/sizeof(char*)+1];
	bool isivars[sizeof(ivarnames)/sizeof(char*)+1];	// flag indicating if variable was read

	// string variable names to read
	char* svarnames[] = {"METEO_DATASET", "METEO_DIR", "MAGNETIC_DECLINATION_DATASET", "MAGNETIC_DECLINATION_DIR", "LNAV_FILE", "VNAV_FILE", "SPEED_MODE", "HDG_MODE", "ENGINE_MODE", "OUTPUT_TRAJ_FILE", "OUTPUT_BTO_BFO_FILE", "OUTPUT_BTO_BFO_IDX",  "DOPPLER_METHOD", "AES_DOPPLER_METHOD", "AES_EARTH_MODEL", "OPT_PING_IDX"};
	char* svars[sizeof(svarnames)/sizeof(char*)+1];

	// floating point variables
	int n = sizeof(fvarnames)/sizeof(char*);
	for (int i=0; i<n; i++)
	{
		isfvars[i] = false;
		int j =0;
		while (j<nLines)
		{
			if (strcmp(fvarnames[i],"T_START")==0 || strcmp(fvarnames[i],"T_END")==0)
			{
				if (parseLine((void*)(&fvars[i]), fvarnames[i], ppLines[j],"time")) {isfvars[i] = true; break;}; // try to read time as floating point or human-readable text
			}
			else
			{
				if (isFPAttribute[i])
				{
					// read FP accompanied by attribute character (either prefix or suffix)
					if (parseLineA(fvars[i], cfvars[i], fvarnames[i], ppLines[j])) {isfvars[i] = true; break;}; // variable was successfully read
				}
				else
				{
					// read FP variable
					if (parseLine((void*)(&fvars[i]), fvarnames[i], ppLines[j],"double")) {isfvars[i] = true; break;}; // variable was successfully read
				}
			}
			j++;
		}
	};


	// integer variables
	n = sizeof(ivarnames)/sizeof(char*);
	for (int i=0; i<n; i++)
	{
		isivars[i] = false;
		int j =0;
		while (j<nLines)
		{
			if (parseLine((void*)(&ivars[i]), ivarnames[i], ppLines[j],"int")) {isivars[i]= true; break;}; // variable was successfully read
			j++;
		}
	};

	// string variables
	n = sizeof(svarnames)/sizeof(char*);
	for (int i=0; i<n; i++)
	{
		svars[i] = NULL;
		int j = 0;
		while (j<nLines)
		{
			if (parseLine((void*)(&svars[i]), svarnames[i], ppLines[j],"string")) break; // variable was successfully read
			j++;
		}
	};

	// ------------------------------------
	// Now we can deallocate input file strings
	if (ppLines)
	{
		for (int j=0; j<nLines; j++) if (ppLines[j]) delete ppLines[j];
		if (ppLines) delete ppLines;
	};

	// -----------------------------------
	// Now check variables and assign return values

	bool iOkFlag = true;

	if (isfvars[0])
	{
		ts = fvars[0];				// start time (s since 2014-03-07 00:00:00 UTC)
	}
	else
	{
		if (outputstream_) fprintf(outputstream_, "Start time must be specified.\n");
		iOkFlag = false;
	}

	//--------
	if (isfvars[1])
	{
		te = fvars[1];				// end time (s since 2014-03-07 00:00:00 UTC)
	}
	else
	{
		if (outputstream_) fprintf(outputstream_, "End time must be specified.\n");
		iOkFlag = false;
	}

	if (isfvars[0] && isfvars[1])
	{
		if (ts>=te)
		{
			if (outputstream_) fprintf(outputstream_, "Error: start time >= end time.\n");
			iOkFlag = false;
		}

		if (ts+1.0>=te)
		{
			if (outputstream_) fprintf(outputstream_, "Warning: duration is less than 1s.\n");
		}

		if (te-ts > 43200.0)
		{
			if (outputstream_) fprintf(outputstream_, "Warning: duration exceeds 12h. Wrong time?\n");
		}

		if (ts<43200.0)
		{
			if (outputstream_) fprintf(outputstream_, "Warning: start time before 2014-03-07 12:00:00?\n");
		}

		if (te>93600.0)
		{
			if (outputstream_) fprintf(outputstream_, "Warning: end time after 2014-03-08 02:00:00?\n");
		}
	}

	//--------
	if (isfvars[2])
	{
		lon0 = fvars[2];	// initial longitude (deg E)
	}
	else
	{
		if (outputstream_) fprintf(outputstream_, "Initial longitude must be specified.\n");
		iOkFlag = false;
	}

	//--------
	if (isfvars[3])
	{
		lat0 = fvars[3];	// initial latitude (deg N)
	}
	else
	{
		if (outputstream_) fprintf(outputstream_, "Initial latitude must be specified.\n");
		iOkFlag = false;
	}

	//--------
	if (isfvars[4])
	{
		alt0 = fvars[4];	//  initial altitude (m)
	}
	else
	{
		if (outputstream_) fprintf(outputstream_, "Initial altitude must be specified.\n");
		iOkFlag = false;
	}

	//--------
	if (isfvars[5])
	{
		wgt0 = fvars[5];	//  initial weight (m)
	}
	else
	{
		if (outputstream_) fprintf(outputstream_, "Initial weight must be specified.\n");
		iOkFlag = false;
	}

	//--------
	// internal heading and track - both hdg and trk can be specified in one setup file, but then depending on specified mode, hd0 = either hdg0__ or trk0__
	double hdg0__;
	double trk0__;
	if (isfvars[6] || isfvars[7])
	{
		if (isfvars[6])
		{
			hdg0__ = fvars[6];	// initial heading (deg from N, clockwise), true or magnetic, depending on mode 

			if (abs(hdg0__)>1.0E+6)
			{
				// make this failsave check 
				if (outputstream_) fprintf(outputstream_, "Initial heading or track is weird (try to make it between 0 and 360 deg).\n");
				iOkFlag = false;
			}
			else
			{
				// simple method to reduce heading to the range [0,360) deg.
				while (hdg0__<0.0) hdg0__ += 360.0; // make sure that 360.0>hdg0__>=0.0
				while (hdg0__>=360.0) hdg0__ -=360.0;
			}

			// heading value type: true, magnetic or gyroscopic
			switch(cfvars[6])
			{
			case 'T':
				hdg0_type= HEADING_MODE_TRUEHDG;
				break;

			case 'M':
				hdg0_type= HEADING_MODE_MAGNHDG;
				break;

			case 'G':
				hdg0_type= HEADING_MODE_GYROHDG;
				break;
			}
		}

		if (isfvars[7])
		{
			trk0__ = fvars[7];	// initial track (deg from N, clockwise), true or magnetic, depending on mode 

			if (abs(trk0__)>1.0E+6)
			{
				// make this failsave check 
				if (outputstream_) fprintf(outputstream_, "Initial heading or track is weird (try to make it between 0 and 360 deg).\n");
				iOkFlag = false;
			}
			else
			{
				// simple method to reduce heading to the range [0,360) deg.
				while (trk0__<0.0) trk0__ += 360.0; // make sure that 360.0>trk0__>=0.0
				while (trk0__>=360.0)trk0__ -=360.0;
			}

			// heading value type: true, magnetic or gyroscopic
			switch(cfvars[6])
			{
			case 'T':
				hdg0_type= HEADING_MODE_TRUEHDG;
				break;

			case 'M':
				hdg0_type= HEADING_MODE_MAGNHDG;
				break;

			case 'G':
				hdg0_type= HEADING_MODE_GYROHDG;
				break;
			}
		}

	}
	else
	{
		if (outputstream_) fprintf(outputstream_, "Initial heading or track must be specified.\n");
		iOkFlag = false;
	}

	// ------------------
	// internal speed and mach - both speed and mach values can be specified in one setup file, but then only one used depending on mode
	double spd0__, mach0__;
	if (isfvars[8] || isfvars[9])
	{
		if (isfvars[8])
		{
			spd0__ = fvars[8];	// initial airspeed (m/s), depending on mode 
			if (spd0__ > AIRCRAFT_VFC*0.514444)
			{
				if (outputstream_) fprintf(outputstream_, "Too high IAS. Resetting to VFC.\n");
				spd0__ = AIRCRAFT_VFC*0.514444;
			}

			if (spd0__ < 50.0)
			{
				if (outputstream_) fprintf(outputstream_, "Too low IAS. Resetting to 50 m/s.\n");
				spd0__ = 50.0;
			}
		}
		else
		{
			mach0__ = fvars[9];	// initial airspeed (mach), depending on mode 
			if (mach0__ > AIRCRAFT_MFC)
			{
				if (outputstream_) fprintf(outputstream_, "Too high Mach. Resetting to MFC\n");
				mach0__ = AIRCRAFT_MFC;
			}
			if (mach0__ < 0.2)
			{
				if (outputstream_) fprintf(outputstream_, "Too low Mach. Resetting to 0.2.\n");
				mach0__ = 0.2;
			}
		}
	}
	else
	{
		if (outputstream_) fprintf(outputstream_, "Initial airspeed (IAS or MACH) must be specified.\n");
		iOkFlag = false;
	}


	//--------
	if (isfvars[10])
	{
		bto_bias_rx1200 = fvars[10];	//  bto bias for rx1200 channel (Hz)
	}
	else
	{
		if (outputstream_) fprintf(outputstream_, "Warning: no bto bias for RX1200 channel was specified. Using default value.\n");
	}


	//--------
	if (isfvars[11])
	{
		bto_bias_rx600 = fvars[11];	//  bto bias for rx600 channel (microseconds)
	}
	else
	{
		if (outputstream_) fprintf(outputstream_, "Warning: no bto bias for RX600 channel was specified. Using default value.\n");
	}

	//--------
	if (isfvars[12])
	{
		bto_bias_tx = fvars[12];	//  bto bias for tx channel (microseconds)
	}
	else
	{
		if (outputstream_) fprintf(outputstream_, "Warning: no bto bias for TX channel was specified. Using default value.\n");
	}

	//--------
	if (isfvars[13])
	{
		bfo_bias_rx = fvars[13];	//  bto bias for rx channels (microseconds)
	}
	else
	{
		if (outputstream_) fprintf(outputstream_, "Warning: no bfo bias for RX channel was specified. Using default value.\n");
	}

	//--------
	if (isfvars[14])
	{
		bfo_bias_tx = fvars[14];	//  bto bias for tx channels (Hz)
	}
	else
	{
		if (outputstream_) fprintf(outputstream_, "Warning: no bfo bias for TX channel was specified. Using default value.\n");
	}

	//--------
	if (isfvars[15])
	{
		bfo_bias = fvars[15];	//  bto bias for other than RX and TX channels (Hz)
	}
	else
	{
		if (outputstream_) fprintf(outputstream_, "Warning: no bfo bias for other than RX and TX channels was specified. Using default value.\n");
	}


	//--------
	if (isfvars[16])
	{
		gyro_hdg_param1 = fvars[16];	//  gyroscopic heading parameter 1
	}

	//--------
	if (isfvars[17])
	{
		gyro_hdg_param2 = fvars[17];	//  gyroscopic heading parameter 2
	}

	//--------
	if (isfvars[18])
	{
		bank_angle = fvars[18];			//  optional - bank angle (deg); if not specified, default (25 deg) is used to connect legs of a trajectory
	}

	//--------
	if (isfvars[19])
	{
		AES_Earth_Radius = fvars[19];	//  optional - Earth radius (m) assumed by AES when spherical model is used to compute Doppler compensation term 
	}

	//--------
	if (isfvars[20])
	{
		AES_GeoSat_Alt = fvars[20];		//  optional - Geosynch. satellite altitude (m) assumed by AES to compute Doppler compensation term
	}

	//--------
	if (isfvars[21])
	{
		AES_Delay = fvars[21];			//  AES delay in position / velocity data to compute Doppler compensation term (can be important during turns); default zero
	}

	// --------------------------------------------------------------------------------------------------------------------------------------

	// integer parameters {"METEO_INTERP_METHOD", "MAGNETIC_DECL_INTERP_METHOD"};
	//--------
	if (isivars[0])
	{
		MeteoInterpMethod = ivars[0];

		if (ivars[0] < 0)
		{
			MeteoInterpMethod = 0;
			if (outputstream_) fprintf(outputstream_, "Interpolation method for meteorology must be 0, 1, 2 or 3. Using default value 0.\n");
		}

		if (ivars[0] > 3)
		{
			MeteoInterpMethod = 3;
			if (outputstream_) fprintf(outputstream_, "Interpolation method for meteorology must be 0, 1, 2 or 3. Using default value 3.\n");
		}
	}
	else
	{
		if (outputstream_) fprintf(outputstream_, "Warning: no interpolation method for meteorology was specified. Using default.\n");
	}


	//--------
	if (isivars[1])
	{
		MagnDeclMethod = ivars[1];

		if (ivars[1] < 0)
		{
			MagnDeclMethod = 0;
			if (outputstream_) fprintf(outputstream_, "Interpolation method for magnetic declination must be 0 or 1. Using default value 0.\n");
		}

		if (ivars[1] > 1)
		{
			MeteoInterpMethod = 1;
			if (outputstream_) fprintf(outputstream_, "Interpolation method for magnetic declination must be 0 or 1. Using default value 1.\n");
		}
	}
	else
	{
		if (outputstream_) fprintf(outputstream_, "Warning: no interpolation method for magnetic declination was specified. Using default.\n");
	}

	//---------------------------------------------------------------------------------------------------------------
	// {"METEO_DATASET", "METEO_DIR", "MAGNETIC_DECLINATION_DATASET", "MAGNETIC_DECLINATION_DIR", "LNAV_FILE", "VNAV_FILE", "SPEED_MODE", "HDG_MODE", "ENGINE_MODE", "OUTPUT_TRAJ_FILE", "OUTPUT_BTO_BFO_FILE", "OUTPUT_BTO_BFO_IDX",  "DOPPLER_METHOD", "AES_DOPPLER_METHOD", "AES_EARTH_MODEL"};

	if (svars[0])
	{
		bool isMeteoName = false;
		if (strcmp(svars[0],"NULL") == 0 || strcmp(svars[0],"null") == 0 || strcmp(svars[0],"Null") == 0 || strcmp(svars[0],"NIL") == 0 || strcmp(svars[0],"Nil") == 0 || strcmp(svars[0],"nil") == 0) {METEO_DATASET = METEO_DATASET_NULL; isMeteoName = true;};
		if (strcmp(svars[0],"GDAS1") == 0 || strcmp(svars[0],"gdas1") == 0 || strcmp(svars[0],"Gdas1") == 0) {METEO_DATASET = METEO_DATASET_GDAS1; isMeteoName = true;};
		if (strcmp(svars[0],"ERA5") == 0 || strcmp(svars[0],"era5") == 0 || strcmp(svars[0],"Era5") == 0) {METEO_DATASET = METEO_DATASET_ERA5; isMeteoName = true;};
		if (!isMeteoName)
		{
			METEO_DATASET = METEO_DATASET_NULL;
			if (outputstream_) fprintf(outputstream_, "Warning: unrecognized meteorological dataset. Using default atmospheric conditions.\n");
		}
	}
	else
	{
		METEO_DATASET = METEO_DATASET_NULL;
		if (outputstream_) fprintf(outputstream_, "Warning: meteorological dataset was not specified. Using default atmospheric conditions.\n");
	}

	// -----
	if (svars[1])
	{
		// create new buffer and copy to avoid potential memory leaks
		pMeteoDir = new char[strlen(svars[1])+1];
		strcpy(pMeteoDir, svars[1]);

		// check if this is valid directory
		bool isDir = false;

		struct _stat info;

		if(_stat(pMeteoDir, &info ) != 0)
		{
			isDir= false;
		}
		else if(info.st_mode & S_IFDIR)
		{
			isDir= true;
		}

		if (!isDir)
		{
			if (outputstream_) fprintf(outputstream_, "Erroneous directory containing meteorological dataset.\n");
			iOkFlag = false;
		}
	}
	else
	{
		if (METEO_DATASET == METEO_DATASET_GDAS1 || METEO_DATASET == METEO_DATASET_ERA5)
		{
			if (outputstream_) fprintf(outputstream_, "Directory containing meteorological dataset must be specified.\n");
			iOkFlag = false;
		}
	}

	// ------
	if (svars[2])
	{
		bool isMagnDeclName = false;
		if (strcmp(svars[2],"NULL") == 0 || strcmp(svars[2],"null") == 0 || strcmp(svars[2],"Null") == 0 || strcmp(svars[2],"NIL") == 0 || strcmp(svars[2],"Nil") == 0 || strcmp(svars[2],"nil") == 0) {MAGN_DECL_DATASET = MAGN_DECL_DATASET_NULL; isMagnDeclName = true;};
		if (strcmp(svars[2],"NOAA") == 0 || strcmp(svars[2],"noaa") == 0 || strcmp(svars[2],"Noaa") == 0) {MAGN_DECL_DATASET = MAGN_DECL_DATASET_NOAA; isMagnDeclName = true;};
		if (!isMagnDeclName)
		{
			MAGN_DECL_DATASET = MAGN_DECL_DATASET_NULL;
			if (outputstream_) fprintf(outputstream_, "Warning: unrecognized magnetic declination dataset. Magnetic declination is set to zero.\n");
		}
	}
	else
	{
		MAGN_DECL_DATASET = MAGN_DECL_DATASET_NULL;
		if (outputstream_) fprintf(outputstream_, "Warning: magnetic declination dataset was not specified. Magnetic declination is set to zero.\n");
	}

	// -----
	if (svars[3])
	{
		// create new buffer and copy to avoid potential memory leaks
		pMagnDeclDir = new char[strlen(svars[3])+1];
		strcpy(pMagnDeclDir, svars[3]);

		// check if this is valid directory
		bool isDir = false;

		struct _stat info;

		if(_stat(pMagnDeclDir, &info ) != 0)
		{
			isDir= false;
		}
		else if(info.st_mode & S_IFDIR)
		{
			isDir= true;
		}

		if (!isDir && (MAGN_DECL_DATASET == MAGN_DECL_DATASET_NOAA))
		{
			if (outputstream_) fprintf(outputstream_, "Erroneous directory containing magnetic declination dataset.\n");
			iOkFlag = false;
		}
	}
	else
	{
		if (MAGN_DECL_DATASET == MAGN_DECL_DATASET_NOAA)
		{
			if (outputstream_) fprintf(outputstream_, "Directory containing magnetic declination dataset must be specified.\n");
			iOkFlag = false;
		}
	}

	// -----
	if (svars[4])
	{
		// create new buffer and copy to avoid potential memory leaks
		if (strlen(svars[4])>0)
		{
			pLNAVfilename = new char[strlen(svars[4])+1];
			strcpy(pLNAVfilename, svars[4]);
		}
	}
	else
	{
		if (outputstream_) fprintf(outputstream_, "Lateral navigation file was not specified. Assume a leg (constant heading or track).\n");
	}

	// -----
	if (svars[5])
	{
		// create new buffer and copy to avoid potential memory leaks
		pVNAVfilename = new char[strlen(svars[5])+1];
		strcpy(pVNAVfilename, svars[5]);
	}
	else
	{
		if (outputstream_) fprintf(outputstream_, "Vertical navigation file was not specified. Assume constant flight level.\n");
	}

	// -----
	if (svars[6])
	{
		bool isSPD = false;
		if (strcmp(svars[6], "IAS") == 0) {SPD_MODE = SPEED_MODE_IAS; isSPD = true;};
		if (strcmp(svars[6], "MACH") == 0) {SPD_MODE = SPEED_MODE_MACH;isSPD = true;};
		if (!isSPD)
		{
			if (isfvars[8])
			{
				if (outputstream_) fprintf(outputstream_, "Unrecognized airspeed mode. Presumed IAS.\n");
			}
			else
			{
				if (isfvars[9])
				{
					if (outputstream_) fprintf(outputstream_, "Unrecognized airspeed mode. Presumed MACH.\n");
				}
				else
				{
					if (outputstream_) fprintf(outputstream_, "Airspeed mode must be specified.\n");
					iOkFlag = false;
				}
			}
		}
	}
	else
	{
		if (isfvars[8])
		{
			SPD_MODE = SPEED_MODE_IAS;
			if (outputstream_) fprintf(outputstream_, "Airspeed mode was not specified, but inferred to be IAS.\n");
		}
		else
		{
			if (isfvars[9])
			{
				SPD_MODE = SPEED_MODE_MACH;
				if (outputstream_) fprintf(outputstream_, "Airspeed mode was not specified, but inferred to be MACH.\n");
			}
			else
			{
				if (outputstream_) fprintf(outputstream_, "Airspeed mode and value must be specified.\n");
				iOkFlag = false;
			}
		}
	}

	if (SPD_MODE == SPEED_MODE_IAS && (!isfvars[8]) && (isfvars[9]))
	{
		if (outputstream_) fprintf(outputstream_, "Error: mismatch of specified airspeed mode (IAS) and value (Mach).\n");
		iOkFlag = false;
	}

	if (SPD_MODE == SPEED_MODE_MACH && (isfvars[8]) && (!isfvars[9]))
	{
		if (outputstream_) fprintf(outputstream_, "Error: mismatch of specified airspeed mode (Mach) and value (IAS).\n");
		iOkFlag = false;
	}

	if (SPD_MODE == SPEED_MODE_IAS) spd0 = spd0__;
	if (SPD_MODE == SPEED_MODE_MACH) spd0 = mach0__;


	// -----
	if (svars[7])
	{
		bool isHDG = false;
		if (strcmp(svars[7], "TRUE_HDG") == 0) {HDG_MODE = HEADING_MODE_TRUEHDG; isHDG = true;};
		if (strcmp(svars[7], "MAGN_HDG") == 0) {HDG_MODE = HEADING_MODE_MAGNHDG; isHDG = true;};
		if (strcmp(svars[7], "MAG_HDG") == 0)  {HDG_MODE = HEADING_MODE_MAGNHDG; isHDG = true;};
		if (strcmp(svars[7], "TRUE_TRK") == 0) {HDG_MODE = HEADING_MODE_TRUETRK; isHDG = true;};
		if (strcmp(svars[7], "MAGN_TRK") == 0) {HDG_MODE = HEADING_MODE_MAGNTRK; isHDG = true;};
		if (strcmp(svars[7], "MAG_TRK") == 0)  {HDG_MODE = HEADING_MODE_MAGNTRK; isHDG = true;};
		if (strcmp(svars[7], "GYRO_HDG") == 0) {HDG_MODE = HEADING_MODE_GYROHDG; isHDG = true;};

		if (!isHDG)
		{
			if (outputstream_) fprintf(outputstream_, "Heading mode %s is not recognized. Must be one of the following: TRUE_HDG, MAGN_HDG, TRUE_TRK, MAGN_TRK, GYRO_HDG.\n", svars[6]);
			iOkFlag = false;
		}
	}
	else
	{
		if (isfvars[6])
		{
			HDG_MODE = HEADING_MODE_TRUEHDG;
			if (outputstream_) fprintf(outputstream_, "Heading mode was not specified, but presumed TRUE_HDG.\n");
		}
		else
		{
			if (isfvars[7])
			{
				HDG_MODE = HEADING_MODE_TRUETRK;
				if (outputstream_) fprintf(outputstream_, "Heading mode was not specified, but presumed TRUE_TRK.\n");
			}
			else
			{
				if (outputstream_) fprintf(outputstream_, "Heading mode must be specified (TRUE_HDG, MAGN_HDG, TRUE_TRK, MAGN_TRK, GYRO_HDG).\n");
				iOkFlag = false;
			}
		}
	}

	if ((HDG_MODE == HEADING_MODE_TRUEHDG || HDG_MODE == HEADING_MODE_MAGNHDG) && (!isfvars[6]) && (isfvars[7]))
	{
		if (outputstream_) fprintf(outputstream_, "Warning: mismatch of heading mode and specified heading value. Specified mode is prioritized.\n");
	}

	if ((HDG_MODE == HEADING_MODE_TRUETRK || HDG_MODE == HEADING_MODE_MAGNTRK) && (isfvars[6]) && (!isfvars[7]))
	{
		if (outputstream_) fprintf(outputstream_, "Warning: mismatch of track mode and specified track value. Specified mode is prioritized.\n");
	}

	if (HDG_MODE == HEADING_MODE_GYROHDG)
	{
		if (!isfvars[16])
		{
			if (outputstream_) fprintf(outputstream_, "Warning: gyroscopic heading was specified, but paramer 1 absent. Using default value.\n");
		}

		if (!isfvars[17])
		{
			if (outputstream_) fprintf(outputstream_, "Warning: gyroscopic heading was specified, but paramer 2 absent. Using default value.\n");
		}
	}

	if (HDG_MODE == HEADING_MODE_TRUEHDG || HDG_MODE == HEADING_MODE_MAGNHDG || HDG_MODE == HEADING_MODE_GYROHDG) hdg0 = hdg0__;
	if (HDG_MODE == HEADING_MODE_TRUETRK || HDG_MODE == HEADING_MODE_MAGNTRK) hdg0 = trk0__;


	// ------------------
	if (svars[8])
	{
		bool isEngMode = false;
		if (strcmp(svars[8],"BOTH") == 0 || strcmp(svars[8],"Both") == 0 || strcmp(svars[8],"both") == 0 || strcmp(svars[8],"DUAL") == 0 || strcmp(svars[8],"Dual") == 0 || strcmp(svars[8],"dual") == 0) {ENG_MODE = ENGINE_MODE_DUAL; isEngMode = true;};
		if (strcmp(svars[8],"SINGLE") == 0 || strcmp(svars[8],"Single") == 0 || strcmp(svars[8],"single") == 0) {ENG_MODE = ENGINE_MODE_SINGLE; isEngMode = true;};
		if (!isEngMode)
		{
			ENG_MODE = ENGINE_MODE_DUAL;
			if (outputstream_) fprintf(outputstream_, "Unrecognized engine mode. Using default mode dual.\n");
		}
	}
	else
	{
			ENG_MODE = ENGINE_MODE_DUAL;
			if (outputstream_) fprintf(outputstream_, "Unspecified engine mode. Using default mode dual.\n");
	}

	// -----
	if (svars[9])
	{
		// create new buffer and copy to avoid potential memory leaks
		pOutTrajFilename = new char[strlen(svars[9])+1];
		strcpy(pOutTrajFilename, svars[9]);
	}
	else
	{
		pOutTrajFilename = NULL;
		if (outputstream_) fprintf(outputstream_, "Warning: output trajectory filename is not specified.\n");
	}


	// -----
	if (svars[10])
	{
		// create new buffer and copy to avoid potential memory leaks
		pOutStatFilename = new char[strlen(svars[10])+1];
		strcpy(pOutStatFilename, svars[10]);
	}
	else
	{
		pOutStatFilename = NULL;
		if (outputstream_) fprintf(outputstream_, "Warning: output statistics filename is not specified.\n");
	}

	// -----
	// List of specific output indices to compare BTOs and BFOs in statistics file
	if (svars[11])
	{
		if (strcmp(svars[11],"ALL") == 0 || strcmp(svars[11],"all") == 0 || strcmp(svars[11],"All") == 0)
		{
			nOutputPingIDX = -1;
		}
		else
		{
			// analize string and convert list of comma-separated values to the list of integers
			char* p0 = svars[11];
			char* p = p0;
			int n = strlen(p);

			while(((*p==' ') || (*p=='	')) && p<(p0+n)) p++;	// skip spaces

			while(p<(p0+n))
			{
				while(((*p==' ') || (*p=='	') || (*p==',') || (*p==';')) && (*p!='}') && p<(p0+n)) p++;	// skip spaces, commas, etc.
				if ((*p=='}') || p==(p0+n)) break;

				char* pp;
				int tmp = strtol(p,&pp,10);
				if (p!=pp)
				{
					p = pp+1;
					// add index to the list
					int* pOutputPingIDXNew = new int[nOutputPingIDX+1];
					if (pOutputPingIDX)
					{
						for (int j=0; j<nOutputPingIDX; j++) pOutputPingIDXNew[j] = pOutputPingIDX[j];
						delete pOutputPingIDX;
					}
					pOutputPingIDX = pOutputPingIDXNew;
					pOutputPingIDX[nOutputPingIDX] = tmp;
					nOutputPingIDX++;
				}
				else
				{
					break;
				}
			}
		}
	}
	else
	{
		pOutputPingIDX = NULL;	// list of indices of pings to be included in statistic file when possible
		nOutputPingIDX = 0;		// number of indices of pings in the list to be included in statistic file
	}


	// -----
	if (svars[12])
	{
		// Doppler method (optional) - default undefined
		if (strcmp(svars[12],"SIM")==0 ||strcmp(svars[12],"Simple")==0 || strcmp(svars[12],"SIMPLE")==0) DopplerMethod = DOPPLER_METHOD_SIMPLE;
		if (strcmp(svars[12],"REL")==0 || strcmp(svars[12],"Relativistic")==0 || strcmp(svars[12],"relativistic")==0) DopplerMethod = DOPPLER_METHOD_REL;
		if (strcmp(svars[12],"ADV")==0 || strcmp(svars[12],"Advanced")==0 || strcmp(svars[12],"advanced")==0) DopplerMethod = DOPPLER_METHOD_ADV;
	}

	// -----
	if (svars[13])
	{
		// AES Doppler method (optional) - default undefined
		if (strcmp(svars[13],"SIM")==0 ||strcmp(svars[13],"Simple")==0 || strcmp(svars[13],"SIMPLE")==0) AES_DopplerMethod = DOPPLER_METHOD_SIMPLE;
		if (strcmp(svars[13],"REL")==0 || strcmp(svars[13],"Relativistic")==0 || strcmp(svars[13],"relativistic")==0) AES_DopplerMethod = DOPPLER_METHOD_REL;
		if (strcmp(svars[13],"ADV")==0 || strcmp(svars[13],"Advanced")==0 || strcmp(svars[13],"advanced")==0) AES_DopplerMethod = DOPPLER_METHOD_ADV;
	}

	// -----
	if (svars[14])
	{
		// AES Earth shape (optional) - spherical or WGS'84 ellipsoid
		if (strcmp(svars[14],"SPH")==0 ||strcmp(svars[14],"Spherical")==0 || strcmp(svars[14],"spherical")==0) AES_EarthModel = AES_EARTH_MODEL_SPH;
		if (strcmp(svars[14],"WGS")==0 || strcmp(svars[14],"Ellipsoid")==0 || strcmp(svars[14],"ellipsoid")==0) AES_EarthModel = AES_EARTH_MODEL_WGS;
	}


	// -----
	// List of ping indices for optimization
	if (svars[15])
	{
		// analize string and convert list of comma-separated values to the list of integers
		char* p0 = svars[15];
		char* p = p0;
		int n = strlen(p);

		while(((*p==' ') || (*p=='	')) && p<(p0+n)) p++;	// skip spaces

		while(p<(p0+n))
		{
			while(((*p==' ') || (*p=='	') || (*p==',') || (*p==';')) && (*p!='}') && p<(p0+n)) p++;	// skip spaces, commas, etc.
			if ((*p=='}') || p==(p0+n)) break;

			char* pp;
			int tmp = strtol(p,&pp,10);
			if (p!=pp)
			{
				p = pp+1;
				// add index to the list
				int* pOptPingsIDXNew = new int[nOptPingsIDX+1];
				if (pOptPingsIDX)
				{
					for (int j=0; j<nOptPingsIDX; j++) pOptPingsIDXNew[j] = pOptPingsIDX[j];
					delete pOptPingsIDX;
				}
				pOptPingsIDX = pOptPingsIDXNew;
				pOptPingsIDX[nOptPingsIDX] = tmp;
				nOptPingsIDX++;
			}
			else
			{
				break;
			}
		}
	}
	else
	{
		pOptPingsIDX = NULL;	// list of indices of pings to be used for optimization
		nOptPingsIDX = 0;		// number of indices of pings to be used for optimization
	}

	// -----------------------------------
	// Deallocate local string varibles
	n = sizeof(svarnames)/sizeof(char*);
	for (int i=0; i<n; i++) if (svars[i]) delete svars[i];

	// check specified ping IDs for optimization and output
	checkPingIds_();

	return iOkFlag;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Save parameters
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool CInpParams::saveParams(char* filename)
{
	if (parserVersion_<1)
	{
		if (outputstream_) fprintf(outputstream_, "Unknown parser version.\n");
		return false;
	}

	FILE* fid = fopen(filename, "wt");
	if (fid == NULL) 
	{
		if (outputstream_) fprintf(outputstream_, "Unable to save file.\n");
		return false;
	}

	fprintf(fid,"##################################################################\n");
	fprintf(fid,"# Main parameters file\n");
	fprintf(fid,"##################################################################\n");
	fprintf(fid,"VERSION = %d\n", parserVersion_);

	// ---------
	fprintf(fid," \n");
	fprintf(fid,"# 4D Meteorological dataset:\n");
	switch (METEO_DATASET)
	{
	case METEO_DATASET_GDAS1:
		fprintf(fid,"METEO_DATASET = GDAS1\n");
		break;
	case METEO_DATASET_ERA5:
		fprintf(fid,"METEO_DATASET = ERA5\n");
		break;
	}
	if (pMeteoDir) fprintf(fid,"METEO_DIR = \"%s\"\n", pMeteoDir);
	fprintf(fid,"METEO_INTERP_METHOD = %d \n", MeteoInterpMethod);

	// ---------
	fprintf(fid," \n");
	fprintf(fid,"# 2D Magnetic declination data:\n");
	switch (MAGN_DECL_DATASET)
	{
	case MAGN_DECL_DATASET_NOAA:
		fprintf(fid,"MAGNETIC_DECLINATION_DATASET = NOAA\n");
		break;
	}
	if (pMagnDeclDir) fprintf(fid,"MAGNETIC_DECLINATION_DIR = \"%s\"\n", pMagnDeclDir);
	fprintf(fid,"MAGNETIC_DECL_INTERP_METHOD = %d \n", MagnDeclMethod);

	// ---------
	if (pLNAVfilename)
	{
		fprintf(fid," \n");
		fprintf(fid,"# Lateral navigation file\n");
		fprintf(fid,"LNAV_FILE = \"%s\"\n", pLNAVfilename);
	}
	else
	{
		fprintf(fid," \n");
		fprintf(fid,"# Lateral navigation file\n");
		fprintf(fid,"LNAV_FILE = \"\"\n");
	}

	// ---------
	if (pVNAVfilename)
	{
		fprintf(fid," \n");
		fprintf(fid,"# Vertical navigation file\n");
		fprintf(fid,"VNAV_FILE = \"%s\"\n", pVNAVfilename);
	}
	else
	{
		fprintf(fid," \n");
		fprintf(fid,"# Vertical navigation file\n");
		fprintf(fid,"VNAV_FILE = \"\"\n");
	}

	// ---------
	fprintf(fid," \n");
	fprintf(fid,"# Time definition:\n");
	//fprintf(fid, "T_START = %11.3f \n", ts);
	//fprintf(fid, "T_END = %11.3f \n", te);
	char buf[256];
	print_time_formatted(buf, ts);
	fprintf(fid, "T_START = %s \n", buf);
	print_time_formatted(buf, te);
	fprintf(fid, "T_END = %s \n", buf);

	// ---------
	fprintf(fid," \n");
	fprintf(fid,"# Speed mode:\n");
	switch (SPD_MODE)
	{
	case SPEED_MODE_IAS:
		fprintf(fid,"SPEED_MODE = IAS \n");
		break;
	case SPEED_MODE_MACH:
		fprintf(fid,"SPEED_MODE = MACH \n");
		break;
	}

	// ---------
	fprintf(fid," \n");
	fprintf(fid,"# Heading mode:\n");
	switch (HDG_MODE)
	{
	case HEADING_MODE_TRUEHDG:
		fprintf(fid,"HDG_MODE = TRUE_HDG\n");
		break;
	case HEADING_MODE_MAGNHDG:
		fprintf(fid,"HDG_MODE = MAGN_HDG\n");
		break;
	case HEADING_MODE_TRUETRK:
		fprintf(fid,"HDG_MODE = TRUE_TRK\n");
		break;
	case HEADING_MODE_MAGNTRK:
		fprintf(fid,"HDG_MODE = MAGN_TRK\n");
		break;
	case HEADING_MODE_GYROHDG:
		fprintf(fid,"HDG_MODE = GYRO_HDG\n");
		break;
	}

	// ---------
	fprintf(fid," \n");
	fprintf(fid,"# Engine mode:\n");
	switch (ENG_MODE)
	{
	case ENGINE_MODE_DUAL:
		fprintf(fid,"ENGINE_MODE = DUAL \n");
		break;
	case ENGINE_MODE_SINGLE:
		fprintf(fid,"ENGINE_MODE = SINGLE \n");
		break;
	}

	// ---------
	fprintf(fid," \n");
	fprintf(fid,"# Initial parameters\n");
	fprintf(fid, "LON = %10.6f \n", lon0);
	fprintf(fid, "LAT = %10.6f \n", lat0);
	fprintf(fid, "ALT = %10.2f \n", alt0);
	fprintf(fid, "WGT = %10.2f \n", wgt0);

	fprintf(fid," \n");
	fprintf(fid,"# Initial heading/track, and airspeed/mach\n");
	switch (HDG_MODE)
	{
	case HEADING_MODE_TRUEHDG:
	case HEADING_MODE_MAGNHDG:
	case HEADING_MODE_GYROHDG:
		switch (hdg0_type)
		{
		case -1:
			fprintf(fid, "HDG = %f \n", hdg0);
			break;

		case HEADING_MODE_TRUEHDG:
			fprintf(fid, "HDG = T%f \n", hdg0);
			break;

		case HEADING_MODE_MAGNHDG:
			fprintf(fid, "HDG = M%f \n", hdg0);
			break;

		case HEADING_MODE_GYROHDG:
			fprintf(fid, "HDG = G%f \n", hdg0);
			break;

		default:
			fprintf(fid, "HDG = %f \n", hdg0);
			break;
		}
		break;

	case HEADING_MODE_TRUETRK:
	case HEADING_MODE_MAGNTRK:
		switch (hdg0_type)
		{
		case -1:
			fprintf(fid, "TRK = %f \n", hdg0);
			break;

		case HEADING_MODE_TRUEHDG:
			fprintf(fid, "TRK = T%f \n", hdg0);
			break;

		case HEADING_MODE_MAGNHDG:
			fprintf(fid, "TRK = M%f \n", hdg0);
			break;

		case HEADING_MODE_GYROHDG:
			fprintf(fid, "TRK = G%f \n", hdg0);
			break;

		default:
			fprintf(fid, "TRK = %f \n", hdg0);
			break;
		}
		break;
	}

	switch (SPD_MODE)
	{
	case SPEED_MODE_IAS:
		fprintf(fid,"SPD = %f \n", spd0);
		break;
	case SPEED_MODE_MACH:
		fprintf(fid,"MACH = %f \n", spd0);
		break;
	}

	// ---------
	fprintf(fid," \n");
	fprintf(fid,"# Gyroscopic heading parameters (deg)\n");
	fprintf(fid,"GYRO_HDG_PARAM1 = %11.6f \n", gyro_hdg_param1);
	fprintf(fid,"GYRO_HDG_PARAM2 = %11.6f \n", gyro_hdg_param2);

	// ---------
	if (bank_angle > 0.0) // optional parameter
	{
		fprintf(fid," \n");
		fprintf(fid,"# Bank angle (deg)\n");
		fprintf(fid,"BANK_ANGLE = %11.6f \n", bank_angle);
	}

	// ---------
	fprintf(fid," \n");
	fprintf(fid,"# BTO Biases (microseconds)\n");
	fprintf(fid,"BTO_BIAS_RX1200 = %11.3f \n", bto_bias_rx1200);
	fprintf(fid,"BTO_BIAS_RX600 = %11.3f \n", bto_bias_rx600);
	fprintf(fid,"BTO_BIAS_TX = %11.3f \n", bto_bias_tx);

	// ---------
	fprintf(fid," \n");
	fprintf(fid,"# BFO Biases (Hz)\n");
	fprintf(fid,"BFO_BIAS_RX = %7.2f \n", bfo_bias_rx);
	fprintf(fid,"BFO_BIAS_TX = %7.2f \n", bfo_bias_tx);
	fprintf(fid,"BFO_BIAS = %7.2f \n", bfo_bias);


	// -----
	switch(DopplerMethod) // Doppler method (optional parameter) - print it if it was defined
	{
	case DOPPLER_METHOD_SIMPLE:
		fprintf(fid," \n");
		fprintf(fid,"# Doppler method\n");
		fprintf(fid,"DOPPLER_METHOD = SIM \n");
		break;

	case DOPPLER_METHOD_REL:
		fprintf(fid," \n");
		fprintf(fid,"# Doppler method\n");
		fprintf(fid,"DOPPLER_METHOD = REL \n");
		break;

	case DOPPLER_METHOD_ADV:
		fprintf(fid," \n");
		fprintf(fid,"# Doppler method\n");
		fprintf(fid,"DOPPLER_METHOD = ADV \n");
		break;
	}

	// -----
	switch(AES_DopplerMethod) // AES Doppler method (optional parameter) - print it if it was defined
	{
	case DOPPLER_METHOD_SIMPLE:
		fprintf(fid," \n");
		fprintf(fid,"# AES Doppler method\n");
		fprintf(fid,"AES_DOPPLER_METHOD = Simple \n");
		break;

	case DOPPLER_METHOD_REL:
		fprintf(fid," \n");
		fprintf(fid,"# AES Doppler method\n");
		fprintf(fid,"AES_DOPPLER_METHOD = REL \n");
		break;

	case DOPPLER_METHOD_ADV:
		fprintf(fid," \n");
		fprintf(fid,"# AES Doppler method\n");
		fprintf(fid,"AES_DOPPLER_METHOD = ADV \n");
		break;
	}


	// -----
	switch(AES_EarthModel) // AES Doppler method (optional parameter) - print it if it was defined
	{
	case AES_EARTH_MODEL_SPH:
		fprintf(fid," \n");
		fprintf(fid,"# AES Earth model (spherical or WGS'84)\n");
		fprintf(fid,"AES_EARTH_MODEL = SPH \n");
		break;

	case AES_EARTH_MODEL_WGS:
		fprintf(fid," \n");
		fprintf(fid,"# AES Earth model (spherical or WGS'84)\n");
		fprintf(fid,"AES_EARTH_MODEL = WGS \n");
		break;
	}

	// ---------
	if (AES_Earth_Radius > 0.0) // optional parameter
	{
		fprintf(fid," \n");
		fprintf(fid,"# Earth radius assumed by AES (m) for Doppler compensation term\n");
		fprintf(fid,"AES_EARTH_RADIUS = %11.6f \n", AES_Earth_Radius);
	}


	// ---------
	if (AES_GeoSat_Alt > 0.0) // optional parameter
	{
		fprintf(fid," \n");
		fprintf(fid,"# Geosynchronous orbit altitude assumed by AES (m) for Doppler compensation term\n");
		fprintf(fid,"AES_GEOSAT_ALT = %11.6f \n", AES_GeoSat_Alt);
	}

	// ---------
	fprintf(fid," \n");
	fprintf(fid,"# AES position/velocity data delay (s) for Doppler compensation term\n");
	fprintf(fid,"AES_DELAY = %11.6f \n", AES_Delay);

	// ---------
	if (pOutTrajFilename)
	{
		fprintf(fid," \n");
		fprintf(fid,"# Output of the Whole trajectory:\n");
		fprintf(fid,"OUTPUT_TRAJ_FILE = \"%s\"\n", pOutTrajFilename);
	}

	// ---------
	if (pOutStatFilename)
	{
		fprintf(fid," \n");
		fprintf(fid,"# Comparison with ping data:\n");
		fprintf(fid,"OUTPUT_BTO_BFO_FILE = \"%s\"\n", pOutStatFilename);
	}

	// ---------
	if (nOutputPingIDX>0 && pOutputPingIDX!=NULL)
	{
		fprintf(fid," \n");
		fprintf(fid,"# Indices of pings for comparison:\n");
		fprintf(fid,"OUTPUT_BTO_BFO_IDX = {");
		for (int i=0; i+1<nOutputPingIDX; i++) fprintf(fid,"%d, ", pOutputPingIDX[i]);
		fprintf(fid,"%d} \n", pOutputPingIDX[nOutputPingIDX-1]);
	}

	// ---------
	if (nOptPingsIDX>0 && pOptPingsIDX!=NULL)
	{
		fprintf(fid," \n");
		fprintf(fid,"# Indices of pings to be used for optimization:\n");
		fprintf(fid,"OPT_PING_IDX = {");
		for (int i=0; i+1<nOptPingsIDX; i++) fprintf(fid,"%d, ", pOptPingsIDX[i]);
		fprintf(fid,"%d} \n", pOptPingsIDX[nOptPingsIDX-1]);
	}

	// ---------
	if (nLNAVOptParams>0)
	{
		fprintf(fid," \n");
		fprintf(fid,"--------\n");
		fprintf(fid,"#Note: LNAV parameters (marked *) found as a result of optimization:\n");
		for (int i=0; i<nLNAVOptParams; i++) fprintf(fid,"%f \n", LNAVOptParams[i]);
		fprintf(fid,"--------\n");
	}

	fclose(fid);

	return true;
}



// check specified pings for output and optimization and remove those outside of the specified simulation time interval
void CInpParams::checkPingIds_()
{
	CMH370InmarsatData* pInmarsatData = new CMH370InmarsatData(); // create temporal Inmarsata data object

	int* pOutputPingIDXnew = new int[nOutputPingIDX+1];	// allocate new list of indices of pings to be included in statistic file when possible
	int nOutputPingIDXnew = 0;							// new number of indices of pings in the list to be included in statistic file

	int* pOptPingsIDXnew = new int[nOptPingsIDX+1];	// allocate new list of indices of pings to be used for path optimization 
	int nOptPingsIDXnew = 0;						// new number of pings used for path optimization


	for (int i=0; i<nOutputPingIDX; i++)
	{
		double tping;
		pInmarsatData->GetPingData(tping,  pOutputPingIDX[i]);
		if (tping >= ts && tping<=te)
		{
			// copy ping ID into the new list if this ping's time is in [ts, te]  
			pOutputPingIDXnew[nOutputPingIDXnew] = pOutputPingIDX[i];
			nOutputPingIDXnew++;
		}
	}


	for (int i=0; i<nOptPingsIDX; i++)
	{
		double tping;
		pInmarsatData->GetPingData(tping,  pOptPingsIDX[i]);
		if (tping >= ts && tping<=te)
		{
			// copy ping ID into the new list if this ping's time is in [ts, te]  
			pOptPingsIDXnew[nOptPingsIDXnew] = pOptPingsIDX[i];
			nOptPingsIDXnew++;
		}
	}

	if (outputstream_)
	{
		if (nOutputPingIDXnew < nOutputPingIDX) fprintf(outputstream_, "Warning: reporting indices outside of the time range were removed.\n");
		if (nOptPingsIDXnew < nOptPingsIDX) fprintf(outputstream_, "Warning: optimization indices outside of the time range were removed.\n");
	}

	// reassign arrays and number of pings
	if (pOutputPingIDX) delete pOutputPingIDX;
	pOutputPingIDX = pOutputPingIDXnew;
	nOutputPingIDX = nOutputPingIDXnew;

	if (pOptPingsIDX) delete pOptPingsIDX;
	pOptPingsIDX = pOptPingsIDXnew;
	nOptPingsIDX = nOptPingsIDXnew;

	// delete arrays if no pings in [ts,te] interval are found
	if (nOutputPingIDX == 0)
	{
		delete pOutputPingIDX;
		pOutputPingIDX = NULL;
	}

	if (nOptPingsIDX == 0)
	{
		delete pOptPingsIDX;
		pOptPingsIDX = NULL;
	}

	delete pInmarsatData;
}