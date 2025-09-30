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

// An auxliary class to simplify reading-in and managing input parameters, similar to struct but
// with internal handling of allocation and deallocation
// It can also save parameters, which is useful for optimization

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#ifndef CINPPARAM_CLASS

#define CINPPARAM_CLASS

class CInpParams
{
public:
	CInpParams();
	~CInpParams();

	// ------------------------------------
	// check parser version (reserved to know what parameters to read)
	int checkVersion(char* filename);

	// ------------------------------------
	// parse main setup file 
	bool parseParams(char* filename);

	// ------------------------------------
	// save parameters 
	bool saveParams(char* filename);

	// ------------------------------------
	// set diagnostic stream (can be NULL)
	void SetVerbose(FILE* outputstream);

	double ts;				// start time (s since 2014-03-07 00:00:00 UTC)
	double te;				// end time (s since 2014-03-07 00:00:00 UTC)
	double lon0;			// initial longitude (deg E)
	double lat0;			// initial latitude (deg N)
	double alt0;			// initial altitude (m)
	double wgt0;			// initial weight (kg)
	double hdg0;			// initial heading (deg from N, clockwise) or track, true or magnetic, depending on mode 
	double spd0;			// initial speed (IAS or MACH)
	double bank_angle;		// default bank angle (degrees) used to connect legs if turn is not explicitly specified in lnav file

	bool force_default_bank_angle; // override bank angles in lnav file, which are not marked with * (subject to individual optimization) with default bank angle

	double gyro_hdg_param1;	// two gyroscopic heading parameters
	double gyro_hdg_param2;

	int hdg0_type;			// type of initial heading: true, magnetic, or gyro; default is assumed to be the same as the heading mode

	double bto_bias_rx1200;	// bto bias for RX1200 channel (Hz)
	double bto_bias_rx600;	// bto bias for RX600 channel (Hz)
	double bto_bias_tx;		// bto bias for TX channel (Hz)

	double bfo_bias_rx;		// bfo bias for RX channel
	double bfo_bias_tx;		// bfo bias for TX channel
	double bfo_bias;		// bfo bias for other channels (CX)

	int SPD_MODE;			// speed mode: SPEED_MODE_IAS or SPEED_MODE_MACH
	int HDG_MODE;			// heading mode: HEADING_MODE_TRUEHDG, HEADING_MODE_MAGNHDG, HEADING_MODE_TRUETRK, HEADING_MODE_MAGNTR, or HEADING_MODE_GYRO
	int ENG_MODE;			// engine mode: BOTH or SINGLE

	int METEO_DATASET;		// currently supported METEO_DATASET_GDAS1 or METEO_DATASET_ERA5
	int MAGN_DECL_DATASET;	// currently supported only MAGN_DECL_DATASET_NOAA

	int MeteoInterpMethod;	// interpolation method for meteorological forcing (0,1,2, or 3)
	int MagnDeclMethod;		// interpolation method for magnetic declination (0 or 1)

	char* pMeteoDir;		// directory path, where meteorological data are stored
	char* pMagnDeclDir;		// directory path, where magnetic declination data are stored

	char* pLNAVfilename;	// lateral navigation profile file
	char* pVNAVfilename;	// vertical navigation profile file

	char* pOutTrajFilename;	// output trajectory filename
	char* pOutStatFilename; // output statistics filename

	int* pOutputPingIDX;	// list of indices of pings to be included in statistic file when possible
	int nOutputPingIDX;		// number of indices of pings in the list to be included in statistic file

	int* pOptPingsIDX;		// list of indices of pings to be used for path optimization 
	int nOptPingsIDX;		// number of pings used for path optimization

	int DopplerMethod;		// Doppler computation method (DOPPLER_METHOD_SIMPLE,  DOPPLER_METHOD_REL,  DOPPLER_METHOD_ADV)
	int AES_DopplerMethod;	// Doppler method used by AES
	int AES_EarthModel;		// AES Earth model: SPH (spherical) or WGS (WGS'84 ellipsoid)
	double AES_Earth_Radius;// AES Earth radius - optional parameter (affects AES doppler compensation term)
	double AES_GeoSat_Alt;	// AES geosynch. sat. orbit altitude - optional (affects AES doppler compensation term)
	double AES_Delay;		// Delay in position and velocity data (affects AES doppler compensation term)

	double LNAVOptParams[32]; // LNAV optimization parameters
	int nLNAVOptParams;		// number of optimization parameters in LNAV profile

private:
	void checkPingIds_();	// check specified pings for output and optimization and remove those outside of the specified time interval

	int parserVersion_;		// version number, reserved for compatibility with future versions (corresponding list of parameters, reading, and saving)
	FILE* outputstream_;	// output file stream

};

#endif