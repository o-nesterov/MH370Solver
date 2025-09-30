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

#ifndef C_INMARSAT_DATA_CLASS

#define C_INMARSAT_DATA_CLASS

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Class to handle Inmarsat data:
// pings, satellite position and velocity, AES compensation, etc.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Channel types: 
// R-Channel RX 1200; 
// R-Channel RX 600; 
// T-Channel RX 1200; 
// C-Channel RX 1200; 
// anamalous data (R-Channel RX 1200)

// #define CH_TYPE_RX1200 0
// #define CH_TYPE_RX600 1
// #define CH_TYPE_T_RX 2
// #define CH_TYPE_C_RX 3
// #define CH_TYPE_ANAMAL 4

#define CH_TYPE_R_RX1200 0
#define CH_TYPE_R_RX600 1
#define CH_TYPE_T_RX1200 2
#define CH_TYPE_C_RX1200 3
#define CH_TYPE_ANAMALOUS 4

// Doppler computation methods:
// Undefined (only for input parsing purposes)
// Simple Doppler shift
// Doppler shift with relativistic effect
// Advanced: regular + relativistic + gravitational
#define DOPPLER_METHOD_UDF -1
#define DOPPLER_METHOD_SIMPLE 0
#define DOPPLER_METHOD_REL 1
#define DOPPLER_METHOD_ADV 2

// Assumed EAS Earth models, a choice of
// - undefined (only for input parsing purposes)
// - Sperical
// - WGS'84 ellipsoid
#define AES_EARTH_MODEL_UDF -1
#define AES_EARTH_MODEL_SPH 0
#define AES_EARTH_MODEL_WGS 1

class CMH370InmarsatData
{
public:
	 CMH370InmarsatData(); // constructor
	~CMH370InmarsatData(); // destructor

	// ------------------------------------------------------------------------------------------
	// Get Inmarsat ping data: ts - time stamp (seconds since 2014-03-07 00:00:00 UTC), dist - distance to aircraft (m), bto (microseconds), bfo (Hz),
	// channel type (0 = R-Channel RX 1200; 1 = R-Channel RX 600;  2 = T-Channel RX; 3 = C-Channel RX; 4 - anamalous data / unknown how to interpret),
	// n - is the sample number (input)
	bool GetPingData(double& ts, double& dist, double& bto, double& bfo, int& channelType, unsigned long n);
	bool GetPingData(double& ts, unsigned long n);

	// Return the total number of available pings
	unsigned long GetNumberOfPings();

	// Inmarsat satellite position {X,Y,Z} in ECEF system at given t (s) since 2014-03-07 00:00:00 UTC
	void GetInmarsatPos(double& X, double& Y, double& Z, double t);

	// Inmarsat satellite velocity {U,V,W} in ECEF system at given t (s) since 2014-03-07 00:00:00 UTC
	void GetInmarsatVel(double& U, double& V, double& W, double t);

	// Interpolated AES frequency compensation term, as recorded/measured by Inmarsat (Hz), including the eclipse effect
	// t (s) since 2014-03-07 00:00:00 UTC
	double GetAFC(double t);

	// Inmarsat uplink frequency (Hz)
	double GetUplinkFrq();

	// Inmarsat downlink frequency (Hz)
	double GetDownlinkFrq();

	// Set desired BTO biases instead of default ones: RX1200, RX600 and TX
	// This automatically re-calculates satellite-to-aircraft distances
	void SetBTOBiases(double BTO_Bias_RX1200, double BTO_Bias_RX600, double BTO_Bias_TX);

	// Set desired BFO biases
	void SetBFOBiases(double BFO_Bias, double BFO_Bias_RX, double BFO_Bias_TX);

	// Set Doppler calculation method
	void SetDopplerMethod(int Doppler_Method);

	// Set AES Earth Model: Sperical or WGS'84 Ellipsoid (it affects only calculation of the Doppler compensation term)
	void SetAESEarthModel(int AES_Earth_Model);

	// Set AES Doppler calculation method (it affects only calculation of the Doppler compensation term)
	void SetAESDopplerMethod(int Doppler_Method);

	// Set Earth radius assumed by AES when the spherical shape is assumed (default 6378137.0 m)
	void SetAESEarthRadius(double AES_Earth_Radius);

	// Set altitude of geosynchronous orbit assumed by AES (e.g., nominal or +422 km according to DSTG)
	void SetAESGeoSatAlt(double AES_GeoSat_Alt);

	// Calculate distance from a point {lon,lat,alt} to the satellite (m) at the given time t
	double CalcDistToSat(double t, double lon, double lat, double alt);

	// Calculate BTO
	double CalcBTO(double t, double lon, double lat, double alt, int channelType);

	// Calculate BFO
	double CalcBFO(double t, double lon, double lat, double alt, double u, double v, double w, int channelType);

	double CalcBFO(double t, double lon, double lat, double alt, double u, double v, double w, 
							 double aes_lon, double aes_lat, double aes_alt, double aes_u, double aes_v, int channelType);

	// returns KLIA ECEF coordinates
	void Get_KLIA_ECEF(double& X, double& Y, double& Z);

	// returns KLIA longitude and latitude
	void Get_KLIA_LonLat(double& Lon, double& Lat);

	// returns Perth ECEF coordinates
	void Get_Perth_ECEF(double& X, double& Y, double& Z);

	// returns Perth GES station longitude and latitude 
	void Get_Perth_LonLat(double& Lon, double& Lat);

	// Returns geostationary satellite ECEF coordinates (depending on set Earth model)
	void Get_AES_GeoSAT_ECEF(double& X, double& Y, double& Z);

	// Returns Burum station ECEF coordinates
	void Get_Burum_ECEF(double& X, double& Y, double& Z);

	// Returns Burum longitude and latitude
	void Get_Burum_LonLat(double& Lon, double& Lat);

	// ------------------------------------
	// set diagnostic stream (can be NULL)
	void SetVerbose(FILE* outputstream);

private:

	// Inmarsat ping data
	unsigned long NPings_;		// total number of pings
	double* pPingTime_;			// Ping times (seconds since 2014-03-07 00:00:00 UTC)
	double* pBTO_;				// BTO (microseconds)
	double* pBFO_;				// BFO (Hz)
	double* pDist_;				// Inferred distance from the satellite to the airplane based on BTO (m)
	int* pChannelType_;			// (0 = R-Channel RX 1200; 1 = R-Channel RX 600;  2 = T-Channel RX; 3 = C-Channel RX; 4 - anamalous data)
	void init_ping_data_();		// initialization

	// Inmarsat-recorded BFO compensation data (see description in ATSB 2014 as well as Ashton et. al)
	void initAFCLookupTable_();			// initialize table 
	double GetGESAFC_(double t);		// calculates Burum to Perth AFC at the given time
	double* pAFC_Lookup_Table_;			// lookup table
	int numAFCLookupTableEntries_;		// number of lookup table records
	double ref_time_AFC_lookup_table_;	// reference (start time) of the first record (s since 2014-03-07 00:00:00)
	double time_step_AFC_lookup_table_; // time step of the table
	int AFC_lookup_table_interp_method_;// AFC table interpolation method: 0 - linear, 1 - cubic 

	// constants:
	double Fr_uplink_;			// Uplink frequency (Hz)
	double Fr_downlink_;		// Downlink frequency (Hz)
	double Speed_of_light_;		// Speed of light (m/s)
	double Perth_LonLat_[2];	// Longitude and latitude of Perth GES (Australia) 
	double Burum_LonLat_[2];	// Longitude and latitude of Burum GES (Netherlands)
	double KLIA_LonLat_[2];		// Longitude and latitude of Kuala Lumpur airport
	double Perth_XYZ_[3];		// Perth GES coordinates in ECEF system
	double Burum_XYZ_[3];		// Burum GES coordinates in ECEF system
	double KLIA_XYZ_[3];		// KLIA coordinates in ECEF system
	double AES_Earth_Radius_;	// Earth radius assumed by AES, when sperical Earth shape is assumed
	double AES_GeoSat_lon_;		// Geostationary satellite longitude assumed by AES
	double AES_GeoSat_lat_;		// Geostationary satellite latitude assumed by AES
	double AES_GeoSat_alt_;		// Geostationary satellite altitude assumed by AES
	double AES_GeoSat_SPH_XYZ_[3]; // Geostationary satellite position in ECEF system assumed by AES using sperical Earth model
	double AES_GeoSat_WGS_XYZ_[3]; // Geostationary satellite position in ECEF system assumed by AES usign WGS'84 sllipsoid mode
	
	// BTO Biases for channels RX1200, RX600, and TX found as a result of calibration against KLIA using this model's geometry
	// BTO_Bias = -495679.0 microseconds according to ATSB
	double BTO_Bias_RX1200_;	// Bias for RX1200 channel
	double BTO_Bias_RX600_;		// Bias for RX600 channel
	double BTO_Bias_TX_;		// Bias for TX channel

	// BFO biases, Hz
	// BFO_Bias = 152.5; according to ATSB; 150.0 according to Ashton et al. (2015); 
	// the average is 153.18 according to calibration in this work, with:
	// 151.04 and 153.73 Hz for RX and TX channels respectively 
	double BFO_Bias_;			// Overall BFO bias, assumed for C-cahnnel as there is no calibration available
	double BFO_Bias_RX_;		// Bias for RX channels
	double BFO_Bias_TX_;		// Bias for TX channels

	// Doppler method, a choice of: 
	// - DOPPLER_METHOD_SIMPLE
	// - DOPPLER_METHOD_REL
	// - DOPPLER_METHOD_ADV
	int Doppler_Method_;
	int AES_Doppler_Method_;

	// Assumed EAS Earth models, a choice of
	// - AES_EARTH_MODEL_SPH (spherical)
	// - AES_EARTH_MODEL_WGS (WGS'84 ellipsoid)
	int AES_Earth_Model_;

	// Verbose diganostic stream
	FILE* fileout_; // message stream for diagnostic output (can be NULL) 
};

#endif // C_INMARSAT_DATA_CLASS
