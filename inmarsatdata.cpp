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
#include "inmarsatdata.h"


extern void LonLatH2XYZ(double &X, double &Y, double &Z, double lon, double lat, double hgt);
extern double doppler(double XYZs[3], double XYZr[3], double UVWs[3], double UVWr[3], double Frq, int Method);
extern void uv2uvw(double& XYZ_U, double& XYZ_V, double& XYZ_W, double u, double v, double lon, double lat);
extern void w2uvw(double& XYZ_U, double& XYZ_V, double& XYZ_W, double w, double lon, double lat);

extern void LonLatH2XYZ_sph(double &X, double &Y, double &Z, double lon, double lat, double hgt);
extern void uv2uvw_sph(double& XYZ_U, double& XYZ_V, double& XYZ_W, double u, double v, double lon, double lat);
extern void w2uvw_sph(double& XYZ_U, double& XYZ_V, double& XYZ_W, double w, double lon, double lat);
////////////////////////////////////////////////////////////////////////////////////////
// Class to handle Inmarsat data (Pings, satellite position, AES compensation etc.
////////////////////////////////////////////////////////////////////////////////////////

// Constructor: various initializations to handle environmental data

CMH370InmarsatData::CMH370InmarsatData()
{
	// ------------------------------------------------------
	// Inmarsat ping ring data
	NPings_ = 0;				// total number of pings
	pPingTime_ = NULL;			// Ping times (seconds since 2014-03-07 00:00:00 UTC)
	pBTO_ = NULL;				// BTO (microseconds)
	pBFO_ = NULL;				// BFO (Hz)
	pDist_ = NULL;				// Inferred distance from the satellite to the aircraft
	pChannelType_ = NULL;		// (0 = R-Channel RX 1200; 1 = R-Channel RX 600;  2 = T-Channel RX; 3 = C-Channel RX; 4 - anamalous data)

	// ------------------------------------------------------
	// AFC table (to calculate AFC term in BFO - see ATSB,2014 and Ashton et. al, 2015)
	numAFCLookupTableEntries_ = 0;		// number of entries in the AFC table
	pAFC_Lookup_Table_ = NULL;			// AFC lookup table
	ref_time_AFC_lookup_table_ = 0.0;	// reference (start time) of the first record (s since 2014-03-07 00:00:00)
	time_step_AFC_lookup_table_ = 0.0;	// time step of the table

	// ------------------------------------------------------
	// Uplink frequency (aircraft to sattelite, Hz) - Ashton et al.; ATSB 2014
	Fr_uplink_ = 1646652500.0;

	// Downlink frequency (sattelite to ground, Hz) - Ashton et al.; ATSB 2014
	Fr_downlink_ = 3615152500.0;

	// Speed of light (m/s)
	Speed_of_light_ = 299792500.0;

	// Perth ground station coordinates 
	Perth_LonLat_[0] = 115.8878; // 115.8869444; longitude, deg - Refined based on GoogleEarth
	Perth_LonLat_[1] = -31.8048; // -31.8044444; latitude, deg - Refined based on GoogleEarth

	// Burum ground station (GES) coordinates:   53°17'4"N   6°12'55"E
	Burum_LonLat_[0] = 6.21527778;
	Burum_LonLat_[1] = 53.2844444;

	// Kuala Lumpur airport terminal coordinates
	KLIA_LonLat_[0] = 101.7058;		// 101.7 (Ashton et al.);
	KLIA_LonLat_[1] =  2.7521;		// 2.7 in Ahston et al., which appears south of the airport
	
	// Geostationary satellite position assumed by AES
	// double geostat_radius = 42164000.0; // radius of geostationary orbit, m 
	// geostat_lonlath = [64.51, 0.0, geostat_radius-0.5*(Earth_R1+Earth_R2)];
	// geostat_lonlath = [64.5237, 0.0, geostat_radius-Earth_R2]; %http://www.satellite-calculations.com/Satellite/Catalog/catalogID.php?23839
	// geostat_lonlath = [64.5, 0.0, geostat_radius-Earth_R2]; %http://www.satellite-calculations.com/Satellite/Catalog/catalogID.php?23839
	// geostat_lonlath = [64.5, 0.0, 35788122.0+422000.0]; % nominal position of the satellite as provided by DSTG; altitude for correction was 422 km higher per DSTG 

	AES_Earth_Radius_ = 6378137.0;			// Assume that AES assumes equatorial radius
	AES_GeoSat_lon_ = 64.5;					// Geosynchronous satellite longitude assumed by AES according to DSTG (2015) "Bayesian Methods in the Search for MH370"
	AES_GeoSat_lat_ = 0.0;					// Geosynchronous satellite latitude assumed by AES
	AES_GeoSat_alt_ = 35788122.0+422000.0;	// Nominal altitude of the satellite assumed by AES according to DSTG; 
											// altitude assumed by AES for correction was 422 km higher than the actual geostationary orbit according to DSTG 

	// Perth GES coordinates in ECEF system (elevation assumed to be zero)
	LonLatH2XYZ(Perth_XYZ_[0], Perth_XYZ_[1], Perth_XYZ_[2], Perth_LonLat_[0], Perth_LonLat_[1], 0.0);

	// Burum GES coordinates in ECEF system (elevation assumed to be zero)
	LonLatH2XYZ(Burum_XYZ_[0], Burum_XYZ_[1], Burum_XYZ_[2], Burum_LonLat_[0], Burum_LonLat_[1], 0.0);

	// KLIA terminal coordinates in ECEF system (elevation assumed to be zero)
	LonLatH2XYZ(KLIA_XYZ_[0], KLIA_XYZ_[1], KLIA_XYZ_[2], KLIA_LonLat_[0], KLIA_LonLat_[1], 0.0);
 
	// Geostationary satellite position in ECEF system assumed by AES
	LonLatH2XYZ(AES_GeoSat_WGS_XYZ_[0], AES_GeoSat_WGS_XYZ_[1], AES_GeoSat_WGS_XYZ_[2], AES_GeoSat_lon_, AES_GeoSat_lat_, AES_GeoSat_alt_);
	LonLatH2XYZ_sph(AES_GeoSat_SPH_XYZ_[0], AES_GeoSat_SPH_XYZ_[1], AES_GeoSat_SPH_XYZ_[2], AES_GeoSat_lon_, AES_GeoSat_lat_, AES_Earth_Radius_ + AES_GeoSat_alt_ - 6378137.0);

	// ------------
	// Default biases found as a result of calibration using this work's model
	// (a single -495679.0 microseconds according to ATSB)
	BTO_Bias_RX1200_ = -495724.681;		// Calibrated R-RX1200 bias, microseconds
	BTO_Bias_RX600_ =  -491144.851;		// Calibrated R-RX600 bias, microseconds
	BTO_Bias_TX_ =     -500700.612;		// Calibrated T-RX1200 bias, microseconds

	// BFO_Bias = 152.5; according to ATSB; 150.0 according to Ashton et al. (2015); 
	// the average is 153.18 according to calibration in this work, with:
	// 151.04 and 153.73 Hz for RX and TX channels respectively 
	BFO_Bias_    = 153.18;
	BFO_Bias_RX_ = 151.04;
	BFO_Bias_TX_ = 153.73;


	// Default Doppler calculation methods (DOPPLER_METHOD_SIMPLE, DOPPLER_METHOD_REL,  DOPPLER_METHOD_ADV)
	Doppler_Method_ = DOPPLER_METHOD_ADV;			// 
	AES_Doppler_Method_ = DOPPLER_METHOD_SIMPLE;	// assume it is used by AES

	// Deafault assumed EAS Earth models, a choice of AES_EARTH_MODEL_SPH (spherical), AES_EARTH_MODEL_WGS (WGS'84 ellipsoid)
	AES_Earth_Model_ = AES_EARTH_MODEL_SPH; //WGS;

	// ------------
	init_ping_data_();			// initialization (pings)
	initAFCLookupTable_();		// measured (digitized) AES compensation term
	AFC_lookup_table_interp_method_ = 0; // the table already pre-interpolated to 10-min intervals, so use linear interpolant



	// ------------
	// Verbose diganostic stream
	fileout_ = NULL;			 // message stream for diagnostic output (can be NULL) 
}


////////////////////////////////////////////////////////////////////////////////////////
// Destructor: deallocate all arrays
////////////////////////////////////////////////////////////////////////////////////////
CMH370InmarsatData::~CMH370InmarsatData()
{
	// deallocate ping data
	if (pPingTime_) delete pPingTime_;
	if (pBTO_) delete pBTO_;
	if (pBFO_) delete pBFO_;
	if (pDist_) delete pDist_;
	if (pChannelType_) delete pChannelType_;

	if (pAFC_Lookup_Table_) delete pAFC_Lookup_Table_;
}



// Inmarsat uplink frequency (aircraft to sattelite, Hz)
double CMH370InmarsatData::GetUplinkFrq()
{
	return Fr_uplink_;
}

// Inmarsat downlink frequency (sattelite to ground, Hz)
double CMH370InmarsatData::GetDownlinkFrq()
{
	return Fr_downlink_;
}

// Returns KLIA ECEF coordinates
void CMH370InmarsatData::Get_KLIA_ECEF(double& X, double& Y, double& Z)
{
	X = KLIA_XYZ_[0];
	Y = KLIA_XYZ_[1];
	Z = KLIA_XYZ_[2];
}

// Returns KLIA longitude and latitude
void CMH370InmarsatData::Get_KLIA_LonLat(double& Lon, double& Lat)
{
	Lon = KLIA_LonLat_[0];
	Lat = KLIA_LonLat_[1];
}

// Returns Perth ECEF coordinates
void CMH370InmarsatData::Get_Perth_ECEF(double& X, double& Y, double& Z)
{
	X = Perth_XYZ_[0];
	Y = Perth_XYZ_[1];
	Z = Perth_XYZ_[2];
}

// Returns Perth's longitude and latitude
void CMH370InmarsatData::Get_Perth_LonLat(double& Lon, double& Lat)
{
	Lon = Perth_LonLat_[0];
	Lat = Perth_LonLat_[1];
}


// Returns geostationary satellite ECEF coordinates
void CMH370InmarsatData::Get_AES_GeoSAT_ECEF(double& X, double& Y, double& Z)
{
	if (AES_Earth_Model_==AES_EARTH_MODEL_WGS)
	{
		X = AES_GeoSat_WGS_XYZ_[0];
		Y = AES_GeoSat_WGS_XYZ_[1];
		Z = AES_GeoSat_WGS_XYZ_[2];
	}
	else
	{
		X = AES_GeoSat_SPH_XYZ_[0];
		Y = AES_GeoSat_SPH_XYZ_[1];
		Z = AES_GeoSat_SPH_XYZ_[2];
	}
}

// Returns Burum station ECEF coordinates
void CMH370InmarsatData::Get_Burum_ECEF(double& X, double& Y, double& Z)
{
	X = Burum_XYZ_[0];
	Y = Burum_XYZ_[1];
	Z = Burum_XYZ_[2];
}

// Returns Burum GES longitude and latitude
void CMH370InmarsatData::Get_Burum_LonLat(double& Lon, double& Lat)
{
	Lon = Burum_LonLat_[0];
	Lat = Burum_LonLat_[1];
}

// Set desired BTO biases instead default ones: RX1200, RX600 and TX
// and re-calculate satellite-to-aircraft distances
void CMH370InmarsatData::SetBTOBiases(double BTO_Bias_RX1200, double BTO_Bias_RX600, double BTO_Bias_TX)
{
	BTO_Bias_RX1200_ = BTO_Bias_RX1200;
	BTO_Bias_RX600_ = BTO_Bias_RX600;
	BTO_Bias_TX_ = BTO_Bias_TX;

	init_ping_data_();
}


// Set desired BFO biases instead default onese
void CMH370InmarsatData::SetBFOBiases(double BFO_Bias, double BFO_Bias_RX, double BFO_Bias_TX)
{
	BFO_Bias_ = BFO_Bias;
	BFO_Bias_RX_ = BFO_Bias_RX;
	BFO_Bias_TX_ = BFO_Bias_TX;
}


// Set Doppler calculation method
void CMH370InmarsatData::SetDopplerMethod(int Doppler_Method)
{
	Doppler_Method_ = Doppler_Method;
}

// Set AES Earth Model: Sperical or WGS'84 Ellipsoid (it affects only calculation of the Doppler compensation term)
void CMH370InmarsatData::SetAESEarthModel(int AES_Earth_Model)
{
	if (AES_Earth_Model != AES_EARTH_MODEL_UDF) AES_Earth_Model_ = AES_Earth_Model;
}

// Set AES Doppler calculation method (it affects only calculation of the Doppler compensation term)
void CMH370InmarsatData::SetAESDopplerMethod(int Doppler_Method)
{
	if (Doppler_Method != DOPPLER_METHOD_UDF) AES_Doppler_Method_ = Doppler_Method;
}

// Set Earth radius assumed by AES when SPH model is used 
void CMH370InmarsatData::SetAESEarthRadius(double AES_Earth_Radius)
{
	if (AES_Earth_Radius>0.0)
	{
		AES_Earth_Radius_ = AES_Earth_Radius;

		// update ECEF coordinates assumed by AES using LonLatH2XYZ_sph, which assumes radius of 6378137.0 m
		LonLatH2XYZ_sph(AES_GeoSat_SPH_XYZ_[0], AES_GeoSat_SPH_XYZ_[1], AES_GeoSat_SPH_XYZ_[2], AES_GeoSat_lon_, AES_GeoSat_lat_, AES_Earth_Radius + AES_GeoSat_alt_ - 6378137.0);
	}
}

// Set geosynchronous altitude assumed by AES 
void CMH370InmarsatData::SetAESGeoSatAlt(double AES_GeoSat_alt)
{
	if (AES_GeoSat_alt>0.0)
	{
		AES_GeoSat_alt_ = AES_GeoSat_alt;

		// update ECEF coordinates assumed by AES using LonLatH2XYZ_sph, which assumes radius of 6378137.0 m
		LonLatH2XYZ_sph(AES_GeoSat_SPH_XYZ_[0], AES_GeoSat_SPH_XYZ_[1], AES_GeoSat_SPH_XYZ_[2], AES_GeoSat_lon_, AES_GeoSat_lat_, AES_Earth_Radius_ + AES_GeoSat_alt - 6378137.0);

		// update ECEF coordinates assumed by AES using LonLatH2XYZ
		LonLatH2XYZ(AES_GeoSat_WGS_XYZ_[0], AES_GeoSat_WGS_XYZ_[1], AES_GeoSat_WGS_XYZ_[2], AES_GeoSat_lon_, AES_GeoSat_lat_, AES_GeoSat_alt);
	}
}



////////////////////////////////////////////////////////////////////////////////////////
//  Verbose for diagnostic (e.g., if data loading fails)
////////////////////////////////////////////////////////////////////////////////////////

void CMH370InmarsatData::SetVerbose(FILE* outputstream)
{
	// set verbose stream (can be NULL - no output)
	fileout_ = outputstream;
}




///////////////////////////////////////////////////////////////////////////////////////////
//  Calculate distance from a point {lon,lat,alt} to the satellite (m) at the given time t
///////////////////////////////////////////////////////////////////////////////////////////

double CMH370InmarsatData::CalcDistToSat(double t, double lon, double lat, double alt)
{
	double Xa, Ya, Za;	// aircraft position in ECEF
	double Xs, Ys, Zs;	// sattelite position in ECEF

	// convert lon,lat,alt to ECEF
	LonLatH2XYZ(Xa, Ya, Za, lon, lat, alt);

	// get satellite position
	GetInmarsatPos(Xs, Ys, Zs, t);

	return sqrt((Xa-Xs)*(Xa-Xs) + (Ya-Ys)*(Ya-Ys) + (Za-Zs)*(Za-Zs));
}



///////////////////////////////////////////////////////////////////////////////////////////
//  Calculate BTO (microseconds)
///////////////////////////////////////////////////////////////////////////////////////////

double CMH370InmarsatData::CalcBTO(double t, double lon, double lat, double alt, int channelType)
{
	// Input: t- time (s since 2014-03-07 00:00:00), lon - longitude (deg E), lat - latitude (deg N), channelType - channel type
	// Output: BTO (microseconds) for R-Channel RX 1200, R-Channel RX 600 and T-Channel RX; 

	double Xa, Ya, Za;	// aircraft position in ECEF
	double Xs, Ys, Zs;	// sattelite position in ECEF


	// convert lon,lat,alt to ECEF
	LonLatH2XYZ(Xa, Ya, Za, lon, lat, alt);

	// get satellite position
	GetInmarsatPos(Xs, Ys, Zs, t);


	// Channel Type:
	// CH_TYPE_R_RX1200 = R-Channel RX 1200; 
	// CH_TYPE_R_RX600  = R-Channel RX 600; 
	// CH_TYPE_T_RX1200 = T-Channel RX; 
	// CH_TYPE_C_RX1200 = C-Channel RX;
	// CH_TYPE_ANAMALOUS = Anamalous (R-Channel RX 1200)

	double bto_bias = 0.0;

	switch(channelType)
	{
	case CH_TYPE_R_RX1200:
		bto_bias = BTO_Bias_RX1200_;
		break;

	case CH_TYPE_R_RX600:
		bto_bias = BTO_Bias_RX600_;
		break;

	case CH_TYPE_T_RX1200:
		bto_bias = BTO_Bias_TX_;
		break;

	case CH_TYPE_ANAMALOUS:
		bto_bias = BTO_Bias_RX1200_; // anamalous pings were transmitted through R-RX Channel
		break;
	}

	// satellite to Perth distance
	double dist_sat_to_Perth = sqrt((Xs-Perth_XYZ_[0])*(Xs-Perth_XYZ_[0])+
									(Ys-Perth_XYZ_[1])*(Ys-Perth_XYZ_[1])+
									(Zs-Perth_XYZ_[2])*(Zs-Perth_XYZ_[2]));

	// satellite to aircraft distance
	double dist_sat_to_air = sqrt((Xs-Xa)*(Xs-Xa) + (Ys-Ya)*(Ys-Ya) + (Zs-Za)*(Zs-Za));

	double bto = (dist_sat_to_air + dist_sat_to_Perth)*2.0/(0.000001*Speed_of_light_) + bto_bias;
	if (channelType==CH_TYPE_C_RX1200) bto = -1.0E+10;  // no BTO was recorded for C-Channel, so BTO bias cannot be defined for C-Channel 

	return bto;

}


////////////////////////////////////////////////////////////////////////////////////////
//  Compute BFO (Hz)
////////////////////////////////////////////////////////////////////////////////////////

double CMH370InmarsatData::CalcBFO(double t, double lon, double lat, double alt, double u, double v, double w, 
											 double aes_lon, double aes_lat, double aes_alt, double aes_u, double aes_v, int channelType)
{
	// Input: t- time (s since 2014-03-07 00:00:00), lon - longitude (deg E), lat - latitude (deg N), alt - altitude (m), 
	//        u - W->E velocity component; v - S->N velocity component, w - vertical velocity component, channelType - channel type, which affects bias
	//        aes_lon, eas_lat, aes_alt, aes_u, aes_v - position and ground velocity data fed as input into AES (this data can be slightly delayed, subject to update frequency, but could be important during turns) 

	// Output: BFO (Hz)

	double aircraft_XYZ[3];					// aircraft position in ECEF, WGS'84
	double aircraft_UVW[3];					// aircraft velocity in ECEF, WGS'84
	double satellite_XYZ[3];				// satellite position in ECEF,  WGS'84
	double satellite_UVW[3];				// satellite velocity in ECEF,  WGS'84
	double aircraft_H_UVW[3];				// aircraft's horizontal velocity in ECEF,  WGS'84
	double aircraft_V_UVW[3];				// aircraft's vertical velocity in ECEF, WGS'84
	double aircraft_AES_XYZ[3];				// AES's assumed position (it assumes location at the ground in ECEF, depending on Earth model)
	double aircraft_AES_H_UVW[3];			// aircraft's horizontal velocity in ECEF, depening on Earth model
	double satellite_geostat_XYZ[3];		// geostationary satellite position in ECEF presumed by AES
	double Zero_UVW[] = {0.0, 0.0, 0.0};	// allocate zero velocity array of the same dimensions as satellite_UVW

	// convert lon,lat,alt to ECEF
	LonLatH2XYZ(aircraft_XYZ[0], aircraft_XYZ[1], aircraft_XYZ[2], lon, lat, alt);

	// convert aircraft's horizontal {u,v,0} to velocity in ECEF
	uv2uvw(aircraft_H_UVW[0], aircraft_H_UVW[1], aircraft_H_UVW[2], u, v, lon, lat);

	// convert aircraft's vertical {0,0,w} to velocity in ECEF
	w2uvw(aircraft_V_UVW[0], aircraft_V_UVW[1], aircraft_V_UVW[2], w, lon, lat);

	// full aircraft's velocity in ECEF
	aircraft_UVW[0] = aircraft_H_UVW[0] + aircraft_V_UVW[0];
	aircraft_UVW[1] = aircraft_H_UVW[1] + aircraft_V_UVW[1];
	aircraft_UVW[2] = aircraft_H_UVW[2] + aircraft_V_UVW[2];

	// convert {lon,lat,0} to ECEF assumed by AES
	// and convert velocity to ECEF from the Earth model's system assumed by AES 
	// (sperical or WGS'84 ellipsoid -  which Earth's geometry is actually used by AES?)
	if (AES_Earth_Model_ == AES_EARTH_MODEL_WGS)
	{
		LonLatH2XYZ(aircraft_AES_XYZ[0], aircraft_AES_XYZ[1], aircraft_AES_XYZ[2], aes_lon, aes_lat, 0.0);
		uv2uvw(aircraft_AES_H_UVW[0], aircraft_AES_H_UVW[1], aircraft_AES_H_UVW[2], aes_u, aes_v, aes_lon, aes_lat);
	}
	else
	{
		LonLatH2XYZ_sph(aircraft_AES_XYZ[0], aircraft_AES_XYZ[1], aircraft_AES_XYZ[2], aes_lon, aes_lat, 0.0);
		uv2uvw_sph(aircraft_AES_H_UVW[0], aircraft_AES_H_UVW[1], aircraft_AES_H_UVW[2], aes_u, aes_v, aes_lon, aes_lat);
	}

	// get satellite position in ECEF
	GetInmarsatPos(satellite_XYZ[0], satellite_XYZ[1], satellite_XYZ[2], t);

	// get satellite velocity in ECEF
	GetInmarsatVel(satellite_UVW[0], satellite_UVW[1], satellite_UVW[2], t);

	// get geostationary satellite position presumed by AES (depening on Earth model)
	Get_AES_GeoSAT_ECEF(satellite_geostat_XYZ[0], satellite_geostat_XYZ[1], satellite_geostat_XYZ[2]);

	// Uplink Doppler
	double dF1 = doppler(aircraft_XYZ, satellite_XYZ, aircraft_UVW, satellite_UVW, Fr_uplink_, Doppler_Method_);

	// Downlink Doppler
	double dF2 = doppler(satellite_XYZ, Perth_XYZ_, satellite_UVW, Zero_UVW, Fr_downlink_, Doppler_Method_);

	// aircraft doppler compensation 
	double dF3 = -doppler(aircraft_AES_XYZ, satellite_geostat_XYZ, aircraft_AES_H_UVW, Zero_UVW, Fr_uplink_, AES_Doppler_Method_);

	// compute AFC term based on Ashton et al.(2015): delta f_sat + delta f_AFC
	double dF45 = GetAFC(t);

	// get bias dF6 based on channel's type, as calibrated in this work
	double dF6 = 0.0;

	switch(channelType)
	{
	case CH_TYPE_R_RX1200:
		dF6 = BFO_Bias_RX_;
		break;

	case CH_TYPE_R_RX600:
		dF6 = BFO_Bias_RX_;
		break;

	case CH_TYPE_T_RX1200:
		dF6 = BFO_Bias_TX_;
		break;

	case CH_TYPE_C_RX1200:
		dF6 = BFO_Bias_;
		break;

	case CH_TYPE_ANAMALOUS:
		dF6 = BFO_Bias_RX_; // anamalous pings were transmitted through R-RX Channel
		break;
	}

	// total:
	double dF = dF1 + dF2 + dF3 + dF45 + dF6;

	return dF;
}


////////////////////////////////////////////////////////////////////////////////////////
//  Compute BFO (Hz) - The same as above, but assumed AES position and velocity to be the same as actual ones (no delay)
////////////////////////////////////////////////////////////////////////////////////////
double CMH370InmarsatData::CalcBFO(double t, double lon, double lat, double alt, double u, double v, double w, int channelType)
{
	// just duplicate lon, lat, alt, u, and v
	return CalcBFO(t, lon, lat, alt, u, v, w, lon, lat, alt, u, v, channelType);
}