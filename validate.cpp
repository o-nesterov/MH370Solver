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

#include "fminsearch.h"
#include "inmarsatdata.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Validation
//

extern double doppler(double XYZs[3], double XYZr[3], double UVWs[3], double UVWr[3], double Frq, bool isGravit=false);
extern void splineInterp(double* pY, double* pX, unsigned long n_out, double* pInpX, double* pInpY, unsigned long n_in);

extern void print_time_formatted(char* buf, double t);
extern bool parseLine(void* var, char* varname, char* pline, char* vartype);

extern bool read_time_formatted(double& t, char* buf);


bool readfile(char* filename, int nvars, int& nrec, double*& pTimeOut, double** ppVarsOut);
bool parseValidationSetupFile(char* setupfilename, 
							  char*& filename_alt, char*& filename_lonlat, char*& filename_uv, char*& outfilename, 
							  int& AES_Earth_Model, int& AES_Doppler_Method,
							  double& AES_GeoSat_Alt, double& AES_Earth_Radius, bool isVerbose);



// --------------------------------------------------------------------------------------------
//                                Flight phase
// --------------------------------------------------------------------------------------------

// filename - validation configuration filename
// it must contain 3 input file names for altitude, lon/lat location, and ground velocity data
// (there could be some difference due to the various sources of data, such as FlightRadar, or FI's ADS-B)
// and output statistic file name
// optionally, some parameters of the AES doppler compensation method, which were not disclosed by Inmarsat

void Validate(char* filename, bool isVerbose)
{
	// Parse input setup file
	char* filename_alt = NULL;
	char* filename_lonlat = NULL;
	char* filename_uv = NULL;
	char* outfilename = NULL;

	// Set default
	int AES_Earth_Model = AES_EARTH_MODEL_UDF;			// AES Earth model (optional, reset if specified)
	int AES_Doppler_Method = DOPPLER_METHOD_UDF;		// AES Doppler method (optional, reset if specified)
	double AES_GeoSat_Alt = -999.0;						// Nominal altitude of the satellite assumed by AES (optional, reset if specified)
	double AES_Earth_Radius = -999.0;					// Earth radius assumed by AES (optional, reset if specified)


	if (!parseValidationSetupFile(filename, 
								  filename_alt, filename_lonlat, filename_uv, outfilename, 
								  AES_Earth_Model, AES_Doppler_Method,
								  AES_GeoSat_Alt, AES_Earth_Radius, 
								  isVerbose))
	{
		if (isVerbose) printf("Error reading validation config file.\n");
		if (filename_alt) delete filename_alt;
		if (filename_lonlat) delete filename_lonlat;
		if (filename_uv) delete filename_uv;
		if (outfilename) delete outfilename;
		return;
	}

	// Load altitude, {lon, lat}, {u,v} data
	// Note that samplings of {u,v}, {lon,lat} and alt do not always match the same time, so each of the files has its own time stamps  

	double* ppVars[3]; // array of pointers on arrays
	bool isReadingError = false;

	// read altitude file
	int nAlt = 0;
	double* pTimeAlt = NULL;
	double* pAlt = NULL;
	if (!readfile(filename_alt, 1, nAlt, pTimeAlt, ppVars))
	{
		if (isVerbose) printf("Error reading altitude data file.\n");
		isReadingError = true;
	}
	pAlt = ppVars[0];

	// read lon lat file
	int nLonLat = 0;
	double* pTimeLonLat = NULL;
	double* pLon = NULL;
	double* pLat = NULL;
	if (!readfile(filename_lonlat, 2, nLonLat, pTimeLonLat, ppVars))
	{
		if (isVerbose) printf("Error reading lon/lat data file.\n");
		isReadingError = true;
	}
	pLon = ppVars[0];
	pLat = ppVars[1];

	// read uv file
	int nUV = 0;
	double* pTimeUV = NULL;
	double* pU = NULL;
	double* pV = NULL;
	if (!readfile(filename_uv, 2, nUV, pTimeUV, ppVars))
	{
		if (isVerbose) printf("Error reading uv data file.\n");
		isReadingError = true;
	}
	pU = ppVars[0];
	pV = ppVars[1];


	if (isReadingError || nAlt<1 || nLonLat<1 || nUV<1)
	{
		if (filename_alt) delete filename_alt;
		if (filename_lonlat) delete filename_alt;
		if (filename_uv) delete filename_alt;
		if (outfilename) delete outfilename;
		if (pTimeAlt) delete pTimeAlt;
		if (pAlt) delete pAlt;
		if (pTimeLonLat) delete pTimeLonLat;
		if (pLon) delete pLon;
		if (pLat) delete pLat;
		if (pTimeUV) delete pTimeUV;
		if (pU) delete pU;
		if (pV) delete pV;
		return;
	}

	// All data are read in

	// Now interpolate all data at regular interval (e.g., 1 s)

	double min_t = ((pTimeAlt[0]>pTimeLonLat[0])?pTimeAlt[0]:pTimeLonLat[0]);
	min_t =  ((min_t > pTimeUV[0])?min_t:pTimeUV[0]);

	double max_t = ((pTimeAlt[nAlt-1]<pTimeLonLat[nLonLat-1])?pTimeAlt[nAlt-1]:pTimeLonLat[nLonLat-1]);
	max_t =  ((max_t < pTimeUV[nUV-1])?max_t:pTimeUV[nUV-1]);

	min_t = ceil(min_t);
	max_t = floor(max_t);

	int n = (int)max_t - (int)min_t + 1;

	if (n<=0)
	{
		printf("Error: too short input interval.\n");
		return;
	}

	// Allocate arrays for interpolated data
	double* pTimeInt = new double[n+1];
	double* pAltInt = new double[n+1];
	double* pLonInt = new double[n+1];
	double* pLatInt = new double[n+1];
	double* pUInt = new double[n+1];
	double* pVInt = new double[n+1];
	double* pWInt = new double[n+1];

	// Define time 
	for (int i=0; i<n; i++) pTimeInt[i] = min_t + i;

	// Interpolate:
	splineInterp(pAltInt, pTimeInt, n, pTimeAlt, pAlt, nAlt);
	splineInterp(pLonInt, pTimeInt, n, pTimeLonLat, pLon, nLonLat);
	splineInterp(pLatInt, pTimeInt, n, pTimeLonLat, pLat, nLonLat);
	splineInterp(pUInt,   pTimeInt, n, pTimeUV, pU, nUV);
	splineInterp(pVInt,   pTimeInt, n, pTimeUV, pV, nUV);

	// compute vertical velocity by differentiating
	// however, since input data are truncated, smooth by taking difference over 'width' s intervals
	// the other method would be to compute differences, and then smooth them by applying some filter, e.g. Gaussian.
	
	const int width = 10;
	
	for (int i=0; i<n; i++) pWInt[i] = 0.0; // initialization

	for (int i=width; i<n-width; i++)
	{
		pWInt[i] = (pAltInt[i+width]-pAltInt[i-width])/(pTimeInt[i+width]-pTimeInt[i-width]);	// average vertical speed (m/s)
	}


// --------------------------------------------------------------------------------------------
//                                Create Inmarsat data object

	CMH370InmarsatData* pInmarsatData = new CMH370InmarsatData();
	if (isVerbose) pInmarsatData -> SetVerbose(stdout);
//	pInmarsatData -> SetBFOBiases(153.14, 152.79+0.84, 153.44+0.84+0.17); // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< experimental
	if (AES_Earth_Model != AES_EARTH_MODEL_UDF) pInmarsatData -> SetAESEarthModel(AES_Earth_Model);
	if (AES_Doppler_Method != DOPPLER_METHOD_UDF) pInmarsatData -> SetAESDopplerMethod(AES_Doppler_Method);
	if (AES_Earth_Radius>0.0) pInmarsatData -> SetAESEarthRadius(AES_Earth_Radius);
	if (AES_GeoSat_Alt>0.0) pInmarsatData -> SetAESGeoSatAlt(AES_GeoSat_Alt);

// --------------------------------------------------------------------------------------------
//                                Interpolate and derive differences

	const int idx_s = 139; // 16:41:52.907 - the first ping index
	const int idx_e = 498; // 17:07:48.907 - the last ping index

	double* pTimePings = new double[idx_e-idx_s+1];
	double* pAltPings  = new double[idx_e-idx_s+1];
	double* pLonPings  = new double[idx_e-idx_s+1];
	double* pLatPings  = new double[idx_e-idx_s+1];
	double* pUPings    = new double[idx_e-idx_s+1];
	double* pVPings    = new double[idx_e-idx_s+1];
	double* pWPings    = new double[idx_e-idx_s+1];

	double* pDistDiff = new double[idx_e-idx_s+1];	// difference arrays: computed minus inferred distance (Sat to airplane)
	double* pBTODiff = new double[idx_e-idx_s+1];	// difference computed bto minus measured bto
	double* pBFODiff = new double[idx_e-idx_s+1];	// difference computed bfo minus measured bfo

	// Get ping time
	for (int i=idx_s; i<=idx_e; i++)
	{
		pInmarsatData -> GetPingData(pTimePings[i-idx_s], i);
	}

	// Interpolate at the time of pings:
	splineInterp(pAltPings, pTimePings, idx_e-idx_s+1, pTimeAlt, pAlt, nAlt);
	splineInterp(pLonPings, pTimePings, idx_e-idx_s+1, pTimeLonLat, pLon, nLonLat);
	splineInterp(pLatPings, pTimePings, idx_e-idx_s+1, pTimeLonLat, pLat, nLonLat);
	splineInterp(pUPings,   pTimePings, idx_e-idx_s+1, pTimeUV, pU, nUV);
	splineInterp(pVPings,   pTimePings, idx_e-idx_s+1, pTimeUV, pV, nUV);
	splineInterp(pWPings,   pTimePings, idx_e-idx_s+1, pTimeInt, pWInt, n);

	char buftimetxt[256];
	char typetextbuf[8];

	// Print differences

	FILE* fid = NULL;
	if (outfilename!=NULL)
	{
		if (strlen(outfilename)>0)
		{
			fid = fopen(outfilename, "wt");
		}
	}

	bool verbose = false;
	if (fid == NULL) verbose = true;
	if (isVerbose) verbose = true;

	if (fid)
	{
		fprintf(fid, "Time, Lon(degE), Lat(degE), Alt(m), RoC(m/s), Diff.dist.(km), Diff.bto(microsecond), Diff.bfo(Hz)\n");
	}

	if (verbose)
	{
		printf("Time, Lon(degE), Lat(degE), Alt(m), RoC(m/s), Diff.dist.(km), Diff.bto(microsecond), Diff.bfo(Hz)\n");
	}


	for (int i=0; i<=idx_e-idx_s; i++)
	{
		double ts, dist, bto, bfo;
		int channelType;


		// Get i-th ping data
		pInmarsatData -> GetPingData(ts, dist, bto, bfo, channelType, i+idx_s);

		// calculate bto
		double bto_computed = pInmarsatData -> CalcBTO(pTimePings[i], pLonPings[i], pLatPings[i], pAltPings[i], channelType);

		// calculate distance
		double dist_computed = pInmarsatData -> CalcDistToSat(pTimePings[i], pLonPings[i], pLatPings[i], pAltPings[i]);

		// calculate bfo
		double bfo_computed = pInmarsatData -> CalcBFO(pTimePings[i], pLonPings[i], pLatPings[i], pAltPings[i], pUPings[i], pVPings[i], pWPings[i], channelType);

		print_time_formatted(buftimetxt, ts);

		switch(channelType)
		{
		case CH_TYPE_R_RX1200:
		case CH_TYPE_R_RX600:
			sprintf(typetextbuf,"R");
			break;

		case CH_TYPE_ANAMALOUS:
			sprintf(typetextbuf,"A");
			break;

		case CH_TYPE_T_RX1200:
			sprintf(typetextbuf,"T");
			break;

		case CH_TYPE_C_RX1200:
			sprintf(typetextbuf,"C");
			break;
		};

		pDistDiff[i] = dist_computed - dist;
		pBTODiff[i] = bto_computed - bto;
		pBFODiff[i] = bfo_computed - bfo;

		if (fid)
		{
			fprintf(fid, "%s  %7.3f %6.3f %8.1f %6.2f %5.2f %7.2f %7.2f  %s\n", &buftimetxt[11], pLonPings[i], pLatPings[i], pAltPings[i], pWPings[i], (dist_computed - dist)*0.001, bto_computed - bto, bfo_computed - bfo, typetextbuf);
		}

		if (verbose)
		{
			printf("%s  %7.3f %6.3f %8.1f %6.2f %5.2f %7.2f %7.2f  %s\n", &buftimetxt[11], pLonPings[i], pLatPings[i], pAltPings[i], pWPings[i], (dist_computed - dist)*0.001, bto_computed - bto, bfo_computed - bfo, typetextbuf);
		}
	}

	// ----------------------------------------------
	// Statistics after the turn, ascending
	int count = 0;
	double dist_avg = 0.0;	// average 
	double dist_dev = 0.0;	// deviation
	double dist_max = 0.0;	// max distance deviation
	double bto_avg = 0.0;	// average 
	double bto_dev = 0.0;	// deviation
	double bto_max = 0.0;	// max bto deviation
	double bfo_avg = 0.0;	// average 
	double bfo_dev = 0.0;	// deviation
	double bfo_max = 0.0;	// max bfo deviation

	for (int i=0; i<=idx_e-idx_s; i++)
	{
		if (pTimePings[i]>60577.73)
		{
			dist_avg += pDistDiff[i];
			dist_dev += pDistDiff[i]*pDistDiff[i];
			if (dist_max < abs(pDistDiff[i])) dist_max = abs(pDistDiff[i]);
			bto_avg += pBTODiff[i];
			bto_dev += pBTODiff[i]*pBTODiff[i];
			if (bto_max < abs(pBTODiff[i])) bto_max = abs(pBTODiff[i]);
			bfo_avg += pBFODiff[i];
			bfo_dev += pBFODiff[i]*pBFODiff[i];
			if (bfo_max < abs(pBFODiff[i])) bfo_max = abs(pBFODiff[i]);
			count++;
		}
	}

	dist_avg = dist_avg/(double)count;
	bto_avg = bto_avg/(double)count;
	bfo_avg = bfo_avg/(double)count;	
	dist_dev = sqrt(dist_dev/(double)count);
	bto_dev = sqrt(bto_dev/(double)count);
	bfo_dev = sqrt(bfo_dev/(double)count);

	if (fid)
	{
		fprintf(fid, "\n");
		fprintf(fid, "****************************************************\n");
		fprintf(fid, "\n");
		fprintf(fid, "AES Doppler correction term assumptions:\n");
		if (AES_Earth_Model == AES_EARTH_MODEL_WGS)
		{
			fprintf(fid, "Earth shape: WGS'84 ellipsoid. \n");
		}
		else
		{
			fprintf(fid, "Earth shape: spherical. \n");
			if (AES_Earth_Radius>0.0)
			{
				fprintf(fid, "Earth radius: %11.0f (m). \n", AES_Earth_Radius);
			}
			else
			{
				fprintf(fid, "Earth radius: default.\n");
			}

			if (AES_GeoSat_Alt>0.0)
			{
				fprintf(fid, "Geosynch. orbit altitude: %11.0f (m). \n", AES_GeoSat_Alt);
			}
			else
			{
				fprintf(fid, "Geosynch. orbit altitude: default.\n");
			}
		}

		fprintf(fid,"\n");
		fprintf(fid,"****************************************************\n");

		fprintf(fid, "\n");
		fprintf(fid, "Distance bias (pings after 16:50): %7.2f km.\n", dist_avg*0.001);
		fprintf(fid, "Distance STD  (pings after 16:50): %7.2f km.\n", dist_dev*0.001);
		fprintf(fid, "Distance maxd (pings after 16:50): %7.2f km.\n", dist_max*0.001);
		fprintf(fid, "\n");
		fprintf(fid, "BTO bias: %5.2f microseconds.\n", bto_avg);
		fprintf(fid, "BTO STD : %5.2f microseconds.\n", bto_dev);
		fprintf(fid, "BTO maxd: %5.2f microseconds.\n", bto_max);
		fprintf(fid, "\n");
		fprintf(fid, "BFO bias: %7.2f Hz.\n", bfo_avg);
		fprintf(fid, "BFO STD : %7.2f Hz.\n", bfo_dev);
		fprintf(fid, "BFO maxd: %7.2f Hz.\n", bfo_max);

	} // fid


	if (verbose)
	{
		printf("\n");
		printf("****************************************************\n");
		printf("\n");
		printf("AES Doppler correction term assumptions:\n");
		if (AES_Earth_Model == AES_EARTH_MODEL_WGS)
		{
			printf("Earth shape: WGS'84 ellipsoid. \n");
		}
		else
		{
			printf("Earth shape: spherical. \n");
			if (AES_Earth_Radius>0.0)
			{
				printf("Earth radius: %11.0f (m). \n", AES_Earth_Radius);
			}
			else
			{
				printf("Earth radius: default.\n");
			}

			if (AES_GeoSat_Alt>0.0)
			{
				printf("Geosynch. orbit altitude: %11.0f (m). \n", AES_GeoSat_Alt);
			}
			else
			{
				printf("Geosynch. orbit altitude: default.\n");
			}
		}

		printf("\n");
		printf("****************************************************\n");

		printf("\n");
		printf("Distance bias (pings after 16:50): %7.2f km.\n", dist_avg*0.001);
		printf("Distance STD  (pings after 16:50): %7.2f km.\n", dist_dev*0.001);
		printf("Distance maxd (pings after 16:50): %7.2f km.\n", dist_max*0.001);
		printf("\n");
		printf("BTO bias: %5.2f microseconds.\n", bto_avg);
		printf("BTO STD : %5.2f microseconds.\n", bto_dev);
		printf("BTO maxd: %5.2f microseconds.\n", bto_max);
		printf("\n");
		printf("BFO bias: %7.2f Hz.\n", bfo_avg);
		printf("BFO STD : %7.2f Hz.\n", bfo_dev);
		printf("BFO maxd: %7.2f Hz.\n", bfo_max);
	}

	// ----------------------------------------------
	// Statistics after reaching FL350
	count = 0;
	dist_avg = 0.0;	// average 
	dist_dev = 0.0;	// deviation
	dist_max = 0.0; // max distance deviation
	bto_avg = 0.0;	// average 
	bto_dev = 0.0;	// deviation
	bto_max = 0.0;	// max bto deviation
	bfo_avg = 0.0;	// average 
	bfo_dev = 0.0;	// deviation
	bfo_max = 0.0;	// max bfo deviation

	for (int i=0; i<=idx_e-idx_s; i++)
	{
		if (pTimePings[i]>61282.76)
		{
			dist_avg += pDistDiff[i];
			dist_dev += pDistDiff[i]*pDistDiff[i];
			if (dist_max < abs(pDistDiff[i])) dist_max = abs(pDistDiff[i]);
			bto_avg += pBTODiff[i];
			bto_dev += pBTODiff[i]*pBTODiff[i];
			if (bto_max < abs(pBTODiff[i])) bto_max = abs(pBTODiff[i]);
			bfo_avg += pBFODiff[i];
			bfo_dev += pBFODiff[i]*pBFODiff[i];
			if (bfo_max < abs(pBFODiff[i])) bfo_max = abs(pBFODiff[i]);
			count++;
		}
	}

	dist_avg = dist_avg/(double)count;
	bto_avg = bto_avg/(double)count;
	bfo_avg = bfo_avg/(double)count;	
	dist_dev = sqrt(dist_dev/(double)count);
	bto_dev = sqrt(bto_dev/(double)count);
	bfo_dev = sqrt(bfo_dev/(double)count);

	if (fid)
	{
		fprintf(fid, "\n");
		fprintf(fid, "****************************************************\n");
		fprintf(fid, "\n");
		fprintf(fid, "Distance bias after reaching FL350: %7.2f km.\n", dist_avg*0.001);
		fprintf(fid, "Distance STD  after reaching FL350: %7.2f km.\n", dist_dev*0.001);
		fprintf(fid, "Distance maxd after reaching FL350: %7.2f km.\n", dist_max*0.001);
		fprintf(fid, "\n");
		fprintf(fid, "BTO bias after reaching FL350: %5.2f microseconds.\n", bto_avg);
		fprintf(fid, "BTO STD  after reaching FL350: %5.2f microseconds.\n", bto_dev);
		fprintf(fid, "BTO maxd after reaching FL350: %5.2f microseconds.\n", bto_max);
		fprintf(fid, "\n");
		fprintf(fid, "BFO bias after reaching FL350: %7.2f Hz.\n", bfo_avg);
		fprintf(fid, "BFO STD  after reaching FL350: %7.2f Hz.\n", bfo_dev);
		fprintf(fid, "BFO maxd after reaching FL350: %7.2f Hz.\n", bfo_max);
	}

	if (verbose)
	{
		printf("\n");
		printf("****************************************************\n");
		printf("\n");
		printf("Distance bias after reaching FL350: %7.2f km.\n", dist_avg*0.001);
		printf("Distance STD  after reaching FL350: %7.2f km.\n", dist_dev*0.001);
		printf("Distance maxd after reaching FL350: %7.2f km.\n", dist_max*0.001);
		printf("\n");
		printf("BTO bias after reaching FL350: %5.2f microseconds.\n", bto_avg);
		printf("BTO STD  after reaching FL350: %5.2f microseconds.\n", bto_dev);
		printf("BTO maxd after reaching FL350: %5.2f microseconds.\n", bto_max);
		printf("\n");
		printf("BFO bias after reaching FL350: %7.2f Hz.\n", bfo_avg);
		printf("BFO STD  after reaching FL350: %7.2f Hz.\n", bfo_dev);
		printf("BFO maxd after reaching FL350: %7.2f Hz.\n", bfo_max);
	}


	if (fid) fclose(fid);

	// ----------------------------------------------
	// Free memory
	delete pDistDiff;
	delete pBTODiff;
	delete pBFODiff;


	delete pTimePings;
	delete pAltPings;
	delete pLonPings;
	delete pLatPings;
	delete pUPings;
	delete pVPings;
	delete pWPings;

	delete pTimeInt;
	delete pAltInt;
	delete pLonInt;
	delete pLatInt;
	delete pUInt;
	delete pVInt;
	delete pWInt;


	if (pTimeAlt) delete pTimeAlt;
	if (pAlt) delete pAlt;
	
	if (pTimeLonLat) delete pTimeLonLat;
	if (pLon) delete pLon;
	if (pLat) delete pLat;

	if (pTimeUV) delete pTimeUV;
	if (pU) delete pU;
	if (pV) delete pV;

	if (filename_alt) delete filename_alt;
	if (filename_lonlat) delete filename_lonlat;
	if (filename_uv) delete filename_uv;
	if (outfilename) delete outfilename;

	delete pInmarsatData;
}






//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Parse validation setup file
bool parseValidationSetupFile(char* setupfilename, 
							  char*& filename_alt, char*& filename_lonlat, char*& filename_uv, char*& outfilename, 
							  int& AES_Earth_Model, int& AES_Doppler_Method,
							  double& AES_GeoSat_Alt, double& AES_Earth_Radius, bool isVerbose)
{
	// Input: 
	// Setup file name

	// Output: 
	// filename_alt - file containing altitude data
	// filename_lonlat - file containing longitude-latitude data
	// filename_uv - file containing u (W->E) and v (S->N) velocity components
	// outfilename - output file name
	// AES_Earth_Model - sperical (SPH) or WGS'84 ellipsoid (WGS) model assumed by AES to calculate compensation term
	// AES_Doppler_Method - AES Doppler calculation method
	// AES_GeoSat_Alt - geosynchronous orbit altitude assumed by AES
	// AES_Earth_Radius - Earth radius assumed by AES when AES_Earth_Model = spherical 

	// Returs true if all data are successfully read

	// read in the whole file, skipping comment lines
	FILE* fid = fopen(setupfilename,"rt");

	if (fid == NULL) return false;

	// initialize output filenames
	filename_alt = NULL;
	filename_lonlat = NULL;
	filename_uv = NULL;
	outfilename = NULL;


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
	char* fvarnames[] = {"AES_GEOSAT_ALT", "AES_EARTH_RADIUS"};
	double fvars[sizeof(fvarnames)/sizeof(char*)+1];
	bool isfvars[sizeof(fvarnames)/sizeof(char*)+1];	// flag indicating if variable was read

	// string variable names to read
	char* svarnames[] = {"AES_EARTH_MODEL", "AES_DOPPLER_METHOD", "FILENAME_ALT", "FILENAME_LONLAT", "FILENAME_UV", "FILENAME_OUT"};
	char* svars[sizeof(svarnames)/sizeof(char*)+1];


	// floating point variables
	int n = sizeof(fvarnames)/sizeof(char*);
	for (int i=0; i<n; i++)
	{
		isfvars[i] = false;
		int j =0;
		while (j<nLines)
		{
			if (parseLine((void*)(&fvars[i]), fvarnames[i], ppLines[j],"double")) {isfvars[i] = true; break;}; // variable was successfully read
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
		AES_GeoSat_Alt = fvars[0];		// geosynchronous orbit altitude assumed by AES - optional parameter
	}


	if (isfvars[0])
	{
		AES_Earth_Radius = fvars[1];	// Earth radius assumed by AES - optional parameter
	}


	// Earth shape assumed by AES
	if (svars[0])
	{
		if (strcmp(svars[0],"ellipsoid")==0 || strcmp(svars[0],"WGS")==0 || strcmp(svars[0],"Ellipsoid")==0)
		{
			AES_Earth_Model = AES_EARTH_MODEL_WGS;		// geosynchronous orbit altitude assumed by AES - optional parameter
		}
		else
		{
			if (strcmp(svarnames[0],"spherical")==0 || strcmp(svarnames[0],"SPH")==0 || strcmp(svarnames[0],"Spherical")==0)
			{
				AES_Earth_Model = AES_EARTH_MODEL_SPH;
			}
			else
			{
				if (isVerbose) printf("The Earth shape to be assumed by AES is unrecognized. Assumed spherical.\n");
			}
		}
	}
	else
	{
		if (isVerbose) printf("The Earth shape to be assumed by AES is not specified. Assumed spherical.\n");
		if (isVerbose) printf("To specify explicitly use AES_EARTH_MODEL = SPH or WGS.\n");
	}


	// AES Doppler calculation method
	if (svars[1])
	{
		if (strcmp(svars[1],"simple")==0 || strcmp(svars[1],"Simple")==0)
		{
			AES_Doppler_Method = DOPPLER_METHOD_SIMPLE; // classic Doppler
		}
		else
		{
			if (strcmp(svars[1],"relativistic")==0 || strcmp(svars[1],"Relativistic")==0)
			{
				AES_Doppler_Method = DOPPLER_METHOD_REL; // Doppler with relativistic effect
			}
			else
			{
				if (strcmp(svars[1],"advanced")==0 || strcmp(svars[1],"Advanced")==0)
				{
					AES_Doppler_Method = DOPPLER_METHOD_ADV; // Doppler with relativistic and gravitational (red shift) effects
				}
				else
				{
					if (isVerbose) printf("Unrecognized Doppler method. Assumed classic.\n");
					AES_Doppler_Method = DOPPLER_METHOD_SIMPLE; // classic Doppler
				}
			}
		}
	}
	else
	{
		if (isVerbose) printf("Unspecified AES Doppler method. Assumed classic.\n");
		AES_Doppler_Method = DOPPLER_METHOD_SIMPLE; // classic Doppler
	}

	// All input file names must be specified
	if (svars[2] && svars[3] && svars[4])
	{
		filename_alt = new char[strlen(svars[2])+1];
		filename_lonlat = new char[strlen(svars[3])+1];
		filename_uv = new char[strlen(svars[4])+1];
		strcpy(filename_alt, svars[2]);
		strcpy(filename_lonlat, svars[3]);
		strcpy(filename_uv, svars[4]);
	}
	else
	{
		if (!svarnames[2])
		{
			if (isVerbose) printf("Error. Altitude file must be specified.\n");
			iOkFlag = false;
		}

		if (!svarnames[3])
		{
			if (isVerbose) printf("Error. Lon/Lat file must be specified.\n");
			iOkFlag = false;
		}

		if (!svarnames[4])
		{
			if (isVerbose) printf("Error. Velocity file must be specified.\n");
			iOkFlag = false;
		}
	}

	// Output file name is optional; if not specified, then output to stdio stream only
	if (svars[5])
	{
		outfilename = new char[strlen(svars[5])+1];
		strcpy(outfilename, svars[5]);
	}

	// -----------------------------------
	// Deallocate local string varibles
	n = sizeof(svarnames)/sizeof(char*);
	for (int i=0; i<n; i++) if (svars[i]) delete svars[i];

	return iOkFlag;
}




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Read file, which contains data:
// header
// time, var1, var2, ...
//
// time must be either in seconds since 2014-03-07 00:00:00 UTC, or formatted yyyy-mm-dd hh:mm:ss.sss

bool readfile(char* filename, int nvars, int& nrec, double*& pTimeOut, double** ppVarsOut)
{
	// Input:
	// filename - full filename name
	// navars - number of variables to read (besides time)
	// ppVarsOut - array of size nvars to contain pointers on arrays of variables

	// Output:
	// nrec - number of read records
	// pTimeOut - newly allocated array of length nrec containing time (seconds since 2014-03-07 00:00:00 UTC)
	// ppVarsOut[...] - newly allocated array of length nrec containing variable
	// note that pTimeOut and ppVarsOut[i] must be deallocated by the caller

	char buf[1024];

	int nr = 0; //actual number of records
	int nsize = 0; // size of array
	const int chunksize = 256; // chunk size to increase array size (to avoid reallocation every time).

	if (nvars<=0) return false;

	// initialize in case of failed reading
	nrec = 0;
	pTimeOut = NULL;
	for (int i=0; i<nvars; i++) ppVarsOut[i] = NULL;

	FILE* fid = fopen(filename, "rt");
	if (!fid) return false;


	fgets(buf, 1000, fid); // skip first line

	double* pTime = NULL;						// time
	double** ppVars = new double*[nvars+1];		// array of pointers on buffers
	double** ppVars_new = new double*[nvars+1];	// array of pointers on new buffers
	for (int i=0; i<nvars; i++) ppVars[i] = NULL;
	for (int i=0; i<nvars; i++) ppVars_new[i] = NULL;
	double* vars_buf = new double[nvars+1];		// buffer to accomodate nvars values read from a single line




	for (int i=0; i<nvars; i++)
	{
		ppVars[i] = NULL;
	};



	while(!feof(fid))
	{
		if (fgets(buf, 1000, fid)!=NULL)
		{
			char* p=&buf[0];
			int n = (int)strlen(p);
			while((*p==' ') || (*p=='	') || (*p==',') || (*p==';') && p<(&buf[n])) p++;	// skip spaces
			if (p==&buf[n] || n<=3) continue;	// skip empty line
			if (!((*p>='0' && *p<='9') || *p=='-' || *p=='+' || *p=='.')) continue;	// skip line as the first symbol is not a digit or sign, or decimal

			// -------------------------
			// try to read time as double
			char* pp;
			double ts = strtod(p, &pp);
			if (p == pp) continue; // unknown date format
			if (*pp=='-') // indicative of date as text string yyyy-.... 
			{
				// try to read time as text string
				if (!read_time_formatted(ts, p)) continue;
				while (*p!=':' && p<(&buf[n])) p++; // hh:mm separator
				p++;
				while (*p!=':' && p<(&buf[n])) p++; // mm:ss separator
				p++;
				strtod(p, &pp);
				p = pp+1;
			}
			else
			{
				p = pp+1;
			}

			// -------------------------
			// try to read variables
			int iv;
			for (iv=0; iv<nvars; iv++)
			{

				while((*p==' ') || (*p=='	') || (*p==',') || (*p==';') && p<(&buf[n])) p++;	// skip spaces etc.
				double var = strtod(p, &pp);
				if (p == pp) break; // unknown format
				p = pp+1;
				vars_buf[iv] = var;
			}
			if (iv < nvars) continue; // something went wrong; skip this line


			// everything is ok here
			if (nr==nsize)
			{
				// need to reallocate arrays and add new data record
				double* pTime_new = new double [nsize+chunksize];
				for (iv=0; iv<nvars; iv++) 
				{
					ppVars_new[iv] = new double [nsize+chunksize];
				}

				// copy all previous records into new arrays
				for (int i=0; i<nr; i++) pTime_new[i] = pTime[i];

				for (iv=0; iv<nvars; iv++)
				{
					for (int i=0; i<nr; i++)
					{
						ppVars_new[iv][i] = ppVars[iv][i];
					}
				}

				// add new record
				pTime_new[nr] = ts;
				for (iv=0; iv<nvars; iv++) ppVars_new[iv][nr] = vars_buf[iv];

				// delete old arrayay
				if (pTime) delete  pTime;
				for (iv=0; iv<nvars; iv++) if (ppVars[iv]) delete ppVars[iv];

				// reassign arrays
				pTime = pTime_new;
				for (iv=0; iv<nvars; iv++) ppVars[iv] = ppVars_new[iv];

				nr++;
				nsize += chunksize;
			}
			else
			{
				// just add new record
				pTime[nr] = ts;
				for (iv=0; iv<nvars; iv++) ppVars[iv][nr] = vars_buf[iv];
				nr++;
			}
		} // fgets
	}; //eof

	delete ppVars_new;
	delete vars_buf;

	// output
	nrec = nr;
	pTimeOut = pTime; // allocated array
	for (int iv=0; iv<nvars; iv++) ppVarsOut[iv] = ppVars[iv]; // allocated arrays

	fclose(fid);

	return true;
}

