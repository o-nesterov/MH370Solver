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
// Calibration / Validation
//

extern double doppler(double XYZs[3], double XYZr[3], double UVWs[3], double UVWr[3], double Frq, bool isGravit=false);

extern void print_time_formatted(char* buf, double t);


// -------------------------------------------------------------------------
//                  Calibration is based on tarmac phase data
// -------------------------------------------------------------------------

// Note that the velocity whilst on the tarmac remained zero, effectively setting the Doppler compensation term of BFO to be zero

void Calibrate(char* outfilename, bool isVerbose)
{
	// Create Inmarsat data object
	CMH370InmarsatData* pInmarsatData = new CMH370InmarsatData();
	if (isVerbose) pInmarsatData -> SetVerbose(stdout);

	// Set BTO and BFO biases to zeros to be able to call this object's BTOCalc and BFOCalc functions to calculate all other terms
	// without a need of duplicating detailed code and constants.
	pInmarsatData -> SetBTOBiases(0.0, 0.0, 0.0);	
	pInmarsatData -> SetBFOBiases(0.0, 0.0, 0.0);

	// Open output file stream, if specified
	FILE* fid = NULL;
	if (strlen(outfilename))
	{
		fid = fopen(outfilename, "wt");
	}
	
	// Make output either to the file or on screen
	bool verbose = false;
	if (fid == NULL) verbose = true;
	if (isVerbose) verbose = true;

	
	// Create list of pings to be used for calibration - finding BTO and BFO biases
	// there are 139 data pings on KLIA tarmac from 15:59:55.413 to 16:29:52.406
	// see pings_data.cpp, or use "-printpings" command to print pings
	const int N_TARMAC_PINGS = 139;

	// --------------------------------------------------------------

	// BTO_Bias are different for different channels 
	// CH_TYPE_R_RX1200 = R-Channel RX 1200; 
	// CH_TYPE_R_RX600  = R-Channel RX 600; 
	// CH_TYPE_T_RX1200 = T-Channel RX; 
	// CH_TYPE_C_RX1200 = C-Channel RX; - no BTO data
	// CH_TYPE_ANAMALOUS - anamalous data not to be counted

	// Note that the distance from the satellite to aircraft (see ATSB, 2014; Ashton et al., 2015) is:
	// distance_satellite_to_aircraft = (BTO - Bto_Bias) x Speed_of_light/2.0 - distance_satellite_to_Perth =>
	// BTO_Bias = BTO - (distance_satellite_to_aircraft + distance_satellite_to_Perth)*2/Speed_of_light

	double KLIA_LONLAT[2];	// KLIA lon lat

	// Get KLIA coordinates
	pInmarsatData -> Get_KLIA_LonLat(KLIA_LONLAT[0], KLIA_LONLAT[1]);


	double BTO_Bias_RX1200 = 0.0;	// Bias for RX1200 channel
	double BTO_Bias_RX600 = 0.0;	// Bias for RX600 channel
	double BTO_Bias_TX = 0.0;		// Bias for TX channel
	int nBTO_Bias_RX1200 = 0;		// number of samples for RX1200 channel
	int nBTO_Bias_RX600 = 0;		// number of samples for RX600 channel
	int nBTO_Bias_TX = 0;			// number of samples for TX channel

	for (int i=0; i<N_TARMAC_PINGS; i++)
	{
		double ts;			// time (s since 2014-03-07 00:00:00 UTC);
		double dist;		// used here only to utilize GetPingData (not used in the following calculations)
		double bto, bfo;	// BTO and BFO data
		int channelType;	// Channel type
		
		// Get i-th ping data
		pInmarsatData -> GetPingData(ts, dist, bto, bfo, channelType, i);

		// Since all BTO biases were set to zeros, use CalcBTO function (channel ID does not matter because biases are zeros) 
		// It will calculate  (dist_sat_to_Perth + dist_sat_to_KLIA)*2.0/(0.000001*Speed_of_light)
		double bto0 = pInmarsatData -> CalcBTO(ts, KLIA_LONLAT[0], KLIA_LONLAT[1], 0.0, CH_TYPE_R_RX1200);

		// Compute difference
		double bto_diff = bto - bto0;

		switch(channelType)
		{
		case CH_TYPE_R_RX1200:
			BTO_Bias_RX1200 += bto_diff;
			nBTO_Bias_RX1200++;
			break;

		case CH_TYPE_R_RX600:
			BTO_Bias_RX600 += bto_diff;
			nBTO_Bias_RX600++;
			break;

		case CH_TYPE_T_RX1200:
			BTO_Bias_TX += bto_diff;
			nBTO_Bias_TX++;
			break;
		};
	};


	// Calculate averages:
	BTO_Bias_RX1200 = BTO_Bias_RX1200 / (double)nBTO_Bias_RX1200;
	BTO_Bias_RX600 = BTO_Bias_RX600 / (double)nBTO_Bias_RX600;
	BTO_Bias_TX = BTO_Bias_TX / (double)nBTO_Bias_TX;

	if (fid)
	{
		fprintf(fid, "***************************************************************\n");
		fprintf(fid, "**************************** BTO ******************************\n");
		fprintf(fid, "***************************************************************\n");
		fprintf(fid, "\n");
		fprintf(fid, "BTO biases:\n");
		fprintf(fid, "R-Channel RX 1200 BTO bias = %11.1f (microsecond)\n", BTO_Bias_RX1200);
		fprintf(fid, "R-Channel RX 600  BTO bias = %11.1f (microsecond)\n", BTO_Bias_RX600);
		fprintf(fid, "T-Channel RX 1200 BTO bias = %11.1f (microsecond)\n", BTO_Bias_TX);
		fprintf(fid, "\n");
	}
	
	if (verbose)
	{
		printf("***************************************************************\n");
		printf("**************************** BTO ******************************\n");
		printf("***************************************************************\n");
		printf("\n");
		printf("BTO biases:\n");
		printf("R-Channel RX 1200 BTO bias = %11.1f (microsecond)\n", BTO_Bias_RX1200);
		printf("R-Channel RX 600  BTO bias = %11.1f (microsecond)\n", BTO_Bias_RX600);
		printf("T-Channel RX 1200 BTO bias = %11.1f (microsecond)\n", BTO_Bias_TX);
		printf("\n");
	}


	// ------------------------------------------------------------------------------------
	// Calculate BTO residuals and STD deviations (for reporting)
	// ------------------------------------------------------------------------------------
	
	char buftimetxt[256]; // time string


	// Reset BTO biases
	pInmarsatData -> SetBTOBiases(BTO_Bias_RX1200, BTO_Bias_RX600, BTO_Bias_TX);


	double STD_BTO_RX1200 = 0.0;	// STD RX1200 channel
	double STD_BTO_RX600 = 0.0;		// STD RX600 channel
	double STD_BTO_TX = 0.0;		// STD TX channel
	double STD_BTO = 0.0;			// Overall STD
	double max_BTO_dev = 0.0;		// maximum deviation

	if (fid)
	{
		fprintf(fid,"Time, Diff.BTO(computed minus measured, microsecond)\n");
	}
	
	if (verbose)
	{
		printf("Time, Diff.BTO(computed minus measured, microsecond)\n");
	}

	for (int i=0; i<N_TARMAC_PINGS; i++)
	{
		double ts;			// time (s since 2014-03-07 00:00:00 UTC);
		double dist;		// used here only to utilize GetPingData (not used in the following calculations)
		double bto, bfo;	// BTO and BFO data
		int channelType;	// Channel type
		
		// Get i-th ping data
		pInmarsatData -> GetPingData(ts, dist, bto, bfo, channelType, i);

		// Compute BTO with new biases 
		double bto_computed = pInmarsatData -> CalcBTO(ts, KLIA_LONLAT[0], KLIA_LONLAT[1], 0.0, channelType);

		// Compute residuals BTO (measured BTO minus calculated BTO)
		double bto_diff = bto_computed - bto;

		switch(channelType)
		{
		case CH_TYPE_R_RX1200:
			STD_BTO_RX1200 += bto_diff*bto_diff;
			if (max_BTO_dev < abs(bto_diff)) max_BTO_dev = abs(bto_diff);
			break;

		case CH_TYPE_R_RX600:
			STD_BTO_RX600 += bto_diff*bto_diff;
			if (max_BTO_dev < abs(bto_diff)) max_BTO_dev = abs(bto_diff);
			break;

		case CH_TYPE_T_RX1200:
			STD_BTO_TX += bto_diff*bto_diff;
			if (max_BTO_dev < abs(bto_diff)) max_BTO_dev = abs(bto_diff);
			break;
		};

		if (channelType == CH_TYPE_R_RX1200 || channelType == CH_TYPE_R_RX600 || CH_TYPE_T_RX1200)
		{
			print_time_formatted(buftimetxt, ts);

			if (fid)
			{
				fprintf(fid, "%s  %7.2f\n", &buftimetxt[11], bto_diff);
			}
	
			if (verbose)
			{
				printf("%s  %7.2f\n", &buftimetxt[11], bto_diff);
			}
		}
	};


	// STD deviations:
	STD_BTO = sqrt((STD_BTO_RX1200 + STD_BTO_RX600 + STD_BTO_TX) / (double)(nBTO_Bias_RX1200 + nBTO_Bias_RX600 + nBTO_Bias_TX));
	STD_BTO_RX1200 = sqrt(STD_BTO_RX1200 / (double)nBTO_Bias_RX1200);
	STD_BTO_RX600 = sqrt(STD_BTO_RX600 / (double)nBTO_Bias_RX600);
	STD_BTO_TX = sqrt(STD_BTO_TX / (double)nBTO_Bias_TX);

	if (fid)
	{
		fprintf(fid, "\n");
		fprintf(fid, "****************************************************\n");
		fprintf(fid, "\n");
		fprintf(fid, "BTO deviations:\n");
		fprintf(fid, "R-Channel RX 1200 BTO STD = %6.2f (microsecond) - %d samples\n", STD_BTO_RX1200, nBTO_Bias_RX1200);
		fprintf(fid, "R-Channel RX 600  BTO STD = %6.2f (microsecond) - %d samples\n", STD_BTO_RX600, nBTO_Bias_RX600);
		fprintf(fid, "T-Channel RX 1200 BTO STD = %6.2f (microsecond) - %d samples\n", STD_BTO_TX, nBTO_Bias_TX);
		fprintf(fid, "Overall BTO STD = %6.2f (microsecond) - %d samples\n", STD_BTO, nBTO_Bias_RX1200+nBTO_Bias_RX600+nBTO_Bias_TX);
		fprintf(fid, "Maximum BTO deviation = %6.2f (microsecond)\n", max_BTO_dev);
		fprintf(fid, "\n");
		fprintf(fid, "\n");
	}
	
	if (verbose)
	{
		printf("\n");
		printf("****************************************************\n");
		printf("\n");
		printf("BTO deviations:\n");
		printf("R-Channel RX 1200 BTO STD = %6.2f (microsecond) - %d samples\n", STD_BTO_RX1200, nBTO_Bias_RX1200);
		printf("R-Channel RX 600  BTO STD = %6.2f (microsecond) - %d samples\n", STD_BTO_RX600, nBTO_Bias_RX600);
		printf("T-Channel RX 1200 BTO STD = %6.2f (microsecond) - %d samples\n", STD_BTO_TX, nBTO_Bias_TX);
		printf("Overall BTO STD = %6.2f (microsecond) - %d samples\n", STD_BTO, nBTO_Bias_RX1200+nBTO_Bias_RX600+nBTO_Bias_TX);
		printf("Maximum BTO deviation = %6.2f (microsecond)\n", max_BTO_dev);
		printf("\n");
		printf("\n");
	}


// -------------------------------------------------------------------------
// BFOs can be grouped by channel types, or time clusters, or both
// There are two clusters with apparently different BFO "offsets", t<=16:00:23.407 (idx = 0..7) and 16:29:17.407 <= t <= 16:29:41.907 (idx = 90...122)
// The first is possibly due to SDU startup (see Holland's paper); 
// the second may presumably be because of the movement for taxiing for take off 16:41:58 (take off cleared 16:40:37) - (FI)
// KLIA runway alighnment is THDG 346 deg (see Fig1.1F "Primary radar targets" and Appendix 1.6B in FI)
// However, experiments to reproduce BFOs 16:29:17.407 - 16:29:41.907 assuming the aircraft was moving HDG 346-180 deg for take-off and 
// assuming the same BFO biases so far were unsuccessful. 


	// split into two groups - count or not count
	int idx_group[N_TARMAC_PINGS+1];
	for (int i=0; i<N_TARMAC_PINGS; i++)
	{
		idx_group[i] = 0;
		if (i>=0 && i<=7) continue;		// skip BFOs 15:59:55.413 to 16:00:23.407, which are ~16 Hz higher than others, possibly due to SDU startup
		if (i>=90 && i<=122) continue;	// skip BFOs 16:29:17.407 to 16:29:41.907, which are ~5-10 Hz higher than the rest, possibly due to taxiing for take off
		idx_group[i] = 1;
	}

	// array of BFOs minus biases
	double dF[N_TARMAC_PINGS+1];

	for (int i=0; i<N_TARMAC_PINGS; i++)
	{
		double ts;

		// Get i-th ping time
		pInmarsatData -> GetPingData(ts, i);

		// Be using that all BFO biases are reset to zeros in pInmarsatData object, it is possible to utilize CalcBFO here
		// to compute the sum of all other terms (dF = dF1+dF2+dF3+dF45) without duplicating the code
		// Channel ID does not matter in this case because all biases are zeros, so let it be CH_TYPE_R_RX1200
		dF[i] = pInmarsatData -> CalcBFO(ts, KLIA_LONLAT[0], KLIA_LONLAT[1], 0.0, 0.0, 0.0, 0.0, CH_TYPE_R_RX1200);

	};


	double BFO_Bias = 0.0;		// overall bias; note, there was no C-Channel messages
	int nBFO_Bias = 0;			// number of samples
	double BFO_Bias_RX = 0.0;	// RX BFO bias; note RX600 occurred only two times in the startup cluster (15:59:55.413 and 15:59:56.413), which is excluded from deriving averages
	int nBFO_Bias_RX = 0;		// number of samples
	double BFO_Bias_TX = 0.0;	// TX BFO bias
	int nBFO_Bias_TX = 0;		// number of samples

	for (int i=0; i<N_TARMAC_PINGS; i++)
	{
		if (idx_group[i]==0) continue;		// skip BFOs 15:59:55.413 to 16:00:23.407, which are ~16 Hz higher than others, possibly due to SDU startup
											// skip BFOs 16:29:17.407 to 16:29:41.907, which are ~5-10 Hz higher than the rest, possibly due to taxiing for take off

		double ts;			// time (s since 2014-03-07 00:00:00 UTC);
		double dist;		// used here only to utilize GetPingData (not used in the following calculations)
		double bto, bfo;	// BTO and BFO data
		int channelType;	// Channel type
		
		// Get i-th ping data
		pInmarsatData -> GetPingData(ts, dist, bto, bfo, channelType, i);
	
		BFO_Bias += (bfo - dF[i]); // overall bias
		nBFO_Bias++;

		switch(channelType)
		{
		case CH_TYPE_R_RX1200:
		case CH_TYPE_R_RX600:
			BFO_Bias_RX += (bfo - dF[i]);  // R-RX bias; RX600 occurred only in the startup cluster
			nBFO_Bias_RX++;
			break;

		case CH_TYPE_T_RX1200:
			BFO_Bias_TX += (bfo - dF[i]);  // T-RX bias
			nBFO_Bias_TX++;
			break;
		};
	};

	// compute biases as averages
	BFO_Bias = BFO_Bias/(double)nBFO_Bias;
	BFO_Bias_RX = BFO_Bias_RX/(double)nBFO_Bias_RX;
	BFO_Bias_TX = BFO_Bias_TX/(double)nBFO_Bias_TX;

	if (fid)
	{
		fprintf(fid, "***************************************************************\n");
		fprintf(fid, "**************************** BFO ******************************\n");
		fprintf(fid, "***************************************************************\n");
		fprintf(fid, "\n");
		fprintf(fid, "BFO biases split into R- and T-channel groups:\n");
		fprintf(fid, "R-Channel RX 1200 BFO bias = %6.2f (Hz)\n", BFO_Bias_RX);
		fprintf(fid, "T-Channel RX 1200 BFO bias = %6.2f (Hz)\n", BFO_Bias_TX);
		fprintf(fid, "Overall BFO bias = %6.2f (Hz)\n", BFO_Bias);
		fprintf(fid, "\n");
	}
	
	if (verbose)
	{
		printf("***************************************************************\n");
		printf("**************************** BFO ******************************\n");
		printf("***************************************************************\n");
		printf("\n");
		printf("BFO biases split into R- and T-channel groups:\n");
		printf("R-Channel RX 1200 BFO bias = %6.2f (Hz)\n", BFO_Bias_RX);
		printf("T-Channel RX 1200 BFO bias = %6.2f (Hz)\n", BFO_Bias_TX);
		printf("Overall BFO bias = %6.2f (Hz)\n", BFO_Bias);
		printf("\n");
	}


	// compute statistics for reporting
	double BFO_STD = 0.0;		// overall bias; note, there was no C-Channel messages
	double BFO_STD_RX = 0.0;	// RX BFO bias; note RX600 occurred only two times in the startup cluster (15:59:55.413 and 15:59:56.413), which is excluded from deriving averages
	double BFO_STD_TX = 0.0;	// TX BFO bias
	double max_BFO_dev = 0.0;	// maximum deviation
	double diff;

	// Print output header
	if (fid)
	{
		fprintf(fid,"Time, Diff.BFO(computed minus measured, Hz)\n");
	}
	
	if (verbose)
	{
		printf("Time, Diff.BFO(computed minus measured, Hz)\n");
	}

	for (int i=0; i<N_TARMAC_PINGS; i++)
	{
		if (idx_group[i]==0) continue;		// skip BFOs 15:59:55.413 to 16:00:23.407, which are ~16 Hz higher than others, possibly due to SDU startup
											// skip BFOs 16:29:17.407 to 16:29:41.907, which are ~5-10 Hz higher than the rest, possibly due to taxiing for take off

		double ts;			// time (s since 2014-03-07 00:00:00 UTC);
		double dist;		// used here only to utilize GetPingData (not used in the following calculations)
		double bto, bfo;	// BTO and BFO data
		int channelType;	// Channel type
		
		// Get i-th ping data
		pInmarsatData -> GetPingData(ts, dist, bto, bfo, channelType, i);
	
		switch(channelType)
		{
		case CH_TYPE_R_RX1200:
		case CH_TYPE_R_RX600:
			diff = dF[i]+BFO_Bias_RX-bfo;
			BFO_STD_RX += (diff*diff); // R-RX STD
			if (max_BFO_dev < abs(diff)) max_BFO_dev = abs(diff);
			break;

		case CH_TYPE_T_RX1200:
			diff = dF[i]+BFO_Bias_TX-bfo;
			BFO_STD_TX += (diff*diff); // T-RX STD
			if (max_BFO_dev < abs(diff)) max_BFO_dev = abs(diff);
			break;
		};

		// print output
		if (channelType == CH_TYPE_R_RX1200 || channelType == CH_TYPE_R_RX600 || channelType == CH_TYPE_T_RX1200)
		{
			print_time_formatted(buftimetxt, ts);

			if (fid)
			{
				fprintf(fid, "%s  %7.2f  %s\n", &buftimetxt[11], diff);
			}
	
			if (verbose)
			{
				printf("%s  %7.2f  %s\n", &buftimetxt[11], diff);
			}
		}
	};

	// compute biases as averages
	BFO_STD = sqrt((BFO_STD_RX+BFO_STD_TX)/(double)(nBFO_Bias_RX+nBFO_Bias_TX));
	BFO_STD_RX = sqrt(BFO_STD_RX/(double)nBFO_Bias_RX);
	BFO_STD_TX = sqrt(BFO_STD_TX/(double)nBFO_Bias_TX);

	if (fid)
	{
		fprintf(fid, "\n");
		fprintf(fid, "****************************************************\n");
		fprintf(fid, "\n");
		fprintf(fid, "BFO deviations:\n");
		fprintf(fid, "R-Channel RX 1200 BFO STD = %6.2f (Hz) - %d samples\n", BFO_STD_RX, nBFO_Bias_RX);
		fprintf(fid, "T-Channel RX 1200 BFO STD = %6.2f (Hz) - %d samples\n", BFO_STD_TX, nBFO_Bias_TX);
		fprintf(fid, "Overall BFO STD = %6.2f (Hz) - %d samples\n", BFO_STD, nBFO_Bias);
		fprintf(fid, "Maximum BFO deviation = %6.2f (Hz)\n", max_BFO_dev);
		fprintf(fid, "\n");
		fprintf(fid, "****************************************************\n");
	}
	
	if (verbose)
	{
		printf("\n");
		printf("****************************************************\n");
		printf("\n");
		printf("BFO deviations:\n");
		printf("R-Channel RX 1200 BFO STD = %6.2f (Hz) - %d samples\n", BFO_STD_RX, nBFO_Bias_RX);
		printf("T-Channel RX 1200 BFO STD = %6.2f (Hz) - %d samples\n", BFO_STD_TX, nBFO_Bias_TX);
		printf("Overall BFO STD = %6.2f (Hz) - %d samples\n", BFO_STD, nBFO_Bias);
		printf("Maximum BFO deviation = %6.2f (Hz)\n", max_BFO_dev);
		printf("\n");
		printf("****************************************************\n");
	}


// -------------------------------------------------------------------------------------------------

	// Method #2 - Analyse BFOs by clusters  (see pings_data.cpp)
	// Time interval                 IDX (from-to)
	// 15:59:55.413 - 16:00:23.407      0 -   7
	// 16:00:27.741 - 16:00:32.406      8 -  14
	// 16:01:16.406 - 16:01:28.906     15 -  25
	// 16:06:34.906 - 16:07:21.908     26 -  57
	// 16:09:37.407 - 16:09:46.907     58 -  68
	// 16:11:03.907 - 16:11:13.408     69 -  76
	// 16:27:59.407 - 16:28:15.909     77 -  89
	// 16:29:17.407 - 16:29:41.907     90 - 122
	// 16:29:46.242 - 16:29:52.406    123 - 138


	// Define groups
	const int ngroups = 9;
	int gstart[] = {0, 8,  15, 26, 58, 69, 77, 90,  123};	// start index
	int gend[] =   {7, 14, 25, 57, 68, 76, 89, 122, 138};	// end index

	for (int ig = 0; ig<9; ig++)
	{
		for (int i=gstart[ig]; i<=gend[ig]; i++) idx_group[i] = ig;
	}


	// accumulators and counters
	double group_BFO_TX_biases[ngroups];
	int group_nBFO_TX_biases[ngroups];
	double group_BFO_RX_biases[ngroups];
	int group_nBFO_RX_biases[ngroups];
	double group_BFO_mean_biases[ngroups];
	int group_nBFO_mean_biases[ngroups];

	// iterate by groups:
	for (int ig = 0; ig<9; ig++)
	{
		group_BFO_TX_biases[ig] = 0.0;
		group_BFO_RX_biases[ig] = 0.0;
		group_BFO_mean_biases[ig] = 0.0;

		group_nBFO_TX_biases[ig] = 0;
		group_nBFO_RX_biases[ig] = 0;
		group_nBFO_mean_biases[ig] = 0;

		for (int i=0; i<N_TARMAC_PINGS; i++)
		{
			if (idx_group[i]!=ig) continue;		// skip if index is not within ig-group

			double ts;			// time (s since 2014-03-07 00:00:00 UTC);
			double dist;		// used here only to utilize GetPingData (not used in the following calculations)
			double bto, bfo;	// BTO and BFO data
			int channelType;	// Channel type
		
			// Get i-th ping data
			pInmarsatData -> GetPingData(ts, dist, bto, bfo, channelType, i);
	
			group_BFO_mean_biases[ig] += (bfo - dF[i]); // overall group bias
			group_nBFO_mean_biases[ig]++;

			switch(channelType)
			{
			case CH_TYPE_R_RX1200:
			case CH_TYPE_R_RX600:
				group_BFO_RX_biases[ig] += (bfo - dF[i]);
				group_nBFO_RX_biases[ig]++;
				break;

			case CH_TYPE_T_RX1200:
				group_BFO_TX_biases[ig] += (bfo - dF[i]);
				group_nBFO_TX_biases[ig]++;
				break;
			};
		};

		// compute biases as averages
		if (group_nBFO_mean_biases[ig]>0) group_BFO_mean_biases[ig] = group_BFO_mean_biases[ig] / (double)group_nBFO_mean_biases[ig];
		if (group_nBFO_RX_biases[ig]>0) group_BFO_RX_biases[ig] = group_BFO_RX_biases[ig] / (double)group_nBFO_RX_biases[ig];
		if (group_nBFO_TX_biases[ig]>0) group_BFO_TX_biases[ig] = group_BFO_TX_biases[ig] / (double)group_nBFO_TX_biases[ig];
	} //ig - group counter


	// Print results
	char textbufs[256];	// start time
	char textbufe[256];	// end time

	if (fid)
	{
		fprintf(fid, "\n");
		fprintf(fid, "Method 2 - Grouping BFOs by clusters\n");
		fprintf(fid, "\n");
		fprintf(fid, "Cluster, Interval, Count-R, Count-T, Count, Bias-R(Hz), Bias-T(Hz), Bias(Hz)\n");
	}
	
	if (verbose)
	{
		printf("\n");
		printf("Method 2 - Grouping BFOs by clusters\n");
		printf("\n");
		printf("Cluster, Interval, Count-R, Count-T, Count, Bias-R(Hz), Bias-T(Hz), Bias(Hz)\n");
	}

	for (int ig=0; ig<ngroups; ig++)
	{
		double ts;
		pInmarsatData -> GetPingData(ts, gstart[ig]);
		print_time_formatted(textbufs, ts);		// cluster start time
		pInmarsatData -> GetPingData(ts, gend[ig]);	// cluster end time
		print_time_formatted(textbufe, ts);

		if (group_nBFO_RX_biases[ig] == 0)
		{
			if (fid)
			{
				fprintf(fid, "%d  %s - %s  %2d  %2d  %2d  ******  %6.2f  %6.2f\n", ig, &textbufs[11], &textbufe[11],
						group_nBFO_RX_biases[ig], group_nBFO_TX_biases[ig], group_nBFO_mean_biases[ig],
						group_BFO_TX_biases[ig], group_BFO_mean_biases[ig]);
			}
	
			if (verbose)
			{
				printf("%d  %s - %s  %2d  %2d  %2d  ******  %6.2f  %6.2f\n", ig, &textbufs[11], &textbufe[11],
						group_nBFO_RX_biases[ig], group_nBFO_TX_biases[ig], group_nBFO_mean_biases[ig],
						group_BFO_TX_biases[ig], group_BFO_mean_biases[ig]);
			}
		}

		if (group_nBFO_TX_biases[ig] == 0)
		{
			if (fid)
			{
				fprintf(fid, "%d  %s - %s  %2d  %2d  %2d  %6.2f  ******  %6.2f\n", ig, &textbufs[11], &textbufe[11],
						group_nBFO_RX_biases[ig], group_nBFO_TX_biases[ig], group_nBFO_mean_biases[ig],
						group_BFO_RX_biases[ig], group_BFO_mean_biases[ig]);
			}
	
			if (verbose)
			{
				printf("%d  %s - %s  %2d  %2d  %2d  %6.2f  ******  %6.2f\n", ig, &textbufs[11], &textbufe[11],
						group_nBFO_RX_biases[ig], group_nBFO_TX_biases[ig], group_nBFO_mean_biases[ig],
						group_BFO_RX_biases[ig], group_BFO_mean_biases[ig]);
			}
		}

		if (group_nBFO_RX_biases[ig] != 0 && group_nBFO_TX_biases[ig] != 0)
		{
			if (fid)
			{
				fprintf(fid, "%d  %s - %s  %2d  %2d  %2d  %6.2f  %6.2f  %6.2f\n", ig, &textbufs[11], &textbufe[11],
						group_nBFO_RX_biases[ig], group_nBFO_TX_biases[ig], group_nBFO_mean_biases[ig],
						group_BFO_RX_biases[ig], group_BFO_TX_biases[ig], group_BFO_mean_biases[ig]);
			}
	
			if (verbose)
			{
				printf("%d  %s - %s  %2d  %2d  %2d  %6.2f  %6.2f  %6.2f\n", ig, &textbufs[11], &textbufe[11],
						group_nBFO_RX_biases[ig], group_nBFO_TX_biases[ig], group_nBFO_mean_biases[ig],
						group_BFO_RX_biases[ig], group_BFO_TX_biases[ig], group_BFO_mean_biases[ig]);
			}
		}
	}

	// average over groups
	double group_BFO_RX_bias_avg = 0.0;
	double group_BFO_TX_bias_avg = 0.0;
	double group_BFO_mean_bias_avg = 0.0;
	int group_nBFO_RX_bias_avg = 0;
	int group_nBFO_TX_bias_avg = 0;
	int group_nBFO_mean_bias_avg = 0;
	for (int ig=0; ig<ngroups; ig++)
	{
		if (ig==0 || ig==7) continue; // skip groups 0 and 7

		if (group_nBFO_RX_biases[ig]>0)
		{
			group_BFO_RX_bias_avg += group_BFO_RX_biases[ig];
			group_nBFO_RX_bias_avg++;
		}

		if (group_nBFO_TX_biases[ig]>0)
		{
			group_BFO_TX_bias_avg += group_BFO_TX_biases[ig];
			group_nBFO_TX_bias_avg++;
		}

		if (group_nBFO_mean_biases[ig]>0)
		{
			group_BFO_mean_bias_avg += group_BFO_mean_biases[ig];
			group_nBFO_mean_bias_avg++;
		}
	}

	group_BFO_RX_bias_avg = group_BFO_RX_bias_avg/(double)group_nBFO_RX_bias_avg;
	group_BFO_TX_bias_avg = group_BFO_TX_bias_avg/(double)group_nBFO_TX_bias_avg;
	group_BFO_mean_bias_avg = group_BFO_mean_bias_avg/(double)group_nBFO_mean_bias_avg;

	if (fid)
	{
		fprintf(fid,"\n");
		fprintf(fid, "Averaged by clusters, except nos. 0 and 7:  %6.2f  %6.2f  %6.2f\n", group_BFO_RX_bias_avg, group_BFO_TX_bias_avg, group_BFO_mean_bias_avg);
	}
	
	if (verbose)
	{
		printf("\n");
		printf("Averaged by clusters, except nos. 0 and 7:  %6.2f  %6.2f  %6.2f\n", group_BFO_RX_bias_avg, group_BFO_TX_bias_avg, group_BFO_mean_bias_avg);
	}

	if (fid) fclose(fid);

	delete pInmarsatData;
}
