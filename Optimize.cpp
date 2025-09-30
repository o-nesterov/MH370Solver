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
#include "params_parser.h"
#include "environdata.h"
#include "inmarsatdata.h"
#include "trajectory.h"

// Parameters to be passed
struct TrajectoryParams
{
	CInpParams* pInpParams;
	CMH370InmarsatData* pInmarsatData;
	CMH370EnvironData* pMH370EnvironData;
	CTrajectory* pTrajectory;
	int* idxPings;
	unsigned long nPings;
	__int32 OptFlag;
	__int32 OptCriteria;
	double opt_weight1;
	double opt_weight2;
	FILE* messagestream;
};

// Define optimization bits in the optimization flag (variables to optimize)
#define OPTIMIZE_WGT 1
#define OPTIMIZE_LON 2
#define OPTIMIZE_LAT 4
#define OPTIMIZE_HDG 8
#define OPTIMIZE_SPD 16
#define OPTIMIZE_ALT 32
#define OPTIMIZE_BNK 64
#define OPTIMIZE_BNF 128
#define OPTIMIZE_GY1 256
#define OPTIMIZE_GY2 512
#define OPTIMIZE_DEL 1024
#define OPTIMIZE_LNV 2048

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Main optimization function called from main(). It minimizes target functional.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Declare external functions
extern void BuildTrajectory(CInpParams* pInpParams, CTrajectory* pTrajectory, CMH370EnvironData* pMH370EnvironData, CMH370InmarsatData* pInmarsatData);
extern bool isPingIndexInList(CInpParams* pInpParams, int idx);
extern void print_time_formatted(char* buf, double t);


// Declare optimization functions
CTrajectory*  optimize(CInpParams* pInpParams, CMH370InmarsatData* pInmarsatData, CMH370EnvironData* pMH370EnvironData, 
			  __int32 OptCriteria, __int32 OptFlag, double opt_weight1, double opt_weight2, FILE* messagestream);
CTrajectory*  optimize_weight(CInpParams* pInpParams, CMH370InmarsatData* pInmarsatData, CMH370EnvironData* pMH370EnvironData, CTrajectory* pTrajectory, FILE* messagestream);



void Optimize(char* optname, char* inpfilename, char* outfilename, char* str_optweight1, char* str_optweight2, bool isVerbose)
{
	// set output message stream
	FILE* messagestream = NULL;
	if (isVerbose) messagestream = stdout;

	// Check type of optimization given by the first symbol in optname: 0, 1, 2, m
	__int32 OptCriteria = -1; // optimization criteria type
	switch (optname[0])
	{
	case 'm':
		OptCriteria = 0; // only weight optimization (-optm was called)
		break;
	case '0':
		OptCriteria = 1; // sum: normalized RMS bto + normalized RMS bfo
		break;
	case '1':
		OptCriteria = 2; // max: (normalized RMS bto, normalized RMS bfo)
		break;
	case '2':
		OptCriteria = 3; // max: (normalized bto, normalized bfo)
		break;
	}

	if (OptCriteria<0)
	{
		if (isVerbose) printf("Error: unrecognized optimization option...\n");
		return;
	}

	// --------------------
	// Build optimization flag
	__int32 OptFlag = 0;

	if (OptCriteria>0)
	{
		if (strstr(optname,"m,") || strstr(optname,",m") || strcmp(optname,"m")==0) OptFlag = OptFlag | OPTIMIZE_WGT;			// optimize weight
		if (strstr(optname,"x,") || strstr(optname,",x") || strcmp(optname,"x")==0) OptFlag = OptFlag | OPTIMIZE_LON;			// optimize initial longitude
		if (strstr(optname,"y,") || strstr(optname,",y") || strcmp(optname,"y")==0) OptFlag = OptFlag | OPTIMIZE_LAT;			// optimize initial latitude
		if (strstr(optname,"h,") || strstr(optname,",h") || strcmp(optname,"h")==0) OptFlag = OptFlag | OPTIMIZE_HDG;			// optimize initial heading
		if (strstr(optname,"s,") || strstr(optname,",s") || strcmp(optname,"s")==0) OptFlag = OptFlag | OPTIMIZE_SPD;			// optimize initial IAS/MACH
		if (strstr(optname,"g1,") || strstr(optname,",g1") || strcmp(optname,"g1")==0) OptFlag = OptFlag | OPTIMIZE_GY1;		// optimize gyroscopic parameter 1 (phase) 
		if (strstr(optname,"g2,") || strstr(optname,",g2") || strcmp(optname,"g2")==0) OptFlag = OptFlag | OPTIMIZE_GY2;		// optimize gyroscopic parameter 2 (amplitude)
		if (strstr(optname,"a,") || strstr(optname,",a") || strcmp(optname,"a")==0) OptFlag = OptFlag | OPTIMIZE_ALT;			// optimize altitude
		if (strstr(optname,"b,") || strstr(optname,",b") || strcmp(optname,"b")==0) OptFlag = OptFlag | OPTIMIZE_BNK;			// optimize default bank angle 
		if (strstr(optname,"bb,") || strstr(optname,",bb") || strcmp(optname,"bb")==0) OptFlag = OptFlag | OPTIMIZE_BNF;		// optimize default bank angle and all explicitly specified bank angles in lnav profile, which are not subjected to individual optimization (marked with *), forcing them to become default bank angles
		if (strstr(optname,"d,") || strstr(optname,",d") || strcmp(optname,"d")==0) OptFlag = OptFlag | OPTIMIZE_DEL;			// optimize AED delay parameter
		if (strstr(optname,"lnav,") || strstr(optname,",lnav") || strcmp(optname,"lnav")==0) OptFlag = OptFlag | OPTIMIZE_LNV;	// optimize all parameters marked with * in lnav profile file
	}

	if (OptCriteria==0) OptFlag = OPTIMIZE_WGT;



	// --------------------
	// Read-in setup file
	CInpParams* pInpParams = new CInpParams();	// input parameters
	if (isVerbose) pInpParams -> SetVerbose(stdout);
	if (!pInpParams -> parseParams(inpfilename))
	{
		if (isVerbose) printf("Error reading starting input file...\n");
		delete pInpParams;
		return;
	}

	if ((!(pInpParams->HDG_MODE == HEADING_MODE_GYROHDG)) && (OptFlag & OPTIMIZE_GY1 || OptFlag & OPTIMIZE_GY2))
	{
		if (isVerbose) printf("Heading mode is incompatible with optimization parameters. Ignoring those.\n");
		// clear gyroscopic optimization
		OptFlag = OptFlag & (~OPTIMIZE_GY1);
		OptFlag = OptFlag & (~OPTIMIZE_GY2);
	}

	if (OptFlag==0)
	{
		if (isVerbose) printf("Nothing to optimize... Wrong optimization option or wrong combination?\n");
		delete pInpParams;
		return;
	}

	if (OptFlag!=OPTIMIZE_WGT && (pInpParams->nOptPingsIDX==0 || pInpParams->pOptPingsIDX==NULL))
	{
		if (isVerbose) printf("Error: ping indices to be used for optimization must be specified: OPT_PING_IDX = {..,..,..}. Use -printpings command to list available ones.\n");
		delete pInpParams;
		return;
	}

	if (isVerbose) printf("Optimizing...\n");


	// ---------------------
	// optimization weight coefficients
	double opt_weight1 = -999.0;
	double opt_weight2 = -999.0;

	char* pp;
	if (strlen(str_optweight1)>0)
	{
		double tmp = strtod(str_optweight1, &pp);
		if (pp!=&str_optweight1[0]) opt_weight1 = tmp;
	}
	if (strlen(str_optweight2)>0)
	{
		double tmp = strtod(str_optweight2, &pp);
		if (pp!=&str_optweight2[0]) opt_weight2 = tmp;
	}

	// treat all bank angles explicitly specified in lnav, which are not marked with *, as default and optimize them too 
	if (OptFlag & OPTIMIZE_BNF)
	{
		pInpParams->force_default_bank_angle = true;
	}
	else
	{
		pInpParams->force_default_bank_angle = false;
	}

	// ---------------------
	// Inmarsat data
	CMH370InmarsatData* pInmarsatData = new CMH370InmarsatData();
	if (isVerbose) pInmarsatData -> SetVerbose(stdout);
	// set parameters, which may affect Doppler calculation and AES Doppler compensation term
	if (pInpParams->DopplerMethod != DOPPLER_METHOD_UDF) pInmarsatData -> SetDopplerMethod(pInpParams->DopplerMethod);
	if (pInpParams->AES_EarthModel == AES_EARTH_MODEL_SPH || pInpParams->AES_EarthModel == AES_EARTH_MODEL_WGS)  pInmarsatData -> SetAESEarthModel(pInpParams->AES_EarthModel);
	if (pInpParams->AES_DopplerMethod != DOPPLER_METHOD_UDF)  pInmarsatData -> SetAESDopplerMethod(pInpParams->AES_DopplerMethod);
	if (pInpParams->AES_Earth_Radius > 0.0) pInmarsatData -> SetAESEarthRadius(pInpParams->AES_Earth_Radius);
	if (pInpParams->AES_GeoSat_Alt > 0.0) pInmarsatData -> SetAESGeoSatAlt(pInpParams->AES_GeoSat_Alt);


	// ---------------------
	// Trajectory
	CTrajectory* pTrajectory = NULL;

	// ---------------------
	// Load environmental data:
	CMH370EnvironData* pMH370EnvironData = new CMH370EnvironData();
	if (isVerbose) pMH370EnvironData -> SetVerbose(stdout,1);
	pMH370EnvironData -> LoadMeteo(pInpParams->METEO_DATASET, pInpParams->pMeteoDir);
	pMH370EnvironData -> SetMeteoInterpMethod(pInpParams->MeteoInterpMethod);
	pMH370EnvironData -> LoadMagnetic(pInpParams->MAGN_DECL_DATASET, pInpParams->pMagnDeclDir);
	pMH370EnvironData -> SetDeclInterpMethod(pInpParams->MagnDeclMethod);

	// ---------------------
	// optimizize trajectory options
	if (OptFlag & (~OPTIMIZE_WGT))
	{
		pTrajectory = optimize(pInpParams, pInmarsatData, pMH370EnvironData, OptCriteria, OptFlag, opt_weight1, opt_weight2, messagestream);
	}

	// ---------------------
	// Optimize initial weight
	if (OptFlag & OPTIMIZE_WGT)
	{
		pTrajectory = optimize_weight(pInpParams, pInmarsatData, pMH370EnvironData, pTrajectory, messagestream);
	}

	if (isVerbose && pTrajectory!=NULL) printf("Optimization completed.\n");

	if (pTrajectory==NULL)
	{
		// if optimization failed for whatever reason
		delete pMH370EnvironData;
		delete pInmarsatData;
		delete pInpParams;
		return;
	}

	// ---------------------
	// Save configuration file and lateral navigation file with the optimized parameters, as needed
	// if lnav optimization was done, alter lateral navigation file name by appending "_opt.txt"

	if (OptFlag & OPTIMIZE_LNV)
	{
		char* pNewlnavfilename = NULL;
		if (pInpParams->pLNAVfilename)
		{
			pNewlnavfilename = new char[strlen(pInpParams->pLNAVfilename)+16];
			strcpy(pNewlnavfilename, pInpParams->pLNAVfilename);

			if (strlen(pInpParams->pLNAVfilename)>0)
			{
				// find last slash or last back slash back slash  (folder name may contain '.', whilst filename may come without extension)
				char* pp = &pNewlnavfilename[strlen(pNewlnavfilename)-1];
				char* pt = pp;
				while (pp>&pNewlnavfilename[0] && *pp!='\\' && *pp!='/') pp--;
				// find the last '.' after the last slash or backslash
				while (pt>pp && *pt!='.') pt--;
				if (pp!=pt)
				{
					sprintf(pt, "_opt.txt");
				}
				else
				{
					sprintf(pNewlnavfilename, "%s_opt.txt", pInpParams->pLNAVfilename);
				}
			}
		}

		// replace old lnav filename with new lnav filename
		if (pInpParams->pLNAVfilename) delete pInpParams->pLNAVfilename;
		pInpParams->pLNAVfilename = pNewlnavfilename;

		// Save configuration file with the optimized parameters
		pInpParams -> saveParams(outfilename);
		
		if (pNewlnavfilename)
		{
			if (strlen(pNewlnavfilename)>0)
			{
				if (!pTrajectory->SaveLNAVProfile(pNewlnavfilename))
				{
					if (isVerbose) printf("Error saving lateral navigation file: %s\n", pNewlnavfilename);
				}
			}
		}
	}
	else
	{
		// Save only configuration file with the optimized parameters
		pInpParams -> saveParams(outfilename);
	}


	// ---------------------
	// Now reset output times instead of optimization times
	unsigned long len = (unsigned long) (pInpParams->nOutputPingIDX);
	int* Pings_idx_subset = pInpParams->pOutputPingIDX;
	double* pSpecialOutputTimes = new double[len+1];
	for (unsigned long i=0; i<len; i++)
	{
		double ts;
		pInmarsatData->GetPingData(ts, Pings_idx_subset[i]);
		pSpecialOutputTimes[i] = ts;
	};
	pTrajectory->SetSpecificOutputTimes(pSpecialOutputTimes, len, pInpParams->AES_Delay);
	delete pSpecialOutputTimes;

	// ------------
	// Build trajectory using optimized parameters and already exisiting objects, so file parsing is not required 
	BuildTrajectory(pInpParams, pTrajectory, pMH370EnvironData, pInmarsatData);

	// ------------
	// delete all objects (free memory)
	delete pTrajectory;
	delete pMH370EnvironData;
	delete pInmarsatData;
	delete pInpParams;

}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Weight optimization by forward integration (can be reverse integration, knowing trajectory and final weight )
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


double func2min_weight(unsigned long n, double* args, void* params)
{
	// Final weight (without fuel)
	const double WGT_Final = 174369.0; // (kg); according to FI; == 223200.0-49200.0 kg

	CInpParams* pInpParams = ((TrajectoryParams*)params) -> pInpParams;
	CMH370InmarsatData* pInmarsatData = ((TrajectoryParams*)params) -> pInmarsatData;
	CMH370EnvironData* pMH370EnvironData = ((TrajectoryParams*)params) -> pMH370EnvironData;
	CTrajectory* pTrajectory = ((TrajectoryParams*)params) -> pTrajectory;
	int* idxPings = ((TrajectoryParams*)params) -> idxPings;
	unsigned long nPings = ((TrajectoryParams*)params) -> nPings;
	FILE* messagestream = ((TrajectoryParams*)params) -> messagestream;

	// Some of the following parameters can be changed during optimization
	double t0 = pInpParams->ts;
	double t1 = pInpParams->te;
	int spd_mode = pInpParams->SPD_MODE;
	int hdg_mode = pInpParams->HDG_MODE;

	double lon0 = pInpParams->lon0;
	double lat0 = pInpParams->lat0;
	double alt0 = pInpParams->alt0;		// altitude (m)
	double M0 = pInpParams->wgt0;		// weight (kg)
	double hdg0 = pInpParams->hdg0;		// deg
	double spd0 = pInpParams->spd0;		// IAS m/s
	int hdg0_type = pInpParams->hdg0_type;	// heading value type: true, magnetic or gyroscopic

	double gyro_hdg_param1 = pInpParams->gyro_hdg_param1;	// gyroscopic heading parameters
	double gyro_hdg_param2 = pInpParams->gyro_hdg_param2;

	// replace optimization parameter - weight
	M0 = args[0];


	// compute trajectory
	pTrajectory->BuildTrajectory(t0, t1, lon0, lat0, alt0, M0, spd_mode, spd0, hdg_mode, hdg0, hdg0_type, pMH370EnvironData, gyro_hdg_param1, gyro_hdg_param2);

//	if (hdg_mode == HEADING_MODE_GYROHDG)
//	{
//		pTrajectory->BuildTrajectory(t0, t1, lon0, lat0, alt0, M0, spd_mode, spd0, hdg_mode, hdg0, pMH370EnvironData, gyro_hdg_param1, gyro_hdg_param2);
//	}
//	else
//	{
//		pTrajectory->BuildTrajectory(t0, t1, lon0, lat0, alt0, M0, spd_mode, spd0, hdg_mode, hdg0, pMH370EnvironData, 0.0, 0.0);
//	}

	double M_final = (pTrajectory->GetMOut())[pTrajectory->GetNOut()-1];

	double val = (M_final-WGT_Final)*(M_final-WGT_Final);
	if (messagestream) fprintf(messagestream, "%f\n",val);
	return val;
};





CTrajectory* optimize_weight(CInpParams* pInpParams, CMH370InmarsatData* pInmarsatData, CMH370EnvironData* pMH370EnvironData, CTrajectory* pTrajectoryIn, FILE* messagestream = NULL)
{

	// Trajectory
	CTrajectory* pTrajectory = NULL;
	
	if (pTrajectoryIn)
	{
		pTrajectory = pTrajectoryIn;
	}
	else
	{
		pTrajectory = new CTrajectory();
		if (messagestream) pTrajectory->SetVerbose(messagestream);

		if(!pTrajectory->CreateLNAVProfile(pInpParams))
		{
			if (messagestream) fprintf(messagestream, "Optimization failed. Missing or wrong LNAV file %s.\n", pInpParams->pLNAVfilename);
			delete pTrajectory;
			return NULL;
		}
		if (pInpParams->nLNAVOptParams>0) pTrajectory->ResetLNAVOptimizationParameters(pInpParams->LNAVOptParams, pInpParams->nLNAVOptParams);	// reset LNAV optimization parameters
		pTrajectory->SetDefaultBankAngle(pInpParams->bank_angle, pInpParams->force_default_bank_angle);
		if (pInpParams->ENG_MODE == ENGINE_MODE_SINGLE) pTrajectory->SetEngineModeSingle();
		if (pInpParams->ENG_MODE == ENGINE_MODE_DUAL) pTrajectory->SetEngineModeDual();

		// ------------------------
		// set special output times
		//	unsigned long Pings_idx_subset[] = {564, 565, 566, 567, 595, 596, 597}; 
		//	unsigned long Pings_idx_subset[] = {563, 564, 565, 566, 567, 595, 596, 597}; 
		//	unsigned long Pings_idx_subset[] = {499, 501, 502, 503, 504, 505, 511, 512, 553, 562, 563};
		//	unsigned long Pings_idx_subset[] = {499, 501, 502, 503, 504, 505, 511, 512, 553, 562, 563, 564, 565, 566, 567, 595, 596, 597};

		//	const unsigned long len = sizeof(Pings_idx_subset)/sizeof(unsigned long);

		unsigned long len = (unsigned long) (pInpParams->nOptPingsIDX);
		int* Pings_idx_subset = pInpParams->pOptPingsIDX;
		double* pSpecialOutputTimes = new double[len+1];
		for (unsigned long i=0; i<len; i++)
		{
			double ts;
			pInmarsatData->GetPingData(ts, Pings_idx_subset[i]);
			pSpecialOutputTimes[i] = ts;
		};

		pTrajectory->SetSpecificOutputTimes(pSpecialOutputTimes, len, pInpParams->AES_Delay);
		delete pSpecialOutputTimes;
	}

	if (pInpParams->ENG_MODE == ENGINE_MODE_SINGLE) pTrajectory->SetEngineModeSingle();
	if (pInpParams->ENG_MODE == ENGINE_MODE_DUAL) pTrajectory->SetEngineModeDual();


	// ------------------------

	TrajectoryParams TrParams;
	TrParams.pInpParams = pInpParams;
	TrParams.pInmarsatData = pInmarsatData;
	TrParams.pMH370EnvironData = pMH370EnvironData;
	TrParams.pTrajectory = pTrajectory;
	TrParams.nPings = (unsigned long) (pInpParams->nOptPingsIDX);
	TrParams.idxPings = &(pInpParams->pOptPingsIDX)[0];
	TrParams.messagestream = messagestream;

	// ------------------------

	// Some of the following parameters can be changed during optimization
	double t0 = pInpParams->ts;
	double t1 = pInpParams->te;
	int spd_mode = pInpParams->SPD_MODE;
	int hdg_mode = pInpParams->HDG_MODE;

	double lon0 = pInpParams->lon0;
	double lat0 = pInpParams->lat0;
	double alt0 = pInpParams->alt0;		// altitude (m)
	double M0 = pInpParams->wgt0;		// weight (kg)
	double hdg0 = pInpParams->hdg0;		// deg
	double spd = pInpParams->spd0;		// IAS m/s


	// ------------------------

	double args_in[1];
	double args_out[1];
	double tol_args[1]; //tolerance;

	FMINSearchOptions options;
	FMINSearchStats stats;

	args_in[0] = M0;

	tol_args[0] = 1.0;

	options.MaxFunEvals = 500;
	options.MaxIter = 500;
	options.outstream = stdout;
	options.pWorkspace = NULL;
	options.TolFun = 1.0;
	options.TolX = &tol_args[0];
	options.outstream = messagestream;

	// ------------------------
	bool r = fminsearch(func2min_weight, 1, args_in, args_out, (void*)(&TrParams), options, stats);

	if (messagestream) fprintf(messagestream, "Initial weight (kg): %f; func = %f\n", args_out[0], stats.fVal);

	// save weight
	pInpParams->wgt0 = args_out[0];

	return pTrajectory;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Final leg (19:41 -...) parameters optimization, TRUE HDG, TRUE TRK, MAG HDG, MAG TRK modes
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


double func2min(unsigned long n, double* args, void* params)
{
	CInpParams* pInpParams = ((TrajectoryParams*)params) -> pInpParams;
	CMH370InmarsatData* pInmarsatData = ((TrajectoryParams*)params) -> pInmarsatData;
	CMH370EnvironData* pMH370EnvironData = ((TrajectoryParams*)params) -> pMH370EnvironData;
	CTrajectory* pTrajectory = ((TrajectoryParams*)params) -> pTrajectory;
	int* idxPings = ((TrajectoryParams*)params) -> idxPings;
	unsigned long nPings = ((TrajectoryParams*)params) -> nPings;
	FILE* messagestream = ((TrajectoryParams*)params) -> messagestream;

	__int32 OptCriteria = ((TrajectoryParams*)params) -> OptCriteria;
	__int32 OptFlag = ((TrajectoryParams*)params) -> OptFlag;
	double opt_weight1 = ((TrajectoryParams*)params) -> opt_weight1;
	double opt_weight2 = ((TrajectoryParams*)params) -> opt_weight2;

	// Some of the following parameters can be changed during optimization
	double t0 = pInpParams->ts;
	double t1 = pInpParams->te;
	int spd_mode = pInpParams->SPD_MODE;
	int hdg_mode = pInpParams->HDG_MODE;

	double lon0 = pInpParams->lon0;		// longitude (deg)
	double lat0 = pInpParams->lat0;		// latitude (deg)
	double alt0 = pInpParams->alt0;		// altitude (m)
	double M0 = pInpParams->wgt0;		// weight (kg)
	double hdg0 = pInpParams->hdg0;		// deg
	double spd0 = pInpParams->spd0;		// IAS m/s
	double bank_angle = pInpParams->bank_angle; // bank angle (deg)
	int hdg0_type = pInpParams->hdg0_type;	// initial heading value type: true, magnetic or gyroscopic

	double gyro_hdg_param1 = pInpParams->gyro_hdg_param1;	// gyroscopic heading parameters
	double gyro_hdg_param2 = pInpParams->gyro_hdg_param2;

	// replace optimization parameters
	unsigned long num_args = 0;
	if (OptFlag & OPTIMIZE_LON) {lon0 = args[num_args]; num_args++;};
	if (OptFlag & OPTIMIZE_LAT) {lat0 = args[num_args]; num_args++;};
	if (OptFlag & OPTIMIZE_HDG) {hdg0 = args[num_args]; num_args++;};
	if (OptFlag & OPTIMIZE_SPD) {spd0 = args[num_args]; num_args++;};
	if (OptFlag & OPTIMIZE_ALT) {alt0 = args[num_args]; num_args++;};
	if ((OptFlag & OPTIMIZE_BNK) || (OptFlag & OPTIMIZE_BNF)) {pTrajectory->SetDefaultBankAngle(args[num_args], pInpParams->force_default_bank_angle); num_args++;};
	if (OptFlag & OPTIMIZE_GY1) {gyro_hdg_param1 = args[num_args]; num_args++;};
	if (OptFlag & OPTIMIZE_GY2) {gyro_hdg_param2 = args[num_args]; num_args++;};
	if (OptFlag & OPTIMIZE_DEL) 
	{
		// Reset special output times 
		double aes_delay = args[num_args]; 
		num_args++;

		double* pSpecialOutputTimes = new double [nPings+1];
		for (unsigned long i=0; i<nPings; i++)
		{
			double ts;
			pInmarsatData->GetPingData(ts, idxPings[i]);
			pSpecialOutputTimes[i] = ts;
		}
		pTrajectory->SetSpecificOutputTimes(pSpecialOutputTimes, nPings, aes_delay);
		delete pSpecialOutputTimes;
	};

	if (OptFlag & OPTIMIZE_LNV) 
	{
		// reset LNAV optimization parameters if those were recalculated in a result of optimization
		if (pInpParams->nLNAVOptParams>0) pTrajectory->ResetLNAVOptimizationParameters(&args[num_args], pInpParams->nLNAVOptParams);
		num_args += pInpParams->nLNAVOptParams;
	};

	// ---------------------------------------------------------
	// compute trajectory

	pTrajectory->BuildTrajectory(t0, t1, lon0, lat0, alt0, M0, spd_mode, spd0, hdg_mode, hdg0, hdg0_type, pMH370EnvironData, gyro_hdg_param1, gyro_hdg_param2);

	//  special output
	double diff_dist = 0.0;	// distance difference metric
	double diff_bfo = 0.0;	// bfo difference meteric

	double diff_dist_max = -1.0;
	double diff_bfo_max = -1.0;


	// measured BFO = 182 Hz at 2014-03-08 00:19:29.416 is unreliable - 2nd SDU reboot, likely descent as fuel ran out
	// measured BTO at 2014-03-08 00:19:29.416 occurs after the plane ran out of fuel. Whilst it is valid, it would not correspond the assumed mode of flight.
	// so compute metric excluding the last ping

	for (unsigned long i=0; i<nPings; i++)
	{
		double lon, lat, alt, M, u, v, w, q;									// actual position and velocity data
		double aes_lon, aes_lat, aes_alt, aes_M, aes_u, aes_v, aes_w, aes_q;	// delayed position and velocity data used for AES compsation term calculation

		pTrajectory->GetSpecificOutput(lon, lat, alt, M, u, v, w, q, i, false);
		pTrajectory->GetSpecificOutput(aes_lon, aes_lat, aes_alt, aes_M, aes_u, aes_v, aes_w, aes_q, i, true);	// get delayed position and velocity data

		double ts, dist, bto, bfo;
		int channelType;
		pInmarsatData->GetPingData(ts, dist, bto, bfo, channelType,  idxPings[i]);

		double dist_computed = pInmarsatData->CalcDistToSat(ts, lon, lat, alt);
//		double bto_computed = pInmarsatData->CalcBTO(ts, lon, lat, alt, channelType);
		double bfo_computed = pInmarsatData->CalcBFO(ts, lon, lat, alt, u, v, w, aes_lon, aes_lat, aes_alt, aes_u, aes_v, channelType);


		if (abs(bto)<1.0E+9) 
		{
			double diff = dist_computed-dist;
			diff_dist += diff*diff;
			if (diff_dist_max < abs(diff)) diff_dist_max = abs(diff);
		};

		if (abs(bfo)<1.0E+9)
		{
			double diff = bfo_computed-bfo;
			diff_bfo += diff*diff;
			if (diff_bfo_max < abs(diff)) diff_bfo_max = abs(diff);
		};
	}


	// Standard deviation for BTO R1200 is 29 microseconds (DSTG, 2015), p.27, which is equivalent to 4.35 km of satellite to airplane distance (BTO counts double distance)
	// Standard deviation for BFO is 4.3-5.5 Hz (DSTG, 2015), p. 31, Table 5.1
	// STD BTO in this work 4.17 km (after reaching FL350), and STD BFO 1.61 Hz (after reaching FL350, using 2 biases)
	// double a = 0.00023;	// normalization (weighting) coefficient, 1/m
	// double b = 0.2;		// normalization (weighting) coefficient, 1/Hz

	double val = 0.0; // minimization functional value
	double val1;

	switch(OptCriteria)
	{
	case 1:
		val = opt_weight1*sqrt(diff_dist) + opt_weight2*sqrt(diff_bfo);
		break;

	case 2:
		val =  opt_weight1*sqrt(diff_dist);
		val1 = opt_weight2*sqrt(diff_bfo);
		if (val1>val) val = val1;
		break;
		
	case 3:
		val  = opt_weight1*diff_dist_max;
		val1 = opt_weight2*diff_bfo_max;
		if (val1>val) val = val1;
		break;
	}

	if(messagestream) fprintf(messagestream,"%f\n",val);
	return val;
};







CTrajectory* optimize(CInpParams* pInpParams, CMH370InmarsatData* pInmarsatData, CMH370EnvironData* pMH370EnvironData, 
					  int OptCriteria, __int32 OptFlag, double opt_weight1, double opt_weight2, FILE* messagestream = NULL)
{

	// Trajectory
	CTrajectory* pTrajectory = new CTrajectory();
	if (messagestream) pTrajectory->SetVerbose(messagestream);

	if(!pTrajectory->CreateLNAVProfile(pInpParams))
	{
		if (messagestream) fprintf(messagestream, "Optimization failed. Missing or wrong LNAV file %s.\n", pInpParams->pLNAVfilename);
		delete pTrajectory;
		return NULL;
	}

	pTrajectory->SetDefaultBankAngle(pInpParams->bank_angle, pInpParams->force_default_bank_angle);
	if (pInpParams->ENG_MODE == ENGINE_MODE_SINGLE) pTrajectory->SetEngineModeSingle();
	if (pInpParams->ENG_MODE == ENGINE_MODE_DUAL) pTrajectory->SetEngineModeDual();


	// ------------------------
	// set special output times
//	unsigned long Pings_idx_subset[] = {564, 565, 566, 567, 595, 596, 597};
//	unsigned long Pings_idx_subset[] = {563, 564, 565, 566, 567, 595, 596, 597};
//	unsigned long Pings_idx_subset[] = {499, 501, 502, 503, 504, 505, 511, 512, 553, 562, 563};
//	unsigned long Pings_idx_subset[] = {499, 501, 502, 503, 504, 505, 511, 512, 553, 562, 563, 564, 565, 566, 567, 595, 596, 597};

	unsigned long len = (unsigned long) (pInpParams->nOptPingsIDX);
	int* Pings_idx_subset = pInpParams->pOptPingsIDX;
//	const unsigned long len = sizeof(Pings_idx_subset)/sizeof(unsigned long);
	double* pSpecialOutputTimes = new double[len+1];
	for (unsigned long i=0; i<len; i++)
	{
		double ts;
		pInmarsatData->GetPingData(ts, Pings_idx_subset[i]);
		pSpecialOutputTimes[i] = ts;
	};
	pTrajectory->SetSpecificOutputTimes(pSpecialOutputTimes, len, pInpParams->AES_Delay);
	delete pSpecialOutputTimes;


	// ------------------------

	TrajectoryParams TrParams;
	TrParams.pInpParams = pInpParams;
	TrParams.pInmarsatData = pInmarsatData;
	TrParams.pMH370EnvironData = pMH370EnvironData;
	TrParams.pTrajectory = pTrajectory;
	TrParams.nPings = len;
	TrParams.idxPings = &Pings_idx_subset[0];
	TrParams.messagestream = messagestream;

	// ------------------------
	TrParams.OptCriteria = OptCriteria;
	TrParams.OptFlag = OptFlag;
	TrParams.opt_weight1 = 0.0001;	// default value (1/m)
	TrParams.opt_weight2 = 0.2;		// default value (1/Hz)
	if (opt_weight1>=0.0) TrParams.opt_weight1 = opt_weight1;
	if (opt_weight2>=0.0) TrParams.opt_weight2 = opt_weight2;


	// ------------------------

	unsigned long num_args = 0;

	double args_in[32];
	double args_out[32];
	double tol_args[32]; //tolerance;

	FMINSearchOptions options;
	FMINSearchStats stats;
	
	int nLNAVparams = 0; // number of LNAV parameters to optimize

	if (OptFlag & OPTIMIZE_LON) {args_in[num_args] = pInpParams->lon0; tol_args[num_args] = 1.0E-4; num_args++;};
	if (OptFlag & OPTIMIZE_LAT) {args_in[num_args] = pInpParams->lat0; tol_args[num_args] = 1.0E-4; num_args++;};
	if (OptFlag & OPTIMIZE_HDG) {args_in[num_args] = pInpParams->hdg0; tol_args[num_args] = 1.0E-4; num_args++;};
	if (OptFlag & OPTIMIZE_SPD) {args_in[num_args] = pInpParams->spd0; tol_args[num_args] = 1.0E-4; num_args++;};
	if (OptFlag & OPTIMIZE_ALT) {args_in[num_args] = pInpParams->alt0; tol_args[num_args] = 1.0E-1; num_args++;};
	if (OptFlag & OPTIMIZE_BNK) {args_in[num_args] = pInpParams->bank_angle; tol_args[num_args] = 1.0E-4; num_args++;};
	if (OptFlag & OPTIMIZE_GY1) {args_in[num_args] = pInpParams->gyro_hdg_param1; tol_args[num_args] = 1.0E-4; num_args++;};
	if (OptFlag & OPTIMIZE_GY2) {args_in[num_args] = pInpParams->gyro_hdg_param2; tol_args[num_args] = 1.0E-4; num_args++;};
	if (OptFlag & OPTIMIZE_DEL) {args_in[num_args] = pInpParams->AES_Delay; tol_args[num_args] = 1.0E-4; num_args++;};
	if (OptFlag & OPTIMIZE_LNV) 
	{
		pTrajectory->GetLNAVOptimizationParameters(&args_in[num_args], nLNAVparams, 16); // 16 is max number of parameters
		pInpParams->nLNAVOptParams = nLNAVparams;
		for (int i=0; i<nLNAVparams; i++) tol_args[num_args+i] = 1.0E-4;
		num_args+=nLNAVparams;
	}

	options.MaxFunEvals = 1000;
	options.MaxIter = 1000;
	options.outstream = stdout;
	options.pWorkspace = NULL;
	options.TolFun = 1.0E-4;
	options.TolX = &tol_args[0];
	options.outstream = messagestream;

	if (num_args==0)
	{
		if (messagestream) 
		{
			fprintf(messagestream, "Nothing to optimize. Erroneous optimization settings?\n");
			if ((OptFlag & OPTIMIZE_LNV) && nLNAVparams==0) fprintf(messagestream, "LNAV optimization was specified, but no optimization parameters are marked in LNAV file.\n");
		}
		return pTrajectory;
	}

	// ------------------------
	bool r = fminsearch(func2min, num_args, args_in, args_out, (void*)(&TrParams), options, stats);

	if (messagestream) 
	{
		for (int i=0; i<num_args; i++) fprintf(messagestream, "%f ", args_out[i]);
		fprintf(messagestream, "fVal: %f\n", stats.fVal);
	}

	// Save new parameters
	num_args = 0;
	if (OptFlag & OPTIMIZE_LON) {pInpParams->lon0 = args_out[num_args]; num_args++;};
	if (OptFlag & OPTIMIZE_LAT) {pInpParams->lat0 = args_out[num_args]; num_args++;};
	if (OptFlag & OPTIMIZE_HDG) {pInpParams->hdg0 = args_out[num_args]; num_args++;};
	if (OptFlag & OPTIMIZE_SPD) {pInpParams->spd0 = args_out[num_args]; num_args++;};
	if (OptFlag & OPTIMIZE_ALT) {pInpParams->alt0 = args_out[num_args]; num_args++;};
	if (OptFlag & OPTIMIZE_BNK) {pInpParams->bank_angle = args_out[num_args]; num_args++;};
	if (OptFlag & OPTIMIZE_GY1) {pInpParams->gyro_hdg_param1 = args_out[num_args]; num_args++;};
	if (OptFlag & OPTIMIZE_GY2) {pInpParams->gyro_hdg_param2 = args_out[num_args]; num_args++;};
	if (OptFlag & OPTIMIZE_DEL) {pInpParams->AES_Delay = args_out[num_args]; num_args++;};
	if (OptFlag & OPTIMIZE_LNV) 
	{
		pInpParams->nLNAVOptParams = nLNAVparams;
		for (int i=0; i<nLNAVparams; i++) pInpParams->LNAVOptParams[i] = args_out[num_args+i];
		num_args += nLNAVparams;
	}

	return pTrajectory;
}




/*

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Distance minimization to create arcs
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct DistMinParams
{
	CMH370InmarsatData* pInmarsatData;
	double ts;
	double ang;
	double dist;
	double alt;
	FILE* messagestream;
};


double func2min_dist(unsigned long n, double* args, void* params)
{
	CMH370InmarsatData* pInmarsatData = (( DistMinParams*)params) -> pInmarsatData;
	double ts = ((DistMinParams*)params) -> ts;	// time
	double ang = ((DistMinParams*)params) -> ang; // fictitious angle parameter
	double alt = ((DistMinParams*)params) -> alt; // altitude
	double dist = ((DistMinParams*)params) -> dist; // distance to satellite
	FILE* messagestream = ((DistMinParams*)params) -> messagestream;

	double lon = args[0];
	double lat = args[1];

	double dist_comp = pInmarsatData->CalcDistToSat(ts,lon,lat,alt);

	double angd = sqrt((lon-64.5)*(lon-64.5)+lat*lat);

	double x0_fictitious = dist*cos((90.0-ang)*3.1415926535898/180.0);
	double y0_fictitious = dist*sin((90.0-ang)*3.1415926535898/180.0);

	double x1_fictitious = dist_comp*(lon-64.5)/angd;
	double y1_fictitious = dist_comp*lat/angd;

	double diff_x = x1_fictitious - x0_fictitious;
	double diff_y = y1_fictitious - y0_fictitious;

	double val = diff_x*diff_x + diff_y*diff_y;
	if (messagestream) fprintf(messagestream, "%f\n",val);
	return val;
};



void createPingRings()
{
	CMH370InmarsatData* pInmarsatData = new CMH370InmarsatData();
	pInmarsatData -> SetVerbose(stdout);

	double alt = 12000.0; //10668; //0.0; // altitude

	double dist_offset = 0.0; //299792500.0 * 1.0e-5; // equivalent to bto 20 microseconds  

	// ------------------------
	// arc subset 
	unsigned long Pings_idx_subset[] = {261, 365, 367, 499, 501, 502, 504, 507, 511, 563, 564, 565, 566, 596, 597}; 
	const unsigned long len = sizeof(Pings_idx_subset)/sizeof(unsigned long);

	double pingLat[36000];
	double pingLon[36000];

	DistMinParams distParams;

	for (unsigned long n=0; n<len; n++)
	{


		double ts,dist,bto,bfo;
		int channelType;
		pInmarsatData->GetPingData(ts,dist,bto,bfo,channelType, Pings_idx_subset[n]);

		distParams.pInmarsatData = pInmarsatData;
		distParams.dist = dist + dist_offset;
		distParams.ts = ts;
		distParams.alt = alt;
		distParams.messagestream = NULL; //stdout; //NULL;

		for (int i=0; i<36000; i++)
		{
			// ------------------------
			distParams.ang = 0.01*i;

			double args_in[2];
			double args_out[2];
			double tol_args[2]; //tolerance;

			FMINSearchOptions options;
			FMINSearchStats stats;

			if (i>0)
			{
				args_in[0] = pingLon[i-1];
				args_in[1] = pingLat[i-1];
			}
			else
			{
				args_in[0] = 100.0;
				args_in[1] = 10.0;
			}
			tol_args[0] = 0.00001;
			tol_args[1] = 0.00001;

			options.MaxFunEvals = 500;
			options.MaxIter = 500;
			options.outstream = stdout;
			options.pWorkspace = NULL;
			options.TolFun = 1.0E-10;
			options.TolX = &tol_args[0];
			options.outstream = NULL; //stdout;

			// ------------------------
			bool r = fminsearch(func2min_dist, 2, args_in, args_out, (void*)(&distParams), options, stats);

			pingLon[i] = args_out[0];
			pingLat[i] = args_out[1];
			printf("Longitude, deg: %12.5f. Latitude, deg: %12.5f. Fmin value: %f\n", args_out[0], args_out[1], stats.fVal);
		}

		char fname[1024];
		sprintf(fname,"C:\\MH370Solver\\PingRings\\arc_%d_alt%5.5d.txt", Pings_idx_subset[n], ((int)floor(alt)));
		FILE* fid = fopen(fname, "wt");

		for (int i=0; i<36000; i++)
		{
			fprintf(fid, "%12.5f %12.5f\n", pingLon[i], pingLat[i]);
		}

		fclose(fid);

	};

	delete pInmarsatData;
}

*/