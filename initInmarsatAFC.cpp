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


// Doppler shift calculation 
extern double doppler(double XYZs[3], double XYZr[3], double UVWs[3], double UVWr[3], double Frq, int Method);




////////////////////////////////////////////////////////////////////////////////////////
//
// Return Inmarsat's F3 satellite velocity in ECEF system
//
////////////////////////////////////////////////////////////////////////////////////////
void CMH370InmarsatData::initAFCLookupTable_()
{
//	double* pAFC_lookup_table_;			// lookup table
//	int numAFCLookupTableEntries_;		// number of lookup table records
//	double ref_time_AFC_lookup_table_;	// reference (start time) of the first record (s since 2014-03-07 00:00:00)
//	double time_step_AFC_lookup_table_; // time step of the table

	// This is digitized Fig. 11 interpolated to 10-min interval from
	// Ashton C., Bruce A.C., Colledge G., Dickinson M. The Search for MH370.
	// THE JOURNAL OF NAVIGATION (2015), 68, 1–22. doi:10.1017/S037346331400068X

	double AFC_LOOKUP_TABLE[] = {
								 15.932, 16.382, 16.592, 16.704, 16.869, 17.218, 17.720, 18.257, 18.717, 19.122,
								 19.369, 19.476, 19.515, 19.483, 19.264, 18.866, 18.551, 18.397, 18.427, 18.602,
								 18.738, 18.760, 18.663, 18.415, 18.015, 17.547, 17.029, 16.361, 15.659, 14.914,
								 14.006, 13.228, 12.537, 11.661, 10.621, 9.561, 8.5731, 7.5284, 6.4872, 5.5041,
								 4.5013, 3.4886, 2.5883, 1.6774, 0.90689, 0.25986, -0.47901, -1.5045, -2.6689, -3.8234,
								 -4.964, -5.8111, -6.4525, -7.0422, -7.9048, -8.8941, -9.6847, -10.478, -11.194, -12.022,
								 -12.818, -13.424, -14.022, -14.646, -15.285, -15.872, -16.390, -16.870, -17.365, -17.891,
								 -18.344, -18.665, -18.910, -19.125, -19.296, -19.404, -19.510, -19.707, -19.985, -20.213,
								 -20.330, -20.356, -20.304, -20.166, -19.964, -19.795, -19.717, -19.701, -19.707, -19.696,
								 -19.628, -19.479, -19.267, -19.016, -18.770, -18.584, -18.419, -18.155, -17.729, -17.209,
								 -16.642, -16.040, -15.413, -14.761, -14.077, -13.374, -12.668, -11.937, -11.184, -10.433,
								 -9.7018, -8.9605, -8.1850, -7.3821, -6.5417, -5.6537, -4.672, -3.1114, -0.75563, 3.2262,
								 6.3628, 9.4108, 11.553, 13.553, 15.413, 14.592, 13.901, 13.369, 12.983, 12.702,
								 12.501, 12.393, 12.386, 12.526, 12.859, 13.272, 13.671, 14.009, 14.277, 14.527,
								 14.820, 15.167, 15.443, 15.539, 15.932, 16.382, 16.592, 16.704, 16.869, 17.218,
								 17.720, 18.257, 18.717, 19.122, 19.369, 19.476, 19.515};

	numAFCLookupTableEntries_ = sizeof(AFC_LOOKUP_TABLE) / sizeof(double);

	pAFC_Lookup_Table_ = new double[numAFCLookupTableEntries_ +1];

	for (int i=0; i<numAFCLookupTableEntries_; i++) pAFC_Lookup_Table_[i] = AFC_LOOKUP_TABLE[i];

	ref_time_AFC_lookup_table_ = 0.0;		// this table is interpolated / extrapolated to cover the whole day 2014-03-07 and beginning of 2014-03-08
	time_step_AFC_lookup_table_ = 600.0;	// time step (s)

}



////////////////////////////////////////////////////////////////////////////////////////
//
// Return measured (interpolated) AFC term (Hz) at the given time t
//
////////////////////////////////////////////////////////////////////////////////////////
double CMH370InmarsatData::GetAFC(double t)
{
	// Input t is the time (s) since 2014-03-07 00:00:00
	// Output is the AFC term (Hz).
	AFC_lookup_table_interp_method_ = 0; // the table already pre-interpolated to 10-min intervals, so it appears to make little sense to apply cubic interpolation, given the 'jitter' in the digitized Ashton et al. (2015) plot.

	//
	int ti = (int)(floor(t/time_step_AFC_lookup_table_));							// Compute index
	double tf = (t-time_step_AFC_lookup_table_*ti)/time_step_AFC_lookup_table_;		// fraction ranging from 0 to 1

	double fi;

	if (AFC_lookup_table_interp_method_ == 0)
	{
		// linear interpolant
		double tl = pAFC_Lookup_Table_[ti];													// lower lookup table value 
		double th = pAFC_Lookup_Table_[ti+1];												// higher lookup table value 
		fi = tl+(th-tl)*tf;
	}
	else
	{
		if (ti>0 || ti+1<numAFCLookupTableEntries_)
		{
			// cubic interpolant
			// Adapted from https://www.paulinternet.nl/?page=bicubic
			double* p = &pAFC_Lookup_Table_[ti-1];											// pointer so that:
			// p[0] = pAFC_Lookup_Table_[ti-1];												// second lower lookup table value 
			// p[1] = pAFC_Lookup_Table_[ti];												// lower lookup table value 
			// p[2] = pAFC_Lookup_Table_[ti+1];												// higher lookup table value 
			// p[3] = pAFC_Lookup_Table_[ti+2];												// second higher lookup table value 
			fi = p[1] + 0.5 * tf*(p[2] - p[0] + tf*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + tf*(3.0*(p[1] - p[2]) + p[3] - p[0])));
		}
		else
		{
			if (ti==0)
			{
				double v1 = pAFC_Lookup_Table_[ti];												// lower lookup table value 
				double v2 = pAFC_Lookup_Table_[ti+1];											// higher lookup table value 
				double v3 = pAFC_Lookup_Table_[ti+2];											// second higher lookup table value 
				double v0 = v1-(v2-v1);															// assume second lower lookup table value 
				fi = v1 + 0.5 * tf*(v2 - v0 + tf*(2.0*v0 - 5.0*v1 + 4.0*v2 - v3 + tf*(3.0*(v1 - v2) + v3 - v0)));
			}
			else
			{
				double v0 = pAFC_Lookup_Table_[ti-1];											// second lower lookup table value 
				double v1 = pAFC_Lookup_Table_[ti];												// lower lookup table value 
				double v2 = pAFC_Lookup_Table_[ti+1];											// higher lookup table value 
				double v3 = v2+(v2-v1);															// assume second higher lookup table value 
				fi = v1 + 0.5 * tf*(v2 - v0 + tf*(2.0*v0 - 5.0*v1 + 4.0*v2 - v3 + tf*(3.0*(v1 - v2) + v3 - v0)));
			}
		}
	}; // cubic interpolation method


	// calculate Doppler GES AFC (uses Burum to Perth station link)
	double GES_AFC = GetGESAFC_(t);

	return (fi - GES_AFC);
}


// GES AFC doppler shift component based on pilot signal shift from Burum GES (Netherlands) to Perth GES (Australia)
double CMH370InmarsatData::GetGESAFC_(double t)
{
	// Input: t - time in (s) since 2014-03-07 00:00:00
	// Output: GES AFC doppler shift component based on pilot signal shift Burum GES to Perth GES
	// This is analogous to Fig. 10 presented in
	// Ashton C., Bruce A.C., Colledge G., Dickinson M. The Search for MH370.
	// THE JOURNAL OF NAVIGATION (2015), 68, 1–22. doi:10.1017/S037346331400068X

	// -----------------------------------------------
	// get Inmarsat satellite position and velocity
	double satellite_XYZ[3];				// satellite position
	double satellite_UVW[3];				// satellite velocity
	double Zero_UVW[3] = {0.0, 0.0, 0.0};	// zero velocity

	GetInmarsatPos(satellite_XYZ[0], satellite_XYZ[1], satellite_XYZ[2], t);	// get satellite position in ECEF at the time t (s since UTC 2014-03-07 00:00:00)
	GetInmarsatVel(satellite_UVW[0], satellite_UVW[1], satellite_UVW[2], t);	// get satellite velocity in ECEF at the time t (s since UTC 2014-03-07 00:00:00)

	// -----------------------------------------------
	// Uplink Doppler ("L Band uplink"- see Ashton et al.)
	double dF_u = doppler(Burum_XYZ_, satellite_XYZ, Zero_UVW, satellite_UVW, Fr_uplink_, Doppler_Method_);

	// Downlink Doppler ("C Band downlink"- see Ashton et al.)
	double dF_d = doppler(satellite_XYZ, Perth_XYZ_, satellite_UVW, Zero_UVW, Fr_downlink_, Doppler_Method_);

	// Total:
	double dF = dF_u + dF_d;

	return dF;
}
