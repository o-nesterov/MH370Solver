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


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Class to handle environmental data:
// 1). Meteorological 
// 2). Magnetic declination
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef MH370ENVIRONDATA_CLASS

#define MH370ENVIRONDATA_CLASS

// meteorological dataset: 1 - GDAS1, 2 - ERA5
#define METEO_DATASET_NULL 0
#define METEO_DATASET_GDAS1 1
#define METEO_DATASET_ERA5 2

// magnetic declination dataset: 1 - NOAA
#define MAGN_DECL_DATASET_NULL 0
#define MAGN_DECL_DATASET_NOAA 1


class CMH370EnvironData
{
public:
	 CMH370EnvironData(); // constructor
	~CMH370EnvironData(); // destructor

	// Load meteorology. Supported meteorological datasets: "GDAS1", "ERA5".
	// If data already loaded, it unloads the previously loaded dataset and re-loads data
	bool LoadMeteo(int dataset, char* foldername);

    // Supported magnetic declination dataset: "NOAA"
	// If data already loaded, it unloads the previously loaded dataset and re-loads data
	bool LoadMagnetic(int dataset, char* foldername);

	// Get interpolated wind velocity {u,v} (m/s) components, air temperature t (deg C), air pressure p (Pa), and relative humidity r (%) 
	// at the spcified longitude lon (deg E), latitutde lat (deg N), altitude alt (m) and time ts (seconds since 2014-03-07 00:00:00 UTC).
	bool GetMeteo(double& u, double& v, double& t, double& p, double& r,
				  double ts, double lon, double lat, double alt);

	// Set default values to be returned if meteorological data cannot be loaded for whatever reason
	void SetDefault(double u, double v, double t, double p, double r);

	// Get magnetic declination decl at the given longitude lon (deg E) and latitude lat (deg N). 
	// Variation over the day 2014-03-07 and altitude are considered neglibible.
	double GetMagneticDeclination(double lon, double lat);


	// Set 4-D interpolation method for meteorology. Currently supported methods:
	// 0 - linear in time, bilinear horizontally, linear vertically
	// 1 - tricubic in time x horizontally, linear vertically
	// 2 - linear in time, bilinear horizontally, cubic spline vertically
	// 3 - tricubic in time x horizontally, cubic spline vertically
	void SetMeteoInterpMethod(int method);

	// Set interpolation method for magnetic declination. Currently supported methods:
	// 0 - bilinear
	// 1 - bicubic
	void SetDeclInterpMethod(int method);


	// Verbose for diagnostic (e.g., if data loading fails)
	void SetVerbose(FILE* outputstream, int verbose); // set verbose level (0 - silent, 1 - make diagnostic output during data load)

private:

	// ---------------------------------------------------------------------------------------
	// Meteorology:
	int meteo_dataset_;			// -1 - not loaded, 0 - GDAS1, 1 - ERA5

	double meteo_Start_Time_;	// dataset start time since 2014-03-07 00:00:00 UTC
	double meteo_dt_;			// time step (seconds)

	// dimensions:
	int meteo_nt_;				// number of time steps 
	int meteo_nx_;				// number of grid nodes in lon-direction
	int meteo_ny_;				// number of grid nodes in lat-direction
	int meteo_nz_;				// number of layers

	// bounding box
	double meteo_LON_MIN_;		// western boundary longitude
	double meteo_LON_MAX_;		// eastern boundary longitude
	double meteo_LAT_MIN_;		// southern boundary latitude
	double meteo_LAT_MAX_;		// northern boundary latitude

	// arrays of pointers on data [NumTs,nz,nx,ny]
	float*** pWIND_U_;			// Wind U-component (W->E, m/s)
	float*** pWIND_V_;			// Wind V-component (S->N, m/s)
	float*** pLAYERH_;			// Geopotential height (m)
	float*** pTEMPER_;			// Air temperature (deg C)
	float*** pRHUMID_;			// Relative humidity (%)

	float default_u_;			// Default wind U-component (m/s)
	float default_v_;			// Default wind V-component (m/s)
	float default_t_;			// Default air temperature (deg C)
	float default_p_;			// Default air pressure (Pa)
	float default_r_;			// Default relative humidity (%)

	double* pPressureLevels_;	// Pressure levels

	int meteo_interp_method_;	// interpolation method

	bool loadmeteo_gdas1_(char* foldername);	// load GDAS1
	bool loadmeteo_era5_(char* foldername);		// load ERA5

	double* p_meteo_workspace_;	// allocated workspace for interpolations (for spline etc.)

	// ------------------------------------------------------
	// Magnetic declination (deg)
	int magn_nx_;				// number of grid nodes in lon-direction
	int magn_ny_;				// number of grid nodes in lat-direction
	float* pMagneticDecl_;		// magnetic declination (deg)

	// bounding box
	double magn_LON_MIN_;		// western boundary longitude
	double magn_LON_MAX_;		// eastern boundary longitude
	double magn_LAT_MIN_;		// southern boundary latitude
	double magn_LAT_MAX_;		// northern boundary latitude

	int decl_interp_method_;	// interpolation method
	bool loadmagnetic_noaa_(char* foldername);	// load 0.25 x 0.25 deg NOAA data


	// ------------------------------------------------------------------------------------
	// Verbose diganostic stream
	FILE* fileout_; // message stream for diagnostic output (can be NULL) 
	int verboseLevel_; // verbose level: 0 - no output, 1 - output only for data loading diagnostic (missing files, etc.), and 2 - output for interpolation errors (e.g., outside of domain)
};


#endif
