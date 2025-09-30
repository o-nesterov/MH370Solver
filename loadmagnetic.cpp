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
#include "environdata.h"

////////////////////////////////////////////////////////////////////////////////////////
//
// Class to handle environmental and Inmarsat data - load magnetic declination data
//
////////////////////////////////////////////////////////////////////////////////////////

bool CMH370EnvironData::LoadMagnetic(int dataset, char* foldername)
{
    // Currently supported magnetic declination dataset: NOAA's calculator 0.25 deg x 0.25 deg
	if (pMagneticDecl_) delete pMagneticDecl_;
	pMagneticDecl_ = NULL;


	// check if data type is supported
	if (!(dataset==MAGN_DECL_DATASET_NULL || dataset==MAGN_DECL_DATASET_NOAA))
	{
		if (fileout_) fprintf(fileout_,"Unsupported magnetic declination data type. Currently available type is NOAA only.\n");
		return false;
	};

	if (dataset == MAGN_DECL_DATASET_NULL)
	{
		return true;	// pMagneticDecl_ = NULL here
	}

	if (dataset == MAGN_DECL_DATASET_NOAA)
	{
		if (loadmagnetic_noaa_(foldername)) return true;
	}

	if (pMagneticDecl_)
	{
		delete pMagneticDecl_; 
		pMagneticDecl_ = NULL;
	};

	return false;
}


////////////////////////////////////////////////////////////////////////////////////////
// Load NOAA 0.25 deg x 0.25 deg magnetic declination data
////////////////////////////////////////////////////////////////////////////////////////


bool CMH370EnvironData::loadmagnetic_noaa_(char* foldername)
{
    // "NOAA": 0.25 x 0.25 deg resolution
    // References: http://maps.ngdc.noaa.gov/viewers/historical_declination/, 
    //             http://www.ngdc.noaa.gov/geomag/WMM/calculators.shtml, 
    //             http://www.ngdc.noaa.gov/geomag-web/#igrfgrid

	magn_nx_ = 1440;
	magn_ny_ = 719;

	magn_LON_MIN_ = -179.75;								// western boundary longitude
	magn_LON_MAX_ = magn_LON_MIN_ + 0.25*(magn_nx_-1);		// eastern boundary longitude
	magn_LAT_MIN_ = -89.75;									// southern boundary latitude
	magn_LAT_MAX_ = magn_LAT_MIN_ + 0.25*(magn_ny_-1);		// northern boundary latitude


	pMagneticDecl_ = new float[magn_nx_*magn_ny_+1];
	if (pMagneticDecl_ == NULL) 
	{
		if (fileout_) fprintf(fileout_, "Failed to allocate memory for magnetic declination data.\n");
		return false;
	}

	// ----------------------------------------
	// read data
	if (fileout_) fprintf(fileout_, "Loading magnetic declination...\n");

	char filename[1024];
	sprintf(filename, "%s\\Declination.bin", foldername);

	size_t nread = 0;
			
	FILE* fid = fopen(filename,"rb");
	if (fid)
	{
		nread = fread(pMagneticDecl_, sizeof(float), magn_nx_*magn_ny_, fid);
		fclose(fid);
		if (nread!=magn_nx_*magn_ny_)
		{
			if (fileout_) fprintf(fileout_, "Error reading file: %s\n", filename);
			return false;
		}
	}
	else
	{
		if (fileout_) fprintf(fileout_, "Error opening file: %s\n", filename);
		return false;
	}

	if (fileout_) fprintf(fileout_, "Magnetic declination loaded.\n");
	return true;
}
