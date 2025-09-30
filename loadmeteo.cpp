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
// Class to handle environmental and Inmarsat data - load meteorology
//
////////////////////////////////////////////////////////////////////////////////////////

bool CMH370EnvironData::LoadMeteo(int dataset, char* foldername)
{
	// deallocate U-arrays
	if (pWIND_U_)
	{
		for (int n=0; n<meteo_nt_; n++)
		{
			if (pWIND_U_[n])
			{
				for (int k=0; k<meteo_nz_; k++)
				{
					if (pWIND_U_[n][k]) delete pWIND_U_[n][k];
				}
				delete pWIND_U_[n];
			}
		}
	}

	// deallocate V-arrays
	if (pWIND_V_)
	{
		for (int n=0; n<meteo_nt_; n++)
		{
			if (pWIND_V_[n])
			{
				for (int k=0; k<meteo_nz_; k++)
				{
					if (pWIND_V_[n][k]) delete pWIND_V_[n][k];
				}
				delete pWIND_V_[n];
			}
		}
	}


	// deallocate geopotential height arrays
	if (pLAYERH_)
	{
		for (int n=0; n<meteo_nt_; n++)
		{
			if (pLAYERH_[n])
			{
				for (int k=0; k<meteo_nz_; k++)
				{
					if (pLAYERH_[n][k]) delete pLAYERH_[n][k];
				}
				delete pLAYERH_[n];
			}
		}
	}

	// deallocate temperature arrays
	if (pTEMPER_)
	{
		for (int n=0; n<meteo_nt_; n++)
		{
			if (pTEMPER_[n])
			{
				for (int k=0; k<meteo_nz_; k++)
				{
					if (pTEMPER_[n][k]) delete pTEMPER_[n][k];
				}
				delete pTEMPER_[n];
			}
		}
	}

	// deallocate humidity arrays
	if (pRHUMID_)
	{
		for (int n=0; n<meteo_nt_; n++)
		{
			if (pRHUMID_[n])
			{
				for (int k=0; k<meteo_nz_; k++)
				{
					if (pRHUMID_[n][k]) delete pRHUMID_[n][k];
				}
				delete pRHUMID_[n];
			}
		}
	}

	if (pPressureLevels_) delete pPressureLevels_;

	if (p_meteo_workspace_) delete p_meteo_workspace_;

	pWIND_U_ = NULL;
	pWIND_V_ = NULL;
	pLAYERH_ = NULL;
	pTEMPER_ = NULL;
	pRHUMID_ = NULL;
	pPressureLevels_ = NULL;
	p_meteo_workspace_ = NULL;
	meteo_dataset_ = METEO_DATASET_NULL; // meteo data not loaded

	// check if data type is supported
	if (!(dataset==METEO_DATASET_NULL || dataset==METEO_DATASET_GDAS1 || dataset==METEO_DATASET_ERA5))
	{
		if (fileout_) fprintf(fileout_,"Unsupported meteorological type. Currently available types: GDAS1 and ERA5.\n");
		return false;
	};

	if (dataset==METEO_DATASET_NULL)
	{
		return true; // nothing to load - default/standard meteorological conditions to be used
	}


	if (dataset==METEO_DATASET_GDAS1)
	{
		if (loadmeteo_gdas1_(foldername)) 
		{
			meteo_dataset_ = METEO_DATASET_GDAS1;
			return true; // GDAS1 was loaded 
		}
	}

	if (dataset==METEO_DATASET_ERA5)
	{

		if (loadmeteo_era5_(foldername)) 
		{
			meteo_dataset_ = METEO_DATASET_ERA5; 
			return true; 		// ERA5 was loaded
		}
	}

	// if we are here, then loading failed

	if (p_meteo_workspace_) delete p_meteo_workspace_; // deallocate workspace buffer
	p_meteo_workspace_ = NULL;

	// deallocate pressure levels
	if (pPressureLevels_) delete pPressureLevels_;
	pPressureLevels_ = NULL;

	// deallocate U-arrays
	if (pWIND_U_)
	{
		for (int n=0; n<meteo_nt_; n++)
		{
			if (pWIND_U_[n])
			{
				for (int k=0; k<meteo_nz_; k++)
				{
					if (pWIND_U_[n][k]) delete pWIND_U_[n][k];
				}
				delete pWIND_U_[n];
			}
		}
	}

	// deallocate V-arrays
	if (pWIND_V_)
	{
		for (int n=0; n<meteo_nt_; n++)
		{
			if (pWIND_V_[n])
			{
				for (int k=0; k<meteo_nz_; k++)
				{
					if (pWIND_V_[n][k]) delete pWIND_V_[n][k];
				}
				delete pWIND_V_[n];
			}
		}
	}


	// deallocate geopotential height arrays
	if (pLAYERH_)
	{
		for (int n=0; n<meteo_nt_; n++)
		{
			if (pLAYERH_[n])
			{
				for (int k=0; k<meteo_nz_; k++)
				{
					if (pLAYERH_[n][k]) delete pLAYERH_[n][k];
				}
				delete pLAYERH_[n];
			}
		}
	}

	// deallocate temperature arrays
	if (pTEMPER_)
	{
		for (int n=0; n<meteo_nt_; n++)
		{
			if (pTEMPER_[n])
			{
				for (int k=0; k<meteo_nz_; k++)
				{
					if (pTEMPER_[n][k]) delete pTEMPER_[n][k];
				}
				delete pTEMPER_[n];
			}
		}
	}

	// deallocate humidity arrays
	if (pRHUMID_)
	{
		for (int n=0; n<meteo_nt_; n++)
		{
			if (pRHUMID_[n])
			{
				for (int k=0; k<meteo_nz_; k++)
				{
					if (pRHUMID_[n][k]) delete pRHUMID_[n][k];
				}
				delete pRHUMID_[n];
			}
		}
	}

	pWIND_U_ = NULL;
	pWIND_V_ = NULL;
	pLAYERH_ = NULL;
	pTEMPER_ = NULL;
	pRHUMID_ = NULL;

	meteo_dataset_ = METEO_DATASET_NULL; // meteo data not loaded

	return false;
}


////////////////////////////////////////////////////////////////////////////////////////
// Load GDAS1 meteorological data
////////////////////////////////////////////////////////////////////////////////////////


bool CMH370EnvironData::loadmeteo_gdas1_(char* foldername)
{
    // "GDAS1": 1.00 x 1.00 deg resolution, 23 vertical presure levels, 3-hourly
    // References: https://www.ready.noaa.gov/gdas1.php, https://www.ready.noaa.gov/data/archives/gdas1/, https://www.ready.noaa.gov/ready/gdas1grid.gif
	// The data was downloaded in the propriatery ARL NOAA format, extracted by time x layer, 
	// The data are on a 360 by 181 latitude-longitude (grid). The lower-left corner (1,1) is (0W,90S). The upper-right corner (360,181) is (1W, 90N) (Ref.: https://www.ready.noaa.gov/gdas1.php)


//	p_meteo_workspace_ = new double[16*meteo_nz_+1]; // allocate workspace buffer
//	if (p_meteo_workspace_==NULL) {if (fileout_) fprintf(fileout_, "Failed to allocate workspace memory.\n"); return false;};


	meteo_Start_Time_ = 12*3600.0;	// start time 12:00 - seconds since 2014-03-07 00:00:00 UTC
	meteo_dt_ = 3*3600.0;			// 3-hourly time steps

	meteo_nt_ = 7;		// number of time steps (updated later): 2014-03-07 12:00, 15:00, 18:00, 21:00, 2014-03-08 00:00, 03:00, 06:00 (the first and last are needed for cubic interpolation)
	meteo_nx_ = 360;	// number of nodes in x-direction
	meteo_ny_ = 181;	// number of cells in y-direction
	meteo_nz_ = 23;		// number of layers

	// bounding box of ERA5 data (as downloaded and extracted)
	meteo_LON_MIN_ =  0.0;
	meteo_LON_MAX_ = 359.0;
	meteo_LAT_MIN_ = -90.0;
	meteo_LAT_MAX_ =  90.0;

	// start date  (not inclusive)
	int start_yyyy = 2014;
	int start_month = 3;
	int start_day = 7;
	int start_hh = 9;

	// end date (inclusive)
	int end_yyyy = 2014;
	int end_month = 3;
	int end_day = 8;
	int end_hh = 6;

	// filenames of data files
	char filename[1024];

	// ------------------------------------------------------
	// allocate all arrays; if allocation fails, return false deallocate all from the calling procedure
	// as data is rather large (5 variables x 7 x 23 x 360 x 181 floats = ~200 Mb) make sure that memory is properly allocated 
	float*** pMeteoArrays[5]; // 'plain' (horizontal) meteo arrays to be allocated and loaded in the same way: pMeteoArrays[0] = pWIND_U_[time][layer][..][..], pWIND_V_, pLAYERH_, pTEMPER_, pRHUMID_

	pLAYERH_ = new float**[meteo_nt_];
	if (pLAYERH_==NULL) {if (fileout_) fprintf(fileout_, "Failed to allocate memory for meteo data (H).\n"); return false;};
	pWIND_U_ = new float**[meteo_nt_];
	if (pWIND_U_==NULL) {if (fileout_) fprintf(fileout_, "Failed to allocate memory for meteo data (U).\n"); return false;};
	pWIND_V_ = new float**[meteo_nt_];
	if (pWIND_V_==NULL) {if (fileout_) fprintf(fileout_, "Failed to allocate memory for meteo data (V).\n"); return false;};
	pTEMPER_ = new float**[meteo_nt_];
	if (pTEMPER_==NULL) {if (fileout_) fprintf(fileout_, "Failed to allocate memory for meteo data (T).\n"); return false;};
	pRHUMID_ = new float**[meteo_nt_];
	if (pRHUMID_==NULL) {if (fileout_) fprintf(fileout_, "Failed to allocate memory for meteo data (R).\n"); return false;};

	p_meteo_workspace_ = new double[16*meteo_nz_+1]; // allocate workspace buffer
	if (p_meteo_workspace_==NULL) {if (fileout_) fprintf(fileout_, "Failed to allocate workspace memory.\n"); return false;};


	pMeteoArrays[0] = pLAYERH_;
	pMeteoArrays[1] = pWIND_U_;
	pMeteoArrays[2] = pWIND_V_;
	pMeteoArrays[3] = pTEMPER_;
	pMeteoArrays[4] = pRHUMID_;

	for (int varId = 0; varId<5; varId++)
	{
		for (int n=0; n<meteo_nt_; n++) pMeteoArrays[varId][n] = NULL; // this is to prevent crash on deallocation in case of the failure to allocate all arrays here
	}

	for (int varId = 0; varId<5; varId++)
	{
		for (int n=0; n<meteo_nt_; n++) 
		{
			pMeteoArrays[varId][n] = new float*[meteo_nz_];
			if (pMeteoArrays[varId][n] == NULL)
			{
				if (fileout_) fprintf(fileout_, "Failed to allocate memory for meteo data.\n"); 
				return false;
			}
			for (int k=0; k<meteo_nz_; k++) pMeteoArrays[varId][n][k] = NULL;  // this is to prevent crash on deallocation in case of the failure to allocate all arrays here
		}

	}


	// read all files and create global arrays
	int yyyy = start_yyyy;
	int mm = start_month;
	int dd = start_day;
	int hh = start_hh;

	int ts = 1; // meteorological time step (in hours)

	if (fileout_) fprintf(fileout_, "Loading GDAS1 meteorological data...\n");

	int dd_in_month[] = {31,28,31,30,31,30,31,31,30,31,30,31};

	while(!(yyyy==end_yyyy && mm==end_month && dd==end_day && hh==end_hh))
	{	
		if ((yyyy & 3) == 0)
		{
			dd_in_month[1] = 29;
		}
		else
		{
			dd_in_month[1] = 28;
		}
		
    
		hh += 3;  // increment
    
		if (hh==24) {hh=0; dd++;};
		if (dd>dd_in_month[mm-1]) {dd=1; mm++;};
		if (mm>12) {mm=1; yyyy++;};

		// Read layers
		if (fileout_) fprintf(fileout_, "Loading meteo... %4.4d-%2.2d-%2.2d-%2.2d\n", yyyy,mm,dd,hh);
    
		for (int layer=1; layer<=meteo_nz_; layer++)
		{
			for (int varId = 0; varId < 5; varId++) // iterations by variables
			{
				// allocate actual array
				pMeteoArrays[varId][ts-1][layer-1] = new float[meteo_nx_*meteo_ny_+1];
				if (pMeteoArrays[varId][ts-1][layer-1] == NULL)
				{
					if (fileout_) fprintf(fileout_, "Failed to allocate memory for meteo data.\n"); 
					return false;
				}

				switch (varId)
				{
				case 0:
					sprintf(filename, "%s\\HGTS\\HGTS-%4.4d-%2.2d-%2.2d-%2.2d-%2.2d.bin",foldername,yyyy,mm,dd,hh,layer); // geopotential height
					break;
				case 1:
					sprintf(filename, "%s\\UWND\\UWND-%4.4d-%2.2d-%2.2d-%2.2d-%2.2d.bin",foldername,yyyy,mm,dd,hh,layer); // u-component
					break;
				case 2:
					sprintf(filename, "%s\\VWND\\VWND-%4.4d-%2.2d-%2.2d-%2.2d-%2.2d.bin",foldername,yyyy,mm,dd,hh,layer); // v-component
					break;
				case 3:
					sprintf(filename, "%s\\TEMP\\TEMP-%4.4d-%2.2d-%2.2d-%2.2d-%2.2d.bin",foldername,yyyy,mm,dd,hh,layer); // temperature
					break;
				case 4:
					sprintf(filename, "%s\\RELH\\RELH-%4.4d-%2.2d-%2.2d-%2.2d-%2.2d.bin",foldername,yyyy,mm,dd,hh,((layer < 21)?layer:21)); // relative humidity; there are only 21 layers of humidity in GDAS1
					break;
				};


				// ----------------------------------------
				// read data
			
				size_t nread = 0;
			
				FILE* fid = fopen(filename,"rb");
				if (fid)
				{
					nread = fread(pMeteoArrays[varId][ts-1][layer-1], sizeof(float), meteo_nx_*meteo_ny_, fid);
					fclose(fid);
					if (nread!=meteo_nx_*meteo_ny_)
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
			} // end variable loop (varId = 0...4)
  

		} // end layer loop (layer = 1...23)
    
		ts = ts + 1; // end time step loop (ts = 1...7)
	}; //end while (time loop)


	// Pressure levels of 3D GDAS1 (as extracted from downloaded netcdf) 
	pPressureLevels_ = new double[meteo_nz_+1];
	if (pPressureLevels_ == NULL) {if (fileout_) fprintf(fileout_, "Error allocating memory for pressure levels.\n"); return false;};
	// Initialize pressure levels
	unsigned short p_column[] = {1000, 975, 950, 925, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100, 50, 20};
	for (unsigned long i = 0; i<meteo_nz_; i++) pPressureLevels_[i] = 100.0* ((double) p_column[i]); // convert to Pa


	if (fileout_) fprintf(fileout_, "Meteo data has been loaded.\n");

	return true;

}




////////////////////////////////////////////////////////////////////////////////////////
// Load ERA5 meteorological data
////////////////////////////////////////////////////////////////////////////////////////


bool CMH370EnvironData::loadmeteo_era5_(char* foldername)
{
    // "ERA5":  0.25 x 0.25 deg resolution, 37 vertical layers, hourly
    // Reference: https://cds.climate.copernicus.eu/
	// The original data was downloaded in NetCDF format in the region [40E,130E]x[60S,40N], then extracted by time x layer, 
	// and then converted and saved into binary float format to avoid using netcdf library


	meteo_Start_Time_ = 15*3600.0;	// start time 15:00 - seconds since 2014-03-07 00:00:00 UTC
	meteo_dt_ = 3600.0;				// 1 hour time steps

	meteo_nt_ = 12;		// number of time steps (updated later): 2014-03-07 15:00, 16:00, 17:00, 18:00, 19:00, 20:00, 21:00, 22:00, 23:00; 2014-03-08 00:00, 01:00, 02:00 
	meteo_nx_ = 361;	// number of nodes in x-direction
	meteo_ny_ = 401;	// number of cells in y-direction
	meteo_nz_ = 37;		// number of layers

	// bounding box of ERA5 data (as downloaded and extracted)
	meteo_LON_MIN_ =  40.0;
	meteo_LON_MAX_ = 130.0;
	meteo_LAT_MIN_ = -60.0;
	meteo_LAT_MAX_ =  40.0;

	// start date  (not inclusive)
	int start_yyyy = 2014;
	int start_month = 3;
	int start_day = 7;
	int start_hh = 14;

	// end date (inclusive)
	int end_yyyy = 2014;
	int end_month = 3;
	int end_day = 8;
	int end_hh = 2;

	// filenames of data files
	char filename[1024];

	// ------------------------------------------------------
	// allocate all arrays; if allocation fails, return false deallocate all from the calling procedure
	// as data is rather large (5 variables x 12 x 37 x 361 x 401 floats = ~1.3 Gb) make sure that memory is properly allocated 
	float*** pMeteoArrays[5]; // 'plain' (horizontal) meteo arrays to be allocated and loaded in the same way: pMeteoArrays[0] = pWIND_U_[time][layer][..][..], pWIND_V_, pLAYERH_, pTEMPER_, pRHUMID_

	pLAYERH_ = new float**[meteo_nt_];
	if (pLAYERH_==NULL) {if (fileout_) fprintf(fileout_, "Failed to allocate memory for meteo data (H).\n"); return false;};
	pWIND_U_ = new float**[meteo_nt_];
	if (pWIND_U_==NULL) {if (fileout_) fprintf(fileout_, "Failed to allocate memory for meteo data (U).\n"); return false;};
	pWIND_V_ = new float**[meteo_nt_];
	if (pWIND_V_==NULL) {if (fileout_) fprintf(fileout_, "Failed to allocate memory for meteo data (V).\n"); return false;};
	pTEMPER_ = new float**[meteo_nt_];
	if (pTEMPER_==NULL) {if (fileout_) fprintf(fileout_, "Failed to allocate memory for meteo data (T).\n"); return false;};
	pRHUMID_ = new float**[meteo_nt_];
	if (pRHUMID_==NULL) {if (fileout_) fprintf(fileout_, "Failed to allocate memory for meteo data (R).\n"); return false;};

	p_meteo_workspace_ = new double[16*meteo_nz_+1]; // allocate workspace buffer
	if (p_meteo_workspace_==NULL) {if (fileout_) fprintf(fileout_, "Failed to allocate workspace memory.\n"); return false;};


	pMeteoArrays[0] = pLAYERH_;
	pMeteoArrays[1] = pWIND_U_;
	pMeteoArrays[2] = pWIND_V_;
	pMeteoArrays[3] = pTEMPER_;
	pMeteoArrays[4] = pRHUMID_;

	for (int varId = 0; varId<5; varId++)
	{
		for (int n=0; n<meteo_nt_; n++) pMeteoArrays[varId][n] = NULL; // this is to prevent crash on deallocation in case of the failure to allocate all arrays here
	}

	for (int varId = 0; varId<5; varId++)
	{
		for (int n=0; n<meteo_nt_; n++) 
		{
			pMeteoArrays[varId][n] = new float*[meteo_nz_];
			if (pMeteoArrays[varId][n] == NULL)
			{
				if (fileout_) fprintf(fileout_, "Failed to allocate memory for meteo data.\n"); 
				return false;
			}
			for (int k=0; k<meteo_nz_; k++) pMeteoArrays[varId][n][k] = NULL;  // this is to prevent crash on deallocation in case of the failure to allocate all arrays here
		}

	}


	// read all files and create global arrays
	int yyyy = start_yyyy;
	int mm = start_month;
	int dd = start_day;
	int hh = start_hh;

	int ts = 1; // meteorological time step (in hours)

	if (fileout_) fprintf(fileout_, "Loading ERA5 meteorological data...\n");

	int dd_in_month[] = {31,28,31,30,31,30,31,31,30,31,30,31};

	while(!(yyyy==end_yyyy && mm==end_month && dd==end_day && hh==end_hh))
	{	
		if ((yyyy & 3) == 0)
		{
			dd_in_month[1] = 29;
		}
		else
		{
			dd_in_month[1] = 28;
		}
		
    
		hh++;  // increment
    
		if (hh==24) {hh=0; dd++;};
		if (dd>dd_in_month[mm-1]) {dd=1; mm++;};
		if (mm>12) {mm=1; yyyy++;};

		// Read layers
		if (fileout_) fprintf(fileout_, "Loading meteo... %4.4d-%2.2d-%2.2d-%2.2d\n", yyyy,mm,dd,hh);
    
		for (int layer=1; layer<=meteo_nz_; layer++)
		{
			for (int varId = 0; varId < 5; varId++) // iterations by variables
			{
				// allocate actual array
				pMeteoArrays[varId][ts-1][layer-1] = new float[meteo_nx_*meteo_ny_+1];
				if (pMeteoArrays[varId][ts-1][layer-1] == NULL)
				{
					if (fileout_) fprintf(fileout_, "Failed to allocate memory for meteo data.\n"); 
					return false;
				}

				switch (varId)
				{
				case 0:
					sprintf(filename, "%s\\HGTS\\HGTS-%4.4d-%2.2d-%2.2d-%2.2d-%2.2d.bin",foldername,yyyy,mm,dd,hh,layer); // geopotential height
					break;
				case 1:
					sprintf(filename, "%s\\UWND\\UWND-%4.4d-%2.2d-%2.2d-%2.2d-%2.2d.bin",foldername,yyyy,mm,dd,hh,layer); // u-component
					break;
				case 2:
					sprintf(filename, "%s\\VWND\\VWND-%4.4d-%2.2d-%2.2d-%2.2d-%2.2d.bin",foldername,yyyy,mm,dd,hh,layer); // v-component
					break;
				case 3:
					sprintf(filename, "%s\\TEMP\\TEMP-%4.4d-%2.2d-%2.2d-%2.2d-%2.2d.bin",foldername,yyyy,mm,dd,hh,layer); // temperature
					break;
				case 4:
					sprintf(filename, "%s\\RELH\\RELH-%4.4d-%2.2d-%2.2d-%2.2d-%2.2d.bin",foldername,yyyy,mm,dd,hh,layer); // relative humidity
					break;
				};


				// ----------------------------------------
				// read data
			
				size_t nread = 0;
			
				FILE* fid = fopen(filename,"rb");
				if (fid)
				{
					nread = fread(pMeteoArrays[varId][ts-1][layer-1], sizeof(float), meteo_nx_*meteo_ny_, fid);
					fclose(fid);
					if (nread!=meteo_nx_*meteo_ny_)
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
			} // end variable loop (varId = 0...4)
  

		} // end layer loop (layer = 1...37)
    
		ts = ts + 1; // end time step loop (ts = 1...12)
	}; //end while (time loop)

	// Pressure levels of 3D ERA5 (as extracted from downloaded netcdf) 
	pPressureLevels_ = new double[meteo_nz_+1];
	if (pPressureLevels_ == NULL) {if (fileout_) fprintf(fileout_, "Error allocating memory for pressure levels.\n"); return false;};
	// Initialize pressure levels
	unsigned short p_column[] = {1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 225, 200, 175, 150, 125, 100, 70, 50, 30, 20, 10, 7, 5, 3, 2, 1};
	for (unsigned long i = 0; i<meteo_nz_; i++) pPressureLevels_[i] = 100.0* ((double) p_column[i]); // convert to Pa

	// convert temperature from K to C
	for (int n=0; n<meteo_nt_; n++) 
	{
		for (int l=0; l<meteo_nz_; l++)
		{
			for (int i=0; i<meteo_nx_*meteo_ny_; i++)
			{
				pTEMPER_[n][l][i] -= (float)273.16;
			}
		}
	};

	if (fileout_) fprintf(fileout_, "Meteo data has been loaded.\n");

	return true;
}
