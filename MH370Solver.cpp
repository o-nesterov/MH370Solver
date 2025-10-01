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
#include "trajectory.h"
#include "params_parser.h"


// Define actions:
// ACT_UDF - Undefined
// ACT_PRT - Print ping information
// ACT_TRJ - build trajectory
// ACT_CAL - run calibration
// ACT_OPT - run optimization
// ACT_MET - extract meteorology
// ACT_COM - compute BTO, satellite-aircraft distance and BFO 
// ACT_VAL - run validation

#define ACT_UDF -1
#define ACT_PRT 0
#define ACT_TRJ 1
#define ACT_CAL 2
#define ACT_OPT 3
#define ACT_MET 4
#define ACT_COM 5
#define ACT_VAL 6


// Actions called from the main()
extern void BuildTrajectory(char* filename, bool isVerbose);
extern void Calibrate(char* outfilename, bool isVerbose);
extern void Validate(char* inpfilename, bool isVerbose);
extern void Optimize(char* optname, char* inpfilename, char* outfilename, char* str_optweight1, char* str_optweight2, bool isVerbose);
extern void PrintPings(char* filename);
extern void GetMeteo(char* datasetname, char* dir, char*str_time, char*str_lon, char*str_lat, char*str_alt, char*inpfilename, char*outfilename, bool isVerbose);
extern void ComputeBTO(char* inpfilename, char* outfilename, 
					   char*str_time, char*str_lon, char*str_lat, char*str_alt);
extern void ComputeDist(char* inpfilename, char* outfilename, 
						char*str_time, char*str_lon, char*str_lat, char*str_alt);
extern void ComputeBFO(char* inpfilename, char* outfilename, 
					   char*str_time, char*str_lon, char*str_lat, char*str_alt, char*str_hdg, char*str_spd, char*str_w, 
					   char* AES_Earth_model, char* AES_Altitude);



void main(int argc, char* argv[])
{
	char inpfilename[1024];	// input configuration file
	char outfilename[1024];	// output configuration file
	char meteodir[1024]; // directory containing meteorological data
	char str_arg1[64];
	char str_arg2[64];
	char str_arg3[64];
	char str_arg4[64];
	char str_arg5[64];
	char str_arg6[64];
	char str_arg7[64];
	char str_arg8[64];
	char str_arg9[64];
	char optname[64];		// option name (optimization, etc.)

	// initialize
	bool isVerbose = false;	// print progress
	int action = ACT_UDF;	// initialization - undefined action

	sprintf(inpfilename,"");
	sprintf(outfilename,"");
	sprintf(meteodir,"");
	sprintf(str_arg1,"");
	sprintf(str_arg2,"");
	sprintf(str_arg3,"");
	sprintf(str_arg4,"");
	sprintf(str_arg5,"");
	sprintf(str_arg6,"");
	sprintf(str_arg7,"");
	sprintf(str_arg8,"");
	sprintf(str_arg9,"");
	sprintf(optname,"");


	if (argc<2)
	{
		printf("Insufficient arguments. Usage:\n");
		printf("-traj configfilename: build trajectory as specified.\n");
		printf("-opt0 a,b,bb,d,x,y,h,s,g1,g2,m,lnav inpconfigfilename outconfigfilename: optimization by sum criteria.\n");
		printf("-opt1 a,b,bb,d,t,x,y,h,s,g1,g2,m,lnav inpconfigfilename outconfigfilename: optimization by max criteria.\n");
		printf("-opt2 a,b,bb,d,t,x,y,h,s,g1,g2,m,lnav inpconfigfilename outconfigfilename: optimization by var max criteria.\n");
		printf("\"a,b,bb,d,t,x,y,h,s,g1,g2,m,lnav\" is the list of possible optimization variables\n");
		printf("-optm inpconfigfilename outconfigfilename: weight optimization.\n");
		printf("-optweight1 val1: optional optimization weight for ping distance. Replaces default value 0.0001 m^-1.\n");
		printf("-optweight2 val2: optional optimization weight for bfo. Replaces default value of 0.2 Hz^-1.\n");
		printf("-calibrate: run calibration.\n");
		printf("-validate filename: run validation using specified setup file.\n");
		printf("-printpings: print Inmarsat ping data.\n");
		printf("-getmeteo dataset dir: interactive extraction of meteorology.\n");
		printf("-getmeteo dataset dir time, lon, lat, alt: extract meteorology at specified location and time.\n");
		printf("-getmeteo dataset dir -file inpfilename outfilename: extract meteorology using list {time,lon,lat,alt}.\n");
		printf("-compute bto t, lon, lat, alt: compute BTO.\n");
		printf("-compute bto -file inpfilename: compute BTO using {time, lon, lat, alt} stored in file.\n");
		printf("-compute bfo t, lon, lat, alt, hdg, spd, w: compute BFO.\n");
		printf("-compute bfo -file inpfilename: compute BFO using {time, lon, lat, alt, u, v, w} stored in file.\n");
		printf("-verbose: print progressive output.\n");
		return;
	}

	int argn = 1;
	while (argn<argc)
	{
		// --------------------
		// build trajectory
		if (strcmp(argv[argn],"-traj")==0 && action==ACT_UDF)
		{
			action = ACT_TRJ;
			if (argn+1 == argc)
			{
				printf("Error. Missing configuration file name.\n");
				return;
			}
			else
			{
				if (strlen(argv[argn+1])>1023)
				{
					printf("Error: too long filename.\n");
					return;
				}
				strcpy(inpfilename, argv[argn+1]);
				argn+=2;
				continue;
			}
		}

		// --------------------
		// run BTO and BFO calibration
		if (strcmp(argv[argn],"-calibrate")==0 && action==ACT_UDF)
		{
			action = ACT_CAL;
			if (argn+1 == argc)
			{
				strcpy(outfilename, "");
				argn++;
				continue;
			}
			else
			{
				if (strlen(argv[argn+1])>1023)
				{
					printf("Error: too long filename.\n");
					return;
				}
				strcpy(outfilename, argv[argn+1]);
				argn+=2;
				continue;
			}
		}

		// --------------------
		// run BTO and BFO validation
		if (strcmp(argv[argn],"-validate")==0 && action==ACT_UDF)
		{
			action = ACT_VAL;
			if (argn+1 == argc)
			{
				printf("Error. Missing validation setup file name.\n");
				return;
			}
			else
			{
				if (strlen(argv[argn+1])>1023)
				{
					printf("Error: too long input filename.\n");
					return;
				}
				else
				{
					strcpy(inpfilename, argv[argn+1]);
				};

				argn+=2;
				continue;
			}
		}


		// --------------------
		// run optimization
		if ((strcmp(argv[argn],"-opt0")==0 || strcmp(argv[argn],"-opt1")==0 || strcmp(argv[argn],"-opt2")==0) && action==ACT_UDF)
		{
			action = ACT_OPT;
			if (argn+3 >= argc)
			{
				printf("Error. Missing arguments: list of optimization variables, input or output file names.\n");
				printf("Usage: -opt0\opt1\opt2 list, input filename, output filename.\n");
				return;
			}
			else
			{
				if (strlen(argv[argn])>62)
				{
					printf("Error: too long list of optimization variables.\n");
					return;
				}

				if (strlen(argv[argn+1])>1023 || strlen(argv[argn+2])>1023)
				{
					printf("Error: too long filename.\n");
					return;
				}

				// copy type of optimization: the first letter 0, 1, 2 or m
				if (strcmp(argv[argn],"-opt0")==0) {optname[0]='0'; optname[1]=',';}
				if (strcmp(argv[argn],"-opt1")==0) {optname[0]='1'; optname[1]=',';}
				if (strcmp(argv[argn],"-opt2")==0) {optname[0]='2'; optname[1]=',';}

				strcpy(&optname[2], argv[argn+1]);
				strcpy(inpfilename, argv[argn+2]);
				strcpy(outfilename, argv[argn+3]);
				argn+=4;
				continue;
			}
		}


		// --------------------
		// run weight optimization
		if (strcmp(argv[argn],"-optm")==0 && action==ACT_UDF)
		{
			action = ACT_OPT;
			if (argn+2 >= argc)
			{
				printf("Error. Missing input or output, or both optimization file names.\n");
				return;
			}
			else
			{
				if (strlen(argv[argn+1])>1023 || strlen(argv[argn+2])>1023)
				{
					printf("Error: too long filename.\n");
					return;
				}
				strcpy(inpfilename, argv[argn+1]);
				strcpy(outfilename, argv[argn+2]);
				strcpy(optname, "m");
				argn+=3;
				continue;
			}
		}

		// --------------------
		// print Inmarsat pings data
		if (strcmp(argv[argn],"-printpings")==0 && action==ACT_UDF)
		{
			action = ACT_PRT;
			if (argn+1 == argc)
			{
				sprintf(outfilename,"");
				isVerbose = true;
				argn++;
				continue;
			}
			else
			{
				if (strlen(argv[argn+1])>1023)
				{
					printf("Error: too long filename.\n");
					return;
				}
				strcpy(outfilename, argv[argn+1]);
				argn+=2;
				continue;
			}
		}



		// --------------------
		// Extract meteorological data at provided location/time
		if (strcmp(argv[argn],"-getmeteo")==0 && action==ACT_UDF)
		{
			action = ACT_MET;

			if (argn+2 >= argc)
			{
					printf("Error: insufficient arguments.\n");
					printf("Usage: -getmeteo ERA5/GDAS1 MeteoDir (for interactive) or\n");
					printf("       -getmeteo ERA5/GDAS1 MeteoDir time lon lat alt.\n");
					printf("       -getmeteo ERA5/GDAS1 MeteoDir -file infilename outfilename\n");
					return;
			}

			if (strcmp(argv[argn+1], "era5")==0 || strcmp(argv[argn+1], "ERA5")==0 || strcmp(argv[argn+1], "Era5")==0 || strcmp(argv[argn+1], "gdas1")==0 || strcmp(argv[argn+1], "GDAS1")==0 || strcmp(argv[argn+1], "Gdas1")==0)
			{
				strcpy(optname, argv[argn+1]);
			}
			else
			{
				printf("Error: unsupported dataset. Must be GDAS1 or ERA5\n");
				return;
			}

			if (strlen(argv[argn+2])>1023)
			{
				printf("Error: too long directory name.\n");
				return;
			}
			else
			{
				strcpy(meteodir, argv[argn+2]);
			}
	
			// -----------------------------
			// now analyze what to do next
			if (argn+3 == argc)
			{
				// This is an interactive mode
				strcpy(inpfilename, "");
				strcpy(outfilename, "");
				strcpy(str_arg1, "");
				strcpy(str_arg2, "");
				strcpy(str_arg3, "");
				strcpy(str_arg4, "");

				argn+=3;
				continue;
			}
			else
			{
				// This is a non-interactive mode - all arguments must be supplied - either filenames or {time, lon, lat, alt}			
				if (strcmp(argv[argn+3], "-file") == 0)
				{

					// this followed by input-output file specifications
					if (argn+5 >= argc)
					{
						printf("Error. Incorrect format or missing arguments.\n");
						printf("Usage: -getmeteo ERA5/GDAS1 MeteoDir (for interactive) or\n");
						printf("       -getmeteo ERA5/GDAS1 MeteoDir time lon lat alt.\n");
						printf("       -getmeteo ERA5/GDAS1 MeteoDir -file infilename outfilename\n");
						return;
					}
					else
					{
						strcpy(inpfilename, argv[argn+4]);
						strcpy(outfilename, argv[argn+5]);
						strcpy(str_arg1,"file");
						argn+=6;
						continue;
					}
				}


				if (argn+6 >= argc)
				{
					printf("Error. Incorrect format or missing arguments.\n");
					printf("Usage: -getmeteo ERA5/GDAS1 MeteoDir (for interactive) or\n");
					printf("       -getmeteo ERA5/GDAS1 MeteoDir time lon lat alt.\n");
					printf("       -getmeteo ERA5/GDAS1 MeteoDir -file infilename outfilename\n");
					return;
				}

				if (strlen(argv[argn+3])>63 || strlen(argv[argn+4])>63  || strlen(argv[argn+5])>63  || strlen(argv[argn+6])>63 )
				{
					printf("Unsupported argument format.\n");
					return;
				}
				strcpy(inpfilename, "");
				strcpy(outfilename, "");
				strcpy(str_arg1, argv[argn+3]);
				strcpy(str_arg2, argv[argn+4]);
				strcpy(str_arg3, argv[argn+5]);
				strcpy(str_arg4, argv[argn+6]);

				argn+=7;
				continue;
			}
		}


		// --------------------
		// Compute BTO, distance and BFO
		if (strcmp(argv[argn],"-compute")==0 && action==ACT_UDF)
		{
			action = ACT_COM;
			if (argn+2 >= argc)
			{
				printf("Error. Missing arguments. Usage: -compute bto/bfo/dist arguments.\n");
				return;
			}
			else
			{
				if (!(strcmp(argv[argn+1], "bto")==0 || strcmp(argv[argn+1], "bfo")==0 || 
					  strcmp(argv[argn+1], "BTO")==0 || strcmp(argv[argn+1], "BFO")==0))
				{
					printf("Error. Usage: -compute bto/bfo/dist arguments\n");
					return;
				}

				// here argv[argn+1] is either bto, or dist or bfo
				strcpy(optname, argv[argn+1]);

				if (strcmp(argv[argn+2], "-file")==0)
				{
					// name of file, which contains data {time, lon, lot, alt} or {time, lon, lot, alt, u, v, w}, and
					// and name of output file must be specified
					if (argn+3 >= argc)
					{
						printf("Error. Input file was not specified.\n");
						printf("Usage: -compute bto/bfo/dist -file inpfilename\n");
						return;
					}
					else
					{
						if (strlen(argv[argn+3])>1024)
						{
							printf("Error: too long filename.\n");
							return;
						}
						strcpy(inpfilename, argv[argn+3]);
						argn+=4;
						strcpy(outfilename,"");
						continue;
					}
				}


				if (strcmp(argv[argn+1], "bto")==0 || strcmp(argv[argn+1], "BTO")==0)
				{
					if (argn+5 >= argc)
					{
						printf("Error. Missing arguments. Usage: -compute bto time lon lat alt.\n");
						return;
					}

					if (strlen(argv[argn+2])>63 || strlen(argv[argn+3])>63  || strlen(argv[argn+4])>63  || strlen(argv[argn+5])>63 )
					{
						printf("Unsupported argument format.\n");
						return;
					}
					
					strcpy(str_arg1, argv[argn+2]);
					strcpy(str_arg2, argv[argn+3]);
					strcpy(str_arg3, argv[argn+4]);
					strcpy(str_arg4, argv[argn+5]);
					strcpy(inpfilename,"");
					strcpy(outfilename,"");
					argn+=6;
					continue;
				}

				if (strcmp(argv[argn+1], "bfo")==0 || strcmp(argv[argn+1], "BFO")==0)
				{
					if (argn+8 >= argc)
					{
						printf("Error. Missing arguments. Usage: -compute bfo time lon lat alt u v w. Optional: -AES_SPH, -AES_WGS, -AES_ALT xxxx.\n");
						return;
					}

					if (strlen(argv[argn+2])>63 || strlen(argv[argn+3])>63  || strlen(argv[argn+4])>63  || strlen(argv[argn+5])>63 || strlen(argv[argn+6])>63 || strlen(argv[argn+7])>63 || strlen(argv[argn+8])>63)
					{
						printf("Unsupported argument format.\n");
						return;
					}
					strcpy(str_arg1, argv[argn+2]);
					strcpy(str_arg2, argv[argn+3]);
					strcpy(str_arg3, argv[argn+4]);
					strcpy(str_arg4, argv[argn+5]);
					strcpy(str_arg5, argv[argn+6]);
					strcpy(str_arg6, argv[argn+7]);
					strcpy(str_arg7, argv[argn+8]);
					strcpy(inpfilename,"");
					strcpy(outfilename,"");
					argn+=9;
					continue;
				}
			}
		}

		// --------------------
		// AES Model - spherical or wgs'84 ellipsoid (optional parameter)
		if (strcmp(argv[argn],"-AES_SPH")==0 || strcmp(argv[argn],"-AES_WGS")==0)
		{
			if (action == ACT_COM)
			{
				strcpy(str_arg8, argv[argn]);
			}
			else
			{
				strcpy(str_arg8, "");
			}
			argn++;
			continue;
		}

		// check if AES altitude is specified (optional parameter)
		if (strstr(argv[argn],"-AES_ALT"))
		{
			if (action == ACT_COM)
			{
				if (strlen(argv[argn])<64)
				{
					strcpy(str_arg9, argv[argn]);
				}
				else
				{
					printf("Ignoring too long argument %s.\n", argv[argn]);
				}
			}
			else
			{
				strcpy(str_arg9, "");
			}
			argn++;
			continue;
		}


		// --------------------
		// Set optimization weight for distance (BTO)
		if (strcmp(argv[argn],"-optweight1")==0)
		{
			if (action == ACT_OPT)
			{
				if (argn+1<argc)
				{
					if (strlen(argv[argn+1])<63)
					{
						strcpy(str_arg1, argv[argn+1]);
						argn+=2;
						continue;
					}
					else
					{
						printf("Parameter is too long - skipping.\n");
						strcpy(str_arg1, "");
						argn+=2;
						continue;
					}
				}
				else
				{
					strcpy(str_arg1, "");
				}
			}
			else
			{
				strcpy(str_arg1, "");
			}
			argn++;
			continue;
		}


		// --------------------
		// Set optimization weight for distance (BFO)
		if (strcmp(argv[argn],"-optweight2")==0)
		{
			if (action == ACT_OPT)
			{
				if (argn+1<argc)
				{
					if (strlen(argv[argn+1])<63)
					{
						strcpy(str_arg2, argv[argn+1]);
						argn+=2;
						continue;
					}
					else
					{
						printf("Parameter is too long - skipping.\n");
						strcpy(str_arg2, "");
						argn+=2;
						continue;
					}
				}
				else
				{
					strcpy(str_arg2, "");
				}
			}
			else
			{
				strcpy(str_arg2, "");
			}
			argn++;
			continue;
		}

		// --------------------
		// verbose or silent
		if (strcmp(argv[argn],"-verbose")==0)
		{
			isVerbose = true;
			argn++;
			continue;
		}

		printf("Skipping unrecognized option: %s\n", argv[argn]);
		argn++;
	}



	// ---------------------------------------
	// Execute action based on command
	switch (action)
	{
	case ACT_UDF:
		printf("Error. Incorrect option. Exiting.\n");
		break;

	case ACT_PRT:
		PrintPings(outfilename);
		break;

	case ACT_TRJ:
		BuildTrajectory(inpfilename, isVerbose);
		break;

	case ACT_CAL:
		Calibrate(outfilename, isVerbose);
		break;

	case ACT_VAL:
		Validate(inpfilename, isVerbose);
		break;

	case ACT_OPT:
		Optimize(optname, inpfilename, outfilename, str_arg1, str_arg2, isVerbose);
		break;

	case ACT_MET:
		GetMeteo(optname, meteodir, str_arg1, str_arg2, str_arg3, str_arg4, inpfilename, outfilename, isVerbose);
		break;

	case ACT_COM:
		if (strcmp(optname,"bto") == 0 || strcmp(optname,"BTO") == 0) ComputeBTO(inpfilename, outfilename, str_arg1, str_arg2, str_arg3, str_arg4);
		if (strcmp(optname,"bfo") == 0 || strcmp(optname,"BFO") == 0) ComputeBFO(inpfilename, outfilename, str_arg1, str_arg2, str_arg3, str_arg4, str_arg5, str_arg6, str_arg7, str_arg8, str_arg9);
		break;

	default:
		printf("Internal error. Incorrect option. Exiting.\n");
		break;
	}
}; // End main


