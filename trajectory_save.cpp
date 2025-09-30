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
#include "trajectory.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Save trajectory; if filename = NULL, then print it into outputstream
//
///////////////////////////////////////////////////////////////////////////////////////////////////

extern void print_time_formatted(char* buf, double t);

bool CTrajectory::SaveTrajectory(char* filename)
{
	unsigned long n = GetNOut();
	double* pTime = GetTimeOut();
	double* pLon =  GetLonOut();
	double* pLat = GetLatOut();
	double* pAlt =  GetAltOut();
	double* pM =    GetMOut();
	double* pu =    GetuOut();
	double* pv =    GetvOut();
	double* pw =   GetwOut();
	double* pq =   GetqOut();

	FILE* fid = NULL;

	if (filename!=NULL)
	{
		fid = fopen(filename,"wt");
		if (!fid) return false;
	}
	else
	{
		fid = fileout_;
	}

	fprintf(fid,"Time(s), Lon(degE), Lat(degN), Alt(m), U(m/s), V(m/s), W(m/s), Weight(kg), FF(kg/h)\n");

//		for (unsigned long i=0; i<n; i++)
//		{
//			fprintf(fid, "%11.3f %7.3f %7.3f %7.1f %9.1f %7.1f %7.1f %6.2f %7.1f\n", pTime[i], pLon[i], pLat[i], pAlt[i], pM[i], pu[i], pv[i], pw[i], pq[i]*3600.0);
//		}

	char buf[256];

	for (unsigned long i=0; i<n; i++)
	{
		print_time_formatted(buf, pTime[i]);
		fprintf(fid, "%s  %7.3f %7.3f %7.1f %7.1f %7.1f %6.2f %9.1f %7.1f\n", buf, pLon[i], pLat[i], pAlt[i], pu[i], pv[i], pw[i], pM[i], pq[i]*3600.0);
	}

	if (filename!=NULL)
	{
		fclose(fid);
	}

	return true;
}


