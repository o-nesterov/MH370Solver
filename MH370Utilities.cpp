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
#include <sys/stat.h>

#include "inmarsatdata.h"
#include "environdata.h"
#include "params_parser.h"


extern void uv2uvw(double& XYZ_U, double& XYZ_V, double& XYZ_W, double u, double v, double lon, double lat);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary function to print human-readable date
//
void print_time_formatted(char* buf, double t)
{
	double tt = t;
	int yyyy=2014;
	int mm=3;
	int dd=7+(int)floor(tt/86400.0);
	tt-=(dd-7)*86400.0;
	int hh=(int)floor(tt/3600.0);
	tt-=hh*3600.0;
	int mn=(int)floor(tt/60.0);
	tt-=mn*60.0;
	int ss = (int)floor(tt);
	tt-=ss;
	int ds = (int)floor(tt*1000.0+0.5);
	if (ds==1000) {ss++; ds=0;};
	if (ss==60) {mn++; ss=0;};
	if (mn==60) {hh++; mn=0;};
	if (hh==24) {dd++; hh=0;};
	sprintf(buf, "%4.4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d.%3.3d", yyyy, mm, dd, hh, mn, ss, ds);
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary function to convert time from human-readable date yyyy-mm-dd hh:mm:ss.sss to double
//
bool read_time_formatted(double& t, char* buf)
{
	char* p = &buf[0];
	unsigned long n = strlen(buf);

	while(((*p==' ') || (*p=='	')) && p<=(&buf[n-1])) p++; // skip spaces

	char* pp;
	unsigned long yyyy = strtoul(p,&pp,10);
	if (p+4!=pp || *pp!='-' || pp>=&buf[n-1]) return false;
	p=pp+1;

	unsigned long mm = strtoul(p,&pp,10);
	if (p+2!=pp || *pp!='-' || pp>=&buf[n-1]) return false;
	p=pp+1;

	unsigned long dd = strtoul(p,&pp,10);
	if (p+2!=pp || *pp!=' ' || pp>=&buf[n-1]) return false;
	p=pp+1;

	while(((*p==' ') || (*p=='	')) && p<(&buf[n-1])) p++; // skip spaces

	unsigned long hh = strtoul(p,&pp,10);
	if (p+2!=pp || *pp!=':' || pp>=&buf[n-1]) return false;
	p=pp+1;

	unsigned long mn = strtoul(p,&pp,10);
	if (p+2!=pp || *pp!=':' || pp>=&buf[n-1]) return false;
	p=pp+1;

	double ss = strtod(p,&pp);
	if (p==pp) return false;

	// here everything is fine; all were read; do some sanity checks
	if (yyyy>9999) return false;
	if (mm>12 || mm==0) return false;
	if (dd>31 || dd==0) return false;
	if (hh>=24) return false;
	if (mn>=60) return false;
	if (ss>=60.0) return false;


	double tt = (hh*60+mn)*60+ss;	// seconds since the beginning of the day yyyy,mm,dd

	// compute number of days since beginning of the year
	// determine if this is a leap year or not
	bool isLeapYear;
	unsigned long yyyy4 = yyyy/4;
	unsigned long yyyy100 = yyyy/100;
	unsigned long yyyy1000 = yyyy/1000;

	if (yyyy!=yyyy4)
	{
		isLeapYear = false;
	}
	else
	{
		if (yyyy != yyyy100*100)
		{
			isLeapYear = true; // divisible by 4, but not by 100
		}
		else
		{
			if (yyyy != yyyy1000*1000)
			{
				isLeapYear = false;	// divisible by 100, but not by 1000
			}
			else
			{
				isLeapYear = true;	// divisible by 1000
			}
		}
	}
	// number of days in a month
	int ndays_in_month[] = {31,28,31,30,31,30,31,31,30,31,30,31};
	if (isLeapYear) ndays_in_month[1] = 29;
	if ((int)dd>ndays_in_month[mm-1]) return false;

	// count number of days since the beginning of the year
	int ndays = 0;
	for (int i=1; i<(int)mm; i++) ndays += ndays_in_month[i-1];
	ndays += (dd-1);

	// add number of days since 0000
	ndays += (365*yyyy + yyyy4 - yyyy100 + yyyy1000); // number of days in a regular year + number of leap years

	// Reference date: 2014-03-07 00:00:00
	int nday20140307 = 365*2014 + 503 - 20 + 2 + (31 + 28 + 6);

	ndays -= nday20140307;

	t = 86400.0*ndays + tt;

	return true;
}



//------------------------------------------------------------------------------------
// scan value of varname in line, depending on vartype, save it in var
bool parseLine(void* var, char* varname, char* pline, char* vartype)
{
	// Input: varname - name of variable
	//		  pline - string
	//		  vartype - type of variable

	// Output: var - scanned value
	//		   return false if variable is not found or scan fails, true if success


	// scan if variable name is present in text line
	char* p = strstr(pline, varname);
	if (p==NULL) return false;

	int n= (int)strlen(pline);
	p += strlen(varname);	// skip varname itself

	while(((*p==' ') || (*p=='	')) && p<(&pline[n])) p++;	// skip spaces
	if (p==&pline[n]) return false;

	if (*p !='=') return false; // varname must be "="
	p++;

	while(((*p==' ') || (*p=='	')) && p<(&pline[n])) p++;	// skip spaces and "
	if (p==&pline[n]) return false;

	// check if quotation marks are present
	char* ps = NULL;
	char* pe = NULL;
	if (*p == '"')
	{
		ps = p+1;
		pe = ps;
		while (*pe!='"' && pe<(&pline[n])) pe++;
		if (pe==ps || pe==&pline[n]) return false;	// the second quotation mark is absent
	}

	// read double fp 
	if (strcmp(vartype,"double")==0)
	{
		char* pp;
		double tmp = strtod(p,&pp);
		if (p!=pp)
		{
			*((double*)var) = tmp;
			return true;
		}
	}

	// read int 
	if (strcmp(vartype,"int")==0)
	{
		char* pp;
		int tmp = strtol(p,&pp,10);
		if (p!=pp)
		{
			*((int*)var) = tmp;
			return true;
		}
	}

	// allocate and read string 
	if (strcmp(vartype,"string")==0)
	{
		int nn = (int)strlen(p);
		char* tmp = new char[nn+1];

		if (ps)
		{
			// copy string, save for quotation marks
			p = ps;
			while (p<pe) {tmp[((unsigned long)p)-((unsigned long)ps)] = *p; p++;}
			tmp[((unsigned long)pe)-((unsigned long)ps)]='\0'; // put of string here
		}
		else
		{
			// check if parentheses are present
			char* pths = NULL;
			char* pthe = NULL;
			if (*p == '{')
			{
				pths = p+1;
				pthe = pths;
				while (*pthe!='}' && pe<(&pline[n])) pthe++;
				if (pths==pthe || pthe==&pline[n]) pthe = NULL;	// the second quotation mark is absent
			}

			// if {} are present, then return content  
			if (pths!=NULL && pthe!=NULL)
			{
				p = pths;
				while (p<pthe) {tmp[((unsigned long)p)-((unsigned long)pths)] = *p; p++;}
				tmp[((unsigned long)pthe)-((unsigned long)pths)]='\0'; // put of string here
			}
			else
			{
				// find end of line
				strcpy(tmp,p);
				char*pp = &tmp[0];
				while(((*pp!=' ' && *pp!='	' && *pp!='#' && *pp!='%' && *pp!='!') && *pp!='\n') && pp<(&tmp[n])) pp++;
				*pp='\0'; // replace ", #, ... with the end of line
			}
		}

		*((char**)var) = tmp;
		return true;
	}

	// time: it can be either float (as s since 2014-03-07 00:00:00 UTC) or as full date
	if (strcmp(vartype,"time")==0)
	{
		if (ps==NULL) // no quotation marks
		{
			// attempt to read time as floating point
			char* pp;
			double tmp = strtod(p,&pp);
			if (p==pp) return false; // nothing was read; unrecognized format
			if (*pp!='-')
			{
				*((double*)var) = tmp; // some FP value was read
				return true;
			}
			else
			{
				bool isDate = read_time_formatted(tmp, p); // try to read it as text date
				if (isDate) {*((double*)var) = tmp;};
				return isDate;
			}
		}
		else
		{
			double tmp;
			bool isDate = read_time_formatted(tmp, p); // try to read it as text date
			if (isDate) {*((double*)var) = tmp;};
			return isDate;
		}
	}

	return false;
};



//------------------------------------------------------------------------------------
// scan FP value of varname in line along with attribute character
bool parseLineA(double& var, char& attr, char* varname, char* pline)
{
	// Input: varname - name of variable
	//		  pline - string

	// Output: var - scanned value
	//		   attr - attribute character
	//		   return false if variable is not found or scan fails, true if success


	// scan if variable name is present in text line
	char* p = strstr(pline, varname);
	if (p==NULL) return false;

	int n= (int)strlen(pline);
	p += strlen(varname);	// skip varname itself

	while(((*p==' ') || (*p=='	')) && p<(&pline[n])) p++;	// skip spaces
	if (p==&pline[n]) return false;

	if (*p !='=') return false; // varname must be "="
	p++;

	while(((*p==' ') || (*p=='	')) && p<(&pline[n])) p++;	// skip spaces and "
	if (p==&pline[n]) return false;

	char* pp;
	var = strtod(p,&pp);
	if (p==pp)
	{
		// try to read prefix character attribute first
		if (((*p)>='A' && (*p)<='Z') || ((*p)>='a' && (*p)<='z'))
		{
			attr = *p;
			p++;
			var = strtod(p,&pp);
			if (p==pp) return false; // incorrect format
			return true;
		}
		else
		{
			return false; // incorrect format
		}
	}
	else
	{
		if (((*pp)>='A' && (*pp)<='Z') || ((*pp)>='a' && (*pp)<='z')) // try to read suffix attribute
		{
			attr = *pp;
			return true;
		}
		else
		{
			attr = 0; // fp variable was read, but there is neither prefix nor suffix
			return true; // still return true 
		}
	}

	return false;
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary function to check if index is in the custom list
//
bool isPingIndexInList(CInpParams* pInpParams, int idx)
{
	bool isIn = false;
	if (pInpParams->nOutputPingIDX==-1 && pInpParams->pOutputPingIDX == NULL) return true; // all suitable pings were specified

	if (pInpParams->pOutputPingIDX != NULL)
	{
		// check if ping index is in the custom list
		for (int j=0; j<pInpParams->nOutputPingIDX; j++) if ((pInpParams->pOutputPingIDX)[j]==idx) isIn = true;
	}
	return isIn;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary function to print full list of Inmarsat ping data
//
void PrintPings(char* filename)
{
	FILE* fid = NULL;
	char buf[256];
	char channel_name[256];

	if (strlen(filename)>0)
	{
		fid = fopen(filename, "wt");
	}

	CMH370InmarsatData* pInmarsatData = new CMH370InmarsatData();

	unsigned long n = pInmarsatData->GetNumberOfPings();

	// Print header
	if (fid)
	{
		fprintf(fid, "Ping No., Time(UTC), Time(s), BTO(microsecond), Distance(km), BFO(Hz), Channel\n");
	}
	else
	{
		printf("Ping No., Time(UTC), Time(s), BTO(microsecond), Distance(km), BFO(Hz), Channel\n");
	}

	// Print all Inmarsat pings (which have BTO/BFO data)
	for (unsigned long i=0; i<n; i++)
	{
		double ts, dist, bto, bfo;
		int channelType;
		pInmarsatData->GetPingData(ts, dist, bto, bfo, channelType, i);
		print_time_formatted(buf, ts);

		switch(channelType)
		{
		case CH_TYPE_R_RX1200:
			sprintf(channel_name, "RX1200-R");
			break;

		case CH_TYPE_R_RX600:
			sprintf(channel_name, "RX600-R");
			break;

		case CH_TYPE_T_RX1200:
			sprintf(channel_name, "RX1200-T");
			break;

		case CH_TYPE_C_RX1200:
			sprintf(channel_name, "RX1200-C");
			break;

		case CH_TYPE_ANAMALOUS:
			sprintf(channel_name, "RX1200-R(anamalous)");
			break;

		default:
			sprintf(channel_name, "");
			break;
		}



		// make output
		if (fid)
		{
			if ((channelType!=CH_TYPE_C_RX1200) && (channelType!=CH_TYPE_ANAMALOUS))
			{
				fprintf(fid, "%3.3d   %s  %6.0f  %9.3f  %4.0f   %s\n", i, buf, bto, dist*0.001, bfo, channel_name);
			}
			else
			{
				if (channelType==CH_TYPE_C_RX1200)
				{
					fprintf(fid, "%3.3d   %s  ******  *********  %4.0f   %s\n", i, buf, bfo, channel_name);
				}
				else
				{
					fprintf(fid, "%3.3d   %s  %6.0f  *********  %4.0f   %s\n", i, buf, bto, dist*0.001, bfo, channel_name);
				}
			}
		}
		else
		{
			if ((channelType!=CH_TYPE_C_RX1200) && (channelType!=CH_TYPE_ANAMALOUS))
			{
				printf("%3.3d   %s  %6.0f  %9.3f  %4.0f   %s\n", i, buf, bto, dist*0.001, bfo, channel_name);
			}
			else
			{
				if (channelType==CH_TYPE_C_RX1200)
				{
					printf("%3.3d   %s  ******  *********  %4.0f   %s\n", i, buf, bfo, channel_name);
				}
				else
				{
					printf("%3.3d   %s  %6.0f  *********  %4.0f   %s\n", i, buf, bto, bfo, channel_name);
				}
			}
		}
	}


	if (fid) fclose(fid);	// close file

	delete pInmarsatData;
}





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Load {time, lon, lat, alt} data from file
bool loadLonLatAlt(int& nOut, double* &pTs, double* &pLon, double* &pLat, double* &pAlt, char* inpfilename)
{
	// initialization
	nOut = 0;
	pTs = NULL;
	pLon = NULL;
	pLat = NULL;
	pAlt = NULL;

	int nr = 0; //actual number of records
	int nsize = 0; // size of array
	const int chunksize = 256; // chunk size to increase array size (to avoid reallocation every time).

	FILE* fid = fopen(inpfilename, "rt");
	if (!fid) return false;
	
	char buf[1024];

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
			// try to read lon
			while((*p==' ') || (*p=='	') || (*p==',') || (*p==';') && p<(&buf[n])) p++;	// skip spaces etc.
			double lon = strtod(p, &pp);
			if (p == pp) continue; // unknown lon format
			p = pp+1;

			// -------------------------
			// try to read lat
			while((*p==' ') || (*p=='	') || (*p==',') || (*p==';') && p<(&buf[n])) p++;	// skip spaces etc.
			double lat = strtod(p, &pp);
			if (p == pp) continue; // unknown lat format
			p = pp+1;

			// -------------------------
			// try to read alt
			while((*p==' ') || (*p=='	') || (*p==',') || (*p==';') && p<(&buf[n])) p++;	// skip spaces etc.
			double alt = strtod(p, &pp);
			if (p == pp) continue; // unknown alt format
			p = pp+1;

			// everything is ok here; {ts, lon, lat, alt} were read
			if (nr==nsize)
			{
				// reallocate arrays and add new data record
				double* pTs_new = new double [nsize+chunksize];
				double* pLon_new = new double [nsize+chunksize];
				double* pLat_new = new double [nsize+chunksize];
				double* pAlt_new = new double [nsize+chunksize];

				if (pTs_new == NULL || pLon_new == NULL || pLat_new == NULL || pAlt_new == NULL)
				{
					if (pTs) delete pTs;
					if (pLon) delete pLon;
					if (pLat) delete pLat;
					if (pAlt) delete pAlt;

					nOut = 0;
					pTs = NULL;
					pLon = NULL;
					pLat = NULL;
					pAlt = NULL;

					fclose(fid);
					return false;
				}

				for (int i=0; i<nr; i++)
				{
					pTs_new[i] = pTs[i];
					pLon_new[i] = pLon[i];
					pLat_new[i] = pLat[i];
					pAlt_new[i] = pAlt[i];
				}

				// add new record
				pTs_new[nr] = ts;
				pLon_new[nr] = lon;
				pLat_new[nr] = lat;
				pAlt_new[nr] = alt;

				// delete old array
				if (pTs) delete pTs;
				if (pLon) delete pLon;
				if (pLat) delete pLat;
				if (pAlt) delete pAlt;

				// reassign arrays
				pTs = pTs_new;
				pLon = pLon_new;
				pLat = pLat_new;
				pAlt = pAlt_new;

				nr++;
				nsize += chunksize;
			}
			else
			{
				// add new record
				pTs[nr] = ts;
				pLon[nr] = lon;
				pLat[nr] = lat;
				pAlt[nr] = alt;
				nr++;
			}
		} // fgets
	}; //eof

	fclose(fid);

	if (nr==0)
	{
		if (pTs) delete pTs;
		if (pLon) delete pLon;
		if (pLat) delete pLat;
		if (pAlt) delete pAlt;

		nOut = 0;
		pTs = NULL;
		pLon = NULL;
		pLat = NULL;
		pAlt = NULL;

		return false;
	}

	nOut = nr;

	return true;
}





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Load {time, lon, lat, alt, u, v, w} data from file
bool loadLonLatAltVel(int& nOut, double* &pTs, double* &pLon, double* &pLat, double* &pAlt, double* &pU, double* &pV, double* &pW, char* inpfilename)
{
	// initialization
	nOut = 0;
	pTs = NULL;
	pLon = NULL;
	pLat = NULL;
	pAlt = NULL;
	pU = NULL;
	pV = NULL;
	pW = NULL;

	int nr = 0; //actual number of records
	int nsize = 0; // size of array
	const int chunksize = 256; // chunk size to increase array size (to avoid reallocation every time).

	FILE* fid = fopen(inpfilename, "rt");
	if (!fid) return false;
	
	char buf[1024];

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
			// try to read lon
			while((*p==' ') || (*p=='	') || (*p==',') || (*p==';') && p<(&buf[n])) p++;	// skip spaces etc.
			double lon = strtod(p, &pp);
			if (p == pp) continue; // unknown lon format
			p = pp+1;

			// -------------------------
			// try to read lat
			while((*p==' ') || (*p=='	') || (*p==',') || (*p==';') && p<(&buf[n])) p++;	// skip spaces etc.
			double lat = strtod(p, &pp);
			if (p == pp) continue; // unknown lat format
			p = pp+1;

			// -------------------------
			// try to read alt
			while((*p==' ') || (*p=='	') || (*p==',') || (*p==';') && p<(&buf[n])) p++;	// skip spaces etc.
			double alt = strtod(p, &pp);
			if (p == pp) continue; // unknown alt format
			p = pp+1;

			// -------------------------
			// try to read u
			while((*p==' ') || (*p=='	') || (*p==',') || (*p==';') && p<(&buf[n])) p++;	// skip spaces etc.
			double vel_u = strtod(p, &pp);
			if (p == pp) continue; // unknown U format
			p = pp+1;

			// -------------------------
			// try to read v
			while((*p==' ') || (*p=='	') || (*p==',') || (*p==';') && p<(&buf[n])) p++;	// skip spaces etc.
			double vel_v = strtod(p, &pp);
			if (p == pp) continue; // unknown V format
			p = pp+1;

			// -------------------------
			// try to read w
			while((*p==' ') || (*p=='	') || (*p==',') || (*p==';') && p<(&buf[n])) p++;	// skip spaces etc.
			double vel_w = strtod(p, &pp);
			if (p == pp) continue; // unknown W format
			p = pp+1;

			// everything is ok here; {ts, lon, lat, alt, vel_u, vel_v, vel_w} were read
			if (nr==nsize)
			{
				// reallocate arrays and add new data record
				double* pTs_new = new double [nsize+chunksize];
				double* pLon_new = new double [nsize+chunksize];
				double* pLat_new = new double [nsize+chunksize];
				double* pAlt_new = new double [nsize+chunksize];
				double* pU_new = new double [nsize+chunksize];
				double* pV_new = new double [nsize+chunksize];
				double* pW_new = new double [nsize+chunksize];

				if (pTs_new == NULL || pLon_new == NULL || pLat_new == NULL || pAlt_new == NULL || pU_new == NULL || pV_new==NULL || pW_new==NULL)
				{
					if (pTs) delete pTs;
					if (pLon) delete pLon;
					if (pLat) delete pLat;
					if (pAlt) delete pAlt;
					if (pU) delete pU;
					if (pV) delete pV;
					if (pW) delete pW;

					nOut = 0;
					pTs = NULL;
					pLon = NULL;
					pLat = NULL;
					pAlt = NULL;
					pU = NULL;
					pV = NULL;
					pW = NULL;

					fclose(fid);
					return false;
				}

				for (int i=0; i<nr; i++)
				{
					pTs_new[i] = pTs[i];
					pLon_new[i] = pLon[i];
					pLat_new[i] = pLat[i];
					pAlt_new[i] = pAlt[i];
					pU_new[i] = pU[i];
					pV_new[i] = pV[i];
					pW_new[i] = pW[i];
				}

				// add new record
				pTs_new[nr] = ts;
				pLon_new[nr] = lon;
				pLat_new[nr] = lat;
				pAlt_new[nr] = alt;
				pU_new[nr] = vel_u;
				pV_new[nr] = vel_v;
				pW_new[nr] = vel_w;

				// delete old array
				if (pTs) delete pTs;
				if (pLon) delete pLon;
				if (pLat) delete pLat;
				if (pAlt) delete pAlt;
				if (pU) delete pU;
				if (pV) delete pV;
				if (pW) delete pW;

				// reassign arrays
				pTs = pTs_new;
				pLon = pLon_new;
				pLat = pLat_new;
				pAlt = pAlt_new;
				pU = pU_new;
				pV = pV_new;
				pW = pW_new;

				nr++;
				nsize += chunksize;
			}
			else
			{
				// add new record
				pTs[nr] = ts;
				pLon[nr] = lon;
				pLat[nr] = lat;
				pAlt[nr] = alt;
				pU[nr] = vel_u;
				pV[nr] = vel_v;
				pW[nr] = vel_w;
				nr++;
			}
		} // fgets
	}; //eof

	fclose(fid);

	if (nr==0)
	{
		if (pTs) delete pTs;
		if (pLon) delete pLon;
		if (pLat) delete pLat;
		if (pAlt) delete pAlt;
		if (pU) delete pU;
		if (pV) delete pV;
		if (pW) delete pW;

		nOut = 0;
		pTs = NULL;
		pLon = NULL;
		pLat = NULL;
		pAlt = NULL;
		pU = NULL;
		pV = NULL;
		pW = NULL;

		return false;
	}

	nOut = nr;

	return true;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Extract meteorological parameters {wind u, wind v, air t, air pressure and relative humidity} 
// at a given 3D location and time
void GetMeteo(char* datasetname, char* pMeteoDir, char*str_time, char*str_lon, char*str_lat, char*str_alt, char* inpfilename, char* outfilename, bool isVerbose)
{
	double ts, lon, lat, alt;

	// initialize file data (in case opted)
	int nTs = -1;
	double* pTs = NULL;
	double* pLon = NULL;
	double* pLat = NULL;
	double* pAlt = NULL;

	char buf[256]; // buffer for time text

	// check if dataset is recognized
	int iMeteoDataset= METEO_DATASET_NULL;
	if (strcmp(datasetname,"GDAS1") == 0 || strcmp(datasetname,"gdas1") == 0 || strcmp(datasetname,"Gdas1") == 0) {iMeteoDataset = METEO_DATASET_GDAS1;};
	if (strcmp(datasetname,"ERA5") == 0 || strcmp(datasetname,"era5") == 0 || strcmp(datasetname,"Era5") == 0) {iMeteoDataset = METEO_DATASET_ERA5;};
	if (iMeteoDataset == METEO_DATASET_NULL)
	{
		printf("Error: unrecognized meteorological dataset.\n");
		return;
	}

	// -----------
	// check dirname is valid directory
	bool isDir = false;

	struct _stat info;

	if(_stat(pMeteoDir, &info ) != 0)
	{
		isDir= false;
	}
	else if(info.st_mode & S_IFDIR)
	{
		isDir= true;
	}

	if (!isDir)
	{
		printf("Erroneous directory containing meteorological dataset.\n");
		return;
	}


	if (strcmp(str_time,"file")==0 || strcmp(str_time,"-file")==0) 
	{
		// data is in file
		if (!loadLonLatAlt(nTs, pTs, pLon, pLat, pAlt, inpfilename))
		{
			if (pTs) delete pTs;
			if (pLon) delete pLon;
			if (pLat) delete pLat;
			if (pAlt) delete pAlt;
			printf("Missing or erroneous input file: %s.\n", inpfilename);
			return;
		}
	}

	// -----------
	if ((nTs==-1) && strlen(str_time)>0)
	{
		// non-interactive mode, and not reading {time, lon, lat, alt} data from file
		char* pp;
		ts = strtod(str_time, &pp);
		if (&str_time[0] == pp)
		{
			printf("Error. Wrong time.\n");
			return;
		}
		if (*pp=='-')
		{
			// time seems to be given as a string yyyy-....
			if (!read_time_formatted(ts, str_time))
			{
				printf("Error. Unrecognized time argument.\n");
				printf("Must be either s since 2014-03-07 00:00:00 UTC, or date yyyy-mm-dd hh:mm:ss.sss\n");
				return;
			}
		}

		lon = strtod(str_lon, &pp);
		if (&str_lon[0] == pp)
		{
			printf("Error. Wrong longitude.\n");
			return;
		}

		lat = strtod(str_lat, &pp);
		if (&str_lat[0] == pp)
		{
			printf("Error. Wrong latitude.\n");
			return;
		}

		alt = strtod(str_alt, &pp);
		if (&str_alt[0] == pp)
		{
			printf("Error. Wrong altitude.\n");
			return;
		}
	}

	// Load meteorological data:
	CMH370EnvironData* pMH370EnvironData = new CMH370EnvironData();
	if (nTs>=0)
	{
		if (isVerbose) pMH370EnvironData -> SetVerbose(stdout,1);
	}
	else
	{
		pMH370EnvironData -> SetVerbose(stdout,1); // assume verbose (interactive mode, and explicit specs)
	}
	pMH370EnvironData -> LoadMeteo(iMeteoDataset, pMeteoDir);
	pMH370EnvironData -> SetMeteoInterpMethod(3);


	// ------------------------------
	// Check if file output is selected
	if (nTs>=0)
	{
		// This is file output
		FILE* fid = fopen(outfilename, "wt");
		if (!fid)
		{
			delete pMH370EnvironData;
			if (pTs) delete pTs;
			if (pLon) delete pLon;
			if (pLat) delete pLat;
			if (pAlt) delete pAlt;
			if (isVerbose) printf("Error opening output file: %s\n", outfilename);
			return;
		}

		// Print header
		fprintf(fid,"Time(UTC), Lon(degE), Lat(degN), Alt(m), WindU(m/s), WindV(m/s), AirT(degC), AirP(Pa), RH(%%)\n");

		for (int i=0; i<nTs; i++)
		{
			double u,v,t,p,r;
			pMH370EnvironData->GetMeteo(u,v,t,p,r, pTs[i], pLon[i], pLat[i], pAlt[i]);
			print_time_formatted(buf, pTs[i]);
			fprintf(fid,"%s  %7.3f  %7.3f  %7.1f  %6.2f  %6.2f  %6.2f  %7.0f  %6.1f\n", buf, pLon[i], pLat[i], pAlt[i], u, v ,t ,p, r); 
		}

		fclose(fid);
		delete pMH370EnvironData;
		if (pTs) delete pTs;
		if (pLon) delete pLon;
		if (pLat) delete pLat;
		if (pAlt) delete pAlt;
		
		if (isVerbose) printf("Extraction completed.\n");

		return;
	}

	// ------------------------------
	// Check non-interactive mode is selected
	if (strlen(str_time)>0)
	{
		// non interactive mode
		double u,v,t,p,r;
		pMH370EnvironData->GetMeteo(u,v,t,p,r, ts,lon,lat,alt);

		printf("WindU: %5.2f (m/s)\n", u);
		printf("WindV: %5.2f (m/s)\n", v);
		printf("AirT: %5.2f (deg C)\n", t);
		printf("Air pressure: %6.0f (Pa)\n", p);
		printf("Relative humidity: %5.1f %%\n", r);
	}
	else
	{
		// Interactive mode is selected
		printf("Interactive mode. To exit type x+Enter or X+Enter anytime.\n");
		printf("To restart input type r+Enter or R+Enter anytime.\n");
		char buf[256];

		bool isInteractive = true;

		while(isInteractive)
		{
			// read time
			printf("Enter time (s since 2014-03-07 00:00:00 UTC, or date yyyy-mm-dd hh:mm:ss.sss):\n");
			while(1)
			{
				scanf("%s", buf, 256);
				char* pp;
				ts = strtod(buf, &pp);
				if (buf[0]=='x' || buf[0]=='X') {isInteractive = false; break;};
				if (buf[0]=='r' || buf[0]=='R') {break;};
				if (&buf[0] == pp)
				{
					printf("Error. Wrong time format.\n");
					printf("Must be either s since 2014-03-07 00:00:00 UTC, or date yyyy-mm-dd hh:mm:ss.sss\n");
					continue;
				}
				if (*pp=='-')
				{
					// time seems to be given as a string yyyy-....
					if (!read_time_formatted(ts, str_time))
					{
						printf("Error. Unrecognized time format.\n");
						printf("Must be either s since 2014-03-07 00:00:00 UTC, or date yyyy-mm-dd hh:mm:ss.sss\n");
						continue;
					}
				}
				if (ts<43200.0) printf("Warning: time before 2014-03-07 12:00:00?\n");
				if (ts>93600.0) printf("Warning: time after 2014-03-08 02:00:00?\n");
				break; // everything is fine here
			}
			if (!isInteractive) break;
			if (buf[0]=='r' || buf[0]=='R') continue;	// restart input


			// read longitude
			printf("Enter longitude (deg E):\n");
			while(1)
			{
				scanf("%s", buf, 256);
				char* pp;
				lon = strtod(buf, &pp);
				if (buf[0]=='x' || buf[0]=='X') {isInteractive = false; break;};
				if (buf[0]=='r' || buf[0]=='R') {break;};
				if (&buf[0] == pp)
				{
					printf("Error. Wrong longitude.\n");
					continue;
				}
				if (lon<0.0 || lon>=360.0)
				{
					printf("Error. Wrong longitude (range 0 to 360 deg).\n");
					continue;
				}
				break; // everything is fine here
			}
			if (!isInteractive) break;
			if (buf[0]=='r' || buf[0]=='R') continue;	// restart input

			// read latitutde
			printf("Enter latitude (deg N):\n");
			while(1)
			{
				scanf("%s", buf, 256);
				char* pp;
				lat = strtod(buf, &pp);
				if (buf[0]=='x' || buf[0]=='X') {isInteractive = false; break;};
				if (buf[0]=='r' || buf[0]=='R') {break;};
				if (&buf[0] == pp)
				{
					printf("Error. Wrong latitude.\n");
					continue;
				}
				if (lat<-90.0 || lat>90.0)
				{
					printf("Error. Wrong latitude (range -90 to 90 deg).\n");
					continue;
				}
				break; // everything is fine here
			}
			if (!isInteractive) break;
			if (buf[0]=='r' || buf[0]=='R') continue;	// restart input

			// read altitude
			printf("Enter altitude (m):\n");
			while(1)
			{
				scanf("%s", buf, 256);
				char* pp;
				alt = strtod(buf, &pp);
				if (buf[0]=='x' || buf[0]=='X') {isInteractive = false; break;};
				if (buf[0]=='r' || buf[0]=='R') {break;};
				if (&buf[0] == pp)
				{
					printf("Error. Wrong altitude.\n");
					continue;
				}
				if (alt<0.0 || alt>50000.0)
				{
					printf("Error. Wrong altitude (range 0 to 50,000 m).\n");
					continue;
				}
				break; // everything is fine here
			}
			if (!isInteractive) break;
			if (buf[0]=='r' || buf[0]=='R') continue;	// restart input


			double u,v,t,p,r;
			pMH370EnvironData->GetMeteo(u,v,t,p,r, ts,lon,lat,alt);

			printf("\n");
			printf("Extracted data:\n");
			printf("WindU: %5.2f (m/s)\n", u);
			printf("WindV: %5.2f (m/s)\n", v);
			printf("AirT: %5.2f (deg C)\n", t);
			printf("Air pressure: %6.0f (Pa)\n", p);
			printf("Relative humidity: %5.1f %%\n", r);
			printf("-----------------------------------\n");
			printf("\n");

		}
	}

	delete pMH370EnvironData;
	if (pTs) delete pTs;
	if (pLon) delete pLon;
	if (pLat) delete pLat;
	if (pAlt) delete pAlt;
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compute BTO
void ComputeBTO(char*inpfilename, char* outfilename, char*str_time, char*str_lon, char*str_lat, char*str_alt)
{
	int nOut = 0;
	double* pTs = NULL;
	double* pLon = NULL;
	double* pLat = NULL;
	double* pAlt = NULL;

	if (inpfilename)
	{
		if (strlen(inpfilename)>0)
		{
			if (!loadLonLatAlt(nOut, pTs, pLon, pLat, pAlt, inpfilename))
			{
				if (pTs) delete pTs;
				if (pLon) delete pLon;
				if (pLat) delete pLat;
				if (pAlt) delete pAlt;
				return;
			}
		}
	}

	if (inpfilename == NULL || strlen(inpfilename)==0)
	{
		char* pp;
		double ts = strtod(str_time, &pp);
		if (&str_time[0] == pp)
		{
			printf("Error. Wrong time.\n");
			return;
		}
		if (*pp=='-')
		{
			// time seems to be given as a string yyyy-....
			if (!read_time_formatted(ts, str_time))
			{
				printf("Error. Unrecognized time argument.\n");
				printf("Must be either s since 2014-03-07 00:00:00 UTC, or date yyyy-mm-dd hh:mm:ss.sss\n");
				return;
			}
		}


		double lon = strtod(str_lon, &pp);
		if (&str_lon[0] == pp)
		{
			printf("Error. Wrong longitude.\n");
			return;
		}

		double lat = strtod(str_lat, &pp);
		if (&str_lat[0] == pp)
		{
			printf("Error. Wrong latitude.\n");
			return;
		}

		double alt = strtod(str_alt, &pp);
		if (&str_alt[0] == pp)
		{
			printf("Error. Wrong altitude.\n");
			return;
		}

		pTs = new double[2];
		pLon = new double[2];
		pLat = new double[2];
		pAlt = new double[2];
		nOut = 1;

		pTs[0] = ts;
		pLon[0] = lon;
		pLat[0] = lat;
		pAlt[0] = alt;
	}

	// Here everything is fine
	CMH370InmarsatData* pInmarsatData = new CMH370InmarsatData();
	char timebuf[256];

	// Try to open output file, if specified
	FILE* fid = NULL;

	if (outfilename)
	{
		if (strlen(outfilename)>0)
		{
			fid = fopen(outfilename, "wt");
		}
	}

	if (fid)
	{
		fprintf(fid, "Time(UTC), BTO-R1200(microsec), BTO-R600(microsec), BTO-T1200(microsec)\n");
	}
	else
	{
		printf("Time(UTC), BTO-R1200(microsec), BTO-R600(microsec), BTO-T1200(microsec)\n");
	}

	for (int i=0; i<nOut; i++)
	{
		double bto_rx1200 = pInmarsatData->CalcBTO(pTs[i], pLon[i], pLat[i], pAlt[i], CH_TYPE_R_RX1200);
		double bto_rx600 = pInmarsatData->CalcBTO(pTs[i], pLon[i], pLat[i], pAlt[i], CH_TYPE_R_RX600);
		double bto_t1200 = pInmarsatData->CalcBTO(pTs[i], pLon[i], pLat[i], pAlt[i], CH_TYPE_T_RX1200);

		print_time_formatted(timebuf, pTs[i]);

		if (fid)
		{
			fprintf(fid, "%s  %11.1f %11.1f %11.1f\n", timebuf, bto_rx1200, bto_rx600, bto_t1200);
		}
		else
		{
			printf("%s  %11.1f %11.1f %11.1f\n", timebuf, bto_rx1200, bto_rx600, bto_t1200);
		}
	}

	if (fid) fclose(fid); // close output file, if it was open

	delete pInmarsatData;
	if (pTs) delete pTs;
	if (pLon) delete pLon;
	if (pLat) delete pLat;
	if (pAlt) delete pAlt;
}





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compute BFO
void ComputeBFO(char*inpfilename, char* outfilename, char*str_time, char*str_lon, char*str_lat, char*str_alt, char*str_u, char*str_v, char*str_w, char* AES_Earth_model, char* AES_Altitude)
{
//	CMH370InmarsatData* pInmarsatData = new CMH370InmarsatData();
//	if (isVerbose) pInmarsatData -> SetVerbose(stdout);
//	pInmarsatData -> SetBFOBiases(153.14, 152.79+0.84, 153.44+0.84+0.17); // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< experimental
//	if (AES_Earth_Model != AES_EARTH_MODEL_UDF) pInmarsatData -> SetAESEarthModel(AES_Earth_Model);
//	if (AES_Doppler_Method != DOPPLER_METHOD_UDF) pInmarsatData -> SetAESDopplerMethod(AES_Doppler_Method);
//	if (AES_Earth_Radius>0.0) pInmarsatData -> SetAESEarthRadius(AES_Earth_Radius);
//	if (AES_GeoSat_Alt>0.0) pInmarsatData -> SetAESGeoSatAlt(AES_GeoSat_Alt);

	int nOut = 0;
	double* pTs = NULL;
	double* pLon = NULL;
	double* pLat = NULL;
	double* pAlt = NULL;
	double* pU = NULL;
	double* pV = NULL;
	double* pW = NULL;

	if (inpfilename)
	{
		if (strlen(inpfilename)>0)
		{
			if (!loadLonLatAltVel(nOut, pTs, pLon, pLat, pAlt, pU, pV, pW, inpfilename))
			{
				if (pTs) delete pTs;
				if (pLon) delete pLon;
				if (pLat) delete pLat;
				if (pAlt) delete pAlt;
				if (pU) delete pU;
				if (pV) delete pV;
				if (pW) delete pW;
				return;
			}
		}
	}

	if (inpfilename == NULL || strlen(inpfilename)==0)
	{
		char* pp;
		double ts = strtod(str_time, &pp);
		if (&str_time[0] == pp)
		{
			printf("Error. Wrong time.\n");
			return;
		}
		if (*pp=='-')
		{
			// time seems to be given as a string yyyy-....
			if (!read_time_formatted(ts, str_time))
			{
				printf("Error. Unrecognized time argument.\n");
				printf("Must be either s since 2014-03-07 00:00:00 UTC, or date yyyy-mm-dd hh:mm:ss.sss\n");
				return;
			}
		}


		double lon = strtod(str_lon, &pp);
		if (&str_lon[0] == pp)
		{
			printf("Error. Wrong longitude.\n");
			return;
		}

		double lat = strtod(str_lat, &pp);
		if (&str_lat[0] == pp)
		{
			printf("Error. Wrong latitude.\n");
			return;
		}

		double alt = strtod(str_alt, &pp);
		if (&str_alt[0] == pp)
		{
			printf("Error. Wrong altitude.\n");
			return;
		}

		double vel_u = strtod(str_u, &pp);
		if (&str_u[0] == pp)
		{
			printf("Error. Wrong u-velocity component.\n");
			return;
		}

		double vel_v = strtod(str_v, &pp);
		if (&str_v[0] == pp)
		{
			printf("Error. Wrong v-velocity component.\n");
			return;
		}

		double vel_w = strtod(str_w, &pp);
		if (&str_w[0] == pp)
		{
			printf("Error. Wrong w-velocity component.\n");
			return;
		}

		pTs = new double[2];
		pLon = new double[2];
		pLat = new double[2];
		pAlt = new double[2];
		pU = new double[2];
		pV = new double[2];
		pW = new double[2];

		nOut = 1;

		pTs[0] = ts;
		pLon[0] = lon;
		pLat[0] = lat;
		pAlt[0] = alt;
		pU[0] = vel_u;
		pV[0] = vel_v;
		pW[0] = vel_w;
	}

	// Here everything is fine
	CMH370InmarsatData* pInmarsatData = new CMH370InmarsatData();
	pInmarsatData -> SetVerbose(stdout);
//	pInmarsatData -> SetBFOBiases(153.14, 152.79+0.84, 153.44+0.84+0.17); // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< experimental
//	pInmarsatData -> SetAESEarthModel(AES_EARTH_MODEL_WGS);
//	pInmarsatData -> SetAESGeoSatAlt(35785863.0);
	if(AES_Earth_model)
	{
		if (strcmp(AES_Earth_model,"-AES_WGS") == 0) pInmarsatData -> SetAESEarthModel(AES_EARTH_MODEL_WGS);
		if (strcmp(AES_Earth_model,"-AES_SPH") == 0) pInmarsatData -> SetAESEarthModel(AES_EARTH_MODEL_SPH);
	}

	if(AES_Altitude)
	{
		char* p = strstr(AES_Altitude, "-AES_ALT");
		if (p)
		{
			p+=strlen("-AES_ALT");
			int len = strlen(AES_Altitude);
			while(*p==' ' && (p < &AES_Altitude[len-1])) p++;
			char* pp;
			double AES_alt = strtod(p, &pp);

			if (p!=pp && AES_alt>0) pInmarsatData -> SetAESGeoSatAlt(AES_alt);
		}
	}

	// buffer for time as a text string
	char timebuf[256];

	// Try to open output file, if specified
	FILE* fid = NULL;

	if (outfilename)
	{
		if (strlen(outfilename)>0)
		{
			fid = fopen(outfilename, "wt");
		}
	}

	if (fid)
	{
		fprintf(fid, "Time(UTC), BFO-R(Hz), BFO-T(Hz), BTO-C(Hz)\n");
	}
	else
	{
		printf("Time(UTC), BFO-R(Hz), BFO-T(Hz), BTO-C(Hz)\n");
	}

	for (int i=0; i<nOut; i++)
	{
		double bfo_r = pInmarsatData->CalcBFO(pTs[i], pLon[i], pLat[i], pAlt[i], pU[i], pV[i], pW[i], CH_TYPE_R_RX1200);
		double bfo_t = pInmarsatData->CalcBFO(pTs[i], pLon[i], pLat[i], pAlt[i], pU[i], pV[i], pW[i], CH_TYPE_T_RX1200);
		double bfo_c = pInmarsatData->CalcBFO(pTs[i], pLon[i], pLat[i], pAlt[i], pU[i], pV[i], pW[i], CH_TYPE_C_RX1200);

		print_time_formatted(timebuf, pTs[i]);

		if (fid)
		{
			fprintf(fid, "%s  %11.1f %11.1f %11.1f\n", timebuf, bfo_r, bfo_t, bfo_c);
		}
		else
		{
			printf("%s  %11.1f %11.1f %11.1f\n", timebuf, bfo_r, bfo_t, bfo_c);
		}
	}

	if (fid) fclose(fid); // close output file, if it was open

	delete pInmarsatData;
	if (pTs) delete pTs;
	if (pLon) delete pLon;
	if (pLat) delete pLat;
	if (pAlt) delete pAlt;
	if (pU) delete pU;
	if (pV) delete pV;
	if (pW) delete pW;
}




////////////////////////////////////////////////////////////////////////////////////////
//
// Heading conversion procedures: True\Magnetic\Gyroscopic
//
////////////////////////////////////////////////////////////////////////////////////////

double True2MagneticHDG(double hdg, double lon, double lat, CMH370EnvironData* pEnvironData)
{
	// Input: hdg - true heading
	//		  lon, lat - longitude and latitude
	// Output: magnetic heading

	double h = 0.0;

	if (pEnvironData) h = pEnvironData->GetMagneticDeclination(lon,lat); // magnetic declination
	h = hdg - h;

	if (h>360.0) h -= 360.0;
	if (h<0.0) h += 360.0;

	return h;
}


double Magnetic2TrueHDG(double hdg, double lon, double lat, CMH370EnvironData* pEnvironData)
{
	// Input: hdg - magnetic heading
	//		  lon, lat - longitude and latitude
	// Output: true heading

	double h = 0.0;

	if (pEnvironData) h = pEnvironData->GetMagneticDeclination(lon,lat); // magnetic declination
	h = hdg + h;

	if (h>360.0) h -= 360.0;
	if (h<0.0) h += 360.0;

	return h;
}


double True2GyroHDG(double hdg, double lon, double lat, double t, double  gyro_heading_param1, double gyro_heading_param2)
{
	// Input: hdg - true heading
	//		  lon, lat - longitude and latitude
	//		  t - time (s) since 2014-03-07 00:00:00
	//		  gyro_heading_param1, gyro_heading_param2 - gyroscopic heading parameters (phase and magnitude)
	// Output: gyroscopic heading

	// Define gyroscopic North vector in the non-spinning reference frame (note ECEF system is spinning). 
	// Consider the two parameters being longitude and latitude in the absolute reference frame. ECEF coincides with this frame on 2014-03-07 00:00:00.
	// Note: this vector remains constant over the integration period
	// Presumably it may match the time when SAARU's north was last set, and is may not necessarily be when the SAARU was started up in the KLIA

		double p0x, p0y, p0z;
		uv2uvw(p0x, p0y, p0z, 0.0, 1.0, gyro_heading_param1, gyro_heading_param2);

		// Now calculate {u,v} direction vectors of the local tangential plane in the non-rotating reference frame
		// use uv2uvw to compute W->E and S->N direction vectors, and then projections of the "gyroscopic north" vector in ECEF (assuming that SAARU thinks it is the true north)
		// on those to determine the deviation of the "gyroscopic reference heading" from the true north
		double n_we_x, n_we_y, n_we_z;
		uv2uvw(n_we_x, n_we_y, n_we_z, 1.0, 0.0, lon-360.0/86400.0*t, lat);	// unit vector west to east in the non-spinning reference frame

		double n_sn_x, n_sn_y, n_sn_z;
		uv2uvw(n_sn_x, n_sn_y, n_sn_z, 0.0, 1.0, lon-360.0/86400.0*t, lat);	// unit vector south to north in the non-spinning reference frame

		// Now we need to present gyroscopic heading vector as the sum of vectors of the local reference frame, which moves with the airplane
		// To do that calculate projections of the gyroscopic N vector on the axes of the tangential plane (vertical component is ignored)
		double n_lon_gyro = p0x*n_we_x + p0y*n_we_y + p0z*n_we_z;
		double n_lat_gyro = p0x*n_sn_x + p0y*n_sn_y + p0z*n_sn_z;

		// When the plane moves along a meridian, this vector becomes not parallel to the tangential plane, so its length becomes < 1.0. 
		// Scale it to have a unit vector in the local system.
		double s = sqrt(n_lon_gyro*n_lon_gyro + n_lat_gyro*n_lat_gyro);
		n_lon_gyro = n_lon_gyro/s;
		n_lat_gyro = n_lat_gyro/s;

		// convert these projections to heading (fron N, CW)
		double h = acos(n_lon_gyro)*180.0/3.1415926535897932384626433832795;
		if (n_lat_gyro<0.0) h = 360.0-h;
		h = 90.0 - h;

		// heading in gyroscopic system
		h = hdg - h;
		while (h>360.0) h -= 360.0;
		while (h<0.0) h += 360.0;
		return h;
}



double Gyro2TrueHDG(double hdg, double lon, double lat, double t, double gyro_heading_param1, double gyro_heading_param2)
{
	// Input: hdg - gyroscopic heading
	//		  lon, lat - longitude and latitude
	//		  t - time (s) since 2014-03-07 00:00:00
	//		  gyro_heading_param1, gyro_heading_param2 - gyroscopic heading parameters (phase and magnitude)
	// Output: true heading

	// Define gyroscopic North vector in the non-spinning reference frame (note ECEF system is spinning). 
	// Consider the two parameters being longitude and latitude in the absolute reference frame. ECEF coincides with this frame on 2014-03-07 00:00:00.
	// Note: this vector remains constant over the integration period
	// Presumably it may match the time when SAARU's north was last set, and is may not necessarily be when the SAARU was started up in the KLIA

		double p0x, p0y, p0z;
		uv2uvw(p0x, p0y, p0z, 0.0, 1.0, gyro_heading_param1, gyro_heading_param2);

		// Now calculate {u,v} direction vectors of the local tangential plane in the non-rotating reference frame
		// use uv2uvw to compute W->E and S->N direction vectors, and then projections of the "gyroscopic north" vector in ECEF (assuming that SAARU thinks it is the true north)
		// on those to determine the deviation of the "gyroscopic reference heading" from the true north
		double n_we_x, n_we_y, n_we_z;
		uv2uvw(n_we_x, n_we_y, n_we_z, 1.0, 0.0, lon-360.0/86400.0*t, lat);	// unit vector west to east in the non-spinning reference frame

		double n_sn_x, n_sn_y, n_sn_z;
		uv2uvw(n_sn_x, n_sn_y, n_sn_z, 0.0, 1.0, lon-360.0/86400.0*t, lat);	// unit vector south to north in the non-spinning reference frame

		// Now we need to present gyroscopic heading vector as the sum of vectors of the local reference frame, which moves with the airplane
		// To do that calculate projections of the gyroscopic N vector on the axes of the tangential plane (vertical component is ignored)
		double n_lon_gyro = p0x*n_we_x + p0y*n_we_y + p0z*n_we_z;
		double n_lat_gyro = p0x*n_sn_x + p0y*n_sn_y + p0z*n_sn_z;

		// When the plane moves along a meridian, this vector becomes not parallel to the tangential plane, so its length becomes < 1.0. 
		// Scale it to have a unit vector in the local system.
		double s = sqrt(n_lon_gyro*n_lon_gyro + n_lat_gyro*n_lat_gyro);
		n_lon_gyro = n_lon_gyro/s;
		n_lat_gyro = n_lat_gyro/s;

		// convert these projections to heading (fron N, CW)
		double h = acos(n_lon_gyro)*180.0/3.1415926535897932384626433832795;
		if (n_lat_gyro<0.0) h = 360.0-h;
		h = 90.0 - h;

		// true heading
		h = hdg + h;
		while (h>360.0) h -= 360.0;
		while (h<0.0) h += 360.0;

		return h;
}


double Magnetic2GyroHDG(double hdg, double lon, double lat, double t, double gyro_heading_param1, double gyro_heading_param2, CMH370EnvironData* pEnvironData)
{
	double h = Magnetic2TrueHDG(hdg, lon, lat, pEnvironData);
	h = True2GyroHDG(h,lon,lat,t,gyro_heading_param1,gyro_heading_param2);
	return h;
}


double Gyro2MagneticHDG(double hdg, double lon, double lat, double t, double gyro_heading_param1, double gyro_heading_param2, CMH370EnvironData* pEnvironData)
{
	double h = Gyro2TrueHDG(hdg,lon,lat,t,gyro_heading_param1,gyro_heading_param2);
	h = True2MagneticHDG(h, lon, lat, pEnvironData);
	return h;
}

