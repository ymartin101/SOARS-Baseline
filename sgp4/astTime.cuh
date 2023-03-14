#ifndef _astTime_h_
#define _astTime_h_
/* --------------------------------------------------------------------
*
*                                astTime.cuh
*
*    this file contains miscallaneous time functions.
*
*    current :
*              30 sep 15  david vallado
*                           fix jd, jdfrac
*    changes :
*               3 nov 14  david vallado
*                           update to msvs2013 c++
*               4 may 09  david vallado
*                           misc updates
*              13 feb 08  david vallado
*                           add getmon
*              15 mar 07  david vallado
*                           3rd edition baseline
*              20 jul 05  david vallado
*                           2nd printing baseline
*              14 may 01  david vallado
*                           2nd edition baseline
*              23 nov 87  david vallado
*                           original baseline
  ---------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#pragma once

#define pi 3.14159265358979323846
#define twopi 2.0 * pi

// global object definitions
typedef enum
    {
       eTo, eFrom
    } edirection;


namespace astTime
{

	__device__ void jday
		(
		int year, int mon, int day, int hr, int minute, double sec,
		double& jd, double& jdFrac
		);

	__host__ __device__ void    days2mdhms
		(
		int year, double days,
		int& mon, int& day, int& hr, int& minute, double& sec
		);

	void    invjday
		(
		double jd, double jdFrac,
		int& year, int& mon, int& day,
		int& hr, int& minute, double& sec
		);

	__host__ __device__ double  gstime
		(
		double jdut1
		);

	__device__ void    hms_sec
		(
		int& hr, int& min, double& sec, edirection direct, double& utsec
		);

};  // namespace

#endif
