/*     ----------------------------------------------------------------
*
*                                asttime.cu
*
*   This file contains fundamental Astrodynamic procedures and functions
*   relating to the time functions. These routines are discussed in Ch 3
*   and Ch 5.
*
*                          Companion code for
*             Fundamentals of Astrodynamics and Applications
*                                  2013
*                            by David Vallado
*
*       (w) 719-573-2600, email dvallado@agi.com, davallado@gmail.com
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
*              21 jul 05  david vallado
*                           2nd printing baseline
*              14 may 01  david vallado
*                           2nd edition baseline
*              23 nov 87  david vallado
*                           original baseline
*       ----------------------------------------------------------------      */

// stdafx.h must be the first in cpp files
//#include "stdafx.h"

#include "astTime.cuh"
#include <stdio.h>      /* printf */
#include <math.h>       /* fmod */

#define pi 3.14159265358979323846
#define twopi 2.0 * pi

namespace astTime
{

	/* -----------------------------------------------------------------------------
	*
	*                           procedure jday
	*
	*  this procedure finds the julian date given the year, month, day, and time.
	*    the julian date is defined by each elapsed day since noon, jan 1, 4713 bc.
	*
	*  algorithm     : calculate the answer in one step for efficiency
	*
	*  author        : david vallado                  719-573-2600    1 mar 2001
	*
	*  inputs          description                    range / units
	*    year        - year                           1900 .. 2100
	*    mon         - month                          1 .. 12
	*    day         - day                            1 .. 28,29,30,31
	*    hr          - universal time hour            0 .. 23
	*    min         - universal time min             0 .. 59
	*    sec         - universal time sec             0.0 .. 59.999
	*
	*  outputs       :
	*    jd          - julian date (days only)           days from 4713 bc
	*    jdFrac      - julian date (fraction of a day)   days from 0 hr of the day
	*
	*  locals        :
	*    none.
	*
	*  coupling      :
	*    none.
	*
	*  references    :
	*    vallado       2013, 183, alg 14, ex 3-4
	*
	* --------------------------------------------------------------------------- */

	__device__ void jday
		(
		int year, int mon, int day, int hr, int minute, double sec,
		double& jd, double& jdFrac
		)
	{
		jd = 367.0 * year -
			floor((7 * (year + floor((mon + 9) / 12.0))) * 0.25) +
			floor(275 * mon / 9.0) +
			day + 1721013.5;  // use - 678987.0 to go to mjd directly
		jdFrac = (sec + minute * 60.0 + hr * 3600.0) / 86400.0;

		// check that the day and fractional day are correct
		if (fabs(jdFrac) >= 1.0)
		{
			double dtt = floor(jdFrac);
			jd = jd + dtt;
			jdFrac = jdFrac - dtt;
		}

		// - 0.5*sgn(100.0*year + mon - 190002.5) + 0.5;
	}

	/* -----------------------------------------------------------------------------
	*
	*                           procedure days2mdhms
	*
	*  this procedure converts the day of the year, days, to the equivalent month
	*    day, hour, minute and second.
	*
	*  algorithm     : set up array for the number of days per month
	*                  find leap year - use 1900 because 2000 is a leap year
	*                  loop through a temp value while the value is < the days
	*                  perform int conversions to the correct day and month
	*                  convert remainder into h m s using type conversions
	*
	*  author        : david vallado                  719-573-2600    1 mar 2001
	*
	*  inputs          description                    range / units
	*    year        - year                           1900 .. 2100
	*    days        - julian day of the year         1.0  .. 366.0
	*
	*  outputs       :
	*    mon         - month                          1 .. 12
	*    day         - day                            1 .. 28,29,30,31
	*    hr          - hour                           0 .. 23
	*    min         - minute                         0 .. 59
	*    sec         - second                         0.0 .. 59.999
	*
	*  locals        :
	*    dayofyr     - day of year
	*    temp        - temporary extended values
	*    inttemp     - temporary int value
	*    i           - index
	*    lmonth[12]  - int array containing the number of days per month
	*
	*  coupling      :
	*    none.
	* --------------------------------------------------------------------------- */

	__host__ __device__ void    days2mdhms
		(
		int year, double days,
		int& mon, int& day, int& hr, int& minute, double& sec
		)
	{
		int i, inttemp, dayofyr;
		double    temp;
		// start indicies at 1
		int lmonth[] = { 0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

		dayofyr = (int)floor(days);
		/* ----------------- find month and day of month ---------------- */
		if ((year % 4) == 0)
			lmonth[2] = 29;

		i = 1;
		inttemp = 0;
		while ((dayofyr > inttemp + lmonth[i]) && (i < 12))
		{
			inttemp = inttemp + lmonth[i];
			i++;
		}
		mon = i;
		day = dayofyr - inttemp;

		/* ----------------- find hours minutes and seconds ------------- */
		temp = (days - dayofyr) * 24.0;
		hr = floor(temp);
		temp = (temp - hr) * 60.0;
		minute = floor(temp);
		sec = (temp - minute) * 60.0;
	}

	/* -----------------------------------------------------------------------------
	*
	*                           procedure invjday
	*
	*  this procedure finds the year, month, day, hour, minute and second
	*  given the julian date. tu can be ut1, tdt, tdb, etc.
	*
	*  algorithm     : set up starting values
	*                  find leap year - use 1900 because 2000 is a leap year
	*                  find the elapsed days through the year in a loop
	*                  call routine to find each individual value
	*
	*  author        : david vallado                  719-573-2600    1 mar 2001
	*
	*  inputs          description                    range / units
	*    jd          - julian date (days only)           days from 4713 bc
	*    jdFrac      - julian date (fraction of a day)   days from 0 hr of the day
	*
	*  outputs       :
	*    year        - year                           1900 .. 2100
	*    mon         - month                          1 .. 12
	*    day         - day                            1 .. 28,29,30,31
	*    hr          - hour                           0 .. 23
	*    min         - minute                         0 .. 59
	*    sec         - second                         0.0 .. 59.999
	*
	*  locals        :
	*    days        - day of year plus fractional
	*                  portion of a day               days
	*    tu          - julian centuries from 0 h
	*                  jan 0, 1900
	*    temp        - temporary double values
	*    leapyrs     - number of leap years from 1900
	*
	*  coupling      :
	*    days2mdhms  - finds month, day, hour, minute and second given days and year
	*
	*  references    :
	*    vallado       2013, 203, alg 22, ex 3-13
	* --------------------------------------------------------------------------- */

	void    invjday
		(
		double jd, double jdfrac,
		int& year, int& mon, int& day,
		int& hr, int& minute, double& sec
		)
	{
		int leapyrs;
		double dt, days, tu, temp;

		// check jdfrac for multiple days
		if (fabs(jdfrac) >= 1.0)
		{
			jd = jd + floor(jdfrac);
			jdfrac = jdfrac - floor(jdfrac);
		}

		// check for fraction of a day included in the jd
		dt = jd - floor(jd) - 0.5;
		if (fabs(dt) > 0.00000001)
		{
			jd = jd - dt;
			jdfrac = jdfrac + dt;
		}

		/* --------------- find year and days of the year --------------- */
		temp = jd - 2415019.5;
		tu = temp / 365.25;
		year = 1900 + (int)floor(tu);
		leapyrs = (int)floor((year - 1901) * 0.25);

		days = floor(temp - ((year - 1900) * 365.0 + leapyrs));

		/* ------------ check for case of beginning of a year ----------- */
		if (days + jdfrac < 1.0)
		{
			year = year - 1;
			leapyrs = (int)floor((year - 1901) * 0.25);
			days = floor(temp - ((year - 1900) * 365.0 + leapyrs));
		}

		/* ----------------- find remaing data  ------------------------- */
		// now add the daily time in to preserve accuracy
		days2mdhms(year, days + jdfrac, mon, day, hr, minute, sec);
	}

	/* -----------------------------------------------------------------------------
	*
	*                           function gstime
	*
	*  this function finds the greenwich sidereal time (iau-82).
	*
	*  author        : david vallado                  719-573-2600    1 mar 2001
	*
	*  inputs          description                    range / units
	*    jdut1       - julian date in ut1             days from 4713 bc
	*
	*  outputs       :
	*    gstime      - greenwich sidereal time        0 to 2pi rad
	*
	*  locals        :
	*    temp        - temporary variable for doubles   rad
	*    tut1        - julian centuries from the
	*                  jan 1, 2000 12 h epoch (ut1)
	*
	*  coupling      :
	*    none
	*
	*  references    :
	*    vallado       2013, 187, eq 3-45
	* --------------------------------------------------------------------------- */

	__host__ __device__ double  gstime
		(
		double jdut1
		)
	{
		double temp, tut1;

		tut1 = (jdut1 - 2451545.0) / 36525.0;
		temp = -6.2e-6* tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1 +
			(876600.0 * 3600 + 8640184.812866) * tut1 + 67310.54841;  // sec
		temp = fmod((temp * (pi / 180.0) / 240.0), 2*pi); //360/86400 = 1/240, to deg, to rad

		// ------------------------ check quadrants ---------------------
		if (temp < 0.0)
			temp = temp + 2*pi;

		return temp;
	}

	/* -----------------------------------------------------------------------------
	*
	*                           procedure hms_sec
	*
	*  this procedure converts hours, minutes and seconds to seconds from the
	*  beginning of the day.
	*
	*  author        : david vallado                  719-573-2600   25 jun 2002
	*
	*  revisions
	*                -
	*
	*  inputs          description                    range / units
	*    utsec       - seconds                        0.0 .. 86400.0
	*
	*  outputs       :
	*    hr          - hours                          0 .. 24
	*    min         - minutes                        0 .. 59
	*    sec         - seconds                        0.0 .. 59.99
	*
	*  locals        :
	*    temp        - temporary variable
	*
	*  coupling      :
	*    none.
	* --------------------------------------------------------------------------- */

	__device__ void    hms_sec
		(
		int& hr, int& min, double& sec, edirection direct, double& utsec
		)
	{
		double temp;

		// ------------------------  implementation   ------------------
		if (direct == eTo)
			utsec = hr * 3600.0 + min * 60.0 + sec;
		else
		{
			temp = utsec / 3600.0;
			hr = (int)floor(temp);
			min = (int)floor((temp - hr)* 60.0);
			sec = (temp - hr - min / 60.0) * 3600.0;
		}
	}

}  // namespace
