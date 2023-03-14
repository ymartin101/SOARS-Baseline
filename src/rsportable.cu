/// rsportable.cpp - All code which is non-standard C++ must go in here, with reasons why
/// Marc Brooker, 25 April 2006
/// Edited by Craig Tong
/// Edited by Yaaseen Martin, 27 August 2019

#include <string.h>
#include <math.h>
#include <stdio.h>
#include <sstream>
#include "rsportable.cuh"
#include "rsdebug.cuh"
#include <boost/thread.hpp>

/// Compare two strings, ignoring case
// There isn no suitable one in the standard library of C or C++
int rsPortable::stricmp(const char *one, const char *two)
{
    return strcasecmp(one, two); //strcasecmp is a GNU extension
}

/// Compute the first order Bessel function of the first kind
rsFloat rsPortable::BesselJ1(rsFloat x)
{
    return j1(x); // J1 is non standard, but found on many platforms
}

/// Round off a floating point number
rsFloat rsPortable::rsRound(rsFloat x)
{
    return round(x);
}

/// Detect the number of CPUs in the machine
int rsPortable::CountProcessors()
{
    int iNHardwareThreads = boost::thread::hardware_concurrency();

    if(!iNHardwareThreads)
    {
        rsDebug::printf(rsDebug::RS_IMPORTANT, "NOTE: Unable to get CPU count, assuming 1.\n");
        return 1;
    }
    else
        return iNHardwareThreads;
}
