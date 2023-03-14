/// xmlimport.cpp - Import simulator world and parameters from a SOARSXML file
/// xmlimport_manual: Used to manually propagate a target
/// Marc Brooker, 26 April 2006
/// Edited by Yaaseen Martin, 17 September 2019

#include <string>
#include <sstream>
#include <stdio.h>
#include <stdexcept>
#include <cmath>
#include "xmlimport.cuh"
#include "rsdebug.cuh"
#include "rsworld.cuh"
#include "rsplatform.cuh"
#include "rstarget.cuh"
#include "rsradar.cuh"
#include "rsportable.cuh"
#include "rsparameters.cuh"
#include "rsantenna.cuh"
#include "rstiming.cuh"
#include "rspython.cuh"
#include "rsmultipath.cuh"
#include "SGP4.cuh"
#include "astTime.cuh"
#include "astMath.cuh"
#include "pugixml.hpp"
#include <sys/time.h>

//struct timeval t3, t4;

// Check for CUDA errors
#define cudaCheckErrors(msg) \
    do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
            fprintf(stderr, "Fatal error: %s (%s at %s:%d)\n", \
                msg, cudaGetErrorString(__err), \
                __FILE__, __LINE__); \
            fprintf(stderr, "*** FAILED - ABORTING\n"); \
            exit(1); \
        } \
    } while (0)

// For processing XML like this:
// <tree>
//   <leaf1>Green</leaf1>
//   <leaf2>Blue</leaf2>
// </tree>
// Pass a handle to tree, and the string "leaf1" to get "Green"

using namespace rs;
using std::string;
using std::vector;

/// XML Parsing Utility Functions

// Exception for reporting an XML parsing error
class XmlImportException: public std::runtime_error {
    public:
        XmlImportException(string error):
            std::runtime_error("ERROR: Error while parsing XML file: "+error)
        {
        }
};

// Function which takes a pugi::xml_node and returns the text contained in its children.
string GetChildText(pugi::xml_node &parent, std::string childname)
{
    string data = parent.child_value(childname.c_str());

    // Return the text
    return data;
}

// Gets child text as a rsFloat; see GetChildText for usage description
rsFloat GetChildRsFloat(pugi::xml_node &parent, string childname)
{
    string data = parent.child_value(childname.c_str());
    if (data == "") // If data is not found in XML
        return -1.0;
    else{
        // Parse the first rsFloat from the data
        rsFloat result;
        std::istringstream iss(data);
        iss >> result;
        return result;
    }
}

// Return the string associated with an attribute or throw an exception on failure
string GetAttributeString(pugi::xml_node &handle, string name, string error, bool optional=false)
{
    string text = handle.attribute(name.c_str()).value();
    if (text != "")
        return text;
    else {
        if (!optional)
            throw XmlImportException(error);
        else
            return text;
    }
}

// Return the bool associated with an attribute
bool GetAttributeBool(pugi::xml_node &handle, string name, string error, bool def, bool optional=true)
{
    string str = GetAttributeString(handle, name, error, optional);
    if (str == "")
        return def;
    return ((str == "true") || (str == "yes"));
}

/// SGP4 functions
// convtime finds the time parameters and Julian century values
// Based on code written by David Vallado
__device__ void convtime(int year, int mon, int day, int hr, int min, rsFloat sec, int timezone, rsFloat dut1, int dat,
    rsFloat& ut1, rsFloat& tut1, rsFloat& jdut1, rsFloat& jdut1Frac, rsFloat& utc, rsFloat& tai,
    rsFloat& tt, rsFloat& ttt, rsFloat& jdtt, rsFloat& jdttFrac, rsFloat& tcg, rsFloat& tdb,
    rsFloat& ttdb, rsFloat& jdtdb, rsFloat& jdtdbFrac, rsFloat& tcb)
{
    rsFloat jd, jdFrac, sectemp;
    int localhr, hrtemp, mintemp;

    // Implementation
    astTime::jday(year, mon, day, 0, 0, 0.0, jd, jdFrac);
    localhr = timezone + hr;
    astTime::hms_sec(localhr, min, sec, eTo, utc);
    ut1 = utc + dut1;
    astTime::hms_sec(hrtemp, mintemp, sectemp, eFrom, ut1);
    astTime::jday(year, mon, day, hrtemp, mintemp, sectemp, jdut1, jdut1Frac);
    jdut1 = jdut1 + jdut1Frac;
    tut1 = (jdut1 - 2451545.0) / 36525.0;
    tai = utc + dat;
    tt = tai + 32.184;   // sec
    astTime::hms_sec(hrtemp, mintemp, sectemp, eFrom, tt);
    astTime::jday(year, mon, day, hrtemp, mintemp, sectemp, jdtt, jdttFrac);
    ttt = (jdtt + jdttFrac - 2451545.0) / 36525.0;
}

// polarm calulates the transformation matrix that accounts for polar motion
// Based on code written by David Vallado
__device__ void polarm(rsFloat xp, rsFloat yp, rsFloat ttt, eOpt opt, rsFloat pm[3][3])
{
    rsFloat convrt, cosxp, cosyp, sinxp, sinyp, sp, cossp, sinsp;

    convrt = 1.0;

    cosxp = cos(xp * convrt);
    sinxp = sin(xp * convrt);
    cosyp = cos(yp * convrt);
    sinyp = sin(yp * convrt);

    if ((opt == e80) | (opt == e96))
    {
        pm[0][0] = cosxp;
        pm[0][1] = 0.0;
        pm[0][2] = -sinxp;
        pm[1][0] = sinxp  *  sinyp;
        pm[1][1] = cosyp;
        pm[1][2] = cosxp  *  sinyp;
        pm[2][0] = sinxp  *  cosyp;
        pm[2][1] = -sinyp;
        pm[2][2] = cosxp  *  cosyp;
    }
    else
    {
        // Approximate sp value in rad
        sp = -47.0e-6 * ttt * pi / (180.0 * 3600.0);
        cossp = cos(sp);
        sinsp = sin(sp);

        // Form the matrix
        pm[0][0] = cosxp * cossp;
        pm[0][1] = -cosyp * sinsp + sinyp * sinxp * cossp;
        pm[0][2] = -sinyp * sinsp - cosyp * sinxp * cossp;
        pm[1][0] = cosxp * sinsp;
        pm[1][1] = cosyp * cossp + sinyp * sinxp * sinsp;
        pm[1][2] = sinyp * cossp - cosyp * sinxp * sinsp;
        pm[2][0] = sinxp;
        pm[2][1] = -sinyp * cosxp;
        pm[2][2] = cosyp * cosxp;

    }
}

// teme_ecef transforms a vector from TEME to ECEF
// Based on code written by David Vallado
__device__ void teme_ecef(rsFloat rteme[3], rsFloat vteme[3], rsFloat ateme[3], edirection direct,
               rsFloat recef[3], rsFloat vecef[3], rsFloat aecef[3],
               rsFloat ttt, rsFloat jdut1, rsFloat lod, rsFloat xp, rsFloat yp, int eqeterms)
{
    rsFloat deg2rad, omega, gmstg, thetasa;
    rsFloat st[3][3];
    rsFloat pm[3][3];
    rsFloat pmp[3][3];
    rsFloat stp[3][3];
    rsFloat omegaearth[3], rpef[3], vpef[3], apef[3], omgxr[3], omgxomgxr[3],
        omgxv[3], tempvec1[3], tempvec[3], gmst;

    deg2rad = pi / 180.0;

    // Find omega from nutation theory
    omega = 125.04452222 + (-6962890.5390 *ttt + 7.455 *ttt*ttt + 0.008 *ttt*ttt*ttt) / 3600.0;
    omega = fmod(omega, 360.0) * deg2rad;

    // Find GMST
    gmst = astTime::gstime(jdut1);

    // TEME does not include the geometric terms here; after 1997, kinematic terms apply
    if ((jdut1 > 2450449.5) && (eqeterms > 0))
    {
        gmstg = gmst
            + 0.00264*pi / (3600.0 * 180.0)*sin(omega)
            + 0.000063*pi / (3600.0 * 180.0)*sin(2.0 *omega);
    }
    else
        gmstg = gmst;
    gmstg = fmod(gmstg, 2.0*pi);

    thetasa = 7.29211514670698e-05 * (1.0 - lod / 86400.0);
    omegaearth[0] = 0.0;
    omegaearth[1] = 0.0;
    omegaearth[2] = thetasa;

    st[0][0] = cos(gmstg);
    st[0][1] = -sin(gmstg);
    st[0][2] = 0.0;
    st[1][0] = sin(gmstg);
    st[1][1] = cos(gmstg);
    st[1][2] = 0.0;
    st[2][0] = 0.0;
    st[2][1] = 0.0;
    st[2][2] = 1.0;

    polarm(xp, yp, ttt, e80, pm);

    if (direct == eTo)
    {
        astMath::mattrans(pm, pmp, 3, 3);
        astMath::mattrans(st, stp, 3, 3);

        astMath::matvecmult(stp, rteme, rpef);
        astMath::matvecmult(pmp, rpef, recef);

        astMath::cross(omegaearth, rpef, omgxr);
        astMath::matvecmult(stp, vteme, tempvec1);
        astMath::addvec(1.0, tempvec1, -1.0, omgxr, vpef);
        astMath::matvecmult(pmp, vpef, vecef);

        astMath::addvec(1.0, tempvec1, -1.0, omgxr, vpef);
        astMath::cross(omegaearth, vpef, omgxv);
        astMath::cross(omegaearth, omgxr, omgxomgxr);
        astMath::matvecmult(stp, ateme, tempvec1);
        astMath::addvec(1.0, tempvec1, -1.0, omgxomgxr, tempvec);
        astMath::addvec(1.0, tempvec, -2.0, omgxv, apef);
        astMath::matvecmult(pmp, apef, aecef);
    }
    else
    {
        astMath::matvecmult(pm, recef, rpef);
        astMath::matvecmult(st, rpef, rteme);

        astMath::matvecmult(pm, vecef, vpef);
        astMath::cross(omegaearth, rpef, omgxr);
        astMath::addvec(1.0, vpef, 1.0, omgxr, tempvec1);
        astMath::matvecmult(st, tempvec1, vteme);

        astMath::matvecmult(pm, aecef, apef);
        astMath::cross(omegaearth, omgxr, omgxomgxr);
        astMath::cross(omegaearth, vpef, omgxv);
        astMath::addvec(1.0, apef, 1.0, omgxomgxr, tempvec);
        astMath::addvec(1.0, tempvec, 2.0, omgxv, tempvec1);
        astMath::matvecmult(st, tempvec1, ateme);
    }
}

// Kernel for SGP4 target propagation
// Based on code written by David Vallado
__global__ void SGP4_kernel(rsFloat* x_arr, rsFloat* y_arr, rsFloat* z_arr, elsetrec* satrec_arr, int targsize, int numTimes, rsFloat* time_arr){

    // Only use threads up until the end of (targsize*numTimes) count
    int tid = int(threadIdx.x + (blockIdx.x * blockDim.x)); // Thread index
    int stride = blockDim.x * gridDim.x;  // Total number of threads spawned

    // Initialise local variables
    rsFloat r[3], v[3], acc[3], r1[3], v1[3], acc1[3];
    rsFloat sec, jd, jdFrac;
    int year, mon, day, hr, minute;

    for (int i = tid; i < (targsize*numTimes); i += stride){ // Will only iterate again if (targsize*numTimes) > (MaxBlocks*MaxThreads); very unlikely

        // Initialise the orbit at sgp4epoch
        int ti = i; // Target index
        while (ti >= targsize){ // Bring ti back to range of targets in array
            ti -= targsize;
        }
        SGP4Funcs::sgp4init(wgs84, 'i', satrec_arr[ti].satnum, (satrec_arr[ti].jdsatepoch + satrec_arr[ti].jdsatepochF) - 2433281.5, satrec_arr[ti].bstar,
                            satrec_arr[ti].ndot, satrec_arr[ti].nddot, satrec_arr[ti].ecco, satrec_arr[ti].argpo, satrec_arr[ti].inclo, satrec_arr[ti].mo,
                            satrec_arr[ti].no_kozai, satrec_arr[ti].nodeo, satrec_arr[ti]);

        // Call the propagator to get the initial state vector value; time in minutes (NO fraction notation)
        int time_idx = floor(float(i/targsize));
        rsFloat time_min = time_arr[time_idx]/60; // Get time from array (based on thread index) and convert to minutes
        SGP4Funcs::sgp4(satrec_arr[ti], time_min, r, v);

        // Define variables for ECEF conversion
        jd = satrec_arr[ti].jdsatepoch;
        jdFrac = satrec_arr[ti].jdsatepochF;
        SGP4Funcs::invjday(jd, jdFrac, year, mon, day, hr, minute, sec);
        rsFloat dut1 = 0.0;
        rsFloat lod  =  0.0;
        int timezone = 0;
        int dat  = 32;
        rsFloat jdut1, ttt, ut1, tut1, jdut1Frac, utc, tai, tt, jdtt, jdttFrac, tcg, tdb, ttdb, jdtdb, jdtdbFrac, tcb;

        // Convert TEME to ECEF
        astTime::days2mdhms(year, satrec_arr[ti].epochdays, mon, day, hr, minute, sec);
        convtime(year, mon, day, hr, minute, sec, timezone, dut1, dat,
                 ut1, tut1, jdut1, jdut1Frac, utc, tai, tt, ttt, jdtt, jdttFrac, tcg, tdb, ttdb, jdtdb, jdtdbFrac, tcb);
        teme_ecef(r, v, acc, eTo, r1, v1, acc1, ttt, jdut1, lod, 0, 0, 0);

        // Convert position to [m]
        x_arr[i] = r1[0]*1000;
        y_arr[i] = r1[1]*1000;
        z_arr[i] = r1[2]*1000;
        // printf("x: %.15lf, y: %.15lf, z: %.15lf\n", x_arr[i], y_arr[i], z_arr[i]);
    }
}

/// Anonymous Namespace
namespace {

    /// Process a Gamma target model entry
    RCSModel* ProcessGammaModel(pugi::xml_node &modelXML)
    {
        rsFloat k = GetChildRsFloat(modelXML, "k");
        return new RCSChiSquare(k);
    }

    /// Process a target XML entry
    void ProcessRCS(pugi::xml_node &rcsXML, Platform *platform, World *world, string name)
    {
        // Create target
        Target *target;

        // Get the RCS value
        if (!rcsXML)
            throw XmlImportException("Target "+name+" does not specify RCS.");
        string rcs_type = GetAttributeString(rcsXML, "type", "RCS attached to target '"+name+"' does not specify type.");

        // Handle the target type (isotropic, etc.)
        if (rcs_type == "isotropic") {
            rsFloat value = GetChildRsFloat(rcsXML, "value");
            target = CreateIsoTarget(platform, name, value);
        }
        else if (rcs_type == "SEM") {
            rsFloat diameter = GetChildRsFloat(rcsXML, "diameter");
            target = CreateSEMTarget(platform, name, diameter);
        }
        else if (rcs_type == "file") {
            string filename = GetAttributeString(rcsXML, "filename", "RCS attached to target '"+name+"' does not specify filename.");
            target = CreateFileTarget(platform, name, filename);
        }
        else {
            throw XmlImportException("RCS type "+rcs_type+" not currently supported.");
        }

        // Handle the target statistical model
        pugi::xml_node modelXML = rcsXML.child("model");
        if (modelXML) {

            // Get the mode type
            string model_type = GetAttributeString(modelXML, "type", "Model attached to target '"+name+"' does not specify type.");
            if (model_type == "constant") {
                RCSConst* model = new RCSConst();
                target->SetFluctuationModel(model);
            }
            else if ((model_type == "chisquare") || (model_type == "gamma")) {
                RCSModel *model = ProcessGammaModel(modelXML);
                target->SetFluctuationModel(model);
            }
            else {
                throw XmlImportException("Target fluctuation model type '"+model_type+"' not recognised.");
            }
        }

        // Add the target to the world
        world->Add(target);
    }

    /// Process a receiver XML entry
    Receiver *ProcessReceiver(pugi::xml_node &recvXML, Platform *platform, World *world)
    {
        rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "Loading receiver ");

        // Get the name of the receiver
        string name = GetAttributeString(recvXML, "name", "Receiver does not specify a name");
        Receiver *receiver = new Receiver(platform, name);

        // Get the name of the antenna
        string ant_name = GetAttributeString(recvXML, "antenna", "Receiver '" + string(name) + "' does not specify an antenna");

        rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "'%s'.\n", receiver->GetName().c_str());

        Antenna *antenna = world->FindAntenna(ant_name);
        if (!antenna)
            throw XmlImportException("Antenna with name '" + ant_name + "' does not exist when processing Receiver " + string(name));

        // Set the receiver's antenna
        receiver->SetAntenna(antenna);

        // Process the noise temperature tag
        try {
            rsFloat temperature;
            temperature = GetChildRsFloat(recvXML, "noise_temp");
            receiver->SetNoiseTemperature(temperature);
        }
            catch (XmlImportException &e) {
        }

        // Process the PRF tag
        rsFloat prf = GetChildRsFloat(recvXML, "prf");
        rsFloat skip = GetChildRsFloat(recvXML, "window_skip");
        rsFloat length = GetChildRsFloat(recvXML, "window_length");
        receiver->SetWindowProperties(length, prf, skip);

        // Get the name of the timing source
        string timing_name = GetAttributeString(recvXML, "timing", "Receiver '"+name+"' does not specify a timing source");
        ClockModelTiming *timing = new ClockModelTiming(timing_name);

        PrototypeTiming *proto = world->FindTiming(timing_name);
        if (!proto)
            throw XmlImportException("Timing source '" + timing_name + "' does not exist when processing receiver '"+name+"'");

        // Initialize the new model from the prototype model
        timing->InitializeModel(proto);

        // Set the receiver's timing source
        receiver->SetTiming(timing);

        // Get the NoDirect flag, which causes direct signals to be ignored
        bool nodirect = GetAttributeBool(recvXML, "nodirect", "NoDirect not specified.", false, true);
        if (nodirect) {
            receiver->SetFlag(rs::Receiver::FLAG_NODIRECT);
            rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "Ignoring direct signals for receiver '%s'\n", receiver->GetName().c_str());
        }

        // Get the NoPropagationLoss flag, which causes propagation loss to be ignored, eg. when propagation loss is calculated with AREPS
        bool noproploss = GetAttributeBool(recvXML, "nopropagationloss", "NoPropagationLoss not specified.", false, true);
        if (noproploss) {
            receiver->SetFlag(rs::Receiver::FLAG_NOPROPLOSS);
            rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "Ignoring propagation losses for receiver '%s'\n", receiver->GetName().c_str());
        }

        // Add the receiver to the world
        world->Add(receiver);

        return receiver;
    }

    /// Process a transmitter XML entry
    Transmitter *ProcessTransmitter(pugi::xml_node &transXML, Platform *platform, World *world)
    {
        // Get the name of the transmitter
        string name = GetAttributeString(transXML, "name", "Transmitter does not specify a name");
        rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "Loading transmitter '%s'.\n", name.c_str());

        // Get the name of the pulse
        string pulse_name = GetAttributeString(transXML, "pulse", "Transmitter '" + name + "' does not specify a pulse");

        // Get the pulse from the table of pulses
        RadarSignal *wave = world->FindSignal(pulse_name);
        if (!wave)
            throw XmlImportException("Pulse with name '" + pulse_name + "' does not exist");

        // Get the Pulse Repetition Frequency
        rsFloat prf = GetChildRsFloat(transXML, "prf");

        // Get the name of the antenna
        string ant_name = GetAttributeString(transXML, "antenna", "Transmitter '" + name + "' does not specify an antenna");
        Antenna *antenna = world->FindAntenna(ant_name);
        if (!antenna)
            throw XmlImportException("Antenna with name '" + ant_name + "' does not exist when processing transmitter " + string(name));

        // Get the name of the timing source
        string timing_name = GetAttributeString(transXML, "timing", "Transmitter '"+name+"' does not specify a timing source");

        ClockModelTiming *timing = new ClockModelTiming(name);

        PrototypeTiming *proto = world->FindTiming(timing_name);
        if (!proto)
            throw XmlImportException("Timing source '" + timing_name + "' does not exist when processing receiver "+name);

        // Initialize the new model from the prototype model
        timing->InitializeModel(proto);

        Transmitter *transmitter = new Transmitter(platform, name, true);

        // Set the transmitter's properties
        transmitter->SetWave(wave);
        transmitter->SetPRF(prf);
        transmitter->SetAntenna(antenna);
        transmitter->SetTiming(timing);

        // Add the transmitter to the world
        world->Add(transmitter);

        return transmitter;
    }

    /// Process a monostatic (Receiver and Transmitter sharing an antenna)
    void ProcessMonostatic(pugi::xml_node &transXML, Platform *platform, World *world)
    {
        Transmitter *trans = ProcessTransmitter(transXML, platform, world);
        Receiver *recv = ProcessReceiver(transXML, platform, world);
        trans->MakeMonostatic(recv);
        recv->MakeMonostatic(trans);
    }

    /// Process a motion path waypoint
    void ProcessWaypoint(pugi::xml_node &handXML, Path *path)
    {
    try {
        rsFloat x, y, z, t;
        x = GetChildRsFloat(handXML, "x");
        y = GetChildRsFloat(handXML, "y");
        z = GetChildRsFloat(handXML, "altitude");
        t = GetChildRsFloat(handXML, "time");
        Coord coord;
        coord.t = t;
        coord.pos = Vec3(x, y, z);
        path->AddCoord(coord);
    }
    catch (XmlImportException &e) {
        rsDebug::printf(rsDebug::RS_VERBOSE, "WARNING: Parse error while importing waypoint. Discarding waypoint.\n");
    }
    }

    /// Process the path's python attributes
    void ProcessPythonPath(pugi::xml_node &pathXML, Path *path)
    {
        // Initialize python, if it isn't done already
        rsPython::InitPython();

        // Get the python path definition
        try {
            pugi::xml_node tmp = pathXML.child("pythonpath");

            // Get the module and function name attributes
            string modname = GetAttributeString(tmp, "module", "Attribute module missing");
            string funcname = GetAttributeString(tmp, "function", "Attribute function missing");

            // Load the Path module
            if (modname != "" & funcname != "") // If a valid value was read
                path->LoadPythonPath(modname, funcname);
        }
        catch (XmlImportException &e) {
            rsDebug::printf(rsDebug::RS_VERBOSE, "%s", e.what());
        }
    }

    /// Process a MotionPath XML entry
    void ProcessSystemPath(pugi::xml_node &mpXML, Platform *platform)
    {
        // Get a pointer to the platform's path
        Path *path = platform->GetMotionPath();

        // Get the interpolation type
        try {
            string rottype = GetAttributeString(mpXML, "interpolation", "No interpolation type specified.");
            if (rottype == "linear")
                path->SetInterp(Path::RS_INTERP_LINEAR);
            else if (rottype == "cubic")
                path->SetInterp(Path::RS_INTERP_CUBIC);
            else if (rottype == "static")
                path->SetInterp(Path::RS_INTERP_STATIC);
            else if (rottype == "python") {
                path->SetInterp(Path::RS_INTERP_PYTHON);
                ProcessPythonPath(mpXML, path);
            }
            else {
                rsDebug::printf(rsDebug::RS_VERBOSE, "WARNING: Unsupported motion path interpolation type for platform '"+platform->GetName()+"'; defaulting to static...\n");
                path->SetInterp(Path::RS_INTERP_STATIC);
            }
        }
        catch (XmlImportException &e) {
            rsDebug::printf(rsDebug::RS_VERBOSE, "Motion path interpolation set to static for platform '"+platform->GetName()+"'...\n");
            path->SetInterp(Path::RS_INTERP_STATIC);
        }

        // Process all the PositionWaypoints
        for (pugi::xml_node tmp = mpXML.child("positionwaypoint"); tmp; tmp = tmp.next_sibling("positionwaypoint")) {
            ProcessWaypoint(tmp, path);
        }

        // Finalise the path after all the waypoints have been loaded
        path->Finalize();
    }

    /// Process an entry for a fixed rotation
    void ProcessRotationConstant(pugi::xml_node &mpXML, Platform* platform)
    {
        RotationPath* path = platform->GetRotationPath();
        try {
            RotationCoord start, rate;
            start.azimuth = GetChildRsFloat(mpXML, "startazimuth") * pi/180;        // Convert to rad
            start.elevation = GetChildRsFloat(mpXML, "startelevation") * pi/180;    // Convert to rad
            rate.azimuth = GetChildRsFloat(mpXML, "azimuthrate") * pi/180;          // Convert to rad/s
            rate.elevation = GetChildRsFloat(mpXML, "elevationrate") * pi/180;      // Convert to rad/s
            path->SetConstantRate(start, rate);
        }
        catch (XmlImportException &e) {
            rsDebug::printf(rsDebug::RS_VERBOSE, "WARNING: Parse error while importing constant rotation.\n");
        }
    }

    /// Process a platform, recursively processing all the elements that are attached to it
    void ProcessSystem(pugi::xml_node &platXML, World *world)
    {
        Platform *platform;

        // Create the platform, using the name from the element
        string name = GetAttributeString(platXML, "name", "ERROR: Platform must specify a name");
        platform = new Platform(string(name));

        // Add the platform to the world
        world->Add(platform);

        // Process all the receivers attached to the platform
        for (pugi::xml_node tmp = platXML.child("receiver"); tmp; tmp = tmp.next_sibling("receiver")) {
            ProcessReceiver(tmp, platform, world);
        }

        // Process all the transmitters attached to the platform
        for (pugi::xml_node tmp = platXML.child("transmitter"); tmp; tmp = tmp.next_sibling("transmitter")) {
            ProcessTransmitter(tmp, platform, world);
        }

        // Process all the monostatics attached to the platform
        for (pugi::xml_node tmp = platXML.child("monostatic"); tmp; tmp = tmp.next_sibling("monostatic")) {
            ProcessMonostatic(tmp, platform, world);
        }

        // Process all the motion paths attached to the platform
        for (pugi::xml_node tmp = platXML.child("motionpath"); tmp; tmp = tmp.next_sibling("motionpath")) {
            ProcessSystemPath(tmp, platform);
        }

        // Process all the rotation paths attached to the platform
        for (pugi::xml_node tmp = platXML.child("fixedrotation"); tmp; tmp = tmp.next_sibling("fixedrotation")) {
            ProcessRotationConstant(tmp, platform);
        }
    }

    /// Process a pulse entry of type rect
    void ProcessAnyPulse(pugi::xml_node &pulseXML, World *world, string name, int manmade, int galactic, int cosmic, int sun, int sky)
    {
        rsFloat carrier = GetChildRsFloat(pulseXML, "carrier");
        rsFloat power = GetChildRsFloat(pulseXML, "power");
        rsFloat bandwidth = GetChildRsFloat(pulseXML, "bandwidth");
        rsFloat length = GetChildRsFloat(pulseXML, "length");

        // If any external noise is included, find the frequency limits
        int i1, i2 = 0;
        rsFloat fc_range [] = {0.01,0.02,0.05,0.08,0.1,0.2,0.5,0.9,1,1.4,2,2.8,5,6.2,10,16,20,24,30,40,50,60,70,80,90,100}; // Define carrier frequency range
        rsFloat carrier_ghz = carrier/1000000000;   // Convert carrier to GHz

        if (manmade == 1 || galactic == 1 || cosmic == 1 || sky == 1 || sun == 1){

            // Find Fc within range
            for (int i = 0; i < (int(sizeof(fc_range)/sizeof(fc_range[0])) - 1); i++){
                if (fc_range[i] <= carrier_ghz && fc_range[i+1] >= carrier_ghz){   // Find the array entry range within which the carrier resides
                    i1 = i;
                    i2 = i + 1;
                }
            }
        }

        // Find the temperature for each noise source using linear interpolation; only apply if both top and bottom of range is non-zero
        rsFloat ant_temp = 0;               // Setup running total of noise temperature

        // Add manmade noise
        if (manmade == 1){
            rsFloat noise_manmade [] = {29000000, 5786260.713, 578626.0713, 145344.2978, 45961.90258, 10289.58829, 3253.853517, 1829.776299,
                                        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            if (noise_manmade[i1] != 0 && noise_manmade[i2] != 0){
                ant_temp += ((carrier_ghz - fc_range[i1])/(fc_range[i2] - fc_range[i1]))*(noise_manmade[i2] - noise_manmade[i1]) + noise_manmade[i1];       // Manmade
            }
        }

        // Add galactic noise for worst-case
        if (galactic == 1){
            rsFloat noise_galcentre [] = {365088.3694, 72844.70651, 11545.10795, 3650.883694, 3253.853517, 728.4470651, 91.70605214,
                                            23.03551881, 18.29776299, 10.28958829, 5.786260713, 4.596190258,
                                            0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            if (noise_galcentre[i1] != 0 && noise_galcentre[i2] != 0){
                ant_temp += ((carrier_ghz - fc_range[i1])/(fc_range[i2] - fc_range[i1]))*(noise_galcentre[i2] - noise_galcentre[i1]) + noise_galcentre[i1]; // Galactic centre
            }
        }

        // Add quiet sun noise
        if (sun == 1){
            rsFloat noise_quiet_sun [] = {0,0,0,0, 728447.0651, 515701.0289, 230355.1881, 163078.9843, 129538.2417, 91706.05214,
                                        64922.91302, 45961.90258, 20530.42775, 16307.89843, 11545.10795, 8173.310501,
                                        7284.470651, 6492.291302, 5786.260713, 5157.010289, 4596.190258, 4596.190258,
                                        4596.190258, 4596.190258, 4596.190258, 4596.190258};
            if (noise_quiet_sun[i1] != 0 && noise_quiet_sun[i2] != 0){
                ant_temp += ((carrier_ghz - fc_range[i1])/(fc_range[i2] - fc_range[i1]))*(noise_quiet_sun[i2] - noise_quiet_sun[i1]) + noise_quiet_sun[i1]; // Quiet sun
            }
        }

        // Add atmospheric noise for worst-case (i.e. horizontal)
        if (sky == 1){
            rsFloat noise_sky_hor [] = {4596190.258, 18297.76299, 0,0,0,0,0, 57.86260713, 64.92291302, 72.84470651, 81.73310501, 102.8958829, 115.4510795,
                                    129.5382417, 145.3442978, 182.9776299, 230.3551881, 290, 290, 290, 290, 290, 290, 290, 290, 290};
            if (noise_sky_hor[i1] != 0 && noise_sky_hor[i2] != 0){
                ant_temp += ((carrier_ghz - fc_range[i1])/(fc_range[i2] - fc_range[i1]))*(noise_sky_hor[i2] - noise_sky_hor[i1]) + noise_sky_hor[i1];       // Sky horizontal
            }
        }

        // Add cosmic noise
        if (sky == 1){
            rsFloat noise_cosmic [] = {0,0,0,0,0,0,0,0,0, 2.7,2.7,2.7,2.7,2.7,2.7,2.7,
                                        0,0,0,0,0,0,0,0,0,0};
            if (noise_cosmic[i1] != 0 && noise_cosmic[i2] != 0){
                ant_temp += ((carrier_ghz - fc_range[i1])/(fc_range[i2] - fc_range[i1]))*(noise_cosmic[i2] - noise_cosmic[i1]) + noise_cosmic[i1];          // Cosmic
            }
        }

        // Add pulse to simulation
        RadarSignal *wave = rsPulseFactory::CreatePulse(name, power, carrier, bandwidth, ant_temp, length, rsParameters::rate());
        world->Add(wave);
    }

    Antenna *ProcessPythonAntenna(pugi::xml_node &antXML, const string &name)
    {
      // Initialize python, if it isn't done already
      rsPython::InitPython();

      // Get the module and function name attributes
      string modname = GetAttributeString(antXML, "module", "Attribute module missing");
      string funcname = GetAttributeString(antXML, "function", "Attribute function missing");

      // Create the antenna
      return rs::CreatePythonAntenna(name, modname, funcname);
    }


    Antenna *ProcessXMLAntenna(pugi::xml_node &antXML, const string &name)
    {
        //Get the module and function name attributes
        string filename = GetAttributeString(antXML, "filename", "Antenna definition must specify a filename");

        //Create the antenna
        return rs::CreateXMLAntenna(name, filename);
    }

    Antenna *ProcessFileAntenna(pugi::xml_node &antXML, const string &name)
    {
        //Get the module and function name attributes
        string filename = GetAttributeString(antXML, "filename", "Antenna definition must specify a filename");

        //Create the antenna
        return rs::CreateFileAntenna(name, filename);
    }

    Antenna *ProcessSincAntenna(pugi::xml_node &antXML, const string &name)
    {
        rsFloat alpha = GetChildRsFloat(antXML, "alpha");
        rsFloat beta = GetChildRsFloat(antXML, "beta");
        rsFloat gamma = GetChildRsFloat(antXML, "gamma");
        return rs::CreateSincAntenna(name, alpha, beta, gamma);
    }

    Antenna *ProcessGaussianAntenna(pugi::xml_node &antXML, const string &name)
    {
        rsFloat azscale = GetChildRsFloat(antXML, "azscale");
        rsFloat elscale = GetChildRsFloat(antXML, "elscale");
        return rs::CreateGaussianAntenna(name, azscale, elscale);
    }

    Antenna *ProcessParabolicAntenna(pugi::xml_node &antXML, const string &name)
    {
        rsFloat diameter = GetChildRsFloat(antXML, "diameter");
        return rs::CreateParabolicAntenna(name, diameter);
    }

    void ProcessAntenna(pugi::xml_node &antXML, World *world)
    {
        // Get the name of the antenna
        string ant_name = GetAttributeString(antXML, "name", "Antennas must specify a name");

        // Get the type of the antenna
        string ant_pattern = GetAttributeString(antXML, "pattern", "Antennas must specify a pattern");
        Antenna *antenna;
        if (ant_pattern == "isotropic")
            antenna = CreateIsotropicAntenna(ant_name);
        else if (ant_pattern == "file")
            antenna = ProcessFileAntenna(antXML, ant_name);
        else if (ant_pattern == "xml")
            antenna = ProcessXMLAntenna(antXML, ant_name);
        else if (ant_pattern == "python")
            antenna = ProcessPythonAntenna(antXML, ant_name);
        else if (ant_pattern == "sinc")
            antenna = ProcessSincAntenna(antXML, ant_name);
        else if (ant_pattern == "gaussian")
            antenna = ProcessGaussianAntenna(antXML, ant_name);
        else if (ant_pattern == "parabolic")
            antenna = ProcessParabolicAntenna(antXML, ant_name);
        else
            throw XmlImportException("Antenna specified unrecognised gain pattern '" + ant_pattern + "'");

        // Notify the debug log
        rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "Loading %s antenna '%s'.\n", ant_pattern.c_str(), ant_name.c_str());

        // Load the efficiency factor
        try {
            rsFloat factor = GetChildRsFloat(antXML, "efficiency");
            antenna->SetEfficiencyFactor(factor);
        } catch (XmlImportException &xe) {
            rsDebug::printf(rsDebug::RS_VERBOSE, "Antenna '%s' does not specify efficiency, assuming unity.\n", ant_name.c_str());
        }

        // Add it to the world
        world->Add(antenna);
    }

    /// Process a multipath surface and add it to the world
    void ProcessMultipath(pugi::xml_node &mpXML, World *world)
    {
        // Get the reflecting factor
        rsFloat factor = GetChildRsFloat(mpXML, "factor");
        rsFloat a = GetChildRsFloat(mpXML, "a");
        rsFloat b = GetChildRsFloat(mpXML, "b");
        rsFloat c = GetChildRsFloat(mpXML, "c");
        rsFloat d = GetChildRsFloat(mpXML, "d");

        // Create the multipath object
        MultipathSurface* mps = new MultipathSurface(a, b, c, d, factor);

        // Add it to the world
        world->AddMultipathSurface(mps);
    }

    /// Process a timing source and add it to the world
    void ProcessTiming(pugi::xml_node &antXML, World *world)
    {
        // Get the name of the antenna
        string name = GetAttributeString(antXML, "name", "Timing sources must specify a name");
        PrototypeTiming *timing = new PrototypeTiming(name);

        // Process all the clock entries
        for (pugi::xml_node plat = antXML.child("noise_entry"); plat; plat = plat.next_sibling("noise_entry")) {
            rsFloat alpha = GetChildRsFloat(plat, "alpha");
            rsFloat weight = GetChildRsFloat(plat, "weight");
            timing->AddAlpha(alpha, weight);
        }

        // Process the frequency offset
        try {
            rsFloat offset = GetChildRsFloat(antXML, "freq_offset");
            if (offset != -1.0) // If a valid value was read
                timing->AddFreqOffset(offset);
        }
        catch (XmlImportException &xe) {
        }
        try {
            rsFloat stdev = GetChildRsFloat(antXML, "random_freq_offset");
            if (stdev != -1.0) // If a valid value was read
                timing->AddRandomFreqOffset(stdev);
        }
        catch (XmlImportException &xe) {
        }

        // Process the phase offset
        try {
            rsFloat offset = GetChildRsFloat(antXML, "phase_offset");
            if (offset != -1.0) // If a valid value was read
                timing->AddPhaseOffset(offset);
        }
        catch (XmlImportException &xe) {
        }
        try {
            rsFloat stdev = GetChildRsFloat(antXML, "random_phase_offset");
            if (stdev != -1.0) // If a valid value was read
                timing->AddRandomPhaseOffset(stdev);
        }
        catch (XmlImportException &xe) {
        }

        // Process the frequency
        try {
            rsFloat freq = GetChildRsFloat(antXML, "frequency");
            timing->SetFrequency(freq);
        }
        catch (XmlImportException &xe) {

        // If there is no frequency, we default to the system sample frequency
        timing->SetFrequency(rsParameters::rate());
        rsDebug::printf(rsDebug::RS_VERBOSE, "Timing source '%s' defaulting to frequency %8.2f Hz...\n", name.c_str(), rsParameters::rate());
        }
        // Process the synconpulse tag
        bool sync = GetAttributeBool(antXML, "synconpulse", "No sync on pulse specified.", false, true);
        if (sync)
            timing->SetSyncOnPulse();

        // Add it to the world
        world->Add(timing);
    }

    /// Process the <parameters> element
    void ProcessParameters(pugi::xml_node &parameters)
    {
        // Get the simulation start and end times
        rsParameters::modify_parms()->SetTime(GetChildRsFloat(parameters, "starttime"), GetChildRsFloat(parameters, "endtime"));

        // Get the propagation speed in air
        try {
            rsFloat c = GetChildRsFloat(parameters, "c");
            rsParameters::modify_parms()->SetC(c);
        }
        catch (XmlImportException &xe)
        {
            rsDebug::printf(rsDebug::RS_VERBOSE, "Using default value of c: %f(m/s)\n", rsParameters::c());
        }

        // Get the export sampling rate
        try {
            rsFloat rate = GetChildRsFloat(parameters, "rate");
            if (rate != -1.0) // If a valid value was read
                rsParameters::modify_parms()->SetRate(rate);
        }
        catch (XmlImportException &xe)
        {
            rsDebug::printf(rsDebug::RS_VERBOSE, "Using default sampling rate...\n");
        }

        // Get the cw Interpolation rate
        try {
            rsFloat cwrate = GetChildRsFloat(parameters, "interprate");
            if (cwrate != -1.0) // If a valid value was read
                rsParameters::modify_parms()->SetCWSampleRate(cwrate);
        }
        catch (XmlImportException &xe)
        {
            //rsDebug::printf(rsDebug::RS_VERBOSE, "Using default value of CW position interpolation rate: %g\n", rsParameters::cw_sample_rate());
        }

        // Get the random seed
        try {
          rsFloat seed = GetChildRsFloat(parameters, "randomseed");
          if (seed != -1.0) // If a valid value was read
              rsParameters::modify_parms()->SetRandomSeed(static_cast<unsigned int>(std::fabs(seed)));
        }
        catch (XmlImportException &xe)
        {
            //rsDebug::printf(rsDebug::RS_VERBOSE, "Using random seed from clock(): %d\n", rsParameters::random_seed());
        }

        // Get the number of ADC bits to simulate
        try {
            rsFloat adc_bits  = GetChildRsFloat(parameters, "adc_bits");
            if (adc_bits != -1.0) // If a valid value was read
                rsParameters::modify_parms()->SetADCBits(static_cast<unsigned int>(std::floor(adc_bits)));
        }
        catch (XmlImportException &xe)
        {
            //rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "Using full precision simulation.\n");
        }

        // Get the oversampling ratio
        try {
            rsFloat ratio  = GetChildRsFloat(parameters, "oversample");
            rsParameters::modify_parms()->SetOversampleRatio(static_cast<unsigned int>(std::floor(ratio)));
        }
        catch (XmlImportException &xe)
        {
            //rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "Oversampling not in use.\n");
        }

        // Process the "export" tag
        pugi::xml_node exporttag = parameters.child("export");
        if (exporttag) {
            bool export_xml = GetAttributeBool(exporttag, "xml", "No XML export.", rsParameters::export_xml(), true);
            bool export_csv = GetAttributeBool(exporttag, "csv", "No CSV export.", rsParameters::export_csv(), true);
            bool export_binary = GetAttributeBool(exporttag, "binary", "No binary export.", rsParameters::export_binary(), true);
            rsParameters::modify_parms()->SetExporters(export_xml, export_csv, export_binary);
        }
    }

    /// Process the XML tree, starting at the root
    void ProcessDocument(pugi::xml_node &root, World *world, bool included, int MaxThreads, int MaxBlocks)
    {
        // Process the parameters
        pugi::xml_node parameters = root.child("parameters");
        ProcessParameters(parameters);
        rsFloat timestep = GetChildRsFloat(parameters, "targetsample"); // Interval at which to sample targets
        int targsize = GetChildRsFloat(parameters, "targsize"); // Number of targets
        rsFloat endtime = GetChildRsFloat(parameters, "endtime");  // Simulation endtime

        // Process the noise sources
        rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "Creating noise sources.\n");
        pugi::xml_node noises = root.child("noise");
        int manmade, galactic, cosmic, sun, sky;
        try {
            manmade = int(GetChildRsFloat(noises, "manmade"));
            galactic = int(GetChildRsFloat(noises, "galactic"));
            cosmic = int(GetChildRsFloat(noises, "cosmic"));
            sun = int(GetChildRsFloat(noises, "sun"));
            sky = int(GetChildRsFloat(noises, "sky"));
        }
        catch (XmlImportException &xe)
        {
            rsDebug::printf(rsDebug::RS_VERBOSE, "No external noise source selection found; ignoring all.\n");
            manmade = 0;
            galactic = 0;
            cosmic = 0;
            sun = 0;
            sky = 0;
        }

        // Process all the pulses
        rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "Reading pulse information.\n");
        for (pugi::xml_node pulse = root.child("pulse"); pulse; pulse = pulse.next_sibling("pulse")) {
            string pulse_name = GetAttributeString(pulse, "name", "Pulses must specify a name.");
            ProcessAnyPulse(pulse, world, pulse_name, manmade, galactic, cosmic, sun, sky);
        }

        // Process all the antennas
        for (pugi::xml_node ant = root.child("antenna"); ant; ant = ant.next_sibling("antenna")) {
            ProcessAntenna(ant, world);
        }

        // Process all the timing sources
        pugi::xml_node timing = root.child("timing");
        ProcessTiming(timing, world);

        // Process all the multipath surfaces
        for (pugi::xml_node mp = root.child("multipath"); mp; mp = mp.next_sibling("multipath")) {
            ProcessMultipath(mp, world);
        }

        // Process all the radar system platforms
        rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "Loading radar systems.\n");
        for (pugi::xml_node sys = root.child("system"); sys; sys = sys.next_sibling("system")) {
            ProcessSystem(sys, world); // Recursively process the platforms
        }

        /// Process all the target platforms
        rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "Loading targets...\n");

        //gettimeofday(&t3, 0);

        vector<elsetrec> satrec_vec(targsize);
        vector<Path*> path_vec(targsize);
        int t = 0;  // Target iterator
        for (pugi::xml_node targ = root.child("target"); targ; targ = targ.next_sibling("target")) {

            // Create the platform, using the name from the element
            Platform *platform;
            string name = GetAttributeString(targ, "name", "ERROR: Platform must specify a name");
            platform = new Platform(string(name));

            // Add the platform to the world
            world->Add(platform);

            // Process target RCS
            pugi::xml_node tmp = targ.child("rcs");
            ProcessRCS(tmp, platform, world, name);

            // Process the path of the target
            tmp = targ.child("motionpath");

            // Get a pointer to the platform's path
            Path *path = platform->GetMotionPath();

            // Get the interpolation type
            try {
                string rottype = GetAttributeString(tmp, "interpolation", "No interpolation type specified.");
                if (rottype == "linear")
                    path->SetInterp(Path::RS_INTERP_LINEAR);
                else if (rottype == "cubic")
                    path->SetInterp(Path::RS_INTERP_CUBIC);
                else if (rottype == "static")
                    path->SetInterp(Path::RS_INTERP_STATIC);
                else if (rottype == "python") {
                    path->SetInterp(Path::RS_INTERP_PYTHON);
                    ProcessPythonPath(tmp, path);
                }
                else {
                    rsDebug::printf(rsDebug::RS_VERBOSE, "WARNING: Unsupported motion path interpolation type for platform '"+platform->GetName()+"'; defaulting to static...\n");
                    path->SetInterp(Path::RS_INTERP_STATIC);
                }
            }
            catch (XmlImportException &e) {
                rsDebug::printf(rsDebug::RS_VERBOSE, "Motion path interpolation set to static for platform '"+platform->GetName()+"'...\n");
                path->SetInterp(Path::RS_INTERP_STATIC);
            }

            // Process all the PositionWaypoints
            if (tmp.child("positionwaypoint")) { // Check if target uses waypoints
                pugi::xml_node tmp1 = tmp.child("positionwaypoint");
                for (tmp1 = tmp.child("positionwaypoint"); tmp1; tmp1 = tmp1.next_sibling("positionwaypoint")) {
                    ProcessWaypoint(tmp1, path);
                }

                // Finalise the path after all the waypoints have been loaded
                path->Finalize();
            }
            // else {  // Otherwise process a TLE
            //     pugi::xml_node tleXML = tmp.child("tle"); // Process TLE entry

            //     // Read TLE set
            //     string TLE_LINE1 = GetChildText(tleXML, "line1");
            //     string TLE_LINE2 = GetChildText(tleXML, "line2");
            //     char* longstr1 = &TLE_LINE1[0];
            //     char* longstr2 = &TLE_LINE2[0];

            //     // Get satrec properties
            //     rsFloat startmfe, stopmfe, deltamin = 0.0;
            //     elsetrec satrec;
            //     SGP4Funcs::twoline2rv(longstr1, longstr2, 'm', 'e', 'i', wgs84, startmfe, stopmfe, deltamin, satrec);

            //     // Push back vectors
            //     satrec_vec[t] = satrec;
            //     path_vec[t] = path;
            // }

            // Process all the rotation paths attached to the platform
            tmp = targ.child("fixedrotation");
            ProcessRotationConstant(tmp, platform);

            // Update target iterator
            t++;
        }
        //
        // // SGP4 timesteps
        // int numTimes = ceil(endtime/timestep) + 1;  // Number of timesteps
        // vector<rsFloat> time_vec(numTimes); // Store all timestep values
        // rsFloat current_time = 0;
        // for (int i = 0; i < numTimes; i++){ // SGP4 for t = 0, t = timestep, t = 2*timestep, ..., t = endtime
        //     time_vec[i] = current_time; // Store time value
        //     current_time = current_time + timestep; // Update current time
        //     if (current_time >= endtime)  // If current_time is beyond the endtime, bring it down
        //         current_time = endtime;   // Update current time to endtime
        // }
        //
        // // Create arrays for GPU
        // vector<rsFloat> x_arr(targsize*numTimes); // Target x array
        // vector<rsFloat> y_arr(targsize*numTimes); // Target y array
        // vector<rsFloat> z_arr(targsize*numTimes); // Target z array
        //
        // // Create device copies of arrays for kernel
        // rsFloat *d_x_arr;
        // rsFloat *d_y_arr;
        // rsFloat *d_z_arr;
        // elsetrec* d_satrec_arr; // Device array of satrec_arr
        // rsFloat* d_time_arr; // Device array of time_arr
        //
        // // Allocate memory for device variable copies
        // cudaMalloc((void **)&d_x_arr, sizeof(rsFloat)*targsize*numTimes);
        // cudaMalloc((void **)&d_y_arr, sizeof(rsFloat)*targsize*numTimes);
        // cudaMalloc((void **)&d_z_arr, sizeof(rsFloat)*targsize*numTimes);
        // cudaMalloc((void **)&d_time_arr, sizeof(rsFloat)*numTimes);
        // cudaMalloc((void **)&d_satrec_arr, sizeof(elsetrec)*targsize);  // Allocate GPU memory
        // cudaCheckErrors("Malloc fail");
        //
        // // Copy arrays to GPU
        // cudaMemcpy(d_satrec_arr, satrec_vec.data(), sizeof(elsetrec)*targsize, cudaMemcpyHostToDevice); // Copy from host to GPU memory
        // cudaMemcpy(d_time_arr, time_vec.data(), sizeof(rsFloat)*numTimes, cudaMemcpyHostToDevice); // Copy from host to GPU memory
        // cudaCheckErrors("Memory fail");
        //
        // // Call kernel
        // if ((targsize*numTimes) <= MaxThreads){ // If there are fewer targets than the number of threads in one block (or an equal number of targets)
        //     SGP4_kernel<<<1, (targsize*numTimes)>>>(d_x_arr, d_y_arr, d_z_arr, d_satrec_arr, targsize, numTimes, d_time_arr);
        // }
        // else if ((targsize*numTimes) > (MaxThreads*MaxBlocks)){   // If there are more targets than the maximum number of parallel threads (across all blocks)
        //     SGP4_kernel<<<MaxBlocks, MaxThreads>>>(d_x_arr, d_y_arr, d_z_arr, d_satrec_arr, targsize, numTimes, d_time_arr);
        // }
        // else{ // If number of targets requires more than 1 block, but not all of them
        //     SGP4_kernel<<<(((targsize*numTimes) + (MaxThreads - 1))/MaxThreads), MaxThreads>>>(d_x_arr, d_y_arr, d_z_arr, d_satrec_arr, targsize, numTimes, d_time_arr);
        // }
        // cudaDeviceSynchronize();
        // cudaCheckErrors("Kernel fail");
        // cudaMemcpy(x_arr.data(), d_x_arr, sizeof(rsFloat)*targsize*numTimes, cudaMemcpyDeviceToHost); // Copy from GPU to host memory
        // cudaMemcpy(y_arr.data(), d_y_arr, sizeof(rsFloat)*targsize*numTimes, cudaMemcpyDeviceToHost); // Copy from GPU to host memory
        // cudaMemcpy(z_arr.data(), d_z_arr, sizeof(rsFloat)*targsize*numTimes, cudaMemcpyDeviceToHost); // Copy from GPU to host memory
        // cudaCheckErrors("Memory fail");
        //
        // // Add coords to path
        // for (int n = 0; n < numTimes; n++){
        //     for (int t = 0; t < targsize; t++){
        //         Coord coord;
        //         coord.t = time_vec[n];
        //         int targ_index = t + (targsize*n);  // Account for timestep variation in arrays
        //         coord.pos = Vec3(x_arr[targ_index], y_arr[targ_index], z_arr[targ_index]);
        //         path_vec[t]->AddCoord(coord);
        //     }
        // }
        //
        // // Free up memory
        // cudaFree(d_x_arr);
        // cudaFree(d_y_arr);
        // cudaFree(d_z_arr);
        // cudaFree(d_satrec_arr);
        // cudaFree(d_time_arr);
        // cudaCheckErrors("Delete fail");
        //
        // // Finalise the paths after all the waypoints have been loaded
        // for (int t = 0; t < targsize; t++){
        //     path_vec[t]->Finalize();
        // }

        //gettimeofday(&t4, 0);
        //rsFloat time = (1000000.0*(t4.tv_sec-t3.tv_sec) + t4.tv_usec-t3.tv_usec)/1000.0;
        //printf("Time to write: %.15lf milliseconds \n", time);
    }

} // End of anonymous namespace

/// Load a XML file into the world with the given filename
void xml::LoadXMLFile(string filename, World *world, int MaxThreads, int MaxBlocks)
{
    pugi::xml_document doc;
    pugi::xml_parse_result document = doc.load_file(filename.c_str());
    if (!document)
        throw std::runtime_error("Cannot open script file.");

    // Process the XML document
    pugi::xml_node root = doc.child("simulation");
    ProcessDocument(root, world, false, MaxThreads, MaxBlocks);

    // Create multipath duals of all objects, if a surface was added
    world->ProcessMultipath();
}
