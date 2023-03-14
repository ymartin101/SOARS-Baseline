/// rsantenna.cpp - Implementation of Antenna Class
/// Marc Brooker, 20 July 2006
/// Edited by Craig Tong
/// Edited by Yaaseen Martin, 27 August 2019

#include <stdexcept>
#include <stdio.h>
#include <algorithm>
#include "rsantenna.cuh"
#include "rsdebug.cuh"
#include <cmath>
#include "rsportable.cuh"
#include "rsinterp.cuh"
#include "rspattern.cuh"
#include "rsradarwaveform.cuh"
#include "pugixml.hpp"

using namespace rs;
using namespace rsAntenna;

// Gets child text as a rsFloat; see GetChildText for usage description
rsFloat GetChildRsDouble(pugi::xml_node &parent, std::string childname)
{
    std::string data = parent.child_value(childname.c_str());

    // Parse the first rsFloat from the data
    rsFloat result;
    std::istringstream iss(data);
    iss >> result;
    return result;
}

namespace {
    // Return sin(x)/x
    rsFloat sinc(rsFloat theta)
    {
        return std::sin(theta)/(theta+std::numeric_limits<rsFloat>::epsilon());
    }

    // Return the first order, first kind bessel function of x, divided by x
    rsFloat j1c(double x)
    {
        if (x == 0)                         // Maximum gain when x = 0
            return 0.5;                     // Maximum gain when (2*j1c(x))^2 = 1, i.e. j1c(x) = 0.5
        return rsPortable::BesselJ1(x)/(x); // If NOT maximum, find J1(x)/x; slightly off boresight --> (J1(x)/x) is just below 0.5
    }
}

// Default constructor for the antenna
Antenna::Antenna(const std::string& name):
    lossFactor(1), // Antenna efficiency default is unity
    name(name)
{
}

// Antenna destructor
Antenna::~Antenna()
{
}

// Set the efficiency factor of the antenna
void Antenna::SetEfficiencyFactor(rsFloat loss)
{
    if (loss > 1)
        rsDebug::printf(rsDebug::RS_IMPORTANT, "Using greater than unity antenna efficiency, results might be inconsistent with reality.\n");
    lossFactor = loss;
}

// Get the efficiency factor
rsFloat Antenna::GetEfficiencyFactor() const
{
    return lossFactor;
}

// Return the name of the antenna
std::string Antenna::GetName() const
{
    return name;
}

// Get the angle (in radians) off boresight
rsFloat Antenna::GetAngle(const SVec3 &angle, const SVec3 &refangle) const
{
    // Get the angle off boresight
    SVec3 normangle(angle);
    normangle.length = 1;
    Vec3 cangle(normangle);
    Vec3 ref(refangle);

    return std::acos(DotProduct(cangle, ref));
}

// / Get the noise temperature of the antenna in a particular direction
rsFloat Antenna::GetNoiseTemperature(const SVec3 &angle) const
{
    return 0; // TODO: Antenna noise temperature calculation
}

/// Isotropic Implementation

// Default constructor
Isotropic::Isotropic(const std::string& name):
    Antenna(name)
{
}

// Default destructor
Isotropic::~Isotropic()
{
}

// Return the gain of the antenna
rsFloat Isotropic::GetGain(const SVec3 &angle, const SVec3 &refangle, rsFloat wavelength) const
{
    return GetEfficiencyFactor();
}

///  Gaussian Implementation

// Constructor
Gaussian::Gaussian(const std::string& name, rsFloat azscale, rsFloat elscale):
    Antenna(name),
    azscale(azscale),
    elscale(elscale)
{
}

// Destructor
Gaussian::~Gaussian()
{
}

// Get the gain at an angle
rsFloat Gaussian::GetGain(const rs::SVec3 &angle, const rs::SVec3 &refangle, rsFloat wavelength) const
{
    SVec3 a = angle - refangle;
    rsFloat azfactor = std::exp(-a.azimuth*a.azimuth*azscale);
    rsFloat elfactor = std::exp(-a.elevation*a.elevation*elscale);
    return azfactor*elfactor;
}

/// Sinc Implemetation

// Constructor
Sinc::Sinc(const std::string& name, rsFloat alpha, rsFloat beta, rsFloat gamma):
  Antenna(name),
  alpha(alpha),
  beta(beta),
  gamma(gamma)
{
}

// Default destructor
Sinc::~Sinc()
{
}

//  Return the gain of the antenna at an angle
rsFloat Sinc::GetGain(const SVec3 &angle, const SVec3 &refangle, rsFloat wavelength) const
{
    // Get the angle off boresight
    rsFloat theta = GetAngle(angle, refangle);

    // FIX 2015/02/18 CA Tong:
    // std::pow<double> returns NaN for negative bases with certain fractional indices as they create an uneven
    // root which is, of course, a complex number. Inputs therefore need to be cast to complex numbers before the
    // calculation as this will return complex number. Then return the magnitude as the beam gain.

    rs::rsComplex complexSinc(::sinc(beta*theta), 0.0);
    rs::rsComplex complexGamma(gamma, 0.0);

    // See "Sinc Pattern" in doc/equations.tex for equation used here
    rsComplex complexGain = alpha * std::pow(complexSinc, complexGamma) * GetEfficiencyFactor();

    return std::abs(complexGain);
}

///  SquareHorn Implementation

// Constructor
SquareHorn::SquareHorn(const std::string& name, rsFloat dimension):
    Antenna(name),
    dimension(dimension)
{
}

// Destructor
SquareHorn::~SquareHorn()
{
}

// Return the gain of the antenna
// See doc/equations.tex for details
rsFloat SquareHorn::GetGain(const SVec3 &angle, const SVec3 &refangle, rsFloat wavelength) const
{
    rsFloat Ge = 4*M_PI*dimension*dimension/(wavelength*wavelength);
    rsFloat x = M_PI*dimension*std::sin(GetAngle(angle, refangle))/wavelength;
    rsFloat gain = Ge*std::pow(::sinc(x), 2);
    return gain*GetEfficiencyFactor();
}

///  Parabolic Dish Antenna

//  Constructor
ParabolicReflector::ParabolicReflector(const std::string& name, rsFloat diameter):
    Antenna(name),
    diameter(diameter)
{
}

// Destructor
ParabolicReflector::~ParabolicReflector()
{
}

// Return the gain of the antenna
// See doc/equations.tex for details
rsFloat ParabolicReflector::GetGain(const SVec3 &angle, const SVec3 &refangle, rsFloat wavelength) const
{
    rsFloat Ge = std::pow(M_PI*diameter/wavelength, 2);                         // Maximum possible gain
    rsFloat x = M_PI*diameter*std::sin(GetAngle(angle, refangle))/wavelength;   // Maximum gain when x = 0, i.e. j1c(x) = 0.5
    rsFloat gain = Ge*std::pow(2*::j1c(x), 2);                                  // Overall gain; maximum when (2*j1c(x))^2 = 1
    return gain*GetEfficiencyFactor();
}

///  FileAntenna implementation

// Constructor
FileAntenna::FileAntenna(const std::string& name, const std::string &filename):
    Antenna(name)
{
    pattern = new Pattern(filename);
}

// Default destructor
FileAntenna::~FileAntenna()
{
    delete pattern;
}

// Get the gain at an angle
rsFloat FileAntenna::GetGain(const rs::SVec3 &angle, const rs::SVec3 &refangle, rsFloat wavelength) const
{
    SVec3 a1 = angle;
    SVec3 a2 = refangle;
    SVec3 in_angle = (a1-a2);

    return pattern->GetGain(in_angle)*GetEfficiencyFactor();
}

///  Antenna with pattern loaded from an XML file

// Constructor
XMLAntenna::XMLAntenna(const std::string& name, const std::string &filename):
    Antenna(name)
{
    // Classes to interpolate across elevation and azimuth
    azi_samples = new InterpSet();
    elev_samples = new InterpSet();

    // Load the XML antenna description data
    LoadAntennaDescription(filename);
}

// Destructor
XMLAntenna::~XMLAntenna()
{
    // Clean up the interpolation classes
    delete azi_samples;
    delete elev_samples;
}

// Return the gain of the antenna; see doc/equations.tex for details
rsFloat XMLAntenna::GetGain(const SVec3 &angle, const SVec3 &refangle, rsFloat wavelength) const
{
    SVec3 t_angle = angle-refangle;
    rsFloat azi_gain = azi_samples->Value(std::fabs(t_angle.azimuth));
    rsFloat elev_gain = elev_samples->Value(std::fabs(t_angle.elevation));
    return azi_gain*elev_gain*max_gain*GetEfficiencyFactor();
}

namespace {

    // Load samples of gain along an axis (not a member of XMLAntenna)
    void LoadAntennaGainAxis(InterpSet *set, pugi::xml_node &axisXML)
    {
        rsFloat angle;
        rsFloat gain;

        // Step through the XML file and load all the gain samples
        for (pugi::xml_node tmp = axisXML.child("gainsample"); tmp; tmp = tmp.next_sibling("gainsample")) {

            // Load the angle of the gain sample
            angle = GetChildRsDouble(tmp, "angle");

            // Load the gain of the gain sample
            gain = GetChildRsDouble(tmp, "gain");

            // Load the values into the interpolation table
            set->InsertSample(angle, gain);
        }
    }
}

// Load the antenna description file
void XMLAntenna::LoadAntennaDescription(const std::string& filename)
{
    pugi::xml_document doc;
    pugi::xml_parse_result document = doc.load_file(filename.c_str());
    if (!document)
        throw std::runtime_error("[ERROR] Could not load antenna description "+filename);

    // Get the XML root node
    pugi::xml_node root = doc.child("filetarget");

    // Load the gain samples along the elevation axis
    pugi::xml_node tmp = root.child("elevation");
    if (!tmp)
        throw std::runtime_error("ERROR: Malformed XML in antenna description: No elevation pattern definition.");
    LoadAntennaGainAxis(elev_samples, tmp);

    // Load the gain samples along the azimuth axis
    tmp = root.child("azimuth");
    if (!tmp)
        throw std::runtime_error("ERROR: Malformed XML in antenna description: No azimuth pattern definition.");
    LoadAntennaGainAxis(azi_samples, tmp);

    // Normalize the antenna patterns and calculate the max gain
    max_gain = std::max(azi_samples->Max(), elev_samples->Max());
    elev_samples->Divide(max_gain);
    azi_samples->Divide(max_gain);

}

///  Antenna with gain pattern calculated by a Python program

// Constructor
PythonAntenna::PythonAntenna(const std::string& name, const std::string &module, const std::string& function):
    Antenna(name),
    py_antenna(module, function)
{
}

// Destructor
PythonAntenna::~PythonAntenna()
{
}

// Return the gain of the antenna
rsFloat PythonAntenna::GetGain(const SVec3 &angle, const SVec3 &refangle, rsFloat wavelength) const
{
    SVec3 angle_bore = angle - refangle; // Calculate the angle off boresight
    rsFloat gain = py_antenna.GetGain(angle_bore);
    return gain*GetEfficiencyFactor();
}

///  Functions to create Antenna objects with a variety of properties

// Create an isotropic antenna with the specified name
Antenna* rs::CreateIsotropicAntenna(const std::string &name)
{
    rsAntenna::Isotropic *iso = new rsAntenna::Isotropic(name);
    return iso;
}

// Create a Sinc pattern antenna with the specified name, alpha and beta
Antenna* rs::CreateSincAntenna(const std::string &name, rsFloat alpha, rsFloat beta, rsFloat gamma)
{
    rsAntenna::Sinc *sinc = new rsAntenna::Sinc(name, alpha, beta, gamma);
    return sinc;
}

// Create a Gaussian pattern antenna
Antenna* rs::CreateGaussianAntenna(const std::string &name, rsFloat azscale, rsFloat elscale)
{
    rsAntenna::Gaussian *gau = new rsAntenna::Gaussian(name, azscale, elscale);
    return gau;
}

// Create a square horn antenna
Antenna* rs::CreateHornAntenna(const std::string &name, rsFloat dimension)
{
    rsAntenna::SquareHorn *sq = new rsAntenna::SquareHorn(name, dimension);
    return sq;
}

// Create a parabolic reflector antenna
Antenna* rs::CreateParabolicAntenna(const std::string &name, rsFloat diameter)
{
    rsAntenna::ParabolicReflector *pd = new rsAntenna::ParabolicReflector(name, diameter);
    return pd;
}

// Create an antenna with it's gain pattern stored in an XML file
Antenna* rs::CreateXMLAntenna(const std::string &name, const std::string &file)
{
    rsAntenna::XMLAntenna *fa = new rsAntenna::XMLAntenna(name, file);
    return fa;
}

// Create an antenna with it's gain pattern stored in an XML file
Antenna* rs::CreateFileAntenna(const std::string &name, const std::string &file)
{
    rsAntenna::FileAntenna *fa = new rsAntenna::FileAntenna(name, file);
    return fa;
}

// Create an antenna with gain pattern described by a Python program
Antenna* rs::CreatePythonAntenna(const std::string &name, const std::string &module, const std::string &function)
{
    rsAntenna::PythonAntenna* pa = new rsAntenna::PythonAntenna(name, module, function);
    return pa;
}
