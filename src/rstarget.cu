/// rstarget.cpp - Classes for targets and target RCS
/// Marc Brooker, 11 June 2007
/// Edited by Yaaseen Martin, 30 August 2019

#include <cmath>
#include <stdio.h>
#include "rstarget.cuh"
#include "rsdebug.cuh"
#include "rsinterp.cuh"
#include "rsnoise.cuh"
#include "pugixml.hpp"

using namespace rs;

// Gets child text as a rsFloat; see GetChildText for usage description
rsFloat GetChildrsFloat(pugi::xml_node &parent, std::string childname)
{
    std::string data = parent.child_value(childname.c_str());

    // Parse the first rsFloat from the data
    rsFloat result;
    std::istringstream iss(data);
    iss >> result;
    return result;
}
/// RCSModel Implementation

// Destructor
RCSModel::~RCSModel()
{
}

/// RCSConst Implementation

// Return a constant RCS
rsFloat RCSConst::SampleModel()
{
    return 1.0;
}

/// Destructor
RCSConst::~RCSConst()
{
}

/// RCSChiSquare Implementation

// Constructor
RCSChiSquare::RCSChiSquare(rsFloat k)
{
    gen = new GammaGenerator(k);
}

// Destructor
RCSChiSquare::~RCSChiSquare()
{
    delete gen;
}

// Get an RCS based on the Swerling II model and the mean RCS
rsFloat RCSChiSquare::SampleModel()
{
    return gen->GetSample();
}

/// Target Implementation

// Default constructor for Target object
Target::Target(Platform *platform, const std::string &name):
    Object(platform, name)
{
    model = 0;
}

// Default destructor for Target object
Target::~Target()
{
    delete model;
}

// Get the target polarization matrix
PSMatrix Target::GetPolarization() const
{
    return psm;
}

// Set the target polarization matrix
void Target::SetPolarization(const PSMatrix &in)
{
    psm = in;
}

// Set the target fluctuation model
void Target::SetFluctuationModel(RCSModel *in)
{
    model = in;
}

/// IsoTarget Implementation

// Constructor
IsoTarget::IsoTarget(Platform *platform, const std::string &name, rsFloat rcs):
    Target(platform, name),
    rcs(rcs)
{
}

// Destructor
IsoTarget::~IsoTarget()
{
}

// Return the RCS at the given angle
rsFloat IsoTarget::GetRCS(SVec3 &inAngle, SVec3 &outAngle, rsFloat wavelength) const
{
    if (model)
        return rcs*model->SampleModel();
    else
        return rcs;
}

/// SEMTarget Implementation

// Constructor
SEMTarget::SEMTarget(Platform *platform, const std::string &name, rsFloat diameter):
    Target(platform, name),
    diameter(diameter)
{
}

// Destructor
SEMTarget::~SEMTarget()
{
}

// Calculate SEM RCS
rsFloat SEMTarget::NASA_SEM_Model(rsFloat wavelength) const
{
    // Set up variables
    rsFloat x = diameter/wavelength; // Set norm diameter
    rsFloat z;                       // Initialise norm RCS
    rsFloat optical_limit = sqrt(20/M_PI);  // z > 5
    rsFloat rayleigh_limit = pow((0.12/(9*pow(M_PI, 5))), (1.0/6.0));  // z < 0.03

    // Set up x and z tables
    // NOTE: z table is incomplete because values do not extend all the way to z = 5; this causes inaccuracies so we add a point at z = 5
    rsFloat x_table [] = {0.10997, 0.11685, 0.12444, 0.13302, 0.14256, 0.15256, 0.16220, 0.17138, 0.18039, 0.18982, 0.20014,  
        0.21237, 0.22902, 0.25574, 0.30537, 0.42028, 0.56287, 0.71108, 0.86714, 1.0529, 1.2790, 1.5661, 1.8975, sqrt(4*5/M_PI)};
    
    rsFloat z_table [] = {
        0.001220, 0.001735, 0.002468, 0.003511, 0.004993, 0.007102, 0.01010, 0.01437, 0.02044, 0.02907, 0.04135, 0.05881,
        0.08365, 0.1190, 0.1692, 0.2407, 0.3424, 0.4870, 0.6927, 0.9852, 1.401, 1.993, 2.835, 5.0};

    // Check regime
    if (x > optical_limit) {              // Optical regime
        z = M_PI*pow(x, 2)/4;
    }
    else if (x < rayleigh_limit) {        // Rayleigh regime
        z = 9*pow(x, 6)*pow(M_PI, 5)/4;
    }
    else {                              // Mie resonance regime
        
        // Find x between two values in the x table
        int i1, i2 = 0;
        for (int i = 0; i < (int(sizeof(x_table)/sizeof(x_table[0])) - 1); i++){    // Use -1 due to +1 in for loop
            if (x_table[i] <= x && x_table[i + 1] >= x){                            // Find the two table entries within which x resides
                i1 = i;
                i2 = i + 1;
                break;  // Break out of loop once the indices are found
            }
        }

        // Interpolate z using x and the tables
        z = ((x - x_table[i1])/(x_table[i2] - x_table[i1]))*(z_table[i2] - z_table[i1]) + z_table[i1];
    }

    // Get final SEM RCS
    rsFloat rcs = z*pow(wavelength, 2);
    return rcs;
}

// Return the RCS (angles not used)
rsFloat SEMTarget::GetRCS(SVec3 &inAngle, SVec3 &outAngle, rsFloat wavelength) const
{
    // Apply SEM model
    rsFloat RCS = NASA_SEM_Model(wavelength);

    // Apply probabilistic model (if any)
    if (model)
        return RCS*model->SampleModel();
    else
        return RCS;
}

/// FileTarget Implementation

// Constructor
FileTarget::FileTarget(Platform *platform, const std::string &name, const std::string &filename):
    Target(platform, name)
{
    //Create the objects for azimuth and elevation interpolation
    azi_samples = new InterpSet();
    elev_samples = new InterpSet();

    //Load the data from the description file into the interpolation objects
    LoadRCSDescription(filename);
}

// Destructor
FileTarget::~FileTarget()
{
    delete azi_samples;
    delete elev_samples;
}

// Return the RCS at the given angle
rsFloat FileTarget::GetRCS(SVec3 &inAngle, SVec3 &outAngle, rsFloat wavelength) const
{
    // Currently uses a half angle approximation, this needs to be improved
    SVec3 t_angle = inAngle + outAngle;
    rsFloat RCS = std::sqrt(azi_samples->Value(t_angle.azimuth/2.0)*elev_samples->Value(t_angle.elevation/2.0));
    if (model)
        return RCS*model->SampleModel();
    else
        return RCS;
}

namespace {

    // Load samples of gain along an axis (not a member of FileAntenna)
    void LoadTargetGainAxis(InterpSet *set, pugi::xml_node &axisXML)
    {
        rsFloat angle;
        rsFloat gain;

        // Step through the XML file and load all the gain samples
        for (pugi::xml_node tmp = axisXML.child("rcssample"); tmp; tmp = tmp.next_sibling("rcssample")) {

            // Load the angle of the gain sample
            angle = GetChildrsFloat(tmp, "angle");

            // Load the gain of the gain sample
            gain = GetChildrsFloat(tmp, "rcs");

            // Load the values into the interpolation table
            set->InsertSample(angle, gain);
        }
    }

} // End of Anon. namespace

/// Load data from the RCS description file
void FileTarget::LoadRCSDescription(const std::string& filename)
{
    pugi::xml_document doc;
    pugi::xml_parse_result document = doc.load_file(filename.c_str());
    if (!document)
        throw std::runtime_error("ERROR: Could not load target description from "+filename);

    // Get the XML root node
    pugi::xml_node root = doc.child("filetarget");

    // Load the gain samples along the elevation axis
    pugi::xml_node tmp = root.child("elevation");
    if (!tmp)
        throw std::runtime_error("ERROR: Malformed XML in target description: No elevation pattern definition");
    LoadTargetGainAxis(elev_samples, tmp);

    // Load the gain samples along the azimuth axis
    tmp = root.child("azimuth");
    if (!tmp)
        throw std::runtime_error("ERROR: Malformed XML in target description: No azimuth pattern definition");
    LoadTargetGainAxis(azi_samples, tmp);
}

/// Functions for creating objects of various target types

// Create an isometric radiator target
Target* rs::CreateIsoTarget(Platform *platform, const std::string &name, rsFloat rcs)
{
    return new IsoTarget(platform, name, rcs);
}

// Create a NASA SEM target
Target* rs::CreateSEMTarget(Platform *platform, const std::string &name, rsFloat diameter)
{
    return new SEMTarget(platform, name, diameter);
}

// Create a target, loading the RCS pattern from a file
Target* rs::CreateFileTarget(Platform *platform, const std::string &name, const std::string &filename)
{
    return new FileTarget(platform, name, filename);
}
