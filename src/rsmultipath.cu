/// rsmultipath.cpp - Implementation of multipath propagation
/// Marc Brooker, 9 September 2007
/// Edited by Yaaseen Martin, 02 September 2019

#include "rsmultipath.cuh"
#include "rsdebug.cuh"

using namespace rs;

/// MultipathSurface implementation

// Constructor
MultipathSurface::MultipathSurface(rsFloat a, rsFloat b, rsFloat c, rsFloat d, rsFloat factor):
    // Member initialiser sets factor to input variable "factor"
    factor(factor)
{
    // a, b, c, and d describe the coefficients of the plane of reflection (ax + by + cz = d)
    rsFloat *mat = reflection.GetData();

    // Fill the reflection matrix
    // Used to calculate the point R' as a reflection of the point R through the plane
    mat[0] = -a*a+b*b+c*c;
    mat[4] = a*a-b*b+c*c;
    mat[8] = a*a+b*b-c*c;
    mat[1] = mat[3] = -2*a*b;
    mat[2] = mat[6] = -2*a*c;
    mat[5] = mat[7] = -2*b*c;

    // The above matrix only works if (a^2 + b^2 + c^2 = 1), so normalise
    norm_factor = 1/(a*a+b*b+c*c);

    // Get the translation vector
    translation_vector = Vec3(2*a*d, 2*b*d, 2*c*d);
}

// Default destructor
MultipathSurface::~MultipathSurface()
{
}

// Return the point R' as a reflection of the point R through the reflecting surface
Vec3 MultipathSurface::ReflectPoint(const Vec3 &r) const
{
    Vec3 rdash = r;

    // Calculate the reflected position of b using norm_factor*(reflection*b - translation_vector)
    rdash *= reflection;
    rdash -= translation_vector;
    rdash *= norm_factor;
    return rdash;
}

// Return the reflection factor
rsFloat MultipathSurface::GetFactor() const
{
    return factor;
}
