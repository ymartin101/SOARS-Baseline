/// rsmultipath.h - Classes and definitions for multipath propagation
/// Marc Brooker, 25 April 2006
/// Edited by Yaaseen Martin, 02 September 2019

#ifndef __RS_MULTIPATH_H
#define __RS_MULTIPATH_H

#include "rsgeometry.cuh"

namespace rs {

    /// Class which manages multipath propagation
    class MultipathSurface {
        public:
            // These parameters define a plane of the form a*x+b*y+c*z+d=0
            // Constructor
            MultipathSurface(rsFloat a, rsFloat b, rsFloat c, rsFloat d, rsFloat factor);

            // Default destructor
            ~MultipathSurface();

            // Reflect a point in the surface
            Vec3 ReflectPoint(const Vec3 &r) const;

            // Get the reflectance factor
            rsFloat GetFactor() const;

        private:
            rsFloat factor; //!< How energy is reflected from the plane
            Matrix3 reflection; //!< Matrix defining reflection in this plane
            rsFloat norm_factor; //!< Factor to normalize length
            Vec3 translation_vector; //!< Vector to translate position of reflected point
    };
}

#endif
