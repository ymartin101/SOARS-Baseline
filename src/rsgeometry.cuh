/// rsgeometry.cuh - Classes for calculating the geometry of the world
/// Marc Brooker, 26 May 2006
/// Edited by Yaaseen Martin, 27 August 2019

#ifndef __RSGEOMETRY_H
#define __RSGEOMETRY_H

#include <config.h>

namespace rsGeometry
{
}

namespace rs {

    class SVec3;

    /// 3x3 Matrix, with simple operations
    class Matrix3 {
        public:
            rsFloat elements[9];

            /// Default constructor
            __host__ __device__ Matrix3();

            /// Default destructor
            __host__ __device__ ~Matrix3();

            /// Get a pointer to the element array
            __host__ __device__ const rsFloat *GetData() const;
            __host__ __device__ rsFloat *GetData();
    };

    /// The Vec3 class is a rectangular 3 vector
    class Vec3
    {
        public:
            rsFloat x, y, z;
            // Default constructor
            __host__ __device__ Vec3();

            // Constructor which sets co-ordinates
            __host__ __device__ Vec3(rsFloat x, rsFloat y, rsFloat z);

            // Constructor with a spherical vector
            __host__ __device__ Vec3(const SVec3 &svec);

            // Default destructor
            __host__ __device__ ~Vec3();

            // Vector operations
            __host__ __device__ Vec3 &operator+=(const Vec3 &b); // Vector Addition
            __host__ __device__ Vec3 &operator-=(const Vec3 &b); // Vector Subtraction
            __host__ __device__ Vec3 &operator*=(const Vec3 &b); // Componentwise multiplication
            __host__ __device__ Vec3 &operator=(const Vec3 &b); // Assignment operator

            // Matrix operations
            __host__ __device__ Vec3 &operator*=(const Matrix3 &m); // Multiplication by a 3x3 matrix

            // Scalar operations
            __host__ __device__ Vec3 &operator*=(const rsFloat b); // Multiplication by a scalar
            __host__ __device__ Vec3 &operator/=(const rsFloat b); // Division by a scalar
            __host__ __device__ Vec3 &operator+=(const rsFloat b); // Addition of a scalar

            // Return the length of the vector
            __host__ __device__ rsFloat Length() const;
    };

    // Vector operations
    __host__ __device__ rsFloat DotProduct(const Vec3 &a, const Vec3 &b); // Vector Inner product
    __host__ __device__ Vec3 CrossProduct(const Vec3 &a, const Vec3 &b); // Vector Cross product

    // Componentwise vector operations
    __host__ __device__ Vec3 operator*(const Vec3 &a, const Vec3 &b); // Componentwise product
    __host__ __device__ Vec3 operator+(const Vec3 &a, const Vec3 &b); // Componentwise add
    __host__ __device__ Vec3 operator-(const Vec3 &a, const Vec3 &b); // Componentwise subtract
    __host__ __device__ Vec3 operator/(const Vec3 &a, const Vec3 &b); // Componentwise divide

    // Scalar operations
    __host__ __device__ Vec3 operator*(const Vec3 &a, const rsFloat b); // Multiply by a scalar
    __host__ __device__ Vec3 operator/(const Vec3 &a, const rsFloat b); // Division by a scalar
    __host__ __device__ Vec3 operator/(const rsFloat a, const Vec3 &b); // Division by a scalar

    /// The SVec3 class is a vector in R^3, stored in spherical co-ordinates
    class SVec3
    {
        public:
            rsFloat length; // The length of the vector
            rsFloat azimuth; // Angle in the x-y plane (radians)
            rsFloat elevation; // Elevation angle above x-y plane (radians)

            // Default constructor
            __host__ __device__ SVec3();

            // Constructor with initializers
            __host__ __device__ SVec3(rsFloat length, rsFloat azimuth, rsFloat elevation);

            // Copy constructor
            __host__ __device__ SVec3(const SVec3 &svec);

            // Constructor with a rectangular vector
            __host__ __device__ SVec3(const Vec3 &vec);

            // Destructor
            __host__ __device__ ~SVec3();

            __host__ __device__ SVec3 &operator*=(rsFloat b); // Multiplication by a scalar
            __host__ __device__ SVec3 &operator/=(rsFloat b); // Division by a scalar
    };

}
#endif
