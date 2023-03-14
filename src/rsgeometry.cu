/// rsgeometry.cpp - Implementation of geometry classes
/// Marc Brooker, 26 May 2006
/// Edited by Yaaseen Martin, 02 September 2019

#include <cmath>
#include "rsgeometry.cuh"
#include "rsdebug.cuh"

using namespace rs;

/// Vec3 Implementation

// Default constructor
__host__ __device__ Vec3::Vec3():
    x(0), y(0), z(0)
{
}

// Data initialization constructor
__host__ __device__ Vec3::Vec3(rsFloat x, rsFloat y, rsFloat z):
    x(x), y(y), z(z)
{
}

// Default destructor
__host__ __device__ Vec3::~Vec3()
{
}

// Constructor with a spherical vector
__host__ __device__ Vec3::Vec3(const SVec3 &svec)
{
    x = svec.length*std::cos(svec.azimuth)*std::cos(svec.elevation);
    y = svec.length*std::sin(svec.azimuth)*std::cos(svec.elevation);
    z = svec.length*std::sin(svec.elevation);
}

// Addition operator
__host__ __device__ Vec3 &Vec3::operator+=(const Vec3 &b)
{
    x += b.x;
    y += b.y;
    z += b.z;
    return *this;
}

// Subtraction operator
__host__ __device__ Vec3 &Vec3::operator-=(const Vec3 &b)
{
    x -= b.x;
    y -= b.y;
    z -= b.z;
    return *this;
}

// Componentwise multiplication; used by interpolation code
__host__ __device__ Vec3& Vec3::operator*=(const Vec3 &b)
{
    x *= b.x;
    y *= b.y;
    z *= b.z;
    return *this;
}

// Assignment operator
__host__ __device__ Vec3& Vec3::operator=(const Vec3 &b)
{
    x = b.x;
    y = b.y;
    z = b.z;
    return *this;
}

// Multiplication by a 3x3 matrix
__host__ __device__ Vec3& Vec3::operator*=(const Matrix3 &m)
{
    const rsFloat *mat = m.GetData();
    Vec3 v(x, y, z);
    x = mat[0]*v.x + mat[1]*v.y + mat[2]*v.z;
    y = mat[3]*v.x + mat[4]*v.y + mat[5]*v.z;
    z = mat[6]*v.x + mat[7]*v.y + mat[8]*v.z;
    return *this;
}

// Length division operator
__host__ __device__ Vec3 &Vec3::operator/=(const rsFloat b)
{
    x /= b;
    y /= b;
    z /= b;
    return *this;
}

// Scaling operator
__host__ __device__ Vec3 &Vec3::operator*=(const rsFloat b)
{
    x *= b;
    y *= b;
    z *= b;
    return *this;
}

// Addition of a scalar
__host__ __device__ Vec3& Vec3::operator+=(const rsFloat b)
{
    x += b;
    y += b;
    z += b;
    return *this;
}

// Return the length of the vector
__host__ __device__ rsFloat Vec3::Length() const
{
    return sqrt(x*x+y*y+z*z);
}

/// Vec3 Non-member functions

// Inner (dot) product operator
__host__ __device__ rsFloat rs::DotProduct(const Vec3 &a, const Vec3 &b) {
    return a.x*b.x+a.y*b.y+a.z*b.z;
}

// Returns a new vector containing the cross (outer) product of two vectors
__host__ __device__ Vec3 rs::CrossProduct(const Vec3 &a, const Vec3 &b) {
    Vec3 newvec;
    newvec.x = a.y*b.z-a.z*b.y;
    newvec.y = a.z*b.x-a.x*b.z;
    newvec.z = a.x*b.y-a.y*b.x;
    return newvec;
}

// Component-wise multiplication
__host__ __device__ Vec3 rs::operator*(const Vec3 &a, const Vec3 &b) {
    Vec3 c(a);
    c *= b;
    return c;
}

// Componentwise addition
__host__ __device__ Vec3 rs::operator+(const Vec3 &a, const Vec3 &b) {
    Vec3 c(a);
    c += b;
    return c;
}

// Componentwise subtraction
__host__ __device__ Vec3 rs::operator-(const Vec3 &a, const Vec3 &b) {
    Vec3 c(a);
    c -= b;
    return c;
}

// Componentwise division
__host__ __device__ Vec3 rs::operator/(const Vec3 &a, const Vec3 &b) {
    Vec3 c(a);
    c.x /= b.x;
    c.y /= b.y;
    c.z /= b.z;
    return c;
}

// Multiplication by a scalar
__host__ __device__ Vec3 rs::operator*(const Vec3 &a, const rsFloat b) {
    Vec3 c(a);
    c *= b;
    return c;
}

// Division by a scalar
__host__ __device__ Vec3 rs::operator/(const Vec3 &a, const rsFloat b) {
    Vec3 c(a);
    c /= b;
    return c;
}

// Division by a scalar
__host__ __device__ Vec3 rs::operator/(const rsFloat a, const Vec3 &b) {
    Vec3 c(b);
    c.x = a/c.x;
    c.y = a/c.y;
    c.z = a/c.z;
    return c;
}

/// SVec3 implementation

__host__ __device__ SVec3::SVec3():
    length(0), azimuth(0), elevation(0)
{
}

// Constructor with initialization
__host__ __device__ SVec3::SVec3(rsFloat length, rsFloat azimuth, rsFloat elevation):
    length(length), azimuth(azimuth), elevation(elevation)
{
}

// Copy constructor
__host__ __device__ SVec3::SVec3(const SVec3 &svec):
    length(svec.length), azimuth(svec.azimuth), elevation(svec.elevation)
{
}

// Destructor
__host__ __device__ SVec3::~SVec3()
{
}

// Constructor with a rectangular vector
__host__ __device__ SVec3::SVec3(const Vec3 &vec)
{
    length = vec.Length();
    if (length != 0) {
        elevation = std::asin(vec.z/length);
        azimuth = std::atan2(vec.y, vec.x);
    }
    else {
        elevation = 0;
        azimuth = 0;
    }
}

// Multiplication by a scalar
__host__ __device__ SVec3& SVec3::operator*=(rsFloat b)
{
    length *= b;
    return *this;
}

// Division by a scalar
__host__ __device__ SVec3 &SVec3::operator/=(rsFloat b)
{
    length /= b;
    return *this;
}

/// Matrix3 Implementation

// Default constructor
__host__ __device__ Matrix3::Matrix3()
{
    for (int i = 0; i < 9; i++)
        elements[i] = 0;
}

// Default destructor
__host__ __device__ Matrix3::~Matrix3()
{
}

/// Get a pointer to the element array
__host__ __device__ const rsFloat* Matrix3::GetData() const
{
    return elements;
}

__host__ __device__ rsFloat* Matrix3::GetData()
{
    return elements;
}
