#pragma once

#include <cmath>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <string>

#if 0
#include "./Euler.h"
#include "./Cylindrical.h"
#include "./Spherical.h"
#include "./Matrix3.h"
#endif
#include "./Matrix4.h"
#include "./Quaternion.h"

#if 0
#include "../cameras/Camera.h"
#include "../core/BufferAttribute.h"
#endif

namespace three {

template <typename T>
class Vector3 {
public:
    constexpr Vector3(T x = 0, T y = 0, T z = 0)
        : x(x)
        , y(y)
        , z(z)
    {}
    constexpr Vector3(const Vector3<T> &copy)
        : x(copy.x)
        , y(copy.y)
        , z(copy.z)
    {}
    constexpr Vector3 &operator =(const Vector3<T> &rhs) {
        this->x = rhs.x;
        this->y = rhs.y;
        this->z = rhs.z;
        return *this;
    }
    // use default
    constexpr Vector3(Vector3<T> &&move) = default;
    constexpr Vector3 &operator =(Vector3<T> &&rhs) = default;

    Vector3<T> &set(T x, T y, T z) {
        this->x = x;
        this->y = y;
        this->z = z;
        return *this;
    }
    Vector3<T> &set(T x, T y) {
        this->x = x;
        this->y = y;
        return *this;
    }

    Vector3<T> &setScalar(T scalar) {
        this->x = scalar;
        this->y = scalar;
        this->z = scalar;
        return *this;
    }

    Vector3<T> &setX(T x) {
        this->x = x;
        return *this;
    }

    Vector3<T> &setY(T y) {
        this->y = y;
        return *this;
    }

    Vector3<T> &setZ(T z) {
        this->z = z;
        return *this;
    }

    Vector3<T> &setComponent(int index, T value) {
        switch (index) {
            case 0: this->x = value; break;
            case 1: this->y = value; break;
            case 2: this->z = value; break;
            default: throw std::runtime_error("index is out of range: " + std::to_string(index));
        }
        return *this;
    }

    T getComponent(int index) const {
        switch (index) {
            case 0: return this->x;
            case 1: return this->y;
            case 2: return this->z;
            default: throw std::runtime_error("index is out of range: " + std::to_string(index));
        }
    }

    Vector3<T> clone() const {
        return Vector3<T>(this->x, this->y, this->z);
    }

    Vector3<T> &copy(const Vector3<T> &v) {
        this->x = v.x;
        this->y = v.y;
        this->z = v.z;
        return *this;
    }

    Vector3<T> &add(const Vector3<T> &v) {
        this->x += v.x;
        this->y += v.y;
        this->z += v.z;
        return *this;
    }

    Vector3<T> &addScalar(T s) {
        this->x += s;
        this->y += s;
        this->z += s;
        return *this;
    }

    Vector3<T> &addVectors(const Vector3<T> &a, const Vector3<T> &b) {
        this->x = a.x + b.x;
        this->y = a.y + b.y;
        this->z = a.z + b.z;
        return *this;
    }

    Vector3<T> &addScaledVector(const Vector3<T> &v, T s) {
        this->x += v.x * s;
        this->y += v.y * s;
        this->z += v.z * s;
        return *this;
    }

    Vector3<T> &sub(const Vector3<T> &v) {
        this->x -= v.x;
        this->y -= v.y;
        this->z -= v.z;
        return *this;
    }

    Vector3<T> &subScalar(T s) {
        this->x -= s;
        this->y -= s;
        this->z -= s;
        return *this;
    }

    Vector3<T> &subVectors(const Vector3<T> &a, const Vector3<T> &b) {
        this->x = a.x - b.x;
        this->y = a.y - b.y;
        this->z = a.z - b.z;
        return *this;
    }

    Vector3<T> &multiply(const Vector3<T> &v) {
        this->x *= v.x;
        this->y *= v.y;
        this->z *= v.z;
        return *this;
    }

    Vector3<T> &multiplyScalar(T scalar) {
        this->x *= scalar;
        this->y *= scalar;
        this->z *= scalar;
        return *this;
    }

    Vector3<T> &multiplyVectors(const Vector3<T> &a, const Vector3<T> &b) {
        this->x = a.x * b.x;
        this->y = a.y * b.y;
        this->z = a.z * b.z;
        return *this;
    }

#if 0
    Vector3<T> &applyEuler(const Euler<T> &euler) {
        return this->applyQuaternion(Quaternion<T>().setFromEuler(euler));
    }
#endif

    Vector3<T> &applyAxisAngle(const Vector3<T> &axis, T angle) {
        return this->applyQuaternion(Quaternion<T>().setFromAxisAngle(axis, angle));
    }

#if 0
    Vector3<T> &applyMatrix3(const Matrix3<T> &m) {
        const auto x = this->x, y = this->y, z = this->z;
        const auto &e = m.elements;

        this->x = e[ 0 ] * x + e[ 3 ] * y + e[ 6 ] * z;
        this->y = e[ 1 ] * x + e[ 4 ] * y + e[ 7 ] * z;
        this->z = e[ 2 ] * x + e[ 5 ] * y + e[ 8 ] * z;

        return *this;
    }

    Vector3<T> &applyNormalMatrix(const Matrix3<T> &m) {
        return this->applyMatrix3(m).normalize();
    }
#endif

    Vector3<T> &applyMatrix4(const Matrix4<T> &m) {
        const auto x = this->x, y = this->y, z = this->z;
        const auto &e = m.elements;

        const auto w = static_cast<T>(1) / ( e[ 3 ] * x + e[ 7 ] * y + e[ 11 ] * z + e[ 15 ]);

        this->x = (e[ 0 ] * x + e[ 4 ] * y + e[ 8 ] * z + e[ 12 ]) * w;
        this->y = (e[ 1 ] * x + e[ 5 ] * y + e[ 9 ] * z + e[ 13 ]) * w;
        this->z = (e[ 2 ] * x + e[ 6 ] * y + e[ 10 ] * z + e[ 14 ]) * w;

        return *this;
    }

    Vector3<T> &applyQuaternion(const Quaternion<T> &q) {
        const auto x = this->x, y = this->y, z = this->z;
        const auto qx = q.x(), qy = q.y(), qz = q.z(), qw = q.w();

        // calculate quat * vector

        const auto ix = qw * x + qy * z - qz * y;
        const auto iy = qw * y + qz * x - qx * z;
        const auto iz = qw * z + qx * y - qy * x;
        const auto iw = - qx * x - qy * y - qz * z;

        // calculate result * inverse quat

        this->x = ix * qw + iw * - qx + iy * - qz - iz * - qy;
        this->y = iy * qw + iw * - qy + iz * - qx - ix * - qz;
        this->z = iz * qw + iw * - qz + ix * - qy - iy * - qx;

        return *this;
    }

#if 0
    Vector3<T> &project(const Camera &camera) {
        return this->applyMatrix4(camera.matrixWorldInverse).applyMatrix4(camera.projectionMatrix);
    }

    Vector3<T> &unproject(const Camera &camera) {
        return this->applyMatrix4(camera.projectionMatrixInverse).applyMatrix4(camera.matrixWorld);
    }
#endif

    Vector3<T> &transformDirection(const Matrix4<T> &m) {
        const auto x = this->x, y = this->y, z = this->z;
        const auto &e = m.elements;

        this->x = e[ 0 ] * x + e[ 4 ] * y + e[ 8 ] * z;
        this->y = e[ 1 ] * x + e[ 5 ] * y + e[ 9 ] * z;
        this->z = e[ 2 ] * x + e[ 6 ] * y + e[ 10 ] * z;

        return this->normalize();
    }

    Vector3<T> &divide(const Vector3<T> &v) {
        this->x /= v.x;
        this->y /= v.y;
        this->z /= v.z;

        return *this;
    }

    Vector3<T> &divideScalar(T scalar) {
        return this->multiplyScalar(static_cast<T>(1) / scalar);
    }

    Vector3<T> &min(const Vector3<T> &v) {
        this->x = std::min(this->x, v.x);
        this->y = std::min(this->y, v.y);
        this->z = std::min(this->z, v.z);
        return *this;
    }

    Vector3<T> &max(const Vector3<T> &v) {
        this->x = std::max(this->x, v.x);
        this->y = std::max(this->y, v.y);
        this->z = std::max(this->z, v.z);
        return *this;
    }

    Vector3<T> &clamp(const Vector3<T> &min, const Vector3<T> &max) {
        this->x = std::max(min.x, std::min(max.x, this->x));
        this->y = std::max(min.y, std::min(max.y, this->y));
        this->z = std::max(min.z, std::min(max.z, this->z));
        return *this;
    }

    Vector3<T> &clampScalar(T minVal, T maxVal) {
        this->x = std::max(minVal, std::min(maxVal, this->x));
        this->y = std::max(minVal, std::min(maxVal, this->y));
        this->z = std::max(minVal, std::min(maxVal, this->z));
        return *this;
    }

    Vector3<T> &clampLength(T min, T max) {
        const auto length = this->length();
        return this->divideScalar(length == 0 ? 1 : length).multiplyScalar(std::max(min, std::min(max, length)));
    }

    Vector3<T> &floor() {
        this->x = std::floor(this->x);
        this->y = std::floor(this->y);
        this->z = std::floor(this->z);
        return *this;
    }

    Vector3<T> &ceil() {
        this->x = std::ceil(this->x);
        this->y = std::ceil(this->y);
        this->z = std::ceil(this->z);
        return *this;
    }

    Vector3<T> &round() {
        this->x = std::round(this->x);
        this->y = std::round(this->y);
        this->z = std::round(this->z);
        return *this;
    }

    Vector3<T> &roundToZero() {
        this->x = (this->x < 0) ? std::ceil(this->x) : std::floor(this->x);
        this->y = (this->y < 0) ? std::ceil(this->y) : std::floor(this->y);
        this->z = (this->z < 0) ? std::ceil(this->z) : std::floor(this->z);
        return *this;
    }

    Vector3<T> &negate() {
        this->x = - this->x;
        this->y = - this->y;
        this->z = - this->z;
        return *this;
    }

    T dot(const Vector3<T> &v) const {
        return this->x * v.x + this->y * v.y + this->z * v.z;
    }

    // TODO lengthSquared?
    T lengthSq() const {
        return this->x * this->x + this->y * this->y + this->z * this->z;
    }

    T length() const {
        return std::sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
    }

    T manhattanLength() const {
        return std::abs(this->x) + std::abs(this->y) + std::abs(this->z);
    }

    Vector3<T> &normalize() {
        const auto length = this->length();
        return this->divideScalar(length == 0 ? 1 : length);
    }

    Vector3<T> &setLength(T length) {
        return this->normalize().multiplyScalar(length);
    }

    Vector3<T> &lerp(const Vector3<T> &v, T alpha) {
        this->x += (v.x - this->x) * alpha;
        this->y += (v.y - this->y) * alpha;
        this->z += (v.z - this->z) * alpha;
        return *this;
    }

    Vector3<T> &lerpVectors(const Vector3<T> &v1, const Vector3<T> &v2, T alpha) {
        this->x = v1.x + (v2.x - v1.x) * alpha;
        this->y = v1.y + (v2.y - v1.y) * alpha;
        this->z = v1.z + (v2.z - v1.z) * alpha;
        return *this;
    }

    Vector3<T> &cross(const Vector3<T> &v) {
        return this->crossVectors(*this, v);
    }

    Vector3<T> &crossVectors(const Vector3<T> &a, const Vector3<T> &b) {
        const auto ax = a.x, ay = a.y, az = a.z;
        const auto bx = b.x, by = b.y, bz = b.z;

        this->x = ay * bz - az * by;
        this->y = az * bx - ax * bz;
        this->z = ax * by - ay * bx;

        return *this;
    }

    Vector3<T> &projectOnVector(const Vector3<T> &v) {
        const auto denominator = v.lengthSq();

        if ( denominator == 0 ) { return this->set(0, 0, 0); }

        const auto scalar = v.dot(*this) / denominator;

        return this->copy(v).multiplyScalar(scalar);
    }

    Vector3<T> &projectOnPlane(const Vector3<T> &planeNormal) {
        auto _vector = this->clone().projectOnVector(planeNormal);

        return this->sub(_vector);
    }

    Vector3<T> &reflect(const Vector3<T> &normal) {
        // reflect incident vector off plane orthogonal to normal
        // normal is assumed to have unit length
        return this->sub(normal.clone().multiplyScalar(2 * this->dot(normal)));
    }

    T angleTo(const Vector3<T> &v) const {
        const auto denominator = std::sqrt(this->lengthSq() * v.lengthSq());

        if ( denominator == 0 ) { return std::numbers::pi_v<float> / 2; }

        const auto theta = this->dot(v) / denominator;

        // clamp, to handle numerical problems

        return std::acos(std::clamp<T>(theta, -1, 1));
    }

    T distanceTo(const Vector3<T> &v) const {
        return std::sqrt(this->distanceToSquared(v));
    }

    T distanceToSquared(const Vector3<T> &v) const {
        const auto dx = this->x - v.x, dy = this->y - v.y, dz = this->z - v.z;

        return dx * dx + dy * dy + dz * dz;
    }

    T manhattanDistanceTo(const Vector3<T> &v) const {
        return std::abs(this->x - v.x) + std::abs(this->y - v.y) + std::abs(this->z - v.z);
    }

#if 0
    Vector3<T> &setFromSpherical(const Spherical<T> &s) {
        return this->setFromSphericalCoords(s.radius, s.phi, s.theta);
    }
#endif

    Vector3<T> &setFromSphericalCoords(T radius, T phi, T theta) {
        const auto sinPhiRadius = std::sin(phi) * radius;

        this->x = sinPhiRadius * std::sin(theta);
        this->y = std::cos(phi) * radius;
        this->z = sinPhiRadius * std::cos(theta);

        return *this;
    }

#if 0
    Vector3<T> &setFromCylindrical(const Cylindrical<T> &c) {
        return this->setFromCylindricalCoords(c.radius, c.theta, c.y);
    }
#endif

    Vector3<T> &setFromCylindricalCoords(T radius, T theta, T y) {
        this->x = radius * std::sin(theta);
        this->y = y;
        this->z = radius * std::cos(theta);

        return *this;
    }

    Vector3<T> &setFromMatrixPosition(const Matrix4<T> &m) {
        const auto &e = m.elements;

        this->x = e[ 12 ];
        this->y = e[ 13 ];
        this->z = e[ 14 ];

        return *this;
    }

    Vector3<T> &setFromMatrixScale(const Matrix4<T> &m) {
        const auto sx = this->setFromMatrixColumn(m, 0).length();
        const auto sy = this->setFromMatrixColumn(m, 1).length();
        const auto sz = this->setFromMatrixColumn(m, 2).length();

        this->x = sx;
        this->y = sy;
        this->z = sz;

        return *this;
    }

    Vector3<T> &setFromMatrixColumn(const Matrix4<T> &m, int index) {
        return this->fromArray(m.elements, index * 4);
    }

#if 0
    Vector3<T> &setFromMatrix3Column(const Matrix3<T> &m, int index) {
        return this->fromArray(m.elements, index * 3);
    }
#endif

    bool equals(const Vector3<T> &v) const {
        return ((v.x == this->x) && (v.y == this->y) && (v.z == this->z));
    }

    template <typename V>
    Vector3<T> &fromArray(const V &array, int offset = 0) {
        auto it = array.begin() + offset;
        this->x = *it++;
        this->y = *it++;
        this->z = *it;
        return *this;
    }

    template <typename V>
    V &toArray(V &array, int offset = 0) const {
        assert(array.size() >= static_cast<std::size_t>(offset + 3));
        array[ offset ] = this->x;
        array[ offset + 1 ] = this->y;
        array[ offset + 2 ] = this->z;
        return array;
    }
    template <typename V>
    V toArray() const {
        V array;
        array.reserve(3);
        array.push_back(this->x);
        array.push_back(this->y);
        array.push_back(this->z);
        return array;
    }

#if 0
    Vector3<T> &fromBufferAttribute(const BufferAttribute<T> &attribute, int index) {
        this->x = attribute.getX(index);
        this->y = attribute.getY(index);
        this->z = attribute.getZ(index);
        return *this;
    }
#endif

    template <typename G>
    Vector3<T> &random(const G &generator) {
        this->x = generator();
        this->y = generator();
        this->z = generator();
        return *this;
    }

public:
    T x, y, z;

};

}  // namespace three
