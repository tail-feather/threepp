#pragma once

#include <limits>

#include "./Matrix4.h"

namespace three {

template <typename T>
class Quaternion {
public:
    Quaternion(T x = 0, T y = 0, T z = 0, T w = 1)
        : _x(x)
        , _y(y)
        , _z(z)
        , _w(w)
        , _onChangeCallback([](){})
    {}

    static void slerpFlat(T *dst, int dstOffset, const T *src0, int srcOffset0, const T *src1, int srcOffset1, T t) {
        // fuzz free, array-based Quaternion SLERP operation

        auto x0 = src0[ srcOffset0 + 0 ],
            y0 = src0[ srcOffset0 + 1 ],
            z0 = src0[ srcOffset0 + 2 ],
            w0 = src0[ srcOffset0 + 3 ];

        auto x1 = src1[ srcOffset1+ 0 ],
            y1 = src1[ srcOffset1 + 1 ],
            z1 = src1[ srcOffset1 + 2 ],
            w1 = src1[ srcOffset1 + 3 ];

        if ( t == 0 ) {
            dst[ dstOffset + 0 ] = x0;
            dst[ dstOffset + 1 ] = y0;
            dst[ dstOffset + 2 ] = z0;
            dst[ dstOffset + 3 ] = w0;
            return;
        }

        if ( t == 1 ) {
            dst[ dstOffset + 0 ] = x1;
            dst[ dstOffset + 1 ] = y1;
            dst[ dstOffset + 2 ] = z1;
            dst[ dstOffset + 3 ] = w1;
            return;
        }

        if ( w0 != w1 || x0 != x1 || y0 != y1 || z0 != z1 ) {
            auto s = 1 - t;
            const auto cos = x0 * x1 + y0 * y1 + z0 * z1 + w0 * w1,
                dir = ( cos >= 0 ? 1 : - 1 ),
                sqrSin = 1 - cos * cos;

            // Skip the Slerp for tiny steps to avoid numeric problems:
            if ( sqrSin > std::numeric_limits<T>::epsilon() ) {

                const auto sin = std::sqrt( sqrSin ),
                    len = std::atan2( sin, cos * dir );

                s = std::sin( s * len ) / sin;
                t = std::sin( t * len ) / sin;

            }

            const auto tDir = t * dir;

            x0 = x0 * s + x1 * tDir;
            y0 = y0 * s + y1 * tDir;
            z0 = z0 * s + z1 * tDir;
            w0 = w0 * s + w1 * tDir;

            // Normalize in case we just did a lerp:
            if ( s == 1 - t ) {

                const auto f = 1 / std::sqrt( x0 * x0 + y0 * y0 + z0 * z0 + w0 * w0 );

                x0 *= f;
                y0 *= f;
                z0 *= f;
                w0 *= f;

            }

        }

        dst[ dstOffset ] = x0;
        dst[ dstOffset + 1 ] = y0;
        dst[ dstOffset + 2 ] = z0;
        dst[ dstOffset + 3 ] = w0;

    }

    static T *multiplyQuaternionsFlat(T *dst, int dstOffset, const T *src0, int srcOffset0, T *src1, int srcOffset1 ) {

        const auto x0 = src0[ srcOffset0 ];
        const auto y0 = src0[ srcOffset0 + 1 ];
        const auto z0 = src0[ srcOffset0 + 2 ];
        const auto w0 = src0[ srcOffset0 + 3 ];

        const auto x1 = src1[ srcOffset1 ];
        const auto y1 = src1[ srcOffset1 + 1 ];
        const auto z1 = src1[ srcOffset1 + 2 ];
        const auto w1 = src1[ srcOffset1 + 3 ];

        dst[ dstOffset ] = x0 * w1 + w0 * x1 + y0 * z1 - z0 * y1;
        dst[ dstOffset + 1 ] = y0 * w1 + w0 * y1 + z0 * x1 - x0 * z1;
        dst[ dstOffset + 2 ] = z0 * w1 + w0 * z1 + x0 * y1 - y0 * x1;
        dst[ dstOffset + 3 ] = w0 * w1 - x0 * x1 - y0 * y1 - z0 * z1;

        return dst;

    }

    T x() const { return this->_x; }
    void setX(T x) {
        this->_x = x;
        this->_onChangeCallback();
    }

    T y() const { return this->_y; }
    void setY(T y) {
        this->_y = y;
        this->_onChangeCallback();
    }

    T z() const { return this->_z; }
    void setZ(T z) {
        this->_z = z;
        this->_onChangeCallback();
    }

    T w() const { return this->_w; }
    void setW(T w) {
        this->_w = w;
        this->_onChangeCallback();
    }

    Quaternion<T> &set( T x, T y, T z, T w ) {

        this->_x = x;
        this->_y = y;
        this->_z = z;
        this->_w = w;

        this->_onChangeCallback();

        return *this;

    }

    Quaternion<T> clone() const {

        return Quaternion<T>( this->_x, this->_y, this->_z, this->_w );

    }

    Quaternion<T> &copy( const Quaternion<T> &quaternion ) {

        this->_x = quaternion._x;
        this->_y = quaternion._y;
        this->_z = quaternion._z;
        this->_w = quaternion._w;

        this->_onChangeCallback();

        return *this;

    }

#if 0
    setFromEuler( euler, update ) {

        if ( ! ( euler && euler.isEuler ) ) {

            throw new Error( 'THREE.Quaternion: .setFromEuler() now expects an Euler rotation rather than a Vector3 and order.' );

        }

        const x = euler._x, y = euler._y, z = euler._z, order = euler._order;

        // http://www.mathworks.com/matlabcentral/fileexchange/
        //     20696-function-to-convert-between-dcm-euler-angles-quaternions-and-euler-vectors/
        //    content/SpinCalc.m

        const cos = Math.cos;
        const sin = Math.sin;

        const c1 = cos( x / 2 );
        const c2 = cos( y / 2 );
        const c3 = cos( z / 2 );

        const s1 = sin( x / 2 );
        const s2 = sin( y / 2 );
        const s3 = sin( z / 2 );

        switch ( order ) {

            case 'XYZ':
                this->_x = s1 * c2 * c3 + c1 * s2 * s3;
                this->_y = c1 * s2 * c3 - s1 * c2 * s3;
                this->_z = c1 * c2 * s3 + s1 * s2 * c3;
                this->_w = c1 * c2 * c3 - s1 * s2 * s3;
                break;

            case 'YXZ':
                this->_x = s1 * c2 * c3 + c1 * s2 * s3;
                this->_y = c1 * s2 * c3 - s1 * c2 * s3;
                this->_z = c1 * c2 * s3 - s1 * s2 * c3;
                this->_w = c1 * c2 * c3 + s1 * s2 * s3;
                break;

            case 'ZXY':
                this->_x = s1 * c2 * c3 - c1 * s2 * s3;
                this->_y = c1 * s2 * c3 + s1 * c2 * s3;
                this->_z = c1 * c2 * s3 + s1 * s2 * c3;
                this->_w = c1 * c2 * c3 - s1 * s2 * s3;
                break;

            case 'ZYX':
                this->_x = s1 * c2 * c3 - c1 * s2 * s3;
                this->_y = c1 * s2 * c3 + s1 * c2 * s3;
                this->_z = c1 * c2 * s3 - s1 * s2 * c3;
                this->_w = c1 * c2 * c3 + s1 * s2 * s3;
                break;

            case 'YZX':
                this->_x = s1 * c2 * c3 + c1 * s2 * s3;
                this->_y = c1 * s2 * c3 + s1 * c2 * s3;
                this->_z = c1 * c2 * s3 - s1 * s2 * c3;
                this->_w = c1 * c2 * c3 - s1 * s2 * s3;
                break;

            case 'XZY':
                this->_x = s1 * c2 * c3 - c1 * s2 * s3;
                this->_y = c1 * s2 * c3 - s1 * c2 * s3;
                this->_z = c1 * c2 * s3 + s1 * s2 * c3;
                this->_w = c1 * c2 * c3 + s1 * s2 * s3;
                break;

            default:
                console.warn( 'THREE.Quaternion: .setFromEuler() encountered an unknown order: ' + order );

        }

        if ( update !== false ) this->_onChangeCallback();

        return *this;

    }
#endif

    template <typename V>
    Quaternion<T> &setFromAxisAngle( const V &axis, T angle ) {

        // http://www.euclideanspace.com/maths/geometry/rotations/conversions/angleToQuaternion/index.htm

        // assumes axis is normalized

        const auto halfAngle = angle / 2, s = std::sin( halfAngle );

        this->_x = axis.x * s;
        this->_y = axis.y * s;
        this->_z = axis.z * s;
        this->_w = std::cos( halfAngle );

        this->_onChangeCallback();

        return *this;

    }

    Quaternion<T> &setFromRotationMatrix( const Matrix4<T> &m ) {

        // http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm

        // assumes the upper 3x3 of m is a pure rotation matrix (i.e, unscaled)

        const auto &te = m.elements;
        const auto
            m11 = te[ 0 ], m12 = te[ 4 ], m13 = te[ 8 ],
            m21 = te[ 1 ], m22 = te[ 5 ], m23 = te[ 9 ],
            m31 = te[ 2 ], m32 = te[ 6 ], m33 = te[ 10 ],

            trace = m11 + m22 + m33;

        if ( trace > 0 ) {

            const auto s = 0.5 / std::sqrt( trace + 1.0 );

            this->_w = 0.25 / s;
            this->_x = ( m32 - m23 ) * s;
            this->_y = ( m13 - m31 ) * s;
            this->_z = ( m21 - m12 ) * s;

        } else if ( m11 > m22 && m11 > m33 ) {

            const auto s = 2.0 * std::sqrt( 1.0 + m11 - m22 - m33 );

            this->_w = ( m32 - m23 ) / s;
            this->_x = 0.25 * s;
            this->_y = ( m12 + m21 ) / s;
            this->_z = ( m13 + m31 ) / s;

        } else if ( m22 > m33 ) {

            const auto s = 2.0 * std::sqrt( 1.0 + m22 - m11 - m33 );

            this->_w = ( m13 - m31 ) / s;
            this->_x = ( m12 + m21 ) / s;
            this->_y = 0.25 * s;
            this->_z = ( m23 + m32 ) / s;

        } else {

            const auto s = 2.0 * std::sqrt( 1.0 + m33 - m11 - m22 );

            this->_w = ( m21 - m12 ) / s;
            this->_x = ( m13 + m31 ) / s;
            this->_y = ( m23 + m32 ) / s;
            this->_z = 0.25 * s;

        }

        this->_onChangeCallback();

        return *this;

    }

    template <typename V>
    Quaternion<T> &setFromUnitVectors( const V &vFrom, const V &vTo ) {

        // assumes direction vectors vFrom and vTo are normalized

        auto r = vFrom.dot( vTo ) + 1;

        if ( r < std::numeric_limits<T>::epsilon() ) {

            // vFrom and vTo point in opposite directions

            r = 0;

            if ( std::abs( vFrom.x ) > std::abs( vFrom.z ) ) {

                this->_x = - vFrom.y;
                this->_y = vFrom.x;
                this->_z = 0;
                this->_w = r;

            } else {

                this->_x = 0;
                this->_y = - vFrom.z;
                this->_z = vFrom.y;
                this->_w = r;

            }

        } else {

            // crossVectors( vFrom, vTo ); // inlined to avoid cyclic dependency on Vector3

            this->_x = vFrom.y * vTo.z - vFrom.z * vTo.y;
            this->_y = vFrom.z * vTo.x - vFrom.x * vTo.z;
            this->_z = vFrom.x * vTo.y - vFrom.y * vTo.x;
            this->_w = r;

        }

        return this->normalize();

    }

    T angleTo( const Quaternion<T> &q ) const {

        return 2 * std::acos( std::abs( std::clamp<T>( this->dot( q ), - 1, 1 ) ) );

    }

    Quaternion<T> &rotateTowards( const Quaternion<T> &q, T step ) {

        const auto angle = this->angleTo( q );

        if ( angle == 0 ) return *this;

        const auto t = std::min<T>( 1, step / angle );

        this->slerp( q, t );

        return *this;

    }

    Quaternion<T> &identity() {

        return this->set( 0, 0, 0, 1 );

    }

    Quaternion<T> &invert() {

        // quaternion is assumed to have unit length

        return this->conjugate();

    }

    Quaternion<T> &conjugate() {

        this->_x *= - 1;
        this->_y *= - 1;
        this->_z *= - 1;

        this->_onChangeCallback();

        return *this;

    }

    T dot( const Quaternion<T> &v ) const {

        return this->_x * v._x + this->_y * v._y + this->_z * v._z + this->_w * v._w;

    }

    T lengthSq() const {

        return this->_x * this->_x + this->_y * this->_y + this->_z * this->_z + this->_w * this->_w;

    }

    T length() const {

        return std::sqrt( this->_x * this->_x + this->_y * this->_y + this->_z * this->_z + this->_w * this->_w );

    }

    Quaternion<T> &normalize() {

        auto l = this->length();

        if ( l == 0 ) {

            this->_x = 0;
            this->_y = 0;
            this->_z = 0;
            this->_w = 1;

        } else {

            l = 1 / l;

            this->_x = this->_x * l;
            this->_y = this->_y * l;
            this->_z = this->_z * l;
            this->_w = this->_w * l;

        }

        this->_onChangeCallback();

        return *this;

    }

    Quaternion<T> &multiply( const Quaternion<T> &q ) {

        return this->multiplyQuaternions( *this, q );

    }

    Quaternion<T> &premultiply( const Quaternion<T> &q ) {

        return this->multiplyQuaternions( q, *this );

    }

    Quaternion<T> &multiplyQuaternions( const Quaternion<T> &a, const Quaternion<T> &b ) {

        // from http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/code/index.htm

        const auto qax = a._x, qay = a._y, qaz = a._z, qaw = a._w;
        const auto qbx = b._x, qby = b._y, qbz = b._z, qbw = b._w;

        this->_x = qax * qbw + qaw * qbx + qay * qbz - qaz * qby;
        this->_y = qay * qbw + qaw * qby + qaz * qbx - qax * qbz;
        this->_z = qaz * qbw + qaw * qbz + qax * qby - qay * qbx;
        this->_w = qaw * qbw - qax * qbx - qay * qby - qaz * qbz;

        this->_onChangeCallback();

        return *this;

    }

    Quaternion<T> &slerp( const Quaternion<T> &qb, T t ) {

        if ( t == 0 ) return *this;
        if ( t == 1 ) return this->copy( qb );

        const auto x = this->_x, y = this->_y, z = this->_z, w = this->_w;

        // http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/slerp/

        auto cosHalfTheta = w * qb._w + x * qb._x + y * qb._y + z * qb._z;

        if ( cosHalfTheta < 0 ) {

            this->_w = - qb._w;
            this->_x = - qb._x;
            this->_y = - qb._y;
            this->_z = - qb._z;

            cosHalfTheta = - cosHalfTheta;

        } else {

            this->copy( qb );

        }

        if ( cosHalfTheta >= 1.0 ) {

            this->_w = w;
            this->_x = x;
            this->_y = y;
            this->_z = z;

            return *this;

        }

        const auto sqrSinHalfTheta = 1.0 - cosHalfTheta * cosHalfTheta;

        if ( sqrSinHalfTheta <= std::numeric_limits<T>::epsilon() ) {

            const auto s = 1 - t;
            this->_w = s * w + t * this->_w;
            this->_x = s * x + t * this->_x;
            this->_y = s * y + t * this->_y;
            this->_z = s * z + t * this->_z;

            this->normalize();
            this->_onChangeCallback();

            return *this;

        }

        const auto sinHalfTheta = std::sqrt( sqrSinHalfTheta );
        const auto halfTheta = std::atan2( sinHalfTheta, cosHalfTheta );
        const auto ratioA = std::sin( ( 1 - t ) * halfTheta ) / sinHalfTheta,
            ratioB = std::sin( t * halfTheta ) / sinHalfTheta;

        this->_w = ( w * ratioA + this->_w * ratioB );
        this->_x = ( x * ratioA + this->_x * ratioB );
        this->_y = ( y * ratioA + this->_y * ratioB );
        this->_z = ( z * ratioA + this->_z * ratioB );

        this->_onChangeCallback();

        return *this;

    }

    Quaternion<T> &slerpQuaternions( const Quaternion<T> &qa, const Quaternion<T> &qb, T t ) {

        return this->copy( qa ).slerp( qb, t );

    }

    bool equals( const Quaternion<T> &quaternion ) const {

        return ( quaternion._x == this->_x ) && ( quaternion._y == this->_y ) && ( quaternion._z == this->_z ) && ( quaternion._w == this->_w );

    }

    template <typename V>
    Quaternion<T> &fromArray( const V &array, int offset = 0 ) {

        auto it = array.begin() + offset;
        this->_x = *it++;
        this->_y = *it++;
        this->_z = *it++;
        this->_w = *it;

        this->_onChangeCallback();

        return this;

    }

    template <typename V>
    V &toArray(V &array, int offset = 0) const {
        assert(array.size() >= offset + 4);
        array[ offset ] = this->_x;
        array[ offset + 1 ] = this->_y;
        array[ offset + 2 ] = this->_z;
        array[ offset + 3 ] = this->_w;
        return array;
    }
    template <typename V>
    V toArray() const {
        V array;
        array.reserve(4);
        array.push_back(this->_x);
        array.push_back(this->_y);
        array.push_back(this->_z);
        array.push_back(this->_w);
        return array;
    }

#if 0
    Quaternion<T> &fromBufferAttribute( const BufferAttribute<T> &attribute, int index ) {

        this->_x = attribute.getX( index );
        this->_y = attribute.getY( index );
        this->_z = attribute.getZ( index );
        this->_w = attribute.getW( index );

        return *this;

    }
#endif

    Quaternion<T> &_onChange( const std::function<void()> &callback ) {

        this->_onChangeCallback = callback;

        return *this;

    }

protected:
    T _x, _y, _z, _w;
    std::function<void()> _onChangeCallback;

};

}  // namespace three
