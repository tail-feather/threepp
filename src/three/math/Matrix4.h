#pragma once

#include <array>
#include <memory>

namespace three {

template <typename T>
class Vector3;

template <typename T>
class Matrix4 {
public:
  Matrix4()
    : elements({
      1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1, 0,
      0, 0, 0, 1,
    })
    , _v1(nullptr)
    , _m1(nullptr)
    , _zero(nullptr)
    , _one(nullptr)
    , _x(nullptr)
    , _y(nullptr)
    , _z(nullptr)
  {}
  ~Matrix4()
  {}
  Matrix4(const Matrix4<T> &copy)
    : elements(copy.elements)
    , _v1()
    , _m1()
    , _zero()
    , _one()
    , _x()
    , _y()
    , _z()
  {}
  Matrix4<T> &operator =(const Matrix4<T> &copy) {
    this->elements = copy.elements;
    return *this;
  }
  Matrix4(Matrix4<T> &&move) = default;
  Matrix4<T> &operator =(Matrix4<T> &move) = default;

	Matrix4<T> &set(T n11, T n12, T n13, T n14, T n21, T n22, T n23, T n24, T n31, T n32, T n33, T n34, T n41, T n42, T n43, T n44 ) {

		auto &te = this->elements;

		te[ 0 ] = n11; te[ 4 ] = n12; te[ 8 ] = n13; te[ 12 ] = n14;
		te[ 1 ] = n21; te[ 5 ] = n22; te[ 9 ] = n23; te[ 13 ] = n24;
		te[ 2 ] = n31; te[ 6 ] = n32; te[ 10 ] = n33; te[ 14 ] = n34;
		te[ 3 ] = n41; te[ 7 ] = n42; te[ 11 ] = n43; te[ 15 ] = n44;

		return *this;

	}

	Matrix4<T> &identity() {

		this->set(

			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1

		);

		return *this;

	}

	Matrix4<T> clone() const {

		return Matrix4().fromArray( this->elements );

	}

	Matrix4<T> &copy( const Matrix4<T> &m ) {

		auto &te = this->elements;
		const auto &me = m.elements;

		te[ 0 ] = me[ 0 ]; te[ 1 ] = me[ 1 ]; te[ 2 ] = me[ 2 ]; te[ 3 ] = me[ 3 ];
		te[ 4 ] = me[ 4 ]; te[ 5 ] = me[ 5 ]; te[ 6 ] = me[ 6 ]; te[ 7 ] = me[ 7 ];
		te[ 8 ] = me[ 8 ]; te[ 9 ] = me[ 9 ]; te[ 10 ] = me[ 10 ]; te[ 11 ] = me[ 11 ];
		te[ 12 ] = me[ 12 ]; te[ 13 ] = me[ 13 ]; te[ 14 ] = me[ 14 ]; te[ 15 ] = me[ 15 ];

		return *this;

	}

	Matrix4<T> &copyPosition( const Matrix4<T> &m ) {

		auto &te = this->elements;
    const auto &me = m.elements;

		te[ 12 ] = me[ 12 ];
		te[ 13 ] = me[ 13 ];
		te[ 14 ] = me[ 14 ];

		return *this;

	}

#if 0
	Matrix4<T> &setFromMatrix3( const Matrix3<T> &m ) {

		const auto &me = m.elements;

		this->set(

			me[ 0 ], me[ 3 ], me[ 6 ], 0,
			me[ 1 ], me[ 4 ], me[ 7 ], 0,
			me[ 2 ], me[ 5 ], me[ 8 ], 0,
			0, 0, 0, 1

		);

		return *this;

	}
#endif

  template <typename V>
	Matrix4<T> &extractBasis( V &xAxis, V &yAxis, V &zAxis ) {

		xAxis.setFromMatrixColumn( *this, 0 );
		yAxis.setFromMatrixColumn( *this, 1 );
		zAxis.setFromMatrixColumn( *this, 2 );

		return *this;

	}

  template <typename V>
	Matrix4<T> &makeBasis( const V &xAxis, const V &yAxis, const V &zAxis ) {

		this->set(
			xAxis.x, yAxis.x, zAxis.x, 0,
			xAxis.y, yAxis.y, zAxis.y, 0,
			xAxis.z, yAxis.z, zAxis.z, 0,
			0, 0, 0, 1
		);

		return *this;

	}

	Matrix4<T> &extractRotation( const Matrix4<T> &m ) {

		// this method does not support reflection matrices

		auto &te = this->elements;
		const auto &me = m.elements;

    _init_cache();

		const auto scaleX = 1 / _v1->setFromMatrixColumn( m, 0 ).length();
		const auto scaleY = 1 / _v1->setFromMatrixColumn( m, 1 ).length();
		const auto scaleZ = 1 / _v1->setFromMatrixColumn( m, 2 ).length();

		te[ 0 ] = me[ 0 ] * scaleX;
		te[ 1 ] = me[ 1 ] * scaleX;
		te[ 2 ] = me[ 2 ] * scaleX;
		te[ 3 ] = 0;

		te[ 4 ] = me[ 4 ] * scaleY;
		te[ 5 ] = me[ 5 ] * scaleY;
		te[ 6 ] = me[ 6 ] * scaleY;
		te[ 7 ] = 0;

		te[ 8 ] = me[ 8 ] * scaleZ;
		te[ 9 ] = me[ 9 ] * scaleZ;
		te[ 10 ] = me[ 10 ] * scaleZ;
		te[ 11 ] = 0;

		te[ 12 ] = 0;
		te[ 13 ] = 0;
		te[ 14 ] = 0;
		te[ 15 ] = 1;

		return *this;

	}

#if 0
	makeRotationFromEuler( euler ) {

		if ( ! ( euler && euler.isEuler ) ) {

			console.error( 'THREE.Matrix4: .makeRotationFromEuler() now expects a Euler rotation rather than a Vector3 and order.' );

		}

		auto &te = this->elements;

		const auto x = euler.x, y = euler.y, z = euler.z;
		const auto a = std::cos( x ), b = std::sin( x );
		const auto c = std::cos( y ), d = std::sin( y );
		const auto e = std::cos( z ), f = std::sin( z );

		if ( euler.order === 'XYZ' ) {

			const ae = a * e, af = a * f, be = b * e, bf = b * f;

			te[ 0 ] = c * e;
			te[ 4 ] = - c * f;
			te[ 8 ] = d;

			te[ 1 ] = af + be * d;
			te[ 5 ] = ae - bf * d;
			te[ 9 ] = - b * c;

			te[ 2 ] = bf - ae * d;
			te[ 6 ] = be + af * d;
			te[ 10 ] = a * c;

		} else if ( euler.order === 'YXZ' ) {

			const ce = c * e, cf = c * f, de = d * e, df = d * f;

			te[ 0 ] = ce + df * b;
			te[ 4 ] = de * b - cf;
			te[ 8 ] = a * d;

			te[ 1 ] = a * f;
			te[ 5 ] = a * e;
			te[ 9 ] = - b;

			te[ 2 ] = cf * b - de;
			te[ 6 ] = df + ce * b;
			te[ 10 ] = a * c;

		} else if ( euler.order === 'ZXY' ) {

			const ce = c * e, cf = c * f, de = d * e, df = d * f;

			te[ 0 ] = ce - df * b;
			te[ 4 ] = - a * f;
			te[ 8 ] = de + cf * b;

			te[ 1 ] = cf + de * b;
			te[ 5 ] = a * e;
			te[ 9 ] = df - ce * b;

			te[ 2 ] = - a * d;
			te[ 6 ] = b;
			te[ 10 ] = a * c;

		} else if ( euler.order === 'ZYX' ) {

			const ae = a * e, af = a * f, be = b * e, bf = b * f;

			te[ 0 ] = c * e;
			te[ 4 ] = be * d - af;
			te[ 8 ] = ae * d + bf;

			te[ 1 ] = c * f;
			te[ 5 ] = bf * d + ae;
			te[ 9 ] = af * d - be;

			te[ 2 ] = - d;
			te[ 6 ] = b * c;
			te[ 10 ] = a * c;

		} else if ( euler.order === 'YZX' ) {

			const ac = a * c, ad = a * d, bc = b * c, bd = b * d;

			te[ 0 ] = c * e;
			te[ 4 ] = bd - ac * f;
			te[ 8 ] = bc * f + ad;

			te[ 1 ] = f;
			te[ 5 ] = a * e;
			te[ 9 ] = - b * e;

			te[ 2 ] = - d * e;
			te[ 6 ] = ad * f + bc;
			te[ 10 ] = ac - bd * f;

		} else if ( euler.order === 'XZY' ) {

			const ac = a * c, ad = a * d, bc = b * c, bd = b * d;

			te[ 0 ] = c * e;
			te[ 4 ] = - f;
			te[ 8 ] = d * e;

			te[ 1 ] = ac * f + bd;
			te[ 5 ] = a * e;
			te[ 9 ] = ad * f - bc;

			te[ 2 ] = bc * f - ad;
			te[ 6 ] = b * e;
			te[ 10 ] = bd * f + ac;

		}

		// bottom row
		te[ 3 ] = 0;
		te[ 7 ] = 0;
		te[ 11 ] = 0;

		// last column
		te[ 12 ] = 0;
		te[ 13 ] = 0;
		te[ 14 ] = 0;
		te[ 15 ] = 1;

		return this;

	}
#endif

  template <typename Q>
	Matrix4<T> &makeRotationFromQuaternion( const Q &q ) {

    _init_cache();

		return this->compose( *_zero, q, *_one );

	}

  template <typename V>
	Matrix4<T> &lookAt( const V &eye, const V &target, const V &up ) {

		auto &te = this->elements;

    _init_cache();

		_z->subVectors( eye, target );

		if ( _z->lengthSq() == 0 ) {

			// eye and target are in the same position

			_z->z = 1;

		}

		_z->normalize();
		_x->crossVectors( up, *_z );

		if ( _x->lengthSq() == 0 ) {

			// up and z are parallel

			if ( std::abs( up.z ) == 1 ) {

				_z->x += 0.0001;

			} else {

				_z->z += 0.0001;

			}

			_z->normalize();
			_x->crossVectors( up, *_z );

		}

		_x->normalize();
		_y->crossVectors( *_z, *_x );

		te[ 0 ] = _x->x; te[ 4 ] = _y->x; te[ 8 ] = _z->x;
		te[ 1 ] = _x->y; te[ 5 ] = _y->y; te[ 9 ] = _z->y;
		te[ 2 ] = _x->z; te[ 6 ] = _y->z; te[ 10 ] = _z->z;

		return *this;

	}

	Matrix4<T> &multiply( const Matrix4<T> &m ) {

		return this->multiplyMatrices( *this, m );

	}

	Matrix4<T> &premultiply( const Matrix4<T> &m ) {

		return this->multiplyMatrices( m, *this );

	}

	Matrix4<T> &multiplyMatrices( const Matrix4<T> &a, const Matrix4<T> &b ) {

		const auto &ae = a.elements;
		const auto &be = b.elements;
		auto &te = this->elements;

		const auto a11 = ae[ 0 ], a12 = ae[ 4 ], a13 = ae[ 8 ], a14 = ae[ 12 ];
		const auto a21 = ae[ 1 ], a22 = ae[ 5 ], a23 = ae[ 9 ], a24 = ae[ 13 ];
		const auto a31 = ae[ 2 ], a32 = ae[ 6 ], a33 = ae[ 10 ], a34 = ae[ 14 ];
		const auto a41 = ae[ 3 ], a42 = ae[ 7 ], a43 = ae[ 11 ], a44 = ae[ 15 ];

		const auto b11 = be[ 0 ], b12 = be[ 4 ], b13 = be[ 8 ], b14 = be[ 12 ];
		const auto b21 = be[ 1 ], b22 = be[ 5 ], b23 = be[ 9 ], b24 = be[ 13 ];
		const auto b31 = be[ 2 ], b32 = be[ 6 ], b33 = be[ 10 ], b34 = be[ 14 ];
		const auto b41 = be[ 3 ], b42 = be[ 7 ], b43 = be[ 11 ], b44 = be[ 15 ];

		te[ 0 ] = a11 * b11 + a12 * b21 + a13 * b31 + a14 * b41;
		te[ 4 ] = a11 * b12 + a12 * b22 + a13 * b32 + a14 * b42;
		te[ 8 ] = a11 * b13 + a12 * b23 + a13 * b33 + a14 * b43;
		te[ 12 ] = a11 * b14 + a12 * b24 + a13 * b34 + a14 * b44;

		te[ 1 ] = a21 * b11 + a22 * b21 + a23 * b31 + a24 * b41;
		te[ 5 ] = a21 * b12 + a22 * b22 + a23 * b32 + a24 * b42;
		te[ 9 ] = a21 * b13 + a22 * b23 + a23 * b33 + a24 * b43;
		te[ 13 ] = a21 * b14 + a22 * b24 + a23 * b34 + a24 * b44;

		te[ 2 ] = a31 * b11 + a32 * b21 + a33 * b31 + a34 * b41;
		te[ 6 ] = a31 * b12 + a32 * b22 + a33 * b32 + a34 * b42;
		te[ 10 ] = a31 * b13 + a32 * b23 + a33 * b33 + a34 * b43;
		te[ 14 ] = a31 * b14 + a32 * b24 + a33 * b34 + a34 * b44;

		te[ 3 ] = a41 * b11 + a42 * b21 + a43 * b31 + a44 * b41;
		te[ 7 ] = a41 * b12 + a42 * b22 + a43 * b32 + a44 * b42;
		te[ 11 ] = a41 * b13 + a42 * b23 + a43 * b33 + a44 * b43;
		te[ 15 ] = a41 * b14 + a42 * b24 + a43 * b34 + a44 * b44;

		return *this;

	}

	Matrix4<T> &multiplyScalar( T s ) {

		auto &te = this->elements;

		te[ 0 ] *= s; te[ 4 ] *= s; te[ 8 ] *= s; te[ 12 ] *= s;
		te[ 1 ] *= s; te[ 5 ] *= s; te[ 9 ] *= s; te[ 13 ] *= s;
		te[ 2 ] *= s; te[ 6 ] *= s; te[ 10 ] *= s; te[ 14 ] *= s;
		te[ 3 ] *= s; te[ 7 ] *= s; te[ 11 ] *= s; te[ 15 ] *= s;

		return *this;

	}

	T determinant() const {

		const auto &te = this->elements;

		const auto n11 = te[ 0 ], n12 = te[ 4 ], n13 = te[ 8 ], n14 = te[ 12 ];
		const auto n21 = te[ 1 ], n22 = te[ 5 ], n23 = te[ 9 ], n24 = te[ 13 ];
		const auto n31 = te[ 2 ], n32 = te[ 6 ], n33 = te[ 10 ], n34 = te[ 14 ];
		const auto n41 = te[ 3 ], n42 = te[ 7 ], n43 = te[ 11 ], n44 = te[ 15 ];

		//TODO: make this more efficient
		//( based on http://www.euclideanspace.com/maths/algebra/matrix/functions/inverse/fourD/index.htm )

		return (
			n41 * (
				+ n14 * n23 * n32
				 - n13 * n24 * n32
				 - n14 * n22 * n33
				 + n12 * n24 * n33
				 + n13 * n22 * n34
				 - n12 * n23 * n34
			) +
			n42 * (
				+ n11 * n23 * n34
				 - n11 * n24 * n33
				 + n14 * n21 * n33
				 - n13 * n21 * n34
				 + n13 * n24 * n31
				 - n14 * n23 * n31
			) +
			n43 * (
				+ n11 * n24 * n32
				 - n11 * n22 * n34
				 - n14 * n21 * n32
				 + n12 * n21 * n34
				 + n14 * n22 * n31
				 - n12 * n24 * n31
			) +
			n44 * (
				- n13 * n22 * n31
				 - n11 * n23 * n32
				 + n11 * n22 * n33
				 + n13 * n21 * n32
				 - n12 * n21 * n33
				 + n12 * n23 * n31
			)

		);

	}

	Matrix4<T> &transpose() {

    auto &te = this->elements;
		T tmp;

		tmp = te[ 1 ]; te[ 1 ] = te[ 4 ]; te[ 4 ] = tmp;
		tmp = te[ 2 ]; te[ 2 ] = te[ 8 ]; te[ 8 ] = tmp;
		tmp = te[ 6 ]; te[ 6 ] = te[ 9 ]; te[ 9 ] = tmp;

		tmp = te[ 3 ]; te[ 3 ] = te[ 12 ]; te[ 12 ] = tmp;
		tmp = te[ 7 ]; te[ 7 ] = te[ 13 ]; te[ 13 ] = tmp;
		tmp = te[ 11 ]; te[ 11 ] = te[ 14 ]; te[ 14 ] = tmp;

		return *this;

	}

	Matrix4<T> &setPosition( T x, T y, T z ) {

		auto &te = this->elements;

		te[ 12 ] = x;
		te[ 13 ] = y;
		te[ 14 ] = z;

		return *this;

	}
  template <typename V>
  Matrix4<T> &setPosition( const V &v ) {

		auto &te = this->elements;

		te[ 12 ] = v.x;
		te[ 13 ] = v.y;
		te[ 14 ] = v.z;

		return *this;

	}

	Matrix4<T> &invert() {

		// based on http://www.euclideanspace.com/maths/algebra/matrix/functions/inverse/fourD/index.htm
		auto &te = this->elements;

    const auto
			n11 = te[ 0 ], n21 = te[ 1 ], n31 = te[ 2 ], n41 = te[ 3 ],
			n12 = te[ 4 ], n22 = te[ 5 ], n32 = te[ 6 ], n42 = te[ 7 ],
			n13 = te[ 8 ], n23 = te[ 9 ], n33 = te[ 10 ], n43 = te[ 11 ],
			n14 = te[ 12 ], n24 = te[ 13 ], n34 = te[ 14 ], n44 = te[ 15 ],

			t11 = n23 * n34 * n42 - n24 * n33 * n42 + n24 * n32 * n43 - n22 * n34 * n43 - n23 * n32 * n44 + n22 * n33 * n44,
			t12 = n14 * n33 * n42 - n13 * n34 * n42 - n14 * n32 * n43 + n12 * n34 * n43 + n13 * n32 * n44 - n12 * n33 * n44,
			t13 = n13 * n24 * n42 - n14 * n23 * n42 + n14 * n22 * n43 - n12 * n24 * n43 - n13 * n22 * n44 + n12 * n23 * n44,
			t14 = n14 * n23 * n32 - n13 * n24 * n32 - n14 * n22 * n33 + n12 * n24 * n33 + n13 * n22 * n34 - n12 * n23 * n34;

		const auto det = n11 * t11 + n21 * t12 + n31 * t13 + n41 * t14;

		if ( det == 0 ) return this->set( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 );

		const auto detInv = 1 / det;

		te[ 0 ] = t11 * detInv;
		te[ 1 ] = ( n24 * n33 * n41 - n23 * n34 * n41 - n24 * n31 * n43 + n21 * n34 * n43 + n23 * n31 * n44 - n21 * n33 * n44 ) * detInv;
		te[ 2 ] = ( n22 * n34 * n41 - n24 * n32 * n41 + n24 * n31 * n42 - n21 * n34 * n42 - n22 * n31 * n44 + n21 * n32 * n44 ) * detInv;
		te[ 3 ] = ( n23 * n32 * n41 - n22 * n33 * n41 - n23 * n31 * n42 + n21 * n33 * n42 + n22 * n31 * n43 - n21 * n32 * n43 ) * detInv;

		te[ 4 ] = t12 * detInv;
		te[ 5 ] = ( n13 * n34 * n41 - n14 * n33 * n41 + n14 * n31 * n43 - n11 * n34 * n43 - n13 * n31 * n44 + n11 * n33 * n44 ) * detInv;
		te[ 6 ] = ( n14 * n32 * n41 - n12 * n34 * n41 - n14 * n31 * n42 + n11 * n34 * n42 + n12 * n31 * n44 - n11 * n32 * n44 ) * detInv;
		te[ 7 ] = ( n12 * n33 * n41 - n13 * n32 * n41 + n13 * n31 * n42 - n11 * n33 * n42 - n12 * n31 * n43 + n11 * n32 * n43 ) * detInv;

		te[ 8 ] = t13 * detInv;
		te[ 9 ] = ( n14 * n23 * n41 - n13 * n24 * n41 - n14 * n21 * n43 + n11 * n24 * n43 + n13 * n21 * n44 - n11 * n23 * n44 ) * detInv;
		te[ 10 ] = ( n12 * n24 * n41 - n14 * n22 * n41 + n14 * n21 * n42 - n11 * n24 * n42 - n12 * n21 * n44 + n11 * n22 * n44 ) * detInv;
		te[ 11 ] = ( n13 * n22 * n41 - n12 * n23 * n41 - n13 * n21 * n42 + n11 * n23 * n42 + n12 * n21 * n43 - n11 * n22 * n43 ) * detInv;

		te[ 12 ] = t14 * detInv;
		te[ 13 ] = ( n13 * n24 * n31 - n14 * n23 * n31 + n14 * n21 * n33 - n11 * n24 * n33 - n13 * n21 * n34 + n11 * n23 * n34 ) * detInv;
		te[ 14 ] = ( n14 * n22 * n31 - n12 * n24 * n31 - n14 * n21 * n32 + n11 * n24 * n32 + n12 * n21 * n34 - n11 * n22 * n34 ) * detInv;
		te[ 15 ] = ( n12 * n23 * n31 - n13 * n22 * n31 + n13 * n21 * n32 - n11 * n23 * n32 - n12 * n21 * n33 + n11 * n22 * n33 ) * detInv;

		return *this;

	}

  template <typename V>
	Matrix4<T> &scale( const V &v ) {

		auto &te = this->elements;
		const auto x = v.x, y = v.y, z = v.z;

		te[ 0 ] *= x; te[ 4 ] *= y; te[ 8 ] *= z;
		te[ 1 ] *= x; te[ 5 ] *= y; te[ 9 ] *= z;
		te[ 2 ] *= x; te[ 6 ] *= y; te[ 10 ] *= z;
		te[ 3 ] *= x; te[ 7 ] *= y; te[ 11 ] *= z;

		return *this;

	}

	T getMaxScaleOnAxis() const {

		const auto &te = this->elements;

		const auto scaleXSq = te[ 0 ] * te[ 0 ] + te[ 1 ] * te[ 1 ] + te[ 2 ] * te[ 2 ];
		const auto scaleYSq = te[ 4 ] * te[ 4 ] + te[ 5 ] * te[ 5 ] + te[ 6 ] * te[ 6 ];
		const auto scaleZSq = te[ 8 ] * te[ 8 ] + te[ 9 ] * te[ 9 ] + te[ 10 ] * te[ 10 ];

		return std::sqrt( std::max({ scaleXSq, scaleYSq, scaleZSq }) );

	}

	Matrix4<T> &makeTranslation( T x, T y, T z ) {

		this->set(

			1, 0, 0, x,
			0, 1, 0, y,
			0, 0, 1, z,
			0, 0, 0, 1

		);

		return *this;

	}

	Matrix4<T> &makeRotationX( T theta ) {

		const auto c = std::cos( theta ), s = std::sin( theta );

		this->set(

			1, 0, 0, 0,
			0, c, - s, 0,
			0, s, c, 0,
			0, 0, 0, 1

		);

		return *this;

	}

	Matrix4<T> &makeRotationY( T theta ) {

		const auto c = std::cos( theta ), s = std::sin( theta );

		this->set(

			 c, 0, s, 0,
			 0, 1, 0, 0,
			- s, 0, c, 0,
			 0, 0, 0, 1

		);

		return *this;

	}

	Matrix4<T> &makeRotationZ( T theta ) {

		const auto c = std::cos( theta ), s = std::sin( theta );

		this->set(

			c, - s, 0, 0,
			s, c, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1

		);

		return *this;

	}

  template <typename V>
	Matrix4<T> &makeRotationAxis( const V &axis, T angle ) {

		// Based on http://www.gamedev.net/reference/articles/article1199.asp

		const auto c = std::cos( angle );
		const auto s = std::sin( angle );
		const auto t = 1 - c;
		const auto x = axis.x, y = axis.y, z = axis.z;
		const auto tx = t * x, ty = t * y;

		this->set(

			tx * x + c, tx * y - s * z, tx * z + s * y, 0,
			tx * y + s * z, ty * y + c, ty * z - s * x, 0,
			tx * z - s * y, ty * z + s * x, t * z * z + c, 0,
			0, 0, 0, 1

		);

		return *this;

	}

	Matrix4<T> &makeScale( T x, T y, T z ) {

		this->set(

			x, 0, 0, 0,
			0, y, 0, 0,
			0, 0, z, 0,
			0, 0, 0, 1

		);

		return *this;

	}

	Matrix4<T> &makeShear( T xy, T xz, T yx, T yz, T zx, T zy ) {

		this->set(

			1, yx, zx, 0,
			xy, 1, zy, 0,
			xz, yz, 1, 0,
			0, 0, 0, 1

		);

		return *this;

	}

  template <typename V, typename Q>
	Matrix4<T> &compose( const V &position, const Q &quaternion, const V &scale ) {

		auto &te = this->elements;

		const auto x = quaternion.x(), y = quaternion.y(), z = quaternion.z(), w = quaternion.w();
		const auto x2 = x + x,	y2 = y + y, z2 = z + z;
		const auto xx = x * x2, xy = x * y2, xz = x * z2;
		const auto yy = y * y2, yz = y * z2, zz = z * z2;
		const auto wx = w * x2, wy = w * y2, wz = w * z2;

		const auto sx = scale.x, sy = scale.y, sz = scale.z;

		te[ 0 ] = ( 1 - ( yy + zz ) ) * sx;
		te[ 1 ] = ( xy + wz ) * sx;
		te[ 2 ] = ( xz - wy ) * sx;
		te[ 3 ] = 0;

		te[ 4 ] = ( xy - wz ) * sy;
		te[ 5 ] = ( 1 - ( xx + zz ) ) * sy;
		te[ 6 ] = ( yz + wx ) * sy;
		te[ 7 ] = 0;

		te[ 8 ] = ( xz + wy ) * sz;
		te[ 9 ] = ( yz - wx ) * sz;
		te[ 10 ] = ( 1 - ( xx + yy ) ) * sz;
		te[ 11 ] = 0;

		te[ 12 ] = position.x;
		te[ 13 ] = position.y;
		te[ 14 ] = position.z;
		te[ 15 ] = 1;

		return *this;

	}

  template <typename V, typename Q>
	Matrix4<T> &decompose( V &position, Q &quaternion, V &scale ) {

		const auto &te = this->elements;

    _init_cache();
		auto sx = _v1->set( te[ 0 ], te[ 1 ], te[ 2 ] ).length();
		const auto sy = _v1->set( te[ 4 ], te[ 5 ], te[ 6 ] ).length();
		const auto sz = _v1->set( te[ 8 ], te[ 9 ], te[ 10 ] ).length();

		// if determine is negative, we need to invert one scale
		const auto det = this->determinant();
		if ( det < 0 ) sx = - sx;

		position.x = te[ 12 ];
		position.y = te[ 13 ];
		position.z = te[ 14 ];

    _init_cache();

		// scale the rotation part
		_m1->copy( *this );

		const auto invSX = 1 / sx;
		const auto invSY = 1 / sy;
		const auto invSZ = 1 / sz;

		_m1->elements[ 0 ] *= invSX;
		_m1->elements[ 1 ] *= invSX;
		_m1->elements[ 2 ] *= invSX;

		_m1->elements[ 4 ] *= invSY;
		_m1->elements[ 5 ] *= invSY;
		_m1->elements[ 6 ] *= invSY;

		_m1->elements[ 8 ] *= invSZ;
		_m1->elements[ 9 ] *= invSZ;
		_m1->elements[ 10 ] *= invSZ;

		quaternion.setFromRotationMatrix( *_m1 );

		scale.x = sx;
		scale.y = sy;
		scale.z = sz;

		return *this;

	}

	Matrix4<T> &makePerspective( T left, T right, T top, T bottom, T near, T far ) {

		auto &te = this->elements;
		const auto x = 2 * near / ( right - left );
		const auto y = 2 * near / ( top - bottom );

		const auto a = ( right + left ) / ( right - left );
		const auto b = ( top + bottom ) / ( top - bottom );
		const auto c = - ( far + near ) / ( far - near );
		const auto d = - 2 * far * near / ( far - near );

		te[ 0 ] = x;	te[ 4 ] = 0;	te[ 8 ] = a;	te[ 12 ] = 0;
		te[ 1 ] = 0;	te[ 5 ] = y;	te[ 9 ] = b;	te[ 13 ] = 0;
		te[ 2 ] = 0;	te[ 6 ] = 0;	te[ 10 ] = c;	te[ 14 ] = d;
		te[ 3 ] = 0;	te[ 7 ] = 0;	te[ 11 ] = - 1;	te[ 15 ] = 0;

		return *this;

	}

	Matrix4<T> &makeOrthographic( T left, T right, T top, T bottom, T near, T far ) {

		auto &te = this->elements;
		const auto w = 1.0 / ( right - left );
		const auto h = 1.0 / ( top - bottom );
		const auto p = 1.0 / ( far - near );

		const auto x = ( right + left ) * w;
		const auto y = ( top + bottom ) * h;
		const auto z = ( far + near ) * p;

		te[ 0 ] = 2 * w;	te[ 4 ] = 0;	te[ 8 ] = 0;	te[ 12 ] = - x;
		te[ 1 ] = 0;	te[ 5 ] = 2 * h;	te[ 9 ] = 0;	te[ 13 ] = - y;
		te[ 2 ] = 0;	te[ 6 ] = 0;	te[ 10 ] = - 2 * p;	te[ 14 ] = - z;
		te[ 3 ] = 0;	te[ 7 ] = 0;	te[ 11 ] = 0;	te[ 15 ] = 1;

		return *this;

	}

	bool equals( Matrix4<T> &matrix ) const {

		const auto &te = this->elements;
		const auto &me = matrix.elements;

		for ( auto i = 0; i < 16; i ++ ) {

			if ( te[ i ] != me[ i ] ) return false;

		}

		return true;

	}

  template <typename V>
	Matrix4<T> &fromArray( const V &array, int offset = 0 ) {

    auto it = array.begin() + offset;
		for ( auto i = 0; i < 16; i ++ ) {

			this->elements[ i ] = *it++;

		}

		return *this;

	}

  template <typename V>
  V &toArray(V &array, int offset = 0) const {
  	const auto &te = this->elements;

    assert(array.size() >= offset + 16);
		array[ offset ] = te[ 0 ];
		array[ offset + 1 ] = te[ 1 ];
		array[ offset + 2 ] = te[ 2 ];
		array[ offset + 3 ] = te[ 3 ];

		array[ offset + 4 ] = te[ 4 ];
		array[ offset + 5 ] = te[ 5 ];
		array[ offset + 6 ] = te[ 6 ];
		array[ offset + 7 ] = te[ 7 ];

		array[ offset + 8 ] = te[ 8 ];
		array[ offset + 9 ] = te[ 9 ];
		array[ offset + 10 ] = te[ 10 ];
		array[ offset + 11 ] = te[ 11 ];

		array[ offset + 12 ] = te[ 12 ];
		array[ offset + 13 ] = te[ 13 ];
		array[ offset + 14 ] = te[ 14 ];
		array[ offset + 15 ] = te[ 15 ];
    return array;
  }
  template <typename V>
  V toArray() const {
  	const auto &te = this->elements;
    V array;
    array.reserve(16);
		array.push_back(te[ 0 ]);
		array.push_back(te[ 1 ]);
		array.push_back(te[ 2 ]);
		array.push_back(te[ 3 ]);

		array.push_back(te[ 4 ]);
		array.push_back(te[ 5 ]);
		array.push_back(te[ 6 ]);
		array.push_back(te[ 7 ]);

		array.push_back(te[ 8 ]);
		array.push_back(te[ 9 ]);
		array.push_back(te[ 10 ]);
		array.push_back(te[ 11 ]);

		array.push_back(te[ 12 ]);
		array.push_back(te[ 13 ]);
		array.push_back(te[ 14 ]);
		array.push_back(te[ 15 ]);
    return array;
  }

public:
  std::array<T, 16> elements;

private:
  inline void _init_cache();
  std::unique_ptr<Vector3<T>> _v1;
  std::unique_ptr<Matrix4<T>> _m1;
  std::unique_ptr<Vector3<T>> _zero;
  std::unique_ptr<Vector3<T>> _one;
  std::unique_ptr<Vector3<T>> _x;
  std::unique_ptr<Vector3<T>> _y;
  std::unique_ptr<Vector3<T>> _z;

};

}  // namespace three

#include "./Vector3.h"

template <typename T>
inline void three::Matrix4<T>::_init_cache() {
  if ( ! _v1 ) {
    _v1 = std::make_unique<three::Vector3<T>>();
  }
  if ( ! _m1 ) {
    _m1 = std::make_unique<three::Matrix4<T>>();
  }
  if ( ! _zero ) {
    _zero = std::make_unique<three::Vector3<T>>(0, 0, 0);
  }
  if ( ! _one ) {
    _one = std::make_unique<three::Vector3<T>>(1, 1, 1);
  }
  if ( ! _x ) {
    _x = std::make_unique<three::Vector3<T>>();
  }
  if ( ! _y ) {
    _y = std::make_unique<three::Vector3<T>>();
  }
  if ( ! _z ) {
    _z = std::make_unique<three::Vector3<T>>();
  }
}