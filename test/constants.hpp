#pragma once

#include <limits>

//#include "three/math/Vector2.h"
#include "three/math/Vector3.h"

namespace constants {

constexpr int x = 2;
constexpr int y = 3;
constexpr int z = 4;
constexpr int w = 5;

constexpr double Infinity = std::numeric_limits<double>::infinity();

#if 0
constexpr three::Vector2<double> negInf2( -Infinity, -Infinity );
constexpr three::Vector2<double> posInf2( Infinity, Infinity );

constexpr three::Vector2<double> negOne2( - 1, - 1 );

constexpr three::Vector2<double> zero2();
constexpr three::Vector2<double> one2( 1, 1 );
constexpr three::Vector2<double> two2( 2, 2 );
#endif

constexpr three::Vector3<double> negInf3( - Infinity, - Infinity, - Infinity );
constexpr three::Vector3<double> posInf3( Infinity, Infinity, Infinity );

constexpr three::Vector3<double> zero3();
constexpr three::Vector3<double> one3( 1, 1, 1 );
constexpr three::Vector3<double> two3( 2, 2, 2 );

constexpr double eps = 0.0001;

}  // namespace constants
