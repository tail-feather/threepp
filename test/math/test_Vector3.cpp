#include "../test.hpp"
#include "../constants.hpp"

#include "three/math/Vector3.h"

BOOST_AUTO_TEST_SUITE(three_Vector3)

BOOST_AUTO_TEST_CASE(Instancing)
{
    auto a = three::Vector3<double>();
    BOOST_CHECK_EQUAL(a.x, 0);
    BOOST_CHECK_EQUAL(a.y, 0);
    BOOST_CHECK_EQUAL(a.z, 0);

    using constants::x;
    using constants::y;
    using constants::z;
    a = three::Vector3<double>(x, y, z);
    BOOST_CHECK_EQUAL(a.x, x);
    BOOST_CHECK_EQUAL(a.y, y);
    BOOST_CHECK_EQUAL(a.z, z);
}

BOOST_AUTO_TEST_CASE(set)
{
    auto a = three::Vector3<double>();
    BOOST_CHECK_EQUAL(a.x, 0);
    BOOST_CHECK_EQUAL(a.y, 0);
    BOOST_CHECK_EQUAL(a.z, 0);

    using constants::x;
    using constants::y;
    using constants::z;
    a.set(x, y, z);
    BOOST_CHECK_EQUAL(a.x, x);
    BOOST_CHECK_EQUAL(a.y, y);
    BOOST_CHECK_EQUAL(a.z, z);
}

BOOST_AUTO_TEST_CASE(copy)
{
    using constants::x;
    using constants::y;
    using constants::z;
    auto a = three::Vector3<double>(x, y, z);
    auto b = three::Vector3<double>().copy(a);
    BOOST_CHECK_EQUAL(b.x, x);
    BOOST_CHECK_EQUAL(b.y, y);
    BOOST_CHECK_EQUAL(b.z, z);

    // ensure that it is a true copy
    a.x = 0;
    a.y = -1;
    a.z = -2;
    BOOST_CHECK_EQUAL(b.x, x);
    BOOST_CHECK_EQUAL(b.y, y);
    BOOST_CHECK_EQUAL(b.z, z);
}

BOOST_AUTO_TEST_CASE(add)
{
    using constants::x;
    using constants::y;
    using constants::z;
    auto a = three::Vector3<double>(x, y, z);
    auto b = three::Vector3<double>(-x, -y, -z);

    a.add(b);
    BOOST_CHECK_EQUAL(a.x, 0);
    BOOST_CHECK_EQUAL(a.y, 0);
    BOOST_CHECK_EQUAL(a.z, 0);

    auto c = three::Vector3<double>().addVectors(b, b);
    BOOST_CHECK_EQUAL(c.x, -2 * x);
    BOOST_CHECK_EQUAL(c.y, -2 * y);
    BOOST_CHECK_EQUAL(c.z, -2 * z);
}

BOOST_AUTO_TEST_CASE(addScaledVector)
{
    using constants::x;
    using constants::y;
    using constants::z;
    auto a = three::Vector3<double>(x, y, z);
    auto b = three::Vector3<double>(2, 3, 4);
    auto s = 3;

    a.addScaledVector(b, 3);
    BOOST_CHECK_EQUAL(a.x, x + b.x * s);
    BOOST_CHECK_EQUAL(a.y, y + b.y * s);
    BOOST_CHECK_EQUAL(a.z, z + b.z * s);
}

BOOST_AUTO_TEST_CASE(sub)
{
    using constants::x;
    using constants::y;
    using constants::z;
    auto a = three::Vector3<double>(x, y, z);
    auto b = three::Vector3<double>(-x, -y, -z);

    a.sub(b);
    BOOST_CHECK_EQUAL(a.x, 2 * x);
    BOOST_CHECK_EQUAL(a.y, 2 * y);
    BOOST_CHECK_EQUAL(a.z, 2 * z);

    auto c = three::Vector3<double>().subVectors(a, a);
    BOOST_CHECK_EQUAL(c.x, 0);
    BOOST_CHECK_EQUAL(c.y, 0);
    BOOST_CHECK_EQUAL(c.z, 0);
}

BOOST_AUTO_TEST_CASE(multiplyVectors)
{
    using constants::x;
    using constants::y;
    using constants::z;
    auto a = three::Vector3<double>(x, y, z);
    auto b = three::Vector3<double>(2, 3, -5);

    auto c = three::Vector3<double>().multiplyVectors(a, b);

    BOOST_CHECK_EQUAL(c.x, x * 2);
    BOOST_CHECK_EQUAL(c.y, y * 3);
    BOOST_CHECK_EQUAL(c.z, z * -5);
}

#if 0
BOOST_AUTO_TEST_CASE(applyEuler)
{
    using constants::x;
    using constants::y;
    using constants::z;
    auto a = three::Vector3<double>(x, y, z);
    auto euler = three::Euler<double>(90, -45, 0);
    auto expected = three::Vector3<double>(-2.352970120501014, -4.7441750936226645, 0.9779234597246458);

    using constants::eps;
    a.applyEuler(euler);
    BOOST_CHECK(a.x - expected.x <= eps);
    BOOST_CHECK(a.y - expected.y <= eps);
    BOOST_CHECK(a.z - expected.z <= eps);
}
#endif

BOOST_AUTO_TEST_CASE(applyAxisAngle)
{
    using constants::x;
    using constants::y;
    using constants::z;
    auto a = three::Vector3<double>(x, y, z);
    auto axis = three::Vector3<double>(0, 1, 0);
    auto angle = std::numbers::pi / 4.0;
    auto expected = three::Vector3<double>(3 * std::sqrt(2.0), 3, std::sqrt(2.0));

    using constants::eps;
    a.applyAxisAngle(axis, angle);
    BOOST_CHECK(a.x - expected.x <= eps);
    BOOST_CHECK(a.y - expected.y <= eps);
    BOOST_CHECK(a.z - expected.z <= eps);
}

#if 0
BOOST_AUTO_TEST_CASE(applyMatrix3)
{
    auto a = three::Vector3<double>(x, y, z);
    auto m = three::Matrix3<double>().set(2, 3, 5, 7, 11, 13, 17, 19, 23);

    a.applyMatrix3(m);
    BOOST_CHECK_EQUAL(a.x, 33);
    BOOST_CHECK_EQUAL(a.y, 99);
    BOOST_CHECK_EQUAL(a.z, 183);
}
#endif

#if 0
BOOST_AUTO_TEST_CASE(applyMatrix4)
{
    auto a = three::Vector3<double>(x, y, z);
    auto b = three::Vector4<double>(x, y, z, 1);

    auto m = three::Matrix4().makeRotationX(std::numbers::pi);
    a.applyMatrix4(m);
    b.applyMatrix4(m);
    BOOST_CHECK_EQUAL(a.x, b.x / b.w);
    BOOST_CHECK_EQUAL(a.y, b.y / b.w);
    BOOST_CHECK_EQUAL(a.z, b.z / b.w);

    m = three::Matrix4().makeTranslation(3, 2, 1);
    a.applyMatrix4(m);
    b.applyMatrix4(m);
    BOOST_CHECK_EQUAL(a.x, b.x / b.w);
    BOOST_CHECK_EQUAL(a.y, b.y / b.w);
    BOOST_CHECK_EQUAL(a.z, b.z / b.w);

    m = three::Matrix4().set(
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 1, 0,
    );
    a.applyMatrix4(m);
    b.applyMatrix4(m);
    BOOST_CHECK_EQUAL(a.x, b.x / b.w);
    BOOST_CHECK_EQUAL(a.y, b.y / b.w);
    BOOST_CHECK_EQUAL(a.z, b.z / b.w);
}
#endif

BOOST_AUTO_TEST_CASE(applyQuaternion)
{
    using constants::x;
    using constants::y;
    using constants::z;
    using constants::w;
    auto a = three::Vector3<double>(x, y, z);

    a.applyQuaternion(three::Quaternion<double>());
    BOOST_CHECK_EQUAL(a.x, x);
    BOOST_CHECK_EQUAL(a.y, y);
    BOOST_CHECK_EQUAL(a.z, z);

    a.applyQuaternion(three::Quaternion<double>(x, y, z, w));
    BOOST_CHECK_EQUAL(a.x, 108);
    BOOST_CHECK_EQUAL(a.y, 162);
    BOOST_CHECK_EQUAL(a.z, 216);
}

BOOST_AUTO_TEST_CASE(transformDirection)
{
    using constants::x;
    using constants::y;
    using constants::z;
    auto a = three::Vector3<double>(x, y, z);
    auto m = three::Matrix4<double>();
    auto transformed = three::Vector3<double>(0.3713906763541037, 0.5570860145311556, 0.7427813527082074);

    using constants::eps;
    a.transformDirection(m);
    //BOOST_CHECK(a.x - transformed.x <= eps);
    //BOOST_CHECK(a.y - transformed.y <= eps);
    //BOOST_CHECK(a.z - transformed.z <= eps);
    BOOST_CHECK_EQUAL(a.x, transformed.x);
    BOOST_CHECK_EQUAL(a.y, transformed.y);
    BOOST_CHECK_EQUAL(a.z, transformed.z);
}

BOOST_AUTO_TEST_CASE(clampScalar)
{
    auto a = three::Vector3<double>(-0.01, 0.5, 1.5);
    auto clamped = three::Vector3<double>(0.1, 0.5, 1.0);

    using constants::eps;
    a.clampScalar(0.1, 1.0);
    BOOST_CHECK(std::abs(a.x - clamped.x) <= eps);
    BOOST_CHECK(std::abs(a.y - clamped.y) <= eps);
    BOOST_CHECK(std::abs(a.z - clamped.z) <= eps);
}

BOOST_AUTO_TEST_CASE(negate)
{
    using constants::x;
    using constants::y;
    using constants::z;
    auto a = three::Vector3<double>(x, y, z);

    a.negate();
    BOOST_CHECK_EQUAL(a.x, -x);
    BOOST_CHECK_EQUAL(a.y, -y);
    BOOST_CHECK_EQUAL(a.z, -z);
}

BOOST_AUTO_TEST_CASE(dot)
{
    using constants::x;
    using constants::y;
    using constants::z;
    auto a = three::Vector3<double>(x, y, z);
    auto b = three::Vector3<double>(-x, -y, -z);
    auto c = three::Vector3<double>();

    auto result = a.dot(b);
    BOOST_CHECK_EQUAL(result, (-x * x - y * y - z * z));

    result = a.dot(c);
    BOOST_CHECK_EQUAL(result, 0);
}

BOOST_AUTO_TEST_CASE(manhattanLength)
{
    using constants::x;
    using constants::y;
    using constants::z;
    auto a = three::Vector3<double>(x, 0, 0);
    auto b = three::Vector3<double>(0, -y, 0);
    auto c = three::Vector3<double>(0, 0, z);
    auto d = three::Vector3<double>();

    BOOST_CHECK_EQUAL(a.manhattanLength(), x);
    BOOST_CHECK_EQUAL(b.manhattanLength(), y);
    BOOST_CHECK_EQUAL(c.manhattanLength(), z);
    BOOST_CHECK_EQUAL(d.manhattanLength(), 0);

    a.set(x, y, z);
    BOOST_CHECK_EQUAL(a.manhattanLength(), std::abs(x) + std::abs(y) + std::abs(z));
}

BOOST_AUTO_TEST_CASE(normalize)
{
    using constants::x;
    using constants::y;
    using constants::z;
    auto a = three::Vector3<double>(x, 0, 0);
    auto b = three::Vector3<double>(0, -y, 0);
    auto c = three::Vector3<double>(0, 0, z);

    a.normalize();
    BOOST_CHECK_EQUAL(a.length(), 1);
    BOOST_CHECK_EQUAL(a.x, 1);

    b.normalize();
    BOOST_CHECK_EQUAL(b.length(), 1);
    BOOST_CHECK_EQUAL(b.y, -1);

    c.normalize();
    BOOST_CHECK_EQUAL(c.length(), 1);
    BOOST_CHECK_EQUAL(c.z, 1);
}

BOOST_AUTO_TEST_CASE(setLength)
{
    using constants::x;
    using constants::y;
    auto a = three::Vector3<double>(x, 0, 0);
    BOOST_CHECK_EQUAL(a.length(), x);
    a.setLength(y);
    BOOST_CHECK_EQUAL(a.length(), y);

    a = three::Vector3<double>(0, 0, 0);
    BOOST_CHECK_EQUAL(a.length(), 0);
    a.setLength(y);
    BOOST_CHECK_EQUAL(a.length(), 0);
    //a.setLength();
    //BOOST_CHECK(std::isnan(a.length()));
}

BOOST_AUTO_TEST_CASE(cross)
{
    using constants::x;
    using constants::y;
    using constants::z;
    auto a = three::Vector3<double>(x, y, z);
    auto b = three::Vector3<double>(2 * x, -y, 0.5 * z);
    auto crossed = three::Vector3<double>(18, 12, -18);

    using constants::eps;
    a.cross(b);
    BOOST_CHECK(a.x - crossed.x <= eps);
    BOOST_CHECK(a.y - crossed.y <= eps);
    BOOST_CHECK(a.z - crossed.z <= eps);
}

BOOST_AUTO_TEST_CASE(crossVectors)
{
    using constants::x;
    using constants::y;
    using constants::z;
    auto a = three::Vector3<double>(x, y, z);
    auto b = three::Vector3<double>(x, -y, z);
    auto c = three::Vector3<double>();
    auto crossed = three::Vector3<double>(24, 0, -12);

    using constants::eps;
    c.crossVectors(a, b);
    BOOST_CHECK(c.x - crossed.x <= eps);
    BOOST_CHECK(c.y - crossed.y <= eps);
    BOOST_CHECK(c.z - crossed.z <= eps);
}

BOOST_AUTO_TEST_CASE(projectOnVector)
{
    auto a = three::Vector3<double>(1, 0, 0);
    auto b = three::Vector3<double>();
    auto normal = three::Vector3<double>(10, 0, 0);

    BOOST_CHECK(b.copy(a).projectOnVector(normal).equals(three::Vector3<double>(1, 0, 0)));

    a.set(0, 1, 0);
    BOOST_CHECK(b.copy(a).projectOnVector(normal).equals(three::Vector3<double>(0, 0, 0)));

    a.set(0, 0, -1);
    BOOST_CHECK(b.copy(a).projectOnVector(normal).equals(three::Vector3<double>(0, 0, 0)));

    a.set(-1, 0, 0);
    BOOST_CHECK(b.copy(a).projectOnVector(normal).equals(three::Vector3<double>(-1, 0, 0)));
}

BOOST_AUTO_TEST_CASE(projectOnPlane)
{
    auto a = three::Vector3<double>(1, 0, 0);
    auto b = three::Vector3<double>();
    auto normal = three::Vector3<double>(1, 0, 0);

    BOOST_CHECK(b.copy(a).projectOnPlane(normal).equals(three::Vector3<double>(0, 0, 0)));

    a.set(0, 1, 0);
    BOOST_CHECK(b.copy(a).projectOnPlane(normal).equals(three::Vector3<double>(0, 1, 0)));

    a.set(0, 0, -1);
    BOOST_CHECK(b.copy(a).projectOnPlane(normal).equals(three::Vector3<double>(0, 0, -1)));

    a.set(-1, 0, 0);
    BOOST_CHECK(b.copy(a).projectOnPlane(normal).equals(three::Vector3<double>(0, 0, 0)));
}

BOOST_AUTO_TEST_CASE(reflect)
{
    auto a = three::Vector3<double>();
    auto normal = three::Vector3<double>(0, 1, 0);
    auto b = three::Vector3<double>();

    a.set(0, -1, 0);
    BOOST_CHECK(b.copy(a).reflect(normal).equals(three::Vector3<double>(0, 1, 0)));

    a.set(1, -1, 0);
    BOOST_CHECK(b.copy(a).reflect(normal).equals(three::Vector3<double>(1, 1, 0)));

    a.set(1, -1, 0);
    normal.set(0, -1, 0);
    BOOST_CHECK(b.copy(a).reflect(normal).equals(three::Vector3<double>(1, 1, 0)));
}

BOOST_AUTO_TEST_CASE(angleTo)
{
    auto a = three::Vector3<double>(0, - 0.18851655680720186, 0.9820700116639124);
    auto b = three::Vector3<double>(0, 0.18851655680720186, - 0.9820700116639124);

    BOOST_CHECK_EQUAL(a.angleTo(a), 0);
    BOOST_CHECK_EQUAL(a.angleTo(b), std::numbers::pi);

    auto x = three::Vector3<double>(1, 0, 0);
    auto y = three::Vector3<double>(0, 1, 0);
    auto z = three::Vector3<double>(0, 0, 1);

    BOOST_CHECK_EQUAL(x.angleTo(y), std::numbers::pi / 2);
    BOOST_CHECK_EQUAL(x.angleTo(z), std::numbers::pi / 2);
    BOOST_CHECK_EQUAL(z.angleTo(x), std::numbers::pi / 2);

    BOOST_CHECK(std::abs(x.angleTo(three::Vector3<double>(1, 1, 0)) - std::numbers::pi / 4) <= 0.0000001);
}

#if 0
BOOST_AUTO_TEST_CASE(setFromSpherical)
{
    auto a = three::Vector3<double>();
    const auto phi = std::acos(-0.5);
    const auto theta = std::sqrt(std::numbers::pi) * phi;
    const auto sph = three::Spherical<double>(10, phi, theta);
    const auto expected = three::Vector3<double>(- 4.677914006701843, - 5, - 7.288149322420796);

    using constants::eps;
    a.setFromSpherical(sph);
    BOOST_CHECK(std::abs(a.x - expected.x) <= eps);
    BOOST_CHECK(std::abs(a.y - expected.y) <= eps);
    BOOST_CHECK(std::abs(a.z - expected.z) <= eps);
}

BOOST_AUTO_TEST_CASE(setFromCylindrical)
{
    auto a = three::Vector3<double>();
    const auto phi = std::acos(-0.5);
    const auto theta = std::sqrt(std::numbers::pi) * phi;
    const auto cyl = three::Cylindrical<double>(10, Math.PI * 0.125, 20);
    const auto expected = three::Vector3<double>(3.826834323650898, 20, 9.238795325112868);

    using constants::eps;
    a.setFromSpherical(cyl);
    BOOST_CHECK(std::abs(a.x - expected.x) <= eps);
    BOOST_CHECK(std::abs(a.y - expected.y) <= eps);
    BOOST_CHECK(std::abs(a.z - expected.z) <= eps);
}
#endif

BOOST_AUTO_TEST_CASE(setFromMatrixPosition)
{
    auto a = three::Vector3<double>();
    auto m = three::Matrix4<double>().set(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53);

    a.setFromMatrixPosition(m);

    BOOST_CHECK_EQUAL(a.x, 7);
    BOOST_CHECK_EQUAL(a.y, 19);
    BOOST_CHECK_EQUAL(a.z, 37);
}

BOOST_AUTO_TEST_CASE(setFromMatrixScale)
{
    auto a = three::Vector3<double>();
    auto m = three::Matrix4<double>().set(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53);
    auto expected = three::Vector3<double>(25.573423705088842, 31.921779399024736, 35.70714214271425);

    using constants::eps;
    a.setFromMatrixScale(m);
    BOOST_CHECK(std::abs(a.x - expected.x) <= eps);
    BOOST_CHECK(std::abs(a.y - expected.y) <= eps);
    BOOST_CHECK(std::abs(a.z - expected.z) <= eps);
}

BOOST_AUTO_TEST_CASE(setFromMatrixColumn)
{
    auto a = three::Vector3<double>();
    auto m = three::Matrix4<double>().set(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53);

    a.setFromMatrixColumn(m, 0);
    BOOST_CHECK_EQUAL(a.x, 2);
    BOOST_CHECK_EQUAL(a.y, 11);
    BOOST_CHECK_EQUAL(a.z, 23);

    a.setFromMatrixColumn(m, 2);
    BOOST_CHECK_EQUAL(a.x, 5);
    BOOST_CHECK_EQUAL(a.y, 17);
    BOOST_CHECK_EQUAL(a.z, 31);
}

BOOST_AUTO_TEST_CASE(equals)
{
    using constants::x;
    using constants::y;
    using constants::z;
    auto a = three::Vector3<double>(x, 0, z);
    auto b = three::Vector3<double>(0, -y, 0);

    BOOST_CHECK(a.x != b.x);
    BOOST_CHECK(a.y != b.y);
    BOOST_CHECK(a.z != b.z);

    BOOST_CHECK( ! a.equals(b) );
    BOOST_CHECK( ! b.equals(a) );

    a.copy(b);
    BOOST_CHECK(a.x == b.x);
    BOOST_CHECK(a.y == b.y);
    BOOST_CHECK(a.z == b.z);

    BOOST_CHECK( a.equals(b) );
    BOOST_CHECK( b.equals(a) );
}

BOOST_AUTO_TEST_CASE(fromArray)
{
    auto a = three::Vector3<double>();
    std::vector<double> array = {
        1, 2, 3, 4, 5, 6,
    };

    a.fromArray(array);
    BOOST_CHECK_EQUAL(a.x, 1);
    BOOST_CHECK_EQUAL(a.y, 2);
    BOOST_CHECK_EQUAL(a.z, 3);

    a.fromArray(array, 3);
    BOOST_CHECK_EQUAL(a.x, 4);
    BOOST_CHECK_EQUAL(a.y, 5);
    BOOST_CHECK_EQUAL(a.z, 6);
}

BOOST_AUTO_TEST_CASE(toArray)
{
    using constants::x;
    using constants::y;
    using constants::z;
    auto a = three::Vector3<double>(x, y, z);

    auto array = a.toArray<std::vector<double>>();
    BOOST_CHECK_EQUAL(array[0], x);
    BOOST_CHECK_EQUAL(array[1], y);
    BOOST_CHECK_EQUAL(array[2], z);

    array.clear();
    array.resize(3, -1);
    a.toArray(array);
    BOOST_CHECK_EQUAL(array[0], x);
    BOOST_CHECK_EQUAL(array[1], y);
    BOOST_CHECK_EQUAL(array[2], z);

    array.clear();
    array.resize(4, -1);
    a.toArray(array, 1);
    BOOST_CHECK_EQUAL(array[0], -1);
    BOOST_CHECK_EQUAL(array[1], x);
    BOOST_CHECK_EQUAL(array[2], y);
    BOOST_CHECK_EQUAL(array[3], z);
}

#if 0
BOOST_AUTO_TEST_CASE(fromBufferAttribute)
{
    auto a = three::Vector3<double>();
    auto attr = three::BufferAttribute<double>({1, 2, 3, 4, 5, 6}, 3);

    a.fromBufferAttribute(attr, 0);
    BOOST_CHECK_EQUAL(a.x, 1);
    BOOST_CHECK_EQUAL(a.y, 2);
    BOOST_CHECK_EQUAL(a.z, 3);

    a.fromBufferAttribute(attr, 1);
    BOOST_CHECK_EQUAL(a.x, 4);
    BOOST_CHECK_EQUAL(a.y, 5);
    BOOST_CHECK_EQUAL(a.z, 6);
}
#endif

BOOST_AUTO_TEST_CASE(setX_setY_setZ)
{
    auto a = three::Vector3<double>();
    BOOST_CHECK_EQUAL(a.x, 0);
    BOOST_CHECK_EQUAL(a.y, 0);
    BOOST_CHECK_EQUAL(a.z, 0);

    using constants::x;
    using constants::y;
    using constants::z;
    a.setX(x);
    a.setY(y);
    a.setZ(z);

    BOOST_CHECK_EQUAL(a.x, x);
    BOOST_CHECK_EQUAL(a.y, y);
    BOOST_CHECK_EQUAL(a.z, z);
}

BOOST_AUTO_TEST_CASE(setComponent_getComponent_exceptions)
{
    auto a = three::Vector3<double>();

    const auto message_check = [](const std::runtime_error &err){
        return err.what() == std::string("index is out of range: 3");
    };
    BOOST_CHECK_EXCEPTION(a.setComponent(3, 0), std::runtime_error, message_check);

    BOOST_CHECK_EXCEPTION(a.getComponent(3), std::runtime_error, message_check);
}

BOOST_AUTO_TEST_CASE(min_max_clamp)
{
    using constants::x;
    using constants::y;
    using constants::z;
    auto a = three::Vector3<double>(x, y, z);
    auto b = three::Vector3<double>(-x, -y, -z);
    auto c = three::Vector3<double>();

    c.copy(a).min(b);
    BOOST_CHECK_EQUAL(c.x, -x);
    BOOST_CHECK_EQUAL(c.y, -y);
    BOOST_CHECK_EQUAL(c.z, -z);

    c.copy(a).max(b);
    BOOST_CHECK_EQUAL(c.x, x);
    BOOST_CHECK_EQUAL(c.y, y);
    BOOST_CHECK_EQUAL(c.z, z);

    c.set(-2 * x, 2 * y, -2 * z);
    c.clamp(b, a);
    BOOST_CHECK_EQUAL(c.x, -x);
    BOOST_CHECK_EQUAL(c.y, y);
    BOOST_CHECK_EQUAL(c.z, -z);
}

BOOST_AUTO_TEST_CASE(distanceTo_distanceToSquared)
{
    using constants::x;
    using constants::y;
    using constants::z;
    auto a = three::Vector3<double>(x, 0, 0);
    auto b = three::Vector3<double>(0, -y, 0);
    auto c = three::Vector3<double>(0, 0, z);
    auto d = three::Vector3<double>();

    BOOST_CHECK_EQUAL(a.distanceTo(d), x);
    BOOST_CHECK_EQUAL(a.distanceToSquared(d), x * x);

    BOOST_CHECK_EQUAL(b.distanceTo(d), y);
    BOOST_CHECK_EQUAL(b.distanceToSquared(d), y * y);

    BOOST_CHECK_EQUAL(c.distanceTo(d), z);
    BOOST_CHECK_EQUAL(c.distanceToSquared(d), z * z);
}

BOOST_AUTO_TEST_CASE(setScalar_addScalar_subScalar)
{
    auto a = three::Vector3<double>();
    const auto s = 3;

    a.setScalar(s);
    BOOST_CHECK_EQUAL(a.x, s);
    BOOST_CHECK_EQUAL(a.y, s);
    BOOST_CHECK_EQUAL(a.z, s);

    a.addScalar(s);
    BOOST_CHECK_EQUAL(a.x, 2 * s);
    BOOST_CHECK_EQUAL(a.y, 2 * s);
    BOOST_CHECK_EQUAL(a.z, 2 * s);

    a.subScalar(2 * s);
    BOOST_CHECK_EQUAL(a.x, 0);
    BOOST_CHECK_EQUAL(a.y, 0);
    BOOST_CHECK_EQUAL(a.z, 0);
}

BOOST_AUTO_TEST_CASE(multiply_divide)
{
    using constants::x;
    using constants::y;
    using constants::z;
    auto a = three::Vector3<double>(x, y, z);
    auto b = three::Vector3<double>(2 * x, 2 * y, 2 * z);
    auto c = three::Vector3<double>(4 * x, 4 * y, 4 * z);

    a.multiply(b);
    BOOST_CHECK_EQUAL(a.x, x * b.x);
    BOOST_CHECK_EQUAL(a.y, y * b.y);
    BOOST_CHECK_EQUAL(a.z, z * b.z);

    using constants::eps;
    b.divide(c);
    BOOST_CHECK(std::abs(b.x - 0.5) <= eps);
    BOOST_CHECK(std::abs(b.y - 0.5) <= eps);
    BOOST_CHECK(std::abs(b.z - 0.5) <= eps);
}

BOOST_AUTO_TEST_CASE(multiplyScalar_divideScalar)
{
    using constants::x;
    using constants::y;
    using constants::z;
    auto a = three::Vector3<double>(x, y, z);
    auto b = three::Vector3<double>(-x, -y, -z);

    a.multiplyScalar(-2);
    BOOST_CHECK_EQUAL(a.x, x * -2);
    BOOST_CHECK_EQUAL(a.y, y * -2);
    BOOST_CHECK_EQUAL(a.z, z * -2);

    b.multiplyScalar(-2);
    BOOST_CHECK_EQUAL(b.x, 2 * x);
    BOOST_CHECK_EQUAL(b.y, 2 * y);
    BOOST_CHECK_EQUAL(b.z, 2 * z);

    a.divideScalar(-2);
    BOOST_CHECK_EQUAL(a.x, x);
    BOOST_CHECK_EQUAL(a.y, y);
    BOOST_CHECK_EQUAL(a.z, z);

    b.divideScalar(-2);
    BOOST_CHECK_EQUAL(b.x, -x);
    BOOST_CHECK_EQUAL(b.y, -y);
    BOOST_CHECK_EQUAL(b.z, -z);
}

#if 0
BOOST_AUTO_TEST_CASE(project_unproject)
{
    using constants::x;
    using constants::y;
    using constants::z;
    auto a = three::Vector3<double>(x, y, z);
    auto camera = three::PerspectiveCamera<double>(75, 16.0 / 9, 0.1, 300.0);
    auto projected = three::Vector3<double>(-0.36653213611158914, -0.9774190296309043, 1.0506835611870624);

    using constants::eps;
    a.project(camera);
    BOOST_CHECK(std::abs(a.x - projected.x) <= eps);
    BOOST_CHECK(std::abs(a.y - projected.y) <= eps);
    BOOST_CHECK(std::abs(a.z - projected.z) <= eps);

    a.unproject(camera);
    BOOST_CHECK(std::abs(a.x - x) <= eps);
    BOOST_CHECK(std::abs(a.y - y) <= eps);
    BOOST_CHECK(std::abs(a.z - z) <= eps);
}
#endif

BOOST_AUTO_TEST_CASE(length_lengthSq)
{
    using constants::x;
    using constants::y;
    using constants::z;
    auto a = three::Vector3<double>(x, 0, 0);
    auto b = three::Vector3<double>(0, -y, 0);
    auto c = three::Vector3<double>(0, 0, z);
    auto d = three::Vector3<double>();

    BOOST_CHECK_EQUAL(a.length(), x);
    BOOST_CHECK_EQUAL(a.lengthSq(), x * x);
    BOOST_CHECK_EQUAL(b.length(), y);
    BOOST_CHECK_EQUAL(b.lengthSq(), y * y);
    BOOST_CHECK_EQUAL(c.length(), z);
    BOOST_CHECK_EQUAL(c.lengthSq(), z * z);
    BOOST_CHECK_EQUAL(d.length(), 0);
    BOOST_CHECK_EQUAL(d.lengthSq(), 0);

    a.set(x, y, z);
    BOOST_CHECK_EQUAL(a.length(), std::sqrt(x * x + y * y + z * z));
    BOOST_CHECK_EQUAL(a.lengthSq(), x * x + y * y + z * z);
}

BOOST_AUTO_TEST_CASE(lerp_clone)
{
    using constants::x;
    using constants::y;
    using constants::z;
    auto a = three::Vector3<double>(x, 0, z);
    auto b = three::Vector3<double>(0, -y, 0);

    BOOST_CHECK(a.lerp(a, 0).equals(a.lerp(a, 0.5)));
    BOOST_CHECK(a.lerp(a, 0).equals(a.lerp(a, 1)));

    BOOST_CHECK(a.clone().lerp(b, 0).equals(a));

    BOOST_CHECK_EQUAL(a.clone().lerp(b, 0.5).x, x * 0.5);
    BOOST_CHECK_EQUAL(a.clone().lerp(b, 0.5).y, -y * 0.5);
    BOOST_CHECK_EQUAL(a.clone().lerp(b, 0.5).z, z * 0.5);

    BOOST_CHECK(a.clone().lerp(b, 1).equals(b));
}

BOOST_AUTO_TEST_SUITE_END()
