#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "./math/Vector3.h"
#include "./math/Matrix4.h"
#include "./math/Quaternion.h"

PYBIND11_MODULE(__three, module) {
    init_math_Vector3<double>(module);
    init_math_Matrix4<double>(module);
    init_math_Quaternion<double>(module);
}
