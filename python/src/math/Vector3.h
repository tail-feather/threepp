#include "three/math/Vector3.h"


template <typename T>
void init_math_Vector3(pybind11::module &module, const std::string &suffix="") {
  pybind11::class_<three::Vector3<T>>(module, ("Vector3" + suffix).c_str())
    .def(pybind11::init<T, T, T>(), pybind11::arg("x") = 0, pybind11::arg("y") = 0, pybind11::arg("z") = 0)
    .def("__repr__", [suffix](three::Vector3<T> &self){
      return "Vector3" + suffix + "(" + std::to_string(self.x) + ", " + std::to_string(self.y) + ", " + std::to_string(self.z) + ")";
    })
    .def("set", pybind11::overload_cast<T, T, T>(&three::Vector3<T>::set), pybind11::return_value_policy::reference_internal)
    .def("set", pybind11::overload_cast<T, T>(&three::Vector3<T>::set), pybind11::return_value_policy::reference_internal)
    .def("setScalar", &three::Vector3<T>::setScalar, pybind11::return_value_policy::reference_internal)
    .def("setX", &three::Vector3<T>::setX, pybind11::return_value_policy::reference_internal)
    .def("setY", &three::Vector3<T>::setY, pybind11::return_value_policy::reference_internal)
    .def("setZ", &three::Vector3<T>::setZ, pybind11::return_value_policy::reference_internal)
    .def("setComponent", &three::Vector3<T>::setComponent)
    .def("clone", &three::Vector3<T>::clone)
    .def("copy", &three::Vector3<T>::copy, pybind11::return_value_policy::reference_internal)
    .def("add", &three::Vector3<T>::add, pybind11::return_value_policy::reference_internal)
    .def("addScalar", &three::Vector3<T>::addScalar, pybind11::return_value_policy::reference_internal)
    .def("addVectors", &three::Vector3<T>::addVectors, pybind11::return_value_policy::reference_internal)
    .def("addScaledVector", &three::Vector3<T>::addScaledVector, pybind11::return_value_policy::reference_internal)
    .def("sub", &three::Vector3<T>::sub, pybind11::return_value_policy::reference_internal)
    .def("subScalar", &three::Vector3<T>::subScalar, pybind11::return_value_policy::reference_internal)
    .def("subVectors", &three::Vector3<T>::subVectors, pybind11::return_value_policy::reference_internal)
    .def("multiply", &three::Vector3<T>::multiply, pybind11::return_value_policy::reference_internal)
    .def("multiplyScalar", &three::Vector3<T>::multiplyScalar, pybind11::return_value_policy::reference_internal)
    .def("multiplyVectors", &three::Vector3<T>::multiplyVectors, pybind11::return_value_policy::reference_internal)
#if 0
    .def("applyEuler", &three::Vector3<T>::applyEuler, pybind11::return_value_policy::reference_internal)
    .def("applyAxisAngle", &three::Vector3<T>::applyAxisAngle, pybind11::return_value_policy::reference_internal)
    .def("applyMatrix3", &three::Vector3<T>::applyMatrix3, pybind11::return_value_policy::reference_internal)
    .def("applyNormalMatrix", &three::Vector3<T>::applyNormalMatrix, pybind11::return_value_policy::reference_internal)
#endif
    .def("applyMatrix4", &three::Vector3<T>::applyMatrix4, pybind11::return_value_policy::reference_internal)
    .def("applyQuaternion", &three::Vector3<T>::applyQuaternion, pybind11::return_value_policy::reference_internal)
#if 0
    .def("project", &three::Vector3<T>::project, pybind11::return_value_policy::reference_internal)
    .def("unproject", &three::Vector3<T>::unproject, pybind11::return_value_policy::reference_internal)
#endif
#if 0
    .def("transformDirection", &three::Vector3<T>::transformDirection, pybind11::return_value_policy::reference_internal)
#endif
    .def("divide", &three::Vector3<T>::divide, pybind11::return_value_policy::reference_internal)
    .def("divideScalar", &three::Vector3<T>::divideScalar, pybind11::return_value_policy::reference_internal)
    .def("min", &three::Vector3<T>::min, pybind11::return_value_policy::reference_internal)
    .def("max", &three::Vector3<T>::max, pybind11::return_value_policy::reference_internal)
    .def("clamp", &three::Vector3<T>::clamp, pybind11::return_value_policy::reference_internal)
    .def("clampScalar", &three::Vector3<T>::clampScalar, pybind11::return_value_policy::reference_internal)
    .def("clampLength", &three::Vector3<T>::clampLength, pybind11::return_value_policy::reference_internal)
    .def("floor", &three::Vector3<T>::floor, pybind11::return_value_policy::reference_internal)
    .def("ceil", &three::Vector3<T>::ceil, pybind11::return_value_policy::reference_internal)
    .def("round", &three::Vector3<T>::round, pybind11::return_value_policy::reference_internal)
    .def("roundToZero", &three::Vector3<T>::roundToZero, pybind11::return_value_policy::reference_internal)
    .def("negate", &three::Vector3<T>::negate, pybind11::return_value_policy::reference_internal)
    .def("dot", &three::Vector3<T>::dot)
    .def("lengthSq", &three::Vector3<T>::lengthSq)
    .def("length", &three::Vector3<T>::length)
    .def("manhattanLength", &three::Vector3<T>::manhattanLength)
    .def("normalize", &three::Vector3<T>::normalize, pybind11::return_value_policy::reference_internal)
    .def("setLength", &three::Vector3<T>::setLength, pybind11::return_value_policy::reference_internal)
    .def("lerp", &three::Vector3<T>::lerp, pybind11::return_value_policy::reference_internal)
    .def("lerpVectors", &three::Vector3<T>::lerpVectors, pybind11::return_value_policy::reference_internal)
    .def("cross", &three::Vector3<T>::cross, pybind11::return_value_policy::reference_internal)
    .def("crossVectors", &three::Vector3<T>::crossVectors, pybind11::return_value_policy::reference_internal)
    .def("projectOnVector", &three::Vector3<T>::projectOnVector, pybind11::return_value_policy::reference_internal)
    .def("projectOnPlane", &three::Vector3<T>::projectOnPlane, pybind11::return_value_policy::reference_internal)
    .def("reflect", &three::Vector3<T>::reflect, pybind11::return_value_policy::reference_internal)
    .def("angleTo", &three::Vector3<T>::angleTo)
    .def("distanceTo", &three::Vector3<T>::distanceTo)
    .def("distanceToSquared", &three::Vector3<T>::distanceToSquared)
    .def("manhattanDistanceTo", &three::Vector3<T>::manhattanDistanceTo)
#if 0
    .def("setFromSpherical", &three::Vector3<T>::setFromSpherical, pybind11::return_value_policy::reference_internal)
#endif
    .def("setFromSphericalCoords", &three::Vector3<T>::setFromSphericalCoords, pybind11::return_value_policy::reference_internal)
#if 0
    .def("setFromCylindrical", &three::Vector3<T>::setFromCylindrical, pybind11::return_value_policy::reference_internal)
#endif
    .def("setFromCylindricalCoords", &three::Vector3<T>::setFromCylindricalCoords, pybind11::return_value_policy::reference_internal)
#if 0
    .def("setFromMatrixPosition", &three::Vector3<T>::setFromMatrixPosition, pybind11::return_value_policy::reference_internal)
    .def("setFromMatrixScale", &three::Vector3<T>::setFromMatrixScale, pybind11::return_value_policy::reference_internal)
    .def("setFromMatrixColumn", &three::Vector3<T>::setFromMatrixColumn, pybind11::return_value_policy::reference_internal)
    .def("setFromMatrix3Column", &three::Vector3<T>::setFromMatrix3Column, pybind11::return_value_policy::reference_internal)
#endif
    .def("equals", &three::Vector3<T>::equals)
    .def("fromArray", [](three::Vector3<T> &self, const pybind11::list &array, int offset) {
      const auto x = pybind11::cast<T>(array[offset]);
      const auto y = pybind11::cast<T>(array[offset + 1]);
      const auto z = pybind11::cast<T>(array[offset + 2]);
      return self.set(x, y, z);
    }, pybind11::return_value_policy::reference_internal, pybind11::arg("array"), pybind11::arg("offset") = 0)
    .def("toArray", [](three::Vector3<T> &self) {
      pybind11::list array;
      array.append(self.x);
      array.append(self.y);
      array.append(self.z);
      return array;
    })
    .def("toArray", [](three::Vector3<T> &self, pybind11::list &array, int offset) {
      array[offset] = self.x;
      array[offset + 1] = self.y;
      array[offset + 2] = self.z;
      return array;
    }, pybind11::return_value_policy::reference_internal, pybind11::arg("array"), pybind11::arg("offset") = 0)
#if 0
    .def("fromBufferAttribute", &three::Vector3<T>::fromBufferAttribute, pybind11::return_value_policy::reference_internal)
#endif
    .def("random", [](three::Vector3<T> &self, const pybind11::object &generator) {
      self.x = pybind11::cast<T>(generator());
      self.y = pybind11::cast<T>(generator());
      self.z = pybind11::cast<T>(generator());
      return self;
    }, pybind11::return_value_policy::reference_internal)
    .def_readwrite("x", &three::Vector3<T>::x)
    .def_readwrite("y", &three::Vector3<T>::y)
    .def_readwrite("z", &three::Vector3<T>::z)
    ;
}
