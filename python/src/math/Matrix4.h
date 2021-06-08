#include "three/math/Quaternion.h"


template <typename T>
void init_math_Matrix4(pybind11::module &module, const std::string &suffix="") {
  pybind11::class_<three::Matrix4<T>>(module, ("Matrix4" + suffix).c_str())
    .def(pybind11::init<>())
    .def("__repr__", [suffix](three::Matrix4<T> &self){
      return "Matrix4" + suffix + "(\n  "
        + std::to_string(self.elements[0]) + ", " + std::to_string(self.elements[1]) + ", " + std::to_string(self.elements[2]) + ", " + std::to_string(self.elements[3]) + ",\n  "
        + std::to_string(self.elements[4]) + ", " + std::to_string(self.elements[5]) + ", " + std::to_string(self.elements[6]) + ", " + std::to_string(self.elements[7]) + ",\n  "
        + std::to_string(self.elements[8]) + ", " + std::to_string(self.elements[9]) + ", " + std::to_string(self.elements[10]) + ", " + std::to_string(self.elements[11]) + ",\n  "
        + std::to_string(self.elements[12]) + ", " + std::to_string(self.elements[13]) + ", " + std::to_string(self.elements[14]) + ", " + std::to_string(self.elements[15]) + ")";
    })
    .def("set", &three::Matrix4<T>::set, pybind11::return_value_policy::reference_internal)
    .def("identity", &three::Matrix4<T>::identity, pybind11::return_value_policy::reference_internal)
    .def("clone", &three::Matrix4<T>::clone)
    .def("copy", &three::Matrix4<T>::copy, pybind11::return_value_policy::reference_internal)
    .def("copyPosition", &three::Matrix4<T>::copyPosition, pybind11::return_value_policy::reference_internal)
#if 0
    .def("setFromMatrix3", &three::Matrix4<T>::setFromMatrix3, pybind11::return_value_policy::reference_internal)
#endif
    .def("copyPosition", &three::Matrix4<T>::copyPosition, pybind11::return_value_policy::reference_internal)
    .def("extractBasis", &three::Matrix4<T>::template extractBasis<three::Vector3<T>>, pybind11::return_value_policy::reference_internal)
    .def("makeBasis", &three::Matrix4<T>::template makeBasis<three::Vector3<T>>, pybind11::return_value_policy::reference_internal)
    .def("extractRotation", &three::Matrix4<T>::extractRotation, pybind11::return_value_policy::reference_internal)
#if 0
    .def("makeRotationFromEuler", &three::Matrix4<T>::makeRotationFromEuler, pybind11::return_value_policy::reference_internal)
#endif
    .def("makeRotationFromQuaternion", &three::Matrix4<T>::template makeRotationFromQuaternion<three::Quaternion<T>>, pybind11::return_value_policy::reference_internal)
    .def("lookAt", &three::Matrix4<T>::template lookAt<three::Vector3<T>>, pybind11::return_value_policy::reference_internal)
    .def("multiply", &three::Matrix4<T>::multiply, pybind11::return_value_policy::reference_internal)
    .def("premultiply", &three::Matrix4<T>::premultiply, pybind11::return_value_policy::reference_internal)
    .def("multiplyMatrices", &three::Matrix4<T>::multiplyMatrices, pybind11::return_value_policy::reference_internal)
    .def("multiplyScalar", &three::Matrix4<T>::multiplyScalar, pybind11::return_value_policy::reference_internal)
    .def("determinant", &three::Matrix4<T>::determinant)
    .def("transpose", &three::Matrix4<T>::transpose, pybind11::return_value_policy::reference_internal)
    //.def("setPosition", pybind11::overload_cast<T, T, T>(&three::Matrix4<T>::setPosition), pybind11::return_value_policy::reference_internal)
    .def("setPosition", static_cast<three::Matrix4<T>&(three::Matrix4<T>::*)(T, T, T)>(&three::Matrix4<T>::setPosition), pybind11::return_value_policy::reference_internal)
    .def("setPosition", static_cast<three::Matrix4<T>&(three::Matrix4<T>::*)(const three::Vector3<T> &)>(&three::Matrix4<T>::setPosition), pybind11::return_value_policy::reference_internal)
    .def("invert", &three::Matrix4<T>::invert, pybind11::return_value_policy::reference_internal)
    .def("scale", &three::Matrix4<T>::template scale<three::Vector3<T>>, pybind11::return_value_policy::reference_internal)
    .def("getMaxScaleOnAxis", &three::Matrix4<T>::getMaxScaleOnAxis)
    .def("makeTranslation", &three::Matrix4<T>::makeTranslation, pybind11::return_value_policy::reference_internal)
    .def("makeRotationX", &three::Matrix4<T>::makeRotationX, pybind11::return_value_policy::reference_internal)
    .def("makeRotationY", &three::Matrix4<T>::makeRotationY, pybind11::return_value_policy::reference_internal)
    .def("makeRotationZ", &three::Matrix4<T>::makeRotationZ, pybind11::return_value_policy::reference_internal)
    .def("makeRotationAxis", &three::Matrix4<T>::template makeRotationAxis<three::Vector3<T>>, pybind11::return_value_policy::reference_internal)
    .def("makeScale", &three::Matrix4<T>::makeScale, pybind11::return_value_policy::reference_internal)
    .def("makeShear", &three::Matrix4<T>::makeShear, pybind11::return_value_policy::reference_internal)
    .def("compose", &three::Matrix4<T>::template compose<three::Vector3<T>, three::Quaternion<T>>, pybind11::return_value_policy::reference_internal)
    .def("decompose", &three::Matrix4<T>::template decompose<three::Vector3<T>, three::Quaternion<T>>, pybind11::return_value_policy::reference_internal)
    .def("makePerspective", &three::Matrix4<T>::makePerspective, pybind11::return_value_policy::reference_internal)
    .def("makeOrthographic", &three::Matrix4<T>::makeOrthographic, pybind11::return_value_policy::reference_internal)
    .def("equals", &three::Matrix4<T>::equals)
    .def("fromArray", [](three::Matrix4<T> &self, const pybind11::list &array, int offset) {
      self.elements[ 0 ] = pybind11::cast<T>(array[offset]);
      self.elements[ 1 ] = pybind11::cast<T>(array[offset + 1]);
      self.elements[ 2 ] = pybind11::cast<T>(array[offset + 2]);
      self.elements[ 3 ] = pybind11::cast<T>(array[offset + 3]);
      self.elements[ 4 ] = pybind11::cast<T>(array[offset + 4]);
      self.elements[ 5 ] = pybind11::cast<T>(array[offset + 5]);
      self.elements[ 6 ] = pybind11::cast<T>(array[offset + 6]);
      self.elements[ 7 ] = pybind11::cast<T>(array[offset + 7]);
      self.elements[ 8 ] = pybind11::cast<T>(array[offset + 8]);
      self.elements[ 9 ] = pybind11::cast<T>(array[offset + 9]);
      self.elements[ 10 ] = pybind11::cast<T>(array[offset + 10]);
      self.elements[ 11 ] = pybind11::cast<T>(array[offset + 11]);
      self.elements[ 12 ] = pybind11::cast<T>(array[offset + 12]);
      self.elements[ 13 ] = pybind11::cast<T>(array[offset + 13]);
      self.elements[ 14 ] = pybind11::cast<T>(array[offset + 14]);
      self.elements[ 15 ] = pybind11::cast<T>(array[offset + 15]);
      return self;
    }, pybind11::return_value_policy::reference_internal, pybind11::arg("array"), pybind11::arg("offset") = 0)
    .def("toArray", [](three::Matrix4<T> &self) {
      pybind11::list array;
      array.append(self.elements[0]);
      array.append(self.elements[1]);
      array.append(self.elements[2]);
      array.append(self.elements[3]);
      array.append(self.elements[4]);
      array.append(self.elements[5]);
      array.append(self.elements[6]);
      array.append(self.elements[7]);
      array.append(self.elements[8]);
      array.append(self.elements[9]);
      array.append(self.elements[10]);
      array.append(self.elements[11]);
      array.append(self.elements[12]);
      array.append(self.elements[13]);
      array.append(self.elements[14]);
      array.append(self.elements[15]);
      return array;
    })
    .def("toArray", [](three::Matrix4<T> &self, pybind11::list &array, int offset) {
      array[offset] = self.elements[0];
      array[offset + 1] = self.elements[1];
      array[offset + 2] = self.elements[2];
      array[offset + 3] = self.elements[3];
      array[offset + 4] = self.elements[4];
      array[offset + 5] = self.elements[5];
      array[offset + 6] = self.elements[6];
      array[offset + 7] = self.elements[7];
      array[offset + 8] = self.elements[8];
      array[offset + 9] = self.elements[9];
      array[offset + 10] = self.elements[10];
      array[offset + 11] = self.elements[11];
      array[offset + 12] = self.elements[12];
      array[offset + 13] = self.elements[13];
      array[offset + 14] = self.elements[14];
      array[offset + 15] = self.elements[15];
      return array;
    }, pybind11::return_value_policy::reference_internal, pybind11::arg("array"), pybind11::arg("offset") = 0)
    ;
}
