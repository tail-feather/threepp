#include "three/math/Quaternion.h"

#include "three/math/Vector3.h"


template <typename T>
void init_math_Quaternion(pybind11::module &module, const std::string &suffix="") {
  pybind11::class_<three::Quaternion<T>>(module, ("Quaternion" + suffix).c_str())
    .def(pybind11::init<T, T, T, T>(), pybind11::arg("x") = 0, pybind11::arg("y") = 0, pybind11::arg("z") = 0, pybind11::arg("w") = 1)
    .def("__repr__", [suffix](three::Quaternion<T> &self){
      return "Quaternion" + suffix + "(" + std::to_string(self.x()) + ", " + std::to_string(self.y()) + ", " + std::to_string(self.z()) + ", " + std::to_string(self.w()) + ")";
    })
    //.def_static("slerpFlat", &three::Quaternion<T>::slerpFlat)
    //.def_static("multiplyQuaternionsFlat", &three::Quaternion<T>::multiplyQuaternionsFlat)
    .def_property("x", &three::Quaternion<T>::x, &three::Quaternion<T>::setX)
    .def_property("y", &three::Quaternion<T>::y, &three::Quaternion<T>::setY)
    .def_property("z", &three::Quaternion<T>::z, &three::Quaternion<T>::setZ)
    .def_property("w", &three::Quaternion<T>::w, &three::Quaternion<T>::setW)
    .def("set", &three::Quaternion<T>::set, pybind11::return_value_policy::reference_internal)
    .def("clone", &three::Quaternion<T>::clone)
    .def("copy", &three::Quaternion<T>::copy, pybind11::return_value_policy::reference_internal)
#if 0
    .def("setFromEuler", &three::Quaternion<T>::setFromEuler, pybind11::return_value_policy::reference_internal)
#endif
    .def("setFromAxisAngle", &three::Quaternion<T>::template setFromAxisAngle<three::Vector3<T>>, pybind11::return_value_policy::reference_internal)
    .def("setFromRotationMatrix", &three::Quaternion<T>::setFromRotationMatrix, pybind11::return_value_policy::reference_internal)
    .def("setFromUnitVectors", &three::Quaternion<T>::template setFromUnitVectors<three::Vector3<T>>, pybind11::return_value_policy::reference_internal)
    .def("angleTo", &three::Quaternion<T>::angleTo)
    .def("rotateTowards", &three::Quaternion<T>::rotateTowards, pybind11::return_value_policy::reference_internal)
    .def("identity", &three::Quaternion<T>::identity, pybind11::return_value_policy::reference_internal)
    .def("invert", &three::Quaternion<T>::invert, pybind11::return_value_policy::reference_internal)
    .def("conjugate", &three::Quaternion<T>::conjugate, pybind11::return_value_policy::reference_internal)
    .def("dot", &three::Quaternion<T>::dot)
    .def("lengthSq", &three::Quaternion<T>::lengthSq)
    .def("length", &three::Quaternion<T>::length)
    .def("normalize", &three::Quaternion<T>::normalize, pybind11::return_value_policy::reference_internal)
    .def("multiply", &three::Quaternion<T>::multiply, pybind11::return_value_policy::reference_internal)
    .def("premultiply", &three::Quaternion<T>::premultiply, pybind11::return_value_policy::reference_internal)
    .def("multiplyQuaternions", &three::Quaternion<T>::multiplyQuaternions, pybind11::return_value_policy::reference_internal)
    .def("slerp", &three::Quaternion<T>::slerp, pybind11::return_value_policy::reference_internal)
    .def("slerpQuaternions", &three::Quaternion<T>::slerpQuaternions, pybind11::return_value_policy::reference_internal)
    .def("equals", &three::Quaternion<T>::equals)
    .def("fromArray", [](three::Quaternion<T> &self, const pybind11::list &array, int offset) {
      const auto x = pybind11::cast<T>(array[offset]);
      const auto y = pybind11::cast<T>(array[offset + 1]);
      const auto z = pybind11::cast<T>(array[offset + 2]);
      const auto w = pybind11::cast<T>(array[offset + 3]);
      return self.set(x, y, z, w);
    }, pybind11::return_value_policy::reference_internal, pybind11::arg("array"), pybind11::arg("offset") = 0)
    .def("toArray", [](three::Quaternion<T> &self) {
      pybind11::list array;
      array.append(self.x());
      array.append(self.y());
      array.append(self.z());
      array.append(self.w());
      return array;
    })
    .def("toArray", [](three::Quaternion<T> &self, pybind11::list &array, int offset) {
      array[offset] = self.x();
      array[offset + 1] = self.y();
      array[offset + 2] = self.z();
      array[offset + 3] = self.w();
      return array;
    }, pybind11::return_value_policy::reference_internal, pybind11::arg("array"), pybind11::arg("offset") = 0)
#if 0
    .def("fromBufferAttribute", &three::Quaternion<T>::fromBufferAttribute, pybind11::return_value_policy::reference_internal)
#endif
    .def("_onChange", [](three::Quaternion<T> &self, const pybind11::object &callback) {
      self._onChange([callback](){
        callback();
      });
      return self;
    })
    ;
}
