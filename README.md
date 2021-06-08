# three++

## Usage

### C++

```cpp
#include <numbers>
#include <three/math/Vector3.h>

int main() {
  three::Vector3<double> v(1, 0, 0);
  auto mtx = three::Matrix4<double>().makeRotationZ(std::numbers::pi);
  v.applyMatrix4(mtx);

  std::cout << "(x: " << v.x << ", y: " << v.y << ", z: " << v.z << ")" << std::endl;
  // (x: -1, y: 1.22465e-16, z: 0)
}
```

### Python

```python3
>>> import math
>>> import three
>>> v = three.Vector3(1, 0, 0)
>>> mtx = three.Matrix4().makeRotationZ(math.pi)
>>> v.applyMatrix4(mtx)
Vector3(-1.000000, 0.000000, 0.000000)
```

#### build

```bash
cd python/
python3 -m venv .venv
source .venv/bin/activate
pip install pybind11
python setup.py install
# or
pip install wheel
python setup.py bdist_wheel
```
