import os
import sys

from setuptools import setup, Extension

try:
    import pybind11
    pybind11_include_path = pybind11.get_include()
except ImportError:
    pybind11_include_path = ""


include_dirs = [
    "../src",
    os.path.expanduser("~/local/include/"),
]
if pybind11_include_path:
    include_dirs.insert(0, pybind11_include_path)

libraries = [
]
library_dirs = [
]
extras = [
    "-std=c++20", "-g",
]

if "GITLAB_CI_INCLUDES" in os.environ:
    include_dirs += os.environ["GITLAB_CI_INCLUDES"].split(":")
if "GITLAB_CI_LIBRARIES" in os.environ:
    library_dirs += os.environ["GITLAB_CI_LIBRARIES"].split(":")


sources = []
for (dirname, _, filenames) in os.walk("src"):
    for filename in filenames:
        if not filename.endswith(".cpp"):
            continue
        sources.append(os.path.join(dirname, filename))


ext = Extension(name="three.__three",
                include_dirs=include_dirs,
                libraries=libraries,
                library_dirs=library_dirs,
                sources=sources,
                extra_compile_args=extras,
                extra_link_args=["-fPIC"])

setup(name="three",
      version="0.0.1",
      description="three++",
      author="AstroArts",
      author_email="mugwort.rc@gmail.com",
      url="https://gitlab.com/mugwort-rc/threepp",
      classifiers=[
          "Development Status :: 3 - Alpha",
          "Intended Audience :: Developers",
          "Programming Language :: C++",
          "Programming Language :: Python",
          "Programming Language :: Python :: 3.7",
          "Programming Language :: Python :: 3.8",
          "Programming Language :: Python :: 3.9",
      ],
      packages=["three"],
      ext_modules=[ext]
)
