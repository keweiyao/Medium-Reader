from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

setup(
    ext_modules = cythonize([
    Extension("medium", 
			  ["medium.pyx"],
			  language="c++",
			  extra_compile_args=["-std=c++11", '-march=native'],
              libraries=["m"])
    ]),
)
