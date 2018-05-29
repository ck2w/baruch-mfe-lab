from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = 'Black Model',
    ext_modules = cythonize("black_model.pyx"),
)
