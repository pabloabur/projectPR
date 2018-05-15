from distutils.core import setup, Extension
from Cython.Build import cythonize

#ext = Extension(name="calculate",
#                sources=["cythonfn.pyx"],
#                include_dirs=[np.get_include()])
setup(ext_modules = cythonize('*.pyx'))
#[Extension("calculate", ["cythonfn.pyx"],include_dirs=[np.get_include()])]
