from distutils.core import setup, Extension
from Cython.Build import cythonize
#from distutils.extension import Extension
#from Cython.Distutils import build_ext

import numpy as np

ext = Extension(name="calculate",
                sources=["cythonfn.pyx"],
                include_dirs=[np.get_include()])
setup( #cmdclass = {'build_ext': build_ext},
      ext_modules = cythonize(ext)#[Extension("calculate", ["cythonfn.pyx"],include_dirs=[np.get_include()])]
      )
