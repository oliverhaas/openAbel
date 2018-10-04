from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize


ecadef = ["-O3", "-Wunused-but-set-variable"]
iddef = ["/usr/local/include/"]#, "./"]
lddef = ["/usr/local/lib/"]
compdir = {'boundscheck': False, 'nonecheck': False, 'wraparound': False, 'cdivision': True, 
           'profile': False, 'infer_types': False, 'binding': True}


extensions = cythonize([
                        Extension('openAbel.abel.base',
                            sources=['openAbel/abel/base.pyx'],
                            extra_compile_args = ecadef,
                            include_dirs = iddef
                            ),
                        Extension('openAbel.abel.hansenLaw',
                            sources=['openAbel/abel/hansenLaw.pyx'],
                            extra_compile_args = ecadef,
                            include_dirs = iddef
                            ),
                        Extension('openAbel.abel.trap',
                            sources=['openAbel/abel/trap.pyx'],
                            extra_compile_args = ecadef,
                            include_dirs = iddef
                            ),
                        Extension('openAbel.abel.fmm',
                            sources=['openAbel/abel/fmm.pyx'],
                            extra_compile_args = ecadef,
                            include_dirs = iddef
                            ),
                        Extension('openAbel.abel.wrap',
                            sources=['openAbel/abel/wrap.pyx'],
                            extra_compile_args = ecadef,
                            include_dirs = iddef
                            ),
                        Extension('openAbel.constants',
                            sources=['openAbel/constants.pyx'],
                            extra_compile_args = ecadef,
                            include_dirs = iddef
                            ),
                        Extension('openAbel.mathFun',
                            sources=['openAbel/mathFun.pyx'],
                            extra_compile_args = ecadef,
                            include_dirs = iddef
                            ),
                        ], 
                        compiler_directives = compdir
                        )


setup(name = 'openAbel',
      version='0.2',
      packages = ['openAbel', 
                  'openAbel.abel'],
      package_data={'openAbel': ['*.pxd'], 
                    'openAbel.abel': ['*.pxd','coeffsData/*']},
      cmdclass = {'build_ext': build_ext},
      ext_modules = extensions
     )

logoArt = """

                             _   _         _ 
        ___ _ __  ___ _ _   /_\ | |__  ___| |
       / _ \ '_ \/ -_) ' \ / _ \| '_ \/ -_) |
       \___/ .__/\___|_||_/_/ \_\_.__/\___|_|
           |_|                               
      
openAbel  Copyright (C) 2016-2018  Oliver Sebastian Haas
                                             
"""
print logoArt
     
     
