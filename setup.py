from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

ecadef = ["-O3", "-Wunused-but-set-variable"]
#ecadef = ["-O0", "-g", "-Wunused-but-set-variable"]
iddef = ["/usr/local/include/", "./"]
lddef = ["/usr/local/lib/"]
compdir = {'boundscheck': False, 'nonecheck': False, 'wraparound': False, 'cdivision': True, 'profile': False, 'infer_types': False}


extensions = cythonize([
                        Extension('openAbel.__init__',
                            sources=['openAbel/__init__.pyx'],
                            extra_compile_args = ecadef,
                            include_dirs = iddef
                            ),
                        Extension('openAbel.abel.__init__',
                            sources=['openAbel/abel/__init__.pyx'],
                            extra_compile_args = ecadef,
                            include_dirs = iddef
                            ),
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
                        Extension('openAbel.abel.desingQuad',
                            sources=['openAbel/abel/desingQuad.pyx'],
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


setup(name = 'openChargeState',
      version='0.1',
      packages = ['openAbel', 
                  'openAbel.abel'],
      package_data={'openAbel': ['*.pxd'], 
                    'openAbel.abel': ['*.pxd','data/*']},
      cmdclass = {'build_ext': build_ext},
      ext_modules = extensions
     )

logoArt = """
                       _   _         _ 
  ___ _ __  ___ _ _   /_\ | |__  ___| |
 / _ \ '_ \/ -_) ' \ / _ \| '_ \/ -_) |
 \___/ .__/\___|_||_/_/ \_\_.__/\___|_|
     |_|                               
                                                   
"""
print logoArt
     
     
