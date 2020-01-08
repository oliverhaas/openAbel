from setuptools import setup, Extension
from Cython.Build import cythonize


ecadef = ["-O3", "-Wunused-but-set-variable"]
compdir = {'boundscheck': False, 'nonecheck': False, 'wraparound': False, 'cdivision': True, 
           'profile': False, 'infer_types': False, 'binding': True, 'language_level' : '3'}


extensions = 	[
		Extension('openAbel.abel.base',
		    sources=['openAbel/abel/base.pyx'],
		    extra_compile_args = ecadef,
		    ),
		Extension('openAbel.abel.hansenLaw',
		    sources=['openAbel/abel/hansenLaw.pyx'],
		    extra_compile_args = ecadef,
		    ),
		Extension('openAbel.abel.trap',
		    sources=['openAbel/abel/trap.pyx'],
		    extra_compile_args = ecadef,
		    ),
		Extension('openAbel.abel.fmm',
		    sources=['openAbel/abel/fmm.pyx'],
		    extra_compile_args = ecadef,
		    ),
		Extension('openAbel.abel.wrap',
		    sources=['openAbel/abel/wrap.pyx'],
		    extra_compile_args = ecadef,
		    ),
		Extension('openAbel.constants',
		    sources=['openAbel/constants.pyx'],
		    extra_compile_args = ecadef,
		    ),
		Extension('openAbel.mathFun',
		    sources=['openAbel/mathFun.pyx'],
		    extra_compile_args = ecadef,
		    ),
		]


vers = '0.5'
setup(name = 'openAbel',
      version = vers,
      packages = ['openAbel', 
                  'openAbel.abel'],
      package_data={'openAbel': ['*.pxd'], 
                    'openAbel.abel': ['*.pxd','coeffsData/*']},
      ext_modules = cythonize(extensions, compiler_directives = compdir)
     )


logoArt = """

                             _   _         _ 
        ___ _ __  ___ _ _   /_\ | |__  ___| |
       / _ \ '_ \/ -_) ' \ / _ \| '_ \/ -_) |
       \___/ .__/\___|_||_/_/ \_\_.__/\___|_|
           |_|                               
      
openAbel """ + vers + """  Copyright (C) 2016-2020  Oliver Sebastian Haas
                                             
"""
print(logoArt)
     
     
