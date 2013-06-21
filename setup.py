try:
    import cython
except:
    raise ImportError('Sorry the required package cython was not found')

try:
    import CoolProp
except ImportError:
    raise ImportError('The required python package CoolProp was not found.  Please go to coolprop.sf.net to obtain a copy')
    
try:
    import psutil
    for proc in psutil.get_process_list():
        cmdline = proc.cmdline
        if cmdline and ''.join(cmdline).find('pycompletionserver.py') > 0:
            proc.terminate()
            break
except ImportError:
    print 'psutil was not found, it is used to kill the python completion server in Eclipse which keeps PDSim from building. psutils can be easy_install-ed or installed using pip'
    
from distutils.core import setup, Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext
from Cython.Distutils.extension import Extension as CyExtension
import sys, shutil, os, glob

version = '2.0pre'

#Modify the __init__ file with this version string
fName = os.path.join('PDSim','__init__.py')
lines = open(fName,'r').readlines()
fp = open(fName,'w')
for line in lines:
    if line.startswith('__version__'):
        line = "__version__ = '" + version + "'"
    fp.write(line+'\n')
fp.close()

if len(sys.argv) == 1:
#    sys.argv += ['build_ext','--inplace','install']
#    sys.argv += ['build','build_ext','install']
#    sys.argv += ['build','install']
    sys.argv += ['clean','build','--force','install']

import Cython

#This will generate HTML to show where there are still pythonic bits hiding out
Cython.Compiler.Options.annotate = True

#Get the numpy include folder
import numpy

#Each of the Pure-Python or PYX files in this list will be compiled to a python module       
pyx_list = [
            "PDSim/core/_bearings.pyx",
            "PDSim/core/_core.pyx",
            "PDSim/core/containers.pyx",
            "PDSim/core/callbacks.pyx",
            "PDSim/core/cycleintegrators.pyx",
            "PDSim/misc/scipylike.pyx",
            "PDSim/flow/flow_models.pyx",
            "PDSim/flow/flow.pyx",
            "PDSim/flow/fanno.pyx",
            "PDSim/scroll/scroll_geo.pyx",
            "PDSim/misc/polymath.pyx",
            "PDSim/misc/stl_utilities.pyx",
            "PDSim/misc/datatypes.pyx",
            "PDSim/recip/_recip.pyx",
            "PDSim/scroll/_scroll.pyx"
            ]

def clean():
    for pyx_file in pyx_list:
        f_root = pyx_file.rsplit('.',1)[0]
        for ending in ['.pyd','.so']:
            fname = f_root+ending
            try:
                os.remove(fname)
                print 'removed',fname
            except OSError:
                pass
            
# Try to remove the generated files in the source tree 
# if you are doing an install to the normal location
#if '--inplace' not in sys.argv or '--clean' in sys.argv:
#    clean()
#    if '--clean' in sys.argv:
#        sys.argv.remove('--clean')

pxd_files = []

ext_module_list = []
for pyx_file in pyx_list:
    sources = [pyx_file]
    #Try to find a PXD backpack if possible
    pxd_file = pyx_file.rsplit('.',1)[0]+'.pxd'
    if os.path.exists(pxd_file):
        sources+=[pxd_file]
    pxd_files.append(pxd_file)
    
    #Build an extension with the sources
    ext_name = pyx_file.rsplit('.',1)[0].replace('/','.')

    ext_module_list.append(CyExtension(ext_name,
                                       sources,
                                       language='c++',
                                       cython_directives=dict(profile = True,
                                                              embed_signature = True)
                                       )
                           )

setup(
  name = 'PDSim',
  version = version,
  author = "Ian Bell",
  author_email='ian.h.bell@gmail.com',
  url='http://pdsim.sourceforge.net',
  description = """A flexible open-source framework for the quasi-steady-state simulation of positive displacement machines including compressors and expanders""",
  packages = ['PDSim','PDSim.core','PDSim.flow','PDSim.plot','PDSim.scroll','PDSim.misc','PDSim.recip'],
  cmdclass={'build_ext': build_ext},
  ext_modules = ext_module_list,
  package_data = {'PDSim':pxd_files,'PDSim.include':glob.glob(os.path.join('CoolProp','*.h'))},
  include_dirs = [numpy.get_include(),CoolProp.get_include_directory()]
)

