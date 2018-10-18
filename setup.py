from __future__ import print_function
try:
    import cython
except:
    raise ImportError('Sorry the required package cython was not found')

try:
    import CoolProp
except ImportError as IE:
    print(IE)
    raise ImportError('The required python package CoolProp was not found or could not be imported.  Please go to coolprop.sf.net to obtain a copy')
    
import warnings, setuptools
from distutils.core import setup
from Cython.Build import cythonize
from Cython.Distutils import build_ext
from Cython.Distutils.extension import Extension
import sys, shutil, os, glob, subprocess

# Get the hash of the git revision
git_hash = '????'
try:
    git_hash = subprocess.check_output(['git', 'rev-parse', 'HEAD']).strip().decode('ascii')
except BaseException as BE:
    print('Error:',)
    print('Unable to extract the git revision, set to placeholder')

# Get the branch of the git revision
git_branch = '????'
try:
    git_branch = subprocess.check_output(['git', 'rev-parse', '--abbrev-ref', 'HEAD']).strip().decode('ascii')
except BaseException as BE:
    print('Error:', BE)
    print('Unable to extract the git branch, set to placeholder')
    
version = '2.10.2'

if len(sys.argv) == 1:
    sys.argv += ['clean','develop']
    #sys.argv += ['clean','install']

# Create the __init__ file with this version string
fName = os.path.join('PDSim','__init__.py')
lines = []
lines.append("__version__ = '" + version + "'")
lines.append("__git_revision__ = '" + git_hash + "'")
lines.append("__git_branch__ = '" + git_branch + "'")
with open(fName,'w') as fp:
    fp.write('\n'.join(lines)+'\n')

# Get the numpy include folder
import numpy

# Each of the Pure-Python or PYX files in this list will be compiled to a python module       
pyx_list = [
            "PDSim/core/_bearings.pyx",
            "PDSim/core/containers.pyx",
            "PDSim/core/callbacks.pyx",
            "PDSim/misc/scipylike.pyx",
            "PDSim/flow/flow_models.pyx",
            "PDSim/flow/flow.pyx",
            "PDSim/flow/fanno.pyx",
            "PDSim/misc/stl_utilities.pyx",
            "PDSim/misc/datatypes.pyx",
            "PDSim/misc/clipper/pyclipper.pyx",
            "PDSim/recip/_recip.pyx",
            "PDSim/scroll/common_scroll_geo.pyx",
            "PDSim/scroll/symm_scroll_geo.pyx",
            "PDSim/scroll/_scroll.pyx"
            ]

def clean():
    for pyx_file in pyx_list:
        f_root = pyx_file.rsplit('.',1)[0]
        for ending in ['.pyd','.so','.h','.cpp']:
            fname = f_root+ending
            try:
                os.remove(fname)
                print('removed',fname)
            except OSError:
                pass

# Try to remove the generated files in the source tree 
# if you are doing an install to the normal location
if '--clean' in sys.argv:
    clean()
    if '--clean' in sys.argv:
        sys.argv.remove('--clean')

for i in range(len(pyx_list)-1,-1,-1):
    if not os.path.exists(pyx_list[i]):
        warnings.warn(pyx_list[i]+' was not found')
        pyx_list.pop(i)

pxd_files = []
ext_module_list = []
package_pxd_files = {}
for pyx_file in pyx_list:
    sources = [pyx_file]
    #Try to find a PXD backpack if possible
    pxd_file = pyx_file.rsplit('.',1)[0]+'.pxd'
    if os.path.exists(pxd_file):
        sources+=[pxd_file]
    pxd_files.append(pxd_file)

    # Add sources for clipper module
    if pyx_file.find('pyclipper') > -1:
        sources += [os.path.join('PDSim','misc','clipper','clipper.cpp')]

    #Build an extension with the sources
    ext_name = pyx_file.rsplit('.',1)[0].replace('/','.')
    
    name = ext_name.rsplit('.',1)[0]
    pxd = pxd_file.split('PDSim/',1)[1].split('/')[1]
    
    if name in package_pxd_files:
        package_pxd_files[name].append(pxd)
    else:
        package_pxd_files[name] = [pxd]

    ext_module_list.append(Extension(ext_name,
                                     sources,
                                     language='c++',
                                     define_macros = [('SWIG',None)], # Workaround to remove dependency on rapidjson in Configuration.h
                                     cython_c_in_temp=True
                                     )
                           )

package_data = package_pxd_files

setup(
  name = 'PDSim',
  version = version,
  author = "Ian Bell",
  author_email='ian.h.bell@gmail.com',
  url='http://pdsim.sourceforge.net',
  requires=['CoolProp','cython','quantities'],
  description = """A flexible open-source framework for the quasi-steady-state simulation of positive displacement machines including compressors and expanders""",
  packages = ['PDSim','PDSim.core','PDSim.flow','PDSim.plot','PDSim.scroll','PDSim.misc','PDSim.recip','PDSim.misc.clipper'],
  cmdclass={'build_ext': build_ext},
  ext_modules = cythonize(ext_module_list, 
                          compiler_directives=dict(profile = True, embedsignature = True, language_level='2'),
                          annotate=True
                          ),
  package_dir = {'PDSim':'PDSim',},
  package_data = package_data,
  include_dirs = [numpy.get_include(), CoolProp.get_include_directory(), "PDSim/misc/clipper", "PDSim/misc/spline"],
  zip_safe = False # no compressed egg; see http://stackoverflow.com/a/29124937/1360263
)
    
