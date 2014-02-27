"""
A convenience script for building CoolProp from the PDSim source directory
"""
import subprocess, os

print subprocess.check_output(['python','setup.py','build','--force','install'], cwd = os.path.join('externals','coolprop','wrappers','Python'))