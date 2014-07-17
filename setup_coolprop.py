"""
A convenience script for building CoolProp from the PDSim source directory
"""
import subprocess, os, sys

subprocess.check_call(['python','setup.py','build','--force','install'], cwd = os.path.join('externals','coolprop','wrappers','Python'), stdout = sys.stdout)