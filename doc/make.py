import subprocess, os, glob, sys

fp = open('log.txt','w')
import PDSim
print(os.path.dirname(PDSim.__file__))
commons = dict(shell=True,stdout = sys.stdout)
subprocess.call('sphinx-apidoc -e -f -o PDSim_apidoc ' + os.path.dirname(PDSim.__file__), **commons) 
#~ print 'building docs... logging to log.txt'
#subprocess.call('make html',stdout = fp,stderr = fp, shell=True)
subprocess.call(['make','latex'], stdout = fp, stderr = fp)
subprocess.call('make html', **commons)