import subprocess, os, glob, sys

fp = open('log.txt','w')
subprocess.call('sphinx-apidoc -e -f -o PDSim_apidoc c:\Miniconda\Lib\site-packages\PDSim',shell=True,stdout = sys.stdout) 
subprocess.call('sphinx-apidoc -e -f -o GUI_apidoc ..\GUI',shell=True,stdout = sys.stdout) 
#~ print 'building docs... logging to log.txt'
#subprocess.call('make html',stdout = fp,stderr = fp, shell=True)
subprocess.call('make html',shell=True, stdout=sys.stdout)
subprocess.call(['make','latex'],stdout = fp,stderr = fp, shell=True)
