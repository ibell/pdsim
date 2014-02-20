import subprocess, os, glob

fp = open('log.txt','w')
subprocess.call('sphinx-apidoc -e -f -o PDSim_apidoc c:\python27\Lib\site-packages\PDSim c:\python27\Lib\site-packages\PDSim\misc\clipper',shell=True) 
print 'building docs... logging to log.txt'
#subprocess.call('make html',stdout = fp,stderr = fp, shell=True)
subprocess.call('make html',shell=True)
#subprocess.call(['make','latex'],stdout = fp,stderr = fp, shell=True)
