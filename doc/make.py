import subprocess, os, glob

fp = open('log.txt','w')
print 'building docs... logging to log.txt'
subprocess.call(['make','html'],stdout = fp,stderr = fp, shell=True)
subprocess.call(['make','latex'],stdout = fp,stderr = fp, shell=True)
