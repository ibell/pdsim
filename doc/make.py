import subprocess, os, glob

print 'Rebuilding pdsim locally i.e. in the PDSim folder rather than to site-packages'
this_location = os.path.abspath('.')
os.chdir('..')
subprocess.call(['python','setup.py','build_ext','--inplace'])
os.chdir(this_location)

#When this module is imported, run sphinx-apidoc for GUI and PDSim.  It is fast.
print 'Use of sphinx-apidoc is disabled'
if 0: #raw_input(r'run sphinx-apidoc (y/[n])? ') == 'y':
    print 'Building API documentation'
    print subprocess.check_output(['sphinx-apidoc','-f','-o','GUI_apidoc','../GUI'])
    import PDSim
    PDSim_path=os.path.join('..','PDSim')
    print subprocess.check_output(['sphinx-apidoc','-f','-o','PDSim_apidoc',PDSim_path],shell=False)
    
##     #Remove the undoc-members in each of the api-doc generated files
##     for file in glob.glob('GUI_apidoc/*.rst')+glob.glob('PDSim_apidoc/*.rst'):
##         lines = open(file,'r').readlines()
##         lines = [line for line in lines if line.find('undoc')<0]
##         open(file,'w').write(''.join(lines))

fp = open('log.txt','w')
subprocess.call(['make','html'],stdout = fp,stderr = fp, shell=True)
