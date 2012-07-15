import subprocess, os, glob


#When this module is imported, run sphinx-apidoc for GUI and PDSim.  It is fast.
if raw_input(r'run sphinx-apidoc (y/[n])? ') == 'y':
    import subprocess
    print 'Building API documentation'
    print subprocess.check_output(['sphinx-apidoc','-f','-o','GUI_apidoc','../GUI'])
    import PDSim
    PDSim_path=os.path.join('..','PDSim')
    #PDSim_path=PDSim.__file__.rsplit(os.sep,1)[0]
    print subprocess.check_output(['sphinx-apidoc','-f','-o','PDSim_apidoc',PDSim_path],shell=False)
    
##     #Remove the undoc-members in each of the api-doc generated files
##     for file in glob.glob('GUI_apidoc/*.rst')+glob.glob('PDSim_apidoc/*.rst'):
##         lines = open(file,'r').readlines()
##         lines = [line for line in lines if line.find('undoc')<0]
##         open(file,'w').write(''.join(lines))

fp = open('log.txt','w')
subprocess.call(['make','html'],stdout = fp,stderr = fp, shell=True)
