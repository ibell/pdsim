import subprocess


#When this module is imported, run sphinx-apidoc for GUI and PDSim.  It is fast.
if raw_input(r'run sphinx-apidoc (y/[n])? ') == 'y':
    import subprocess
    print 'Building API documentation'
    print subprocess.check_output(['sphinx-apidoc','-f','-o','GUI_apidoc','../GUI'])
    import PDSim
    PDSim_path=PDSim.__file__.rsplit(os.sep,1)[0]
    print subprocess.check_output(['sphinx-apidoc','-f','-o','PDSim_apidoc',PDSim_path],shell=False)

fp = open('log.txt','w')
subprocess.call(['make','html'],stdout = fp,shell=True)
