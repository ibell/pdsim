import json, os, glob, sys

def get_defaults():
    this_folder = os.path.dirname(__file__)
    file_path = os.path.join(this_folder, 'cx_freeze_defaults.json')
    return json.load(open(file_path, 'r'))
    
def pack(options, argv = None):
    from cx_Freeze import setup, Executable
    import os, sys, shutil

    if len(sys.argv) == 1 and argv is None:
        sys.argv += ['build_exe','--build-exe=PDSimGUI']
    else:
        if '--more_options' in sys.argv:
            i = sys.argv.index('--more_options')
            sys.argv.pop(i)
            
            # Merge options with the dictionary coming from JSON file
            js = json.load(open(sys.argv.pop(i),'r'))
            for k in options:
                if k in js:
                    options[k] += js[k]
    
    includes = options['includes']

    this_folder = os.path.abspath(os.path.dirname(__file__))
    
    include_files = options['include_files']
    for path, search in [('plugins','*.py'),('panels','*.py'),('ico','*.png'), ('imgs','*.png'), ('configs','*.cfg'), ('families','*.py')]:
        for f in glob.glob(os.path.join(this_folder, path, search)):
            prior_abs_path = os.path.join(this_folder, f) # Abs path to the file to be included
            new_rel_path = f.replace(this_folder+os.path.sep, '') # Strip off the common path, leaving only the relative path to the file (relative to this file)
            include_files += [(prior_abs_path, new_rel_path)]

    # The setup for cx_Freeze is different from py2exe. Here I am going to
    # use the Python class Executable from cx_Freeze

    base = None
    if sys.platform == "win32":
        base = "Win32GUI"

    GUI2Exe_Target_1 = Executable(
        # what to build
        script = os.path.join(this_folder, "PDSimGUI.py"),
        init_script = None,
        base = base,
        target_name = "PDSimGUI.exe",
        icon = None
        )

    import os, glob2, numpy, PDSim
    # explore_dirs = [os.path.dirname(numpy.__file__), os.path.dirname(PDSim.__file__)]

    # for directory in explore_dirs:
    #     # Recursively find all .pyd files
    #     files = glob2.glob( os.path.join(directory, '**', '*.pyd') )

    #     # Now we have a list of .pyd files; iterate to build a list of tuples into 
    #     # include files containing the source path and the basename
    #     for f in files:
    #         packages_folder = os.path.normpath(os.path.join(directory,'..')) + os.sep
    #         fn = f.split(packages_folder, 1)[1].replace('\\', '.').split('.pyd', 1)[0]
    #         includes.append(fn)

    # That's serious now: we have all (or almost all) the options cx_Freeze
    # supports. I put them all even if some of them are usually defaulted
    # and not used. Some of them I didn't even know about.
    
    setup(
        
        version = PDSim.__version__,
        description = "No Description",
        author = "Ian Bell et al.",
        name = "PDSim GUI",
        
        options = {"build_exe": {"include_files": include_files,
                                 "includes": includes,
                                 "excludes": options['excludes'],
                                 "packages": options['packages'],
                                 "path": options['path']
                                 }
                   },
                               
        executables = [GUI2Exe_Target_1]
        )

    # This is a place where any post-compile code may go.
    # You can add as much code as you want, which can be used, for example,
    # to clean up your folders or to do some particular post-compilation
    # actions.

    import subprocess
    # Further windows packaging things
    if sys.platform.startswith('win'):
        if not os.path.exists(os.path.join('PDSimGUI','configs')):
            os.mkdir(os.path.join('PDSimGUI','configs'))
        # Compress the dll with upx
        subprocess.call('upx PDSimGUI/*.*', stdout = sys.stdout, stderr = sys.stderr)
        # Make an installer using InnoSetup
        subprocess.call(['C:\Program Files (x86)\Inno Setup 6\Compil32.exe','/cc','package_gui.iss'])
        # Rename the installer to include the PDSim version
        old_name = os.path.join('Output','SetupPDSimGUI.exe')
        import PDSim
        new_name = os.path.join('Output','SetupPDSimGUI_version-'+PDSim.__version__+'.exe')
        if os.path.exists(new_name):
            os.remove(new_name)
        os.rename(old_name, new_name)
        
    # And we are done. That's a setup script :-D

if __name__=='__main__':
    pack(get_defaults())
