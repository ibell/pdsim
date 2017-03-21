import subprocess, sys, glob, os

def convert_notebooks():
    template = '''
    ***************************
    Converted IPython notebooks
    ***************************

    .. note:: These static pages were auto-generated from the IPython notebooks residing `here on github <https://github.com/ibell/pdsim/tree/master/doc/notebooks>`_ They can be also be `viewed with nbviewer <http://nbviewer.ipython.org/github/ibell/pdsim/tree/master/doc/notebooks/>`_
    '''

    thisdir = os.path.abspath(os.path.dirname(__file__))

    # Convert the notebooks
    subprocess.check_call('jupyter nbconvert --to rst *.ipynb', shell = True, stdout = sys.stdout, stderr = sys.stderr, cwd = thisdir)

    # Write the hidden table to make sphinx happy
    with open(thisdir+'/index.rst', 'w') as fp:
        fp.write(template)
        fp.write('.. toctree::\n\n')
        for file in glob.glob(thisdir+'/*.rst'):
            with open(file, 'r') as fin:
                file_contents= fin.read()
            file_contents = file_contents.replace('%5C','/')
            with open(file, 'w') as fout:
                fout.write(file_contents)
            if 'index.' in file: continue
            else:
                fp.write('    ' + file + '\n')

if __name__=='__main__':
    convert_notebooks()