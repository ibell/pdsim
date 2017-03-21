import subprocess, sys, glob, os

template = '''
***************************
Converted IPython notebooks
***************************

.. note:: These static pages were auto-generated from the IPython notebooks residing `here on github <https://github.com/ibell/pdsim/tree/master/doc/notebooks>`_ They can be also be `viewed with nbviewer <http://nbviewer.ipython.org/github/ibell/pdsim/tree/master/doc/notebooks/>`_

'''

def convert_notebooks():
    
    thisdir = os.path.abspath(os.path.dirname(__file__))

    # Convert the notebooks
    subprocess.check_call('jupyter nbconvert --to rst *.ipynb', shell = True, stdout = sys.stdout, stderr = sys.stderr, cwd = thisdir)

    # Write the hidden table to make sphinx happy
    with open(thisdir+'/index.rst', 'w') as fp:
        fp.write(template)
        fp.write('.. toctree::\n\n')
        for file in glob.glob(thisdir+'/*.rst'):
            # Read in the converted RST file
            with open(file, 'r') as fin:
                file_contents= fin.read()
            # Convert a character
            file_contents = file_contents.replace('%5C','/')
            # Write back to file
            with open(file, 'w') as fout:
                fout.write(file_contents)
            # Write the TOC into the index file again
            if 'index.' in file: continue
            else:
                fp.write('    ' + os.path.split(file)[1].rsplit('.rst',1)[0] + '\n')

if __name__=='__main__':
    convert_notebooks()