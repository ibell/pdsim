import subprocess, sys, glob

template = '''
***************************
Converted IPython notebooks
***************************

.. note:: These static pages were auto-generated from the IPython notebooks residing `here on github <https://github.com/ibell/pdsim/tree/master/doc/notebooks>`_ They can be also be `viewed with nbviewer <http://nbviewer.ipython.org/github/ibell/pdsim/tree/master/doc/notebooks/>`_
'''

# Convert the notebooks
subprocess.check_call('jupyter nbconvert --to rst *.ipynb', shell = True, stdout = sys.stdout, stderr = sys.stderr)

# Write the hidden table to make sphinx happy
with open('index.rst', 'w') as fp:
    fp.write(template)
    fp.write('.. toctree::\n\n')
    for file in glob.glob('*.rst'):
        with open(file, 'r') as fin:
            file_contents= fin.read()
        file_contents = file_contents.replace('%5C','/')
        with open(file, 'w') as fout:
            fout.write(file_contents)
        if 'index.' in file: continue
        else:
            fp.write('    ' + file + '\n')