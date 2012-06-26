Building and Using PDSim
**************************

Get PDSim
===========

The svn repo can be checked out by running a command like::

    svn checkout svn://svn.code.sf.net/p/pdsim/code/trunk pdsim-code
    
Or if you have read/write abilities::

    svn checkout --username=ibell svn+ssh://ibell@svn.code.sf.net/p/pdsim/code/trunk pdsim-code
    
where you should replace both ``ibell`` with your username.

Install PDSim
===============

For Windows development, the use of the `Python(x,y) <http://www.pythonxy.com/>`_ distribution is highly recommended.  It includes the following packages that are required for the use of ``PDSim``::

    Packages:
        cython
        scipy
        numpy
        matplotlib
        nose
    Other packages:
        swig
        mingw
    
A final pre-requisite for the use of ``PDSim`` is `CoolProp <http://coolprop.sf.net>`_.  Development of CoolProp is also ongoing, so the releases of CoolProp are likely not recent enough for the use of ``PDSim``.  Therefore it is recommended to build CoolProp from source.  This is pretty easy.  Use TortoiseSVN to checkout the Coolprop sources, or at the command line::

    svn co https://coolprop.svn.sourceforge.net/svnroot/coolprop coolprop
    
Then ``cd`` into the folder and run::

    python setup.py install
    
Profiling
---------

I use a couple of different tools for profiling.  One is RunSnakeRun.  Install RunSnakeRun using::

    easy_install SquareMap RunSnakeRun
    
and line_profiler, which can be obtained from http://pypi.python.org/pypi/line_profiler .  Open the zip file, go to the folder and run::

    python setup.py install

Use PDSim
===========
It is recommended to use `Eclipse <http://www.eclipse.org/downloads/>`_ (pick the Eclipse IDE for C/C++ development because it is the smallest) to do the development.  Once Eclipse is installed, you will want the Pydev plugin.  Once Eclipse is open, go to the menu Help-->Install New Software... Click *Add...* and add http://pydev.org/updates to the sources.  Then go back and you should be able to install pydev.  Don't install mylyn integration.

SciTE is also nice for doing python development.  Here are the user options I use and recommend::

    tabsize=4
    indent.size=4
    use.tabs=0
    wrap=1
    minimize.to.tray=0
    open.dialog.in.file.directory=1
    buffers=40
    statusbar.visible=1
    title.full.path=1
    
Once PDSim is set up, run the setup.py in the ``trunk`` folder.  This will compile some of the python files using Cython in order to get large improvements in speed.  There is stil more work to be done on this front. 
