Building and Using PDSim modules
********************************

Step-By-Step on Windows
=======================

#. Run the web installer for Microsoft Visual Studio 2008 Express from `VS2008Express <http://go.microsoft.com/?linkid=7729279>`_

#. Install the most recent version of MSYSGit from `full MSYSGit installs <http://code.google.com/p/msysgit/downloads/list?can=2&q=%22Full+installer+for+official+Git+for+Windows%22>`_.  Accept all the defaults.

#. Install TortoseGit from `TortoiseGit installs <http://code.google.com/p/tortoisegit/wiki/Download>`_ (pick the right version for your operating system)

#. Install python using Python (x,y) python installer.  Download the full installer from `Python(x,y) <https://code.google.com/p/pythonxy/wiki/Downloads?tm=2>`_ website.  In the installer, make sure to check Cython in the python section

#. Check that when you go to a command prompt, you get output like::

    C:\Users\Belli>python
    Python 2.7.2 (default, Jun 12 2011, 15:08:59) [MSC v.1500 32 bit (Intel)] on win32
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import scipy
    >>> scipy.__version__
    '0.11.0'
    >>>

#. Ok, python has been successfully installed.

#. Install the remaining dependencies of PDSim using the command::

    pip install -U 
    
You have two steps remaining, building 



Get PDSim
===========

The git repo can be cloned by running a command like::

    git clone git://git.code.sf.net/p/pdsim/git pdsim-git
    
Or if you have read/write abilities::

    git clone ssh://ibell@git.code.sf.net/p/pdsim/git pdsim-git
    
where you should replace ``ibell`` with your sourceforge username.

.. _install-PDSim:

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
    
A final pre-requisite for the use of ``PDSim`` is `CoolProp <http://coolprop.sf.net>`_.  Development of CoolProp is also ongoing, so the releases of CoolProp are likely not recent enough for the use of ``PDSim``.  

If you are on windows and are using 32-bit python 2.7, you can use nightly CoolProp executables from `nightly CoolProp build <https://sourceforge.net/projects/coolprop/files/CoolProp/Nightly/>`_.  Beware, these are development builds, but are likely stable enough.

Failing that, you will have to build from source.  This is pretty easy.  Use TortoiseSVN to checkout the Coolprop sources, or at the command line::

    svn co https://coolprop.svn.sourceforge.net/svnroot/coolprop coolprop
    
Then ``cd`` into the folder and run::

    python setup.py install

Profiling
---------

I use a couple of different tools for profiling.  One is RunSnakeRun.  Install RunSnakeRun using::

    easy_install SquareMap RunSnakeRun
    
and line_profiler, which can be obtained from http://pypi.python.org/pypi/line_profiler .  Open the zip file, go to the folder and run::

    python setup.py install

.. _Use-PDSim:

Use PDSim
=========
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
