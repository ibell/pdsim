Installation and Configuration
******************************

Windows
=======

#. Run the web installer for Microsoft Visual Studio 2008 Express from `VS2008Express <http://go.microsoft.com/?linkid=7729279>`_

#. Install the most recent version of MSYSGit from `full MSYSGit installs <http://code.google.com/p/msysgit/downloads/list?can=2&q=%22Full+installer+for+official+Git+for+Windows%22>`_.  Accept all the defaults.

#. Install TortoiseGit from `TortoiseGit installs <http://code.google.com/p/tortoisegit/wiki/Download>`_ (pick the right version for your operating system)

#. Download python using Anaconda package: `Download link <http://continuum.io/downloads>`_.  Download the python 2.7 64-bit version of anaconda - the python 2.7 32-bit version should be fine too.  Do not use a python 3.x version!!  Run the installer.

#. Check that when you go to a command prompt, you get output something like::

    C:\Users\Belli>python
    Python 2.7.2 (default, Jun 12 2011, 15:08:59) [MSC v.1500 32 bit (Intel)] on win32
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import scipy
    >>> scipy.__version__
    '0.11.0'
    >>>

#. Ok, python has been successfully installed. 
    
#. Now we are going to collect the source code for PDSim and CoolProp.  In Windows Explorer, go to a folder where you want to put the source code for PDSim and CoolProp.  Right-click and select "Git Clone..."

#. Use the URL ``https://github.com/ibell/pdsim``, and select the recursive option (which will force it to also install CoolProp and quantities as submodules of the PDSim repository)

#. Go into the folder you just created.  Run the script setup_coolprop.py (double-click on it), or at the command prompt::

    python setup_coolprop.py
    
#. Install PDSim using::

    python setup.py install
    
#. There are a few other dependencies that must be met before the GUI can be used that do not come standard with Anaconda distribution.  At the command line do::

    conda install wxpython h5py pyyaml
    
#. If you start a command prompt and run the command ``python``, you can try to import each of the required python packages using the code below.  This will allow you to check whether all the necessary files are included, and where they have been installed to.  If none of the following commands give you errors, you are ready to proceed.

.. ipython::

    In [0]: import CoolProp,matplotlib,Cython,PDSim,wx,numpy,scipy,yaml
    
    In [0]: print CoolProp.__file__; print CoolProp.__version__
    
    In [0]: print matplotlib.__file__; print matplotlib.__version__
    
    In [0]: print Cython.__file__; print Cython.__version__
    
    In [0]: print PDSim.__file__; print PDSim.__version__
    
    In [0]: print wx.__file__; print wx.version()
    
    In [0]: print numpy.__file__; print numpy.__version__
    
    In [0]: print scipy.__file__; print scipy.__version__
    
    In [0]: print yaml.__file__; print yaml.__version__

#. Now go into the examples folder, and try to run the file simple_example.py, you should get output like::

    C:\Users\Belli\Documents\Code\pdsim-git\examples>python simple_example.py
    Number of steps taken 7000
    Elapsed time for cycle is 0.4413 s
    Mass flow difference 2.03524833559  %
    ===========
    || # 000 ||
    ===========
    ||................|......@..........................|| energy balance  -0.0313473302944 [399.876899321797]
    ||................|.....@...........................|| discharge state 0.0225653554386 445.300857982
    ||................|...........@.....................|| cycle-cycle     0.145381939657
    Number of steps taken 7000
    Elapsed time for cycle is 0.437784 s
    Mass flow difference 0.324805392183  %
    ===========
    || # 001 ||
    ===========
    ||................|......@..........................|| energy balance  -0.0312974294669 [361.27298241958306]
    ||................|@................................|| discharge state 0.00386685157006 453.064795878
    ||................|......@..........................|| cycle-cycle     0.024111213808
    Number of steps taken 7000
    Elapsed time for cycle is 0.439191 s
    Mass flow difference -0.0570907915056  %
    ===========
    || # 002 ||
    ===========
    ||................|....@............................|| energy balance  -0.0156487147334 [341.97102396847896]
    ||............@...|.................................|| discharge state -0.000716737604005 451.684013021
    ||................|@................................|| cycle-cycle     0.00426059998124
    Number of steps taken 7000
    Elapsed time for cycle is 0.499563 s
    Mass flow difference -0.0664255877057  %
    ===========
    || # 003 ||
    ===========
    ||................|..@..............................|| energy balance  -0.00782435736672 [332.3200447429269]
    ||............@...|.................................|| discharge state -0.000807748849267 450.151035073
    ||................|.@...............................|| cycle-cycle     0.0048139260336
    Number of steps taken 7000
    Elapsed time for cycle is 0.52113 s
    Mass flow difference -0.0374442760086  %
    ===========
    || # 004 ||
    ===========
    ||................|@................................|| energy balance  -0.00391217868336 [327.49455513015084]
    ||...........@....|.................................|| discharge state -0.000497721396576 449.251208618
    ||................@.................................|| cycle-cycle     0.00279279716226
    Number of steps taken 7000
    Elapsed time for cycle is 0.49191 s
    Mass flow difference -0.0184814204442  %
    ===========
    || # 005 ||
    ===========
    ||...............@|.................................|| energy balance  -0.00195608934168 [325.08181032376285]
    ||.........@......|.................................|| discharge state -0.000286132835949 448.779806416
    ||..............@.|.................................|| cycle-cycle     0.0014388268967
    Ntheta is 7001
    mdot*(h2-h1),P-v,Qamb 0.0874497858735 0.07491487626 -0.0119560893417
    Mass flow rate is 0.503773487596 g/s
    Volumetric efficiency is 88.6166506787 %

#. Ok good, PDSim is working!
    
#. Go to the GUI folder from the root of the source.  Double-click on PDSimGUI.py to start

#. Press F5 to run the default scroll compressor

#. Wait about 80 seconds for it to finish

Linux and OSX
=============

The procedure is nearly identical on linux and OSX, apart from the fact that you do not need to install git or Microsoft Visual Studio.  Use the anaconda installer to get python 2.7 64-bit, follow the windows instructions otherwise

Update source code
==================
When the code for PDSim is updated, all that you have to do to update the code is to go to the top-level folder where you installed PDSim, right-click and select "TortoiseGit->Pull..." and then OK

If you want to update CoolProp to the most recent code, you can use TortoiseGit.  Go to the root folder of the cloned git files, right click and select "TortoiseGit->Submodule update...". Make sure the following options are selected::

* Initialize submodules(--init)
* Merge
* Remote tracking branch

Press Ok button.

Alternatively, you can do the update from the command line::

    git.exe submodule update --init --merge --remote "externals/coolprop"
    
See also `StackOverflow question <http://stackoverflow.com/questions/16058917/pulling-git-submodules-with-tortoisegit>`_

Profiling
---------

I use a couple of different tools for profiling.  One is RunSnakeRun.  Install RunSnakeRun using::

    easy_install SquareMap RunSnakeRun
    
and line_profiler, which can be obtained from http://pypi.python.org/pypi/line_profiler .  Open the zip file, go to the folder and run::

    python setup.py install

Uninstallation
==============

To uninstall PDSim, go to the site-packages folder corresponding to the installation of python (probably c:\\Python27\\Lib\\site-packages), delete the folder PDSim.  You might want to also delete any files like ``PDSim-x.x.x-py2.7.egg-info`` where ``x`` are numbers.  For a thorough uninstallation, you might also want to remove the ``build`` folder in the directory where you cloned the git files

Using conda environments
========================

If you have multiple versions of python or PDSim floating around, it can be useful to use conda to create virtual environments that encapsulate the desired versions of each of the pieces.  This can be easily carried out at the command line.  For instance, we might create an environment (named ``pdsim_stable``) with the most up to date version of PDSim and CoolProp version 4.2.5.  This can be achieved using a command like::
    
    C:\Users\XXXX>c:\Miniconda32bit\Scripts\conda.exe create -n pdsim_stable python=2.7 matplotlib numpy scipy h5py cython pip wxpython
    Fetching package metadata: ..
    Solving package specifications: .............
    Package plan for installation in environment c:\Miniconda32bit\envs\pdsim_stable:

    The following packages will be downloaded:

        package                    |            build
        ---------------------------|-----------------
        cython-0.21                |           py27_0         1.6 MB
        h5py-2.3.1                 |       np19py27_0         1.2 MB
        matplotlib-1.4.0           |       np19py27_0        41.7 MB
        numpy-1.9.0                |           py27_0        14.2 MB
        pytz-2014.7                |           py27_0         169 KB
        scipy-0.14.0               |       np19py27_0        33.1 MB
        setuptools-5.8             |           py27_0         729 KB
        six-1.8.0                  |           py27_0          15 KB
        ------------------------------------------------------------
                                               Total:        92.8 MB

    The following packages will be linked:

        package                    |            build
        ---------------------------|-----------------
        cython-0.21                |           py27_0   hard-link
        dateutil-2.1               |           py27_2   hard-link
        h5py-2.3.1                 |       np19py27_0   hard-link
        matplotlib-1.4.0           |       np19py27_0   hard-link
        numpy-1.9.0                |           py27_0   hard-link
        pip-1.5.6                  |           py27_0   hard-link
        pyparsing-2.0.1            |           py27_0   hard-link
        pyqt-4.10.4                |           py27_0   hard-link
        python-2.7.8               |                0   hard-link
        pytz-2014.7                |           py27_0   hard-link
        scipy-0.14.0               |       np19py27_0   hard-link
        setuptools-5.8             |           py27_0   hard-link
        six-1.8.0                  |           py27_0   hard-link
        wxpython-3.0               |           py27_0   hard-link

    Proceed ([y]/n)?
    
when you say yes, miniconda will fetch the required versions of the software packages, as in::

    Fetching packages ...
    cython-0.21-py 100% |###############################| Time: 0:00:02 654.50 kB/s
    h5py-2.3.1-np1 100% |###############################| Time: 0:00:01   1.27 MB/s
    matplotlib-1.4 100% |###############################| Time: 0:00:29   1.49 MB/s
    numpy-1.9.0-py 100% |###############################| Time: 0:00:10   1.47 MB/s
    pytz-2014.7-py 100% |###############################| Time: 0:00:00 362.89 kB/s
    scipy-0.14.0-n 100% |###############################| Time: 0:00:21   1.59 MB/s
    setuptools-5.8 100% |###############################| Time: 0:00:01 738.79 kB/s
    six-1.8.0-py27 100% |###############################| Time: 0:00:00 181.98 kB/s
    Extracting packages ...
    [      COMPLETE      ] |#################################################| 100%
    Linking packages ...
    [      COMPLETE      ] |#################################################| 100%
    #
    # To activate this environment, use:
    # > activate pdsim_stable
    #

To activate this new environment, you do::

    C:\Users\XXXX>c:\Miniconda32bit\Scripts\activate pdsim_stable
    Activating environment "pdsim_stable"...

    [pdsim_stable] C:\Users\XXXX>
    
Normally you would just do ``activate pdsim_stable``, but on my machine, the default miniconda is 64-bit and it gets all confused if you don't call the activate script directly.  Once the environment has been populated, you can pull in the remaining packages using pip::

    [pdsim_stable] C:\Users\Belli>pip install CoolProp==4.2.5 cx_Freeze glob2
    Downloading/unpacking CoolProp==4.2.5
    Downloading/unpacking cx-Freeze
    Downloading/unpacking glob2
      Downloading glob2-0.4.1.tar.gz
      Running setup.py (path:c:\users\belli\appdata\local\temp\pip_build_Belli\glob2\setup.py) egg_info for package glob2

    Installing collected packages: CoolProp, cx-Freeze, glob2
      Running setup.py install for glob2

    Successfully installed CoolProp cx-Freeze glob2
    Cleaning up...
    
And you can check that the right things are setup by doing::

    [pdsim_stable] C:\Users\Belli>python
    Python 2.7.8 |Continuum Analytics, Inc.| (default, Jul  2 2014, 15:13:35) [MSC v.1500 32 bit (Intel)] on win32
    Type "help", "copyright", "credits" or "license" for more information.
    Anaconda is brought to you by Continuum Analytics.
    Please check out: http://continuum.io/thanks and https://binstar.org
    >>> import CoolProp
    >>> CoolProp.__file__
    'c:\\Miniconda32bit\\envs\\pdsim_stable\\lib\\site-packages\\CoolProp\\__init__.pyc'
    
The path should be to a file in your ``envs`` folder of the miniconda installation.

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
    split.vertical=0
    title.full.path=1
    # one instance of SciTE only
    check.if.already.open=1
    are.you.sure.on.reload=1

    font.base=$(font.monospace)
    font.small=$(font.monospace)
    font.comment=$(font.monospace)
    font.text=$(font.monospace)
    font.text.comment=$(font.monospace)
    font.embedded.base=$(font.monospace)
    font.embedded.comment=$(font.monospace)
    font.vbs=$(font.monospace) 
