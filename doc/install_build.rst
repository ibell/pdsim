Building and Using PDSim modules
********************************

Step-By-Step on Windows
=======================

#. Run the web installer for Microsoft Visual Studio 2008 Express from `VS2008Express <http://go.microsoft.com/?linkid=7729279>`_

#. Install the most recent version of MSYSGit from `full MSYSGit installs <http://code.google.com/p/msysgit/downloads/list?can=2&q=%22Full+installer+for+official+Git+for+Windows%22>`_.  Accept all the defaults.

#. Install TortoiseGit from `TortoiseGit installs <http://code.google.com/p/tortoisegit/wiki/Download>`_ (pick the right version for your operating system)

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
    
#. Now we are going to collect the source code for PDSim and CoolProp.  In Windows Explorer, go to a folder where you want to put the source code for PDSim and CoolProp.  Right-click and select "Git Clone..."

#. Use the URL ``https://github.com/ibell/pdsim``, and select the recursive option (which will force it to also install CoolProp as a submodule of the PDSim repository)

#. Go into the folder you just created.  Run the script setup_coolprop.py (double-click on it, or at the command prompt::

    python setup_coolprop.py
    
#. Install PDSim using::

    python setup.py install

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
    
#. Install PyYAML using the installer for win32-py2.7 at `PyYAML installers <http://pyyaml.org/wiki/PyYAML>`_
    
#. Go to the GUI folder from the root of the source.  Double-click on PDSimGUI.py to start

#. Press F5 to run the default scroll compressor

#. Wait about 80 seconds for it to finish

Update source code
===========
When the code for PDSim is updated, all that you have to do to update the code is to go to the top-level folder where you installed PDSim, right-click and select "TortoiseGit->Pull..." and then OK

If you want to update CoolProp to the most recent code (XXXXXXXXXXXXXXXXXXX)

Profiling
---------

I use a couple of different tools for profiling.  One is RunSnakeRun.  Install RunSnakeRun using::

    easy_install SquareMap RunSnakeRun
    
and line_profiler, which can be obtained from http://pypi.python.org/pypi/line_profiler .  Open the zip file, go to the folder and run::

    python setup.py install

Uninstallation
==============

To uninstall PDSim, go to the site-packages folder corrresponding to the installation of python (probably c:\\Python27\\Lib\\site-packages), delete the folder PDSim.  You might want to also delete any files like ```PDSim-x.x.x-py2.7.egg-info``` where ``x`` are numbers.  For a thorough uninstallation, you might also want to remove the ``build`` folder in the directory where you cloned the git files

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
