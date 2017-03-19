Installation and Configuration
******************************

Windows
=======

#. Run the web installer for Microsoft Visual Studio 2008 Express from `VS2008Express <http://go.microsoft.com/?linkid=7729279>`_

#. Install the most recent version of MSYSGit from `full MSYSGit installs <http://code.google.com/p/msysgit/downloads/list?can=2&q=%22Full+installer+for+official+Git+for+Windows%22>`_.  Accept all the defaults.

#. Install TortoiseGit from `TortoiseGit installs <http://code.google.com/p/tortoisegit/wiki/Download>`_ (pick the right version for your operating system)

#. Download python using miniconda package: `Download link <http://conda.pydata.org/miniconda.html>`_.  Download the python 2.7 32-bit version of miniconda - the python 2.7 64-bit version should be fine too, but to compile you will require the professional version of Visual Studio.  So just stick with the 32-bit version of miniconda.  Do not use a python 3.x version!!  

#. Run the installer. In the setup, if you are doing a clean install, it is safe to select the option "Add Anaconda to the system PATH environmental variable".  Otherwise, selecting this option will make this the default conda installation on your computer, which may or may not be what you want.  If you want Miniconda to live peaceably with an existing miniconda (64-bit?), make sure this option is unselected.

#. Populate the python installation with necessary packages.  At the command prompt, do::

    conda install matplotlib numpy scipy h5py cython pip wxpython pyyaml
  
  If you installed Miniconda somewhere else (and/or Miniconda/Scripts is not on the PATH), you might need to give the full path to ``conda``, which would be something like ``c:\Miniconda32bit\Scripts\conda`` on my machine

#. Check that when you go to a command prompt, you get output something like::

    C:\Users\Belli>python
    Python 2.7.2 (default, Jun 12 2011, 15:08:59) [MSC v.1500 32 bit (Intel)] on win32
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import scipy
    >>> scipy.__version__
    '0.11.0'

#. Ok, python has been successfully installed.
    
#. Now we are going to collect the source code for PDSim and CoolProp.  In Windows Explorer, go to a folder where you want to put the source code for PDSim and CoolProp.  Right-click and select "Git Clone..."

#. Use the URL ``https://github.com/ibell/pdsim``, and select the recursive option (which will force it to also install CoolProp and quantities as submodules of the PDSim repository)

#. Go into the folder you just created.  Run the script setup_coolprop.py (double-click on it), or at the command prompt::

    python setup_coolprop.py
    
#. Install PDSim using::

    python setup.py install
    
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

#. Now go into the examples folder, and start ipython, you should get output something like::

.. ipython::

    In [0]: %run 'examples/simple_example.py'

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

If you have multiple versions of python or PDSim floating around, it can be useful to use conda to create conda environments that encapsulate the desired versions of each of the pieces.  This can be easily carried out at the command line.  For instance, we might create an environment (named ``pdsim_stable``) with the most up to date version of PDSim and its dependencies.  There is a file called ``RTDenvironment.yml`` in the root of the repository.  This file can be used to create a conda environment with::
    
    conda env create --name pdsim_stable --file RTDenvironment.yml

To activate this new environment on windows, you do::

    activate pdsim_stable

If on linux/OSX, you need to do::

    source activate pdsim_stable

To remove the environment you created without confirmation (thanks to ``-y``), do::

    conda env remove -y -n pdsim_stable

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
