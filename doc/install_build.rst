Installation and Configuration
******************************

.. warning::

	Having trouble with installation? Open an issue here and we'll get you going: https://github.com/ibell/pdsim/issues/new

Prerequisites
=============

git
---

**WINDOWS** :

#. Install the most recent version for windows from `git for windows <https://git-scm.com/download>`_.  Accept all the defaults.

#. (optional, but recommended) Install TortoiseGit from `TortoiseGit installs <http://code.google.com/p/tortoisegit/wiki/Download>`_ (pick the right version for your operating system).  Other solid graphical git interface options are VSCode or GitKraken

**LINUX/OSX** :

Your operating system should already come with git.  You should not need to do anything extra.  You might need to install g++ via your package manager (e.g., on debian-based distros: ``sudo apt-get install g++``)

Anaconda/Miniconda
------------------

#. Download python using miniconda package: `Download link <http://conda.pydata.org/miniconda.html>`_. Do not use a python 2.7 version!!

#. Run the installer. In the setup, if you are doing a clean install, it is safe to select the option "Add Anaconda to the system PATH environmental variable".  Otherwise, selecting this option will make this the default conda installation on your computer, which may or may not be what you want.  If you want Miniconda to live peaceably with an existing miniconda (64-bit?), make sure this option is unselected.

#. If your shell doesn't understand the command ``conda``, you might want to supercharge your shell by doing ``conda init`` to enable the direct use of the conda command without supplying the complete path.

Compiler
--------

**WINDOWS** :

Install the version of Visual Studio that matches your version of Python you would like to use, or install `the build tools <https://www.visualstudio.com/downloads/#build-tools-for-visual-studio-2019>`_

`More information here <https://wiki.python.org/moin/WindowsCompilers#Which_Microsoft_Visual_C.2B-.2B-_compiler_to_use_with_a_specific_Python_version_.3F>`_

**LINUX/OSX** :

Your operating system might already come with a compiler.  In that case, you should not need to do anything else.  On debian-base linux distros, you probably need to do something like::

    sudo apt install g++

On OSX, install Xcode, or install ``gcc`` via homebrew.

Windows
=======

Conda environments are useful to ensure that all the dependencies needed for a particular task (here to run code with PDSim).  They are recommended for most use cases with PDSim.

#. If you have multiple versions of python or PDSim floating around, it can be useful to use conda to create conda environments that encapsulate the desired versions of each of the pieces.  This can be easily carried out at the command line.  For instance, we might create an environment (named ``pdsim_stable``) with the most up to date version of PDSim and its dependencies.  There is a file called ``conda_environment.yml`` in the root of the repository.  This file can be used to create a conda environment with::
    
    conda env create --name pdsim_stable --file conda_environment.yml

To activate this new environment on windows, you do::

    activate pdsim_stable

If on linux/OSX, you need to do::

    source activate pdsim_stable

To remove the environment you created without confirmation (thanks to ``-y``), do::

    conda env remove -y -n pdsim_stable

#. Check that when you go to a command prompt, you get output something like::

    C:\Users\ian>python
    Python 3.7.6 | packaged by conda-forge | (default, Mar 23 2020, 22:22:21) [MSC v.1916 64 bit (AMD64)] on win32
    Type "help", "copyright", "credits" or "license" for more information.
    >>> quit()

#. Ok, python has been successfully installed.
    
#. Now we are going to collect the source code for PDSim.  

    Option A: If you installed TortoiseGit, In Windows Explorer, go to a folder where you want to put the source code for PDSim.  Right-click and select "Git Clone...". Use the URL ``https://github.com/ibell/pdsim``.
    Option B: At a command prompt, do: ``git clone --recursive https://github.com/ibell/pdsim``
    Option C: Use your git client of choice to clone pdsim following its instructions with the url ``https://github.com/ibell/pdsim``
    
#. Move into the pdsim folder, then install PDSim using::

    python setup.py install
    
#. If you start a command prompt and run the command ``python``, you can try to import each of the required python packages using the code below.  This will allow you to check whether all the necessary files are included, and where they have been installed to.  If none of the following commands give you errors, you are ready to proceed.

.. ipython::

    In [0]: import CoolProp,matplotlib,Cython,PDSim,wx,numpy,scipy,yaml
    
    In [0]: CoolProp.__version__, CoolProp.__file__
    
    In [0]: matplotlib.__version__, matplotlib.__file__
    
    In [0]: Cython.__version__, Cython.__file__
    
    In [0]: PDSim.__version__, PDSim.__file__
    
    In [0]: wx.version(), wx.__file__
    
    In [0]: numpy.__version__, numpy.__file__
    
    In [0]: scipy.__version__, scipy.__file__
    
    In [0]: yaml.__version__, yaml.__file__

#. As an alternative to the previous step, you can do: ``conda list`` which will list all the installed packages and their versions. The list can be rather long....

#. Now go into the examples folder and you you shold be able to run::

    python simple_example.py

#. Ok good, PDSim is working!
    
#. Go to the GUI folder from the root of the source. In the file system explorer, double-click on PDSimGUI.py to start, or at the command prompt, move to the GUI folder, and type: ``python PDSimGUI.py``

#. Press F5 to run the default scroll compressor

#. Wait about 80 seconds for it to finish

Update source code
==================
When the code for PDSim is updated, all that you have to do to update the code is to go to the top-level folder where you installed PDSim, right-click and select "TortoiseGit->Pull..." and then OK

If you want to update CoolProp to the most recent code, you can use TortoiseGit.  Go to the root folder of the cloned git files, right click and select "TortoiseGit->Submodule update...". Make sure the following options are selected::

* Initialize submodules(--init)
* Merge
* Remote tracking branch

Press Ok button.

Alternatively, you can do the update from the command line::

    git submodule update --init --merge --remote "externals/coolprop"
    
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

Manual installation without conda environment
========================

!! This approach is not recommended.  Better to use the environment file. !!

Populate the root conda environment installation the with necessary packages (or see below about using conda environments).  At the command prompt, do::

    conda install matplotlib numpy scipy h5py cython pip wxpython pyyaml
  
If you installed Miniconda somewhere else (and/or Miniconda/Scripts is not on the PATH), you might need to give the full path to ``conda``, which would be something like ``c:\Miniconda3\Scripts\conda`` on my machine

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
