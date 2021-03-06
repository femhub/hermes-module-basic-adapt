Hermes Module Basic Adapt
=========================

This is an adaptive version of the module Basic. 

Equation
--------

 -div(c1 \nabla u) + (c2, c3) \cdot \nabla u + c4 u = c5 

Here:

 * c1 ... equation parameter, element-wise constant
 * c2 ... equation parameter, element-wise constant
 * c3 ... equation parameter, element-wise constant
 * c4 ... equation parameter, element-wise constant
 * c5 ... equation parameter, element-wise constant

Boundary Conditions
-------------------

Dirichlet with piecewise-constant values (u = const)
Neumann with piecewise-constant normal derivatives (du/dn = const)
Newton with piecewise-constant parameters (const_1 u + du/dn = const_2)

Build the Module
----------------

For this, install the Hermes library first - only the 2D version 
without tutorial, examples and benchmarks is enough. Change dir to 
your local Hermes repository, for example /home/pavel/repos/hermes/. 
Add to your CMake.vars file the lines

set(WITH_H1D NO)
set(WITH_H3D NO)
set(H2D_WITH_TUTORIAL NO)
set(H2D_WITH_EXAMPLES NO)
set(H2D_WITH_BENCHMARKS NO)

Remove CMakeChache.txt and rerun "cmake ."

Then, either type "sudo make install" to have the Hermes library 
installed system-wide (typically into /usr/local/), or add into your 
CMake.vars file a line::

   set(CMAKE_INSTALL_PREFIX    /home/pavel/build/hermes)

to have Hermes installed into /home/pavel/build/hermes (adjust
the directory name to your needs). If using the latter option,
you do not need to use sudo. More details about the installation 
process can be found in the INSTALLATION section of the file 
CMake.vars.example in the Hermes repository.

Next, change dir to the directory with your module (this directory),
and in the file CMake.vars set the location of the Hermes library
by adding the line::

    set(HERMES_ROOT /home/pavel/build/hermes)

Then type::

    cmake .
    make

Please note that this module is dependent at the source level on module Basic.
Source code of module Basic is located in the directory src/hermes_basic and 
should be kept up-to-date with the development of that module.

Run the Module on C++ Level
---------------------------

C++ sources are located in the directory src/. Change dir to the directory 
src/ and run the module using::

    ./module-basicadapt model.cfg

The file model.cfg is a text file that emulates input from a GUI. You can 
change the parameters there at your will.


Run the Module on Python Level
------------------------------

If you have uncommented the line 

#  set(PYTHON_INSTALL_PATH   USE_SYSTEM_PYTHON_DIRECTORY)

in the file CMake.vars in your local Hermes repository, 
then you may delete the lines
::

    import sys
    sys.path.append("/home/pavel/build/hermes/lib/python")

in the file module-basic.py. Otherwise, adjust the line to point to the 
location of Hermes Python wrappers which must be compatible with your 
definition of HERMES_ROOT above.

Then run the module by typing::

    python module-basicadapt.py

The file module-basicadapt.py contains a set of parameters analogous to those
which on C++ level are in the file model.cfg. The user can change these
parameters arbitrarily. 

Run the Module in the Online Lab SDK
------------------------------------

To be completed.
