.. _install:

============
Installation
============

-------------
Prerequisites
-------------

In order to compile OpenMC, you will need to install Portable, Extensible
Toolkit for Scientific Computation, PETSC_.  CMFD is based on Fortran 95/2003
and therefore it is recommended to you the latest version of whatever compiler
you choose.  Currently, only gfortran_ has been tested.  It is recommended that
you use version 4.5.0 or above.

If you are using Debian or a Debian derivative such as Ubuntu_, you can install
the gfortran compiler using the following command::

    sudo apt-get install gfortran

.. _PETSC: http://www.mcs.anl.gov/petsc/petsc-as
.. _gfortran: http://gcc.gnu.org/wiki/GFortran
.. _Ubuntu: http://www.ubuntu.com

-------------------
PETSC Configuration
-------------------

Discuss petsc configuration (along with PETSC_DIR AND PETSC_ARCH env vars

------------------
CMFD Configuration
------------------

List CMFD configuration here

---------
Compiling
---------

To compile the code, run the following commands from within the root directory 
for CMFD::

     cd CMFD/src
     make all

This will build an exectable named ``cmfd``.
