#-------- PROJECT -------------------------------------------------------------#

PROJECT( CMFD )
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 FATAL_ERROR )
SET(CMFD_DIR ${CMAKE_CURRENT_SOURCE_DIR})
ENABLE_LANGUAGE(Fortran)

#-------- CMAKE MODULES -------------------------------------------------------#

# CMake files included with CMFD 
# This is the tool that finds modules and located in CMFD/tools/CMake
LIST( APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/tools/CMake" )

#-------- EXTERNAL LIBRARIES --------------------------------------------------#

# Find PACAKAGE calls the Petsc find routine
FIND_PACKAGE( PETSc REQUIRED )
INCLUDE_DIRECTORIES( ${PETSC_INCLUDES} )
ADD_DEFINITIONS( ${PETSC_DEFINITIONS} )

# Add XML library
# -- include directory needed since there are .mod files to be used
INCLUDE_DIRECTORIES( "${PROJECT_SOURCE_DIR}/lib/gfortran" )
# define actual library
SET( XMLPARSE "/home/bherman/Documents/PhD/CMFD/lib/gfortran/libxmlparse.a" )

# ------------ DEFINE SOURCE CODE --------------------------------------------#
SET( CMFD_SRC
    src/cmfd_header.f90
    src/cmfd_execute.F90
    src/cmfd_utils.f90
    src/constants.f90
    src/global.f90
    src/main.f90
    src/timing.f90
    src/xml_cmfd.f90
)

# ------------ DEFINE EXECUTABLE ---------------------------------------------#
ADD_EXECUTABLE(  src/cmfd ${CMFD_SRC} )
TARGET_LINK_LIBRARIES( src/cmfd
                       ${XMLPARSE}
                       ${PETSC_LIBRARIES} )
