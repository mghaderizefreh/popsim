# CMake project file for Random Regression

##################################################
# Define the project and the depencies that it has
##################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)
PROJECT(PopSim Fortran)

# Set the version
SET(VERSION 1.1.2)

# Add our local modlues to the module path
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Uncomment if it is required that Fortran 90 is supported
#IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
#    MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
#ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

# Set some options the user may choose
# Uncomment the below if you want the user to choose a parallelization library
#OPTION(USE_MPI "Use the MPI library for parallelization" OFF)
#OPTION(USE_OPENMP "Use OpenMP for parallelization" OFF)

# This INCLUDE statement executes code that sets the compile flags for DEBUG,
# RELEASE, and TESTING.  You should  review this file and make sure the flags 
# are to your liking.
INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake) 
# Locate and set parallelization libraries.  There are some CMake peculiarities
# taken care of here, such as the fact that the FindOpenMP routine doesn't know
# about Fortran.
INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)
# Setup the LAPACK libraries.  This also takes care of peculiarities, such as
# the fact the searching for MKL requires a C compiler, and that the results
# are not stored in the cache. 
INCLUDE(${CMAKE_MODULE_PATH}/SetUpLAPACK.cmake)

# There is an error in CMAKE with this flag for pgf90.  Unset it
GET_FILENAME_COMPONENT(FCNAME ${CMAKE_Fortran_COMPILER} NAME)
IF(FCNAME STREQUAL "pgf90")
    UNSET(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS)
ENDIF(FCNAME STREQUAL "pgf90")

############################################################
# Define the actual files and folders that make up the build
############################################################

# Define the executable name
SET(SIM simulate.exe)

# Define the library name
SET(UTILLIB util)
SET(BLUPLIB blup)
SET(RNLIB rnmodel)
SET(SORTLIB sort)
SET(RNGLIB rng)
SET(EVOLLIB evolution)
SET(MATHLIB math)

# Define some directories
SET(SRC ${CMAKE_SOURCE_DIR}/src)
SET(LIB1 ${CMAKE_SOURCE_DIR}/lib/util)
SET(LIB2 ${CMAKE_SOURCE_DIR}/lib/sort)
SET(LIB3 ${CMAKE_SOURCE_DIR}/lib/blup)
SET(LIB4 ${CMAKE_SOURCE_DIR}/lib/rnmodel)
SET(LIB5 ${CMAKE_SOURCE_DIR}/lib/rng)
SET(LIB6 ${CMAKE_SOURCE_DIR}/lib/evolution)
SET(LIB7 ${CMAKE_SOURCE_DIR}/lib/math)
SET(BIN ${CMAKE_SOURCE_DIR}/bin)
SET(TEST ${CMAKE_SOURCE_DIR}/tests)
SET(SRCUTIL ${SRC}/utils)
SET(SRCSORT ${SRC}/quickSort)
SET(SRCBLUP ${SRC}/blup)
SET(SRCRNMODEL ${SRC}/rnmodel)
SET(SRCRNG ${SRC}/rng)
SET(SRCEVOL ${SRC}/evolution)
SET(SRCMATH ${SRC}/math)
SET(MAIN ${SRC}/main)

# Have the .mod files placed in the lib folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${LIB1} ${LIB2} ${LIB3} ${LIB5} ${LIB6})

# The source for the UTILS library and have it placed in the lib folder
ADD_SUBDIRECTORY(${SRCUTIL} ${LIB1})
ADD_SUBDIRECTORY(${SRCSORT} ${LIB2})
ADD_SUBDIRECTORY(${SRCBLUP} ${LIB3})
ADD_SUBDIRECTORY(${SRCRNMODEL} ${LIB4})
ADD_SUBDIRECTORY(${SRCRNG} ${LIB5})
ADD_SUBDIRECTORY(${SRCEVOL} ${LIB6})
ADD_SUBDIRECTORY(${SRCMATH} ${LIB7})

# The source for the REML binary and have it placed in the bin folder
ADD_SUBDIRECTORY(${MAIN} ${BIN})

# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean 
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/distclean.cmake
)

enable_testing()

include(FindUnixCommands)

set(INP1 "${TEST}/input1")
#if (BASH)
#    add_test(
#        NAME correlated_RR_Reml
#        COMMAND ${BASH} -c "${BIN}/${analysis} < ${INP1};" WORKING_DIRECTORY ${TEST}
#    )
#else()
#    if (WIN32)
#        add_test(
#            NAME correlated_RR_Reml
#            COMMAND ${CMAKE_COMMAND} -E chdir ${RR} $ENV{ComSpec} /c "${RR} < ${INP1}" WORKING_DIRECTORY ${TEST}
#        )
#        add_test(
#            NAME uncorrelated_RR_Reml
#            COMMAND ${CMAKE_COMMAND} -E chdir ${RR} $ENV{ComSpec} /c "${RR} < ${INP2}" WORKING_DIRECTORY ${TEST}
#        )
#        add_test(
#            NAME single
#            COMMAND ${CMAKE_COMMAND} -E chdir ${RR} /c "${ST} < ${INP3}" WORKING_DIRECTORY ${TEST}
#        )
#    else()
#        message(FATAL_ERROR "Unknown shell command for ${CMAKE_HOST_SYSTEM_NAME}")
#    endif()
#endif()
