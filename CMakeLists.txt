CMAKE_MINIMUM_REQUIRED(VERSION 3.8)

PROJECT(soars LANGUAGES CUDA CXX C)

SET(CMAKE_VERBOSE_MAKEFILE ON)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${CMAKE_BINARY_DIR}/config.h )

# Check for build configuration
IF(NOT CMAKE_BUILD_TYPE)
    MESSAGE(STATUS "No build configuration specified, defaulting to Release")
    SET(CMAKE_BUILD_TYPE "Release")
ENDIF()

include(CheckLanguage)
check_language(CUDA)
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_35 -cudart=shared -rdc=true")

# Set compiler flags
IF (CMAKE_COMPILER_IS_GNUCXX)
    MESSAGE(STATUS "Setting general compiler flags for detected compiler: gnu-g++")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -pthread -ffast-math -I/usr/include/hdf5/serial -L/usr/lib/x86_64-linux-gnu -lhdf5 -lhdf5_hl")

    IF(CMAKE_BUILD_TYPE MATCHES "Release")

        MESSAGE(STATUS "Setting g++ flags for Release configuration")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")   ## Optimize
        SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -s")  ## Strip binary

        ELSEIF(CMAKE_BUILD_TYPE MATCHES "Debug")
            MESSAGE(STATUS "Setting g++ flags for Debug configuration")

    ENDIF()

ELSE()
    MESSAGE(STATUS "Setting general compiler flags")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -pthread -ffast-math -I/usr/local/include/")

    IF(CMAKE_BUILD_TYPE MATCHES "Release")

        MESSAGE(STATUS "Setting flags for Release configuration")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")   ## Optimize
        SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -s")  ## Strip binary

        ELSEIF(CMAKE_BUILD_TYPE MATCHES "Debug")
            MESSAGE(STATUS "Setting flags for Debug configuration")

    ENDIF()

ENDIF()

# SOARS is multithreaded, boost needs to be informed
ADD_DEFINITIONS(-D_REENTRANT)

ADD_SUBDIRECTORY( fftwcpp )
ADD_SUBDIRECTORY( sgp4 )
ADD_SUBDIRECTORY( pugixml )
ADD_SUBDIRECTORY( src )
ADD_SUBDIRECTORY( utilities )
ADD_SUBDIRECTORY( test )
