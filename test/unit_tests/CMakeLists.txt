INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} )
INCLUDE_DIRECTORIES( ${CMAKE_BINARY_DIR} )
INCLUDE_DIRECTORIES( ${soars_SOURCE_DIR}/src )
INCLUDE_DIRECTORIES( ${soars_SOURCE_DIR}/fftwcpp )
INCLUDE_DIRECTORIES( ${soars_SOURCE_DIR}/sgp4 )

########### noisegen_test ###############

SET(noisegen_test_SRCS
   noisegen_test.cu
   ${soars_SOURCE_DIR}/src/rsnoise.cu
   ${soars_SOURCE_DIR}/src/rsnoise.cuh
   ${soars_SOURCE_DIR}/src/rsdsp.cu
   ${soars_SOURCE_DIR}/src/rsdsp.cuh
   ${soars_SOURCE_DIR}/src/rspython.cu
   ${soars_SOURCE_DIR}/src/rspython.cuh
   ${soars_SOURCE_DIR}/src/rsparameters.cu
   ${soars_SOURCE_DIR}/src/rsparameters.cuh
   ${soars_SOURCE_DIR}/src/rsgeometry.cu
   ${soars_SOURCE_DIR}/src/rsgeometry.cuh
   ${soars_SOURCE_DIR}/src/rsdebug.cu
   ${soars_SOURCE_DIR}/src/rsdebug.cuh
   )

ADD_EXECUTABLE(noisegen_test ${noisegen_test_SRCS})

# Link to fftwcpp libraries
TARGET_LINK_LIBRARIES(noisegen_test fftwcpp)

# See if the FFTW3 libraries are available
FIND_LIBRARY(soars_LIB_FFTW3 NAMES fftw3)

IF (soars_LIB_FFTW3)
   TARGET_LINK_LIBRARIES(noisegen_test ${soars_LIB_FFTW3})
ENDIF (soars_LIB_FFTW3)

# Check if we are building on a multithreaded FFTW
FIND_LIBRARY(soars_LIB_FFTW3_THREADS NAMES fftw3_threads)

IF (soars_LIB_FFTW3_THREADS)
   TARGET_LINK_LIBRARIES(noisegen_test ${soars_LIB_FFTW3_THREADS})
ENDIF (soars_LIB_FFTW3_THREADS)

# Include code to check for python libs

INCLUDE(FindPythonLibs)

IF (PYTHONLIBS_FOUND)
   TARGET_LINK_LIBRARIES(noisegen_test ${PYTHON_LIBRARIES})
   INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
ENDIF (PYTHONLIBS_FOUND)

# Link to Boost Thread and System

IF (WIN32)
   SET(Boost_USE_STATIC_LIBS ON) ## Link Boost statically
ELSE ()
   SET(Boost_USE_STATIC_LIBS OFF) ## Link Boost with shared libs
ENDIF ()
FIND_PACKAGE(Boost COMPONENTS system thread REQUIRED)

MESSAGE(STATUS "Linking libraries for Boost: " ${Boost_LIBRARIES})

IF (Boost_LIBRARIES)
   TARGET_LINK_LIBRARIES(noisegen_test ${Boost_LIBRARIES} )
ENDIF ()

########### timing_test ###############

SET(timing_test_SRCS
   timing_test.cu
   ${soars_SOURCE_DIR}/src/rsnoise.cu
   ${soars_SOURCE_DIR}/src/rsnoise.cuh
   ${soars_SOURCE_DIR}/src/rsdsp.cu
   ${soars_SOURCE_DIR}/src/rsdsp.cuh
   ${soars_SOURCE_DIR}/src/rspython.cu
   ${soars_SOURCE_DIR}/src/rspython.cuh
   ${soars_SOURCE_DIR}/src/rsparameters.cu
   ${soars_SOURCE_DIR}/src/rsparameters.cuh
   ${soars_SOURCE_DIR}/src/rsgeometry.cu
   ${soars_SOURCE_DIR}/src/rsgeometry.cuh
   ${soars_SOURCE_DIR}/src/rsdebug.cu
   ${soars_SOURCE_DIR}/src/rsdebug.cuh
   ${soars_SOURCE_DIR}/src/rstiming.cu
   ${soars_SOURCE_DIR}/src/rstiming.cuh
   )

ADD_EXECUTABLE(timing_test ${timing_test_SRCS})

# Link to fftwcpp libraries
TARGET_LINK_LIBRARIES(timing_test fftwcpp)

# See if the FFTW3 libraries are available
FIND_LIBRARY(soars_LIB_FFTW3 NAMES fftw3)

IF (soars_LIB_FFTW3)
   TARGET_LINK_LIBRARIES(timing_test ${soars_LIB_FFTW3})
ENDIF (soars_LIB_FFTW3)

# Check if we are building on a multithreaded FFTW
FIND_LIBRARY(soars_LIB_FFTW3_THREADS NAMES fftw3_threads)

IF (soars_LIB_FFTW3_THREADS)
   TARGET_LINK_LIBRARIES(timing_test ${soars_LIB_FFTW3_THREADS})
ENDIF (soars_LIB_FFTW3_THREADS)

# Include code to check for python libs

INCLUDE(FindPythonLibs)

IF (PYTHONLIBS_FOUND)
   TARGET_LINK_LIBRARIES(timing_test ${PYTHON_LIBRARIES})
   INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
ENDIF (PYTHONLIBS_FOUND)

# Link to Boost Thread and System
IF (Boost_LIBRARIES)
   TARGET_LINK_LIBRARIES(timing_test ${Boost_LIBRARIES} )
ENDIF ()

########### clock_test ###############

SET(clock_test_SRCS
   clock_test.cu
   ${soars_SOURCE_DIR}/src/rsnoise.cu
   ${soars_SOURCE_DIR}/src/rsnoise.cuh
   ${soars_SOURCE_DIR}/src/rsdsp.cu
   ${soars_SOURCE_DIR}/src/rsdsp.cuh
   ${soars_SOURCE_DIR}/src/rspython.cu
   ${soars_SOURCE_DIR}/src/rspython.cuh
   ${soars_SOURCE_DIR}/src/rsparameters.cu
   ${soars_SOURCE_DIR}/src/rsparameters.cuh
   ${soars_SOURCE_DIR}/src/rsgeometry.cu
   ${soars_SOURCE_DIR}/src/rsgeometry.cuh
   ${soars_SOURCE_DIR}/src/rsdebug.cu
   ${soars_SOURCE_DIR}/src/rsdebug.cuh
   ${soars_SOURCE_DIR}/src/rstiming.cu
   ${soars_SOURCE_DIR}/src/rstiming.cuh
   )

ADD_EXECUTABLE(clock_test ${clock_test_SRCS})

# Link to fftwcpp libraries
TARGET_LINK_LIBRARIES(clock_test fftwcpp)

# See if the FFTW3 libraries are available
IF (soars_LIB_FFTW3)
   TARGET_LINK_LIBRARIES(clock_test ${soars_LIB_FFTW3})
ENDIF (soars_LIB_FFTW3)

# Check if we are building on a multithreaded FFTW
IF (soars_LIB_FFTW3_THREADS)
   TARGET_LINK_LIBRARIES(clock_test ${soars_LIB_FFTW3_THREADS})
ENDIF (soars_LIB_FFTW3_THREADS)

IF (PYTHONLIBS_FOUND)
   TARGET_LINK_LIBRARIES(clock_test ${PYTHON_LIBRARIES})
ENDIF (PYTHONLIBS_FOUND)

# Link to Boost thread and system
IF (Boost_LIBRARIES)
   TARGET_LINK_LIBRARIES(clock_test ${Boost_LIBRARIES} )
ENDIF ()


########### filter_test ###############

SET(filter_test_SRCS
   filter_test.cu
   ${soars_SOURCE_DIR}/src/rsnoise.cu
   ${soars_SOURCE_DIR}/src/rsnoise.cuh
   ${soars_SOURCE_DIR}/src/rsdsp.cu
   ${soars_SOURCE_DIR}/src/rsdsp.cuh
   ${soars_SOURCE_DIR}/src/rspython.cu
   ${soars_SOURCE_DIR}/src/rspython.cuh
   ${soars_SOURCE_DIR}/src/rsparameters.cu
   ${soars_SOURCE_DIR}/src/rsparameters.cuh
   ${soars_SOURCE_DIR}/src/rsgeometry.cu
   ${soars_SOURCE_DIR}/src/rsgeometry.cuh
   ${soars_SOURCE_DIR}/src/rsdebug.cu
   ${soars_SOURCE_DIR}/src/rsdebug.cuh
   )

ADD_EXECUTABLE(filter_test ${filter_test_SRCS})

# Link to fftwcpp libraries
TARGET_LINK_LIBRARIES(filter_test fftwcpp)

# See if the FFTW3 libraries are available
FIND_LIBRARY(soars_LIB_FFTW3 NAMES fftw3)

IF (soars_LIB_FFTW3)
   TARGET_LINK_LIBRARIES(filter_test ${soars_LIB_FFTW3})
ENDIF (soars_LIB_FFTW3)

# Check if we are building on a multithreaded FFTW
FIND_LIBRARY(soars_LIB_FFTW3_THREADS NAMES fftw3_threads)

IF (soars_LIB_FFTW3_THREADS)
   TARGET_LINK_LIBRARIES(filter_test ${soars_LIB_FFTW3_THREADS})
ENDIF (soars_LIB_FFTW3_THREADS)

# Include code to check for python libs

INCLUDE(FindPythonLibs)

IF (PYTHONLIBS_FOUND)
   TARGET_LINK_LIBRARIES(filter_test ${PYTHON_LIBRARIES})
   INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
ENDIF (PYTHONLIBS_FOUND)

# Link to Boost Thread and System

IF (Boost_LIBRARIES)
   TARGET_LINK_LIBRARIES(filter_test ${Boost_LIBRARIES} )
ENDIF ()

########### install files ###############
