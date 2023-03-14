INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} )
INCLUDE_DIRECTORIES( ${CMAKE_BINARY_DIR} )

########### next target ###############

SET(fftwcpp_STAT_SRCS
   fftwcpp.cpp
   fftwcpp.h
   )

ADD_LIBRARY(fftwcpp STATIC ${fftwcpp_STAT_SRCS})

#Check if we are building on a multithreaded FFTW
FIND_LIBRARY(SOARS_LIB_FFTW3_THREADS NAMES fftw3_threads)

IF (SOARS_LIB_FFTW3_THREADS)
   ADD_DEFINITIONS(-DFFTW_WITH_THREADS)
ENDIF (SOARS_LIB_FFTW3_THREADS)


########### install files ###############
