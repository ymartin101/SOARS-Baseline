INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} )
INCLUDE_DIRECTORIES( ${CMAKE_BINARY_DIR} )
INCLUDE_DIRECTORIES( ${soars_SOURCE_DIR}/fftwcpp )
INCLUDE_DIRECTORIES( ${soars_SOURCE_DIR}/sgp4 )
INCLUDE_DIRECTORIES( ${soars_SOURCE_DIR}/pugixml )

########### next target ###############

SET(soars_SRCS
   rsdebug.cu
   rsmain.cu
   rsplatform.cu
   rsradar.cu
   rstarget.cu
   rsworld.cu
   xmlimport.cu
   rsgeometry.cu
   rspath.cu
   rssim.cu
   rsthreadedsim.cu
   rspulserender.cu
   rsradarwaveform.cu
   rsparameters.cu
   rsobject.cu
   rsantenna.cu
   rsportable.cu
   rsresponse.cu
   rsnoise.cu
   rssignal.cu
   rstiming.cu
   rshdf5.cu
   rspython.cu
   rsinterp.cu
   rsdsp.cu
   rsmultipath.cu
   rspattern.cu
   rspolarize.cu
   rsdebug.cuh
   rspath.cuh
   rsportable.cuh
   rsradar.cuh
   rssim.cuh
   rsworld.cuh
   rsgeometry.cuh
   rsplatform.cuh
   rstarget.cuh
   xmlimport.cuh
   rsradarwaveform.cuh
   rsparameters.cuh
   rssignal.cuh
   rsthreadedsim.cuh
   rsobject.cuh
   rsantenna.cuh
   rsresponse.cuh
   rsnoise.cuh
   rstiming.cuh
   rshdf5.cuh
   rspython.cuh
   rsinterp.cuh
   rsdsp.cuh
   rsmultipath.cuh
   rspattern.cuh
   rspolarize.cuh
   )

 function(my_add_executable TargetName)
   set(Files ${ARGV})
   list(REMOVE_AT Files 0)
   add_executable(${TargetName} ${Files})
   set_target_properties(${TargetName} PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON
                             RUNTIME_OUTPUT_DIRECTORY
                                 "${CMAKE_SOURCE_DIR}/build")
 endfunction()

my_add_executable(soars ${soars_SRCS})

#Link to the HDF5 libraries
FIND_LIBRARY(SOARS_LIB_HDF5 NAMES hdf5 )

IF (SOARS_LIB_HDF5)
   TARGET_LINK_LIBRARIES(soars ${SOARS_LIB_HDF5})
   ADD_DEFINITIONS(-DHAVE_LIBHDF5)
ENDIF (SOARS_LIB_HDF5)

FIND_LIBRARY(SOARS_LIB_HDF5_HL NAMES hdf5_hl )

IF (SOARS_LIB_HDF5_HL)
   TARGET_LINK_LIBRARIES(soars ${SOARS_LIB_HDF5_HL})
ENDIF (SOARS_LIB_HDF5_HL)

#Link to tinyxml libraries
FIND_LIBRARY(SOARS_LIB_TINYXML NAMES tinyxml)

IF (SOARS_LIB_TINYXML)
   TARGET_LINK_LIBRARIES(soars ${SOARS_LIB_TINYXML})
ENDIF (SOARS_LIB_TINYXML)

#Link to fftwcpp libraries
TARGET_LINK_LIBRARIES(soars fftwcpp)

#Link to sgp4 libraries
TARGET_LINK_LIBRARIES(soars sgp4)

#Link to pugixml libraries
TARGET_LINK_LIBRARIES(soars pugixml)

#See if the FFTW3 libraries are available
FIND_LIBRARY(SOARS_LIB_FFTW3 NAMES fftw3)

IF (SOARS_LIB_FFTW3)
   TARGET_LINK_LIBRARIES(soars ${SOARS_LIB_FFTW3})
ENDIF (SOARS_LIB_FFTW3)

#Check if we are building on a multithreaded FFTW
FIND_LIBRARY(SOARS_LIB_FFTW3_THREADS NAMES fftw3_threads)

IF (SOARS_LIB_FFTW3_THREADS)
   TARGET_LINK_LIBRARIES(soars ${SOARS_LIB_FFTW3_THREADS})
ENDIF (SOARS_LIB_FFTW3_THREADS)

#Include code to check for python libs
#We must use Python2
SET(Python_ADDITIONAL_VERSIONS 2.7 2.6 2.5)
INCLUDE(FindPythonLibs)

IF (PYTHONLIBS_FOUND)
   TARGET_LINK_LIBRARIES(soars ${PYTHON_LIBRARIES})
   INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
ENDIF (PYTHONLIBS_FOUND)

#Link to Boost Thread and System
IF(WIN32)
   SET(Boost_USE_STATIC_LIBS ON) ##Link Boost statically
ELSE ()
   SET(Boost_USE_STATIC_LIBS OFF) ##Link Boost with shared libs
ENDIF ()
FIND_PACKAGE(Boost COMPONENTS system thread REQUIRED)

MESSAGE(STATUS "Linking libraries for Boost: " ${Boost_LIBRARIES})

IF (Boost_LIBRARIES)
   TARGET_LINK_LIBRARIES(soars ${Boost_LIBRARIES} )
ENDIF ()

########### install files ###############

INSTALL( TARGETS soars DESTINATION bin)
