INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} )
INCLUDE_DIRECTORIES( ${CMAKE_BINARY_DIR} )

########### next target ###############

SET(sgp4_SRCS
   astMath.cu
   astMath.cuh
   astTime.cu
   astTime.cuh
   SGP4.cu
   SGP4.cuh
)

ADD_LIBRARY(sgp4 STATIC ${sgp4_SRCS})
