INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} )
INCLUDE_DIRECTORIES( ${CMAKE_BINARY_DIR} )

########### next target ###############

SET(pugixml_SRCS
   pugixml.cpp
   pugixml.hpp
   pugiconfig.hpp
)

ADD_LIBRARY(pugixml STATIC ${pugixml_SRCS})
