/// xmlimport.h - Import a simulator world and simulation parameters from an XML file
/// Marc Brooker, 26 April 2006
/// Edited by Yaaseen Martin, 27 August 2019

#ifndef __XMLIMPORT_H
#define __XMLIMPORT_H

#include <config.h>
#include "rsworld.cuh"

typedef enum
{
    e80,
    e96,
    e00a,
    e00b
}   eOpt;

namespace xml {

    // Load an XML file containing a simulation description
    void LoadXMLFile(std::string filename, rs::World *world, int MaxThreads, int MaxBlocks); // Puts the structure of the world into the World object in the process

    class XMLException
    {
    };

}

#endif
