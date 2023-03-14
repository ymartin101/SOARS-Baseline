/// rsthreadedsim.cuh - Definitions for threaded simulation code
/// Marc Brooker, 19 July 2006
/// Edited by Yaaseen Martin, 27 August 2019

#ifndef __RSTHREADEDSIM_H
#define __RSTHREADEDSIM_H

#include "rsworld.cuh"

namespace rs {
    void RunThreadedSim(rs::World *world, int MaxThreads, int MaxBlocks, int thread_limit);
}

#endif
