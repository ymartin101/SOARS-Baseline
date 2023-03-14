/// rssim.cuh - Declarations and definitions for rssim.cpp and rsthreadedsim.cpp
/// Marc Brooker, 30 May 2006
/// Edited by Yaaseen Martin, 27 August 2019

#ifndef __RSSIM_H
#define __RSSIM_H

#include <config.h>
#include "rsworld.cuh"
#include "rsradar.cuh"

namespace rs {
    // Run the radar simulation specified by world
    // Limit the number of concurrent threads to thread_limit
    void RunThread(int thread_limit, World *world);

    // Run a simulation on the specified receiver/transmitter pair
    void SimulatePair(rs::Transmitter *trans, rs::Receiver *recv, rs::World *world, int MaxBlocks, int MaxThreads);
}

#endif
