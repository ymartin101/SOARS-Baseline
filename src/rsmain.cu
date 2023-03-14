/// rsmain.cu - Contains the main function and some support code for SOARS (now with CUDA)
/// Marc Brooker, 25 April 2006
/// Edited by Yaaseen Martin, 19 March 2020

#include <sys/time.h>
#include <stdexcept>
#include "rsworld.cuh"
#include "xmlimport.cuh"
#include "rsdebug.cuh"
#include "rsthreadedsim.cuh"
#include "rsnoise.cuh"
#include "fftwcpp.h"
#include "rsparameters.cuh"
#include "rsportable.cuh"
#include <cstring>

/// SOARS main function - runs on host only
int main(int argc, char *argv[])
{
    // Timer for full runtime
    struct timeval timer;
    gettimeofday(&timer, NULL);
    double StartTime = timer.tv_sec + (timer.tv_usec/1000000.0);

    // SOARS start-up sequence
    rsDebug::printf(rsDebug::RS_CRITICAL, "\n/--------------------------------------------------------\\\n");
    rsDebug::printf(rsDebug::RS_CRITICAL, "| SOARS - Space Object Astrodynamics and Radar Simulator |\n");
    rsDebug::printf(rsDebug::RS_CRITICAL, "| Version 0.20 (feat. CUDA)                              |\n");
    rsDebug::printf(rsDebug::RS_CRITICAL, "\\--------------------------------------------------------/\n\n");

    if (argc != 2 || !strncmp(argv[1], "--help", 6))
    {
        rsDebug::printf(rsDebug::RS_CRITICAL, "Use '%s <SOARSXML FILE>' to run the simulation defined by the file!\n\n", argv[0]);
        return 2;
    }

    try
    {
        // Set the number of CPU boost threads
        rs::rsParameters::modify_parms()->SetThreads(rsPortable::CountProcessors());

        // Set the GPU thread and block limits (assumes CUDA Compute Capability of 3.0 or higher)
        int MaxThreads = 256;      // 1024 for Compute Capability 2.x and higher
        int MaxBlocks = 40000000; // Round number; Compute Capability 3.0 and higher

        // Create the world container
        rs::World *world = new rs::World();

        // Initialise the RNG code
        rsNoise::InitializeNoise();

        // Start the simulation
        rsDebug::printf(rsDebug::RS_CRITICAL, "---------------------------------------------------------\n");
        rsDebug::printf(rsDebug::RS_VERBOSE, "SIMULATION: '%s'\n", argv[1]);
        rsDebug::printf(rsDebug::RS_CRITICAL, "---------------------------------------------------------\n\n");
        rsDebug::printf(rsDebug::RS_VERBOSE, "Loading SOARSXML file.\n");

        // Load the script file
        xml::LoadXMLFile(argv[1], world, MaxThreads, MaxBlocks);

        // Start the threaded simulation (CUDA used for threaded simulation)
        rs::RunThreadedSim(world, MaxThreads, MaxBlocks, rsPortable::CountProcessors());
        rsDebug::printf(rsDebug::RS_VERBOSE, "Finishing simulation...\n");

        // Clean up the world model
        delete world;

        // Clean up singleton objects
        rsNoise::CleanUpNoise();

        // Complete the simulation
        rsDebug::printf(rsDebug::RS_CRITICAL, "\n---------------------------------------------------------\n");
        rsDebug::printf(rsDebug::RS_CRITICAL, "Success! The SOARS simulation was completed.\n");
        rsDebug::printf(rsDebug::RS_CRITICAL, "---------------------------------------------------------\n\n");

        // Timer for full runtime
        gettimeofday(&timer, NULL);
        double total_time = timer.tv_sec + (timer.tv_usec/1000000.0) - StartTime;
        rsDebug::printf(rsDebug::RS_CRITICAL, "Full simulation took %lf seconds.\n", total_time);
        // printf("Full simulation took %lf seconds.\n", total_time);

        return 0;
    }
    catch (std::exception &ex)
    {
        rsDebug::printf(rsDebug::RS_CRITICAL, "\n---------------------------------------------------------\n");
        rsDebug::printf(rsDebug::RS_CRITICAL, "Fail! SOARS has encountered the following error:\n%s\n", ex.what());
        rsDebug::printf(rsDebug::RS_CRITICAL, "---------------------------------------------------------\n\n");
        return 1;
    }
}
