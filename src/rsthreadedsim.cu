/// rsthreadedsim.cpp - Thread management for the simulator
/// Marc Brooker, 19 July 2006
/// Edited by Yaaseen Martin, 27 August 2019

// One of the goals for SOARS is to support multiple processors using multithreading
// One simulation is performed in for each transmitter-receiver pair
// Several simulations are run in parallel based on the system's number of CPUs (or cores)

#include <sys/time.h>
#include <vector>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/xtime.hpp>
#include <boost/version.hpp>

#if BOOST_VERSION < 105000
#define TIME_UTC_ TIME_UTC
#endif

#include "rssim.cuh"
#include "rsworld.cuh"
#include "rsradar.cuh"
#include "rsdebug.cuh"
#include "rsthreadedsim.cuh"

using namespace rs;

// Counter of currently running threads
volatile int threads;
boost::mutex threadsMutex; // Mutex to protect it

// Flag to set if a thread encounters an error
volatile int error;
boost::mutex errorMutex;

// Method to decrease the running thread count
static void decThreads()
{
    boost::mutex::scoped_lock lock(threadsMutex);
    threads--;
}

// Method to flag if a thread experienced an error
static void setError()
{
    boost::mutex::scoped_lock lock(errorMutex);
    error = 1;
}

// SimThread contains a simulator thread
class SimThread {
    public:
        // Constructor
        SimThread(Transmitter *transmitter, Receiver *receiver, World *world, int MaxThreads, int MaxBlocks):
        trans(transmitter), recv(receiver), world(world), MaxThreads(MaxThreads), MaxBlocks(MaxBlocks)
        {
        }

        // "operator ()" is executed when the thread is created
        void operator()() {
            if (trans->IsMonostatic() && (trans->GetAttached() == recv))
                rsDebug::printf(rsDebug::RS_VERBOSE, "Creating simulation thread for radar '%s'.\n", trans->GetName().c_str());
            else
                rsDebug::printf(rsDebug::RS_VERBOSE, "Creating simulation thread for transmitter '%s' and receiver '%s'.\n", trans->GetName().c_str(), recv->GetName().c_str());
            try {
                SimulatePair(trans, recv, world, MaxThreads, MaxBlocks);
            }
            catch (std::exception &ex) {
                rsDebug::printf(rsDebug::RS_CRITICAL, "ERROR: First pass thread terminated with unexpected error:\n\t%s\nSimulator will terminate\n", ex.what());
                setError();
            }
            decThreads();
        }

    protected:
        // The transmitter-receiver pair to simulate
        Transmitter *trans;
        Receiver *recv;
        World *world;
        int MaxThreads;
        int MaxBlocks;

};

// RenderThread contains a thread which performs the second phase of the simulation
class RenderThread {
    public:
        // Constructor
        RenderThread(Receiver *recv):
        recv(recv)
        {
        }

        // "operator ()" is executed when we create the thread
        void operator()() {
            if (recv->IsMonostatic())
                rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "Creating render thread for radar '%s'...\n", recv->GetName().c_str());
            else
                rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "Creating render thread for receiver '%s'...\n", recv->GetName().c_str());
            try {
                recv->Render();
            }
            catch (std::exception &ex) {
                rsDebug::printf(rsDebug::RS_CRITICAL, "ERROR: Render thread terminated with unexpected error:\n\t%s\nSimulator will terminate\n", ex.what());
                setError();
            }
            decThreads();
        }
    protected:
        Receiver *recv; // The receiver to render

};

// Increase the count of running threads
static void IncThreads()
{
    boost::mutex::scoped_lock lock(threadsMutex);
    threads++;
}

// Run the SOARS simulation
void rs::RunThreadedSim(World *world, int MaxThreads, int MaxBlocks, int thread_limit) {

    // Timer for overall radar time
    struct timeval timer;
    gettimeofday(&timer, NULL);
    double radar_StartTime = timer.tv_sec + (timer.tv_usec/1000000.0);

    // GPU to be used for target processing for each Tx-Rx pair
    std::vector<boost::thread *> running;
    std::vector<Receiver*>::iterator ri;
    std::vector<Transmitter*>::const_iterator ti;
    boost::thread mainthrd();
    rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "Starting CUDA simulation...\n");

    // Phase 1: Do first pass of simulator - go through the lists for transmitters and receivers
    for (ri = world->receivers.begin(); ri != world->receivers.end(); ri++)
    {
        for (ti = world->transmitters.begin(); ti != world->transmitters.end(); ti++)
        {
            IncThreads();
            SimThread sim(*ti, *ri, world, MaxThreads, MaxBlocks);
            boost::thread *thrd = new boost::thread(sim);

            // Delay until a thread is terminated, if we have reached the limit
            while (threads >= thread_limit)
                boost::thread::yield();

            // If a thread ended in error, abort the simulation
            if (error)
                throw std::runtime_error("Thread terminated with error. Aborting simulation.");

            // Add the thread pointers to a vector to be freed later
            running.push_back(thrd);
        }
    }

    // Wait for all the first pass threads to finish before continuing
    while (threads)
    {
        boost::thread::yield();
        if (error)
            throw std::runtime_error("Thread terminated with error. Aborting simulation.");
    }

    // Clean all the thread pointers
    for (std::vector<boost::thread*>::iterator i = running.begin(); i != running.end(); i++)
        delete *i;
    running.clear();

    // Timer for overall radar time
    gettimeofday(&timer, NULL);
    double radar_time = timer.tv_sec + (timer.tv_usec/1000000.0) - radar_StartTime;
    rsDebug::printf(rsDebug::RS_CRITICAL, "Radar computation took %lf seconds.\n", radar_time);

    // // Phase 2: Do render pass of simulation - go through the lists of receivers and set each to render XML/CSV file(s)
    // for (ri = world->receivers.begin(); ri != world->receivers.end(); ri++)
    // {
    //     // Print info about Rx rendering
    //     if ((*ri)->IsMonostatic())
    //         rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "Rendering %d responses for radar '%s'.\n", (*ri)->CountResponses(), (*ri)->GetName().c_str());
    //     else
    //         rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "Rendering %d responses for receiver '%s'.\n", (*ri)->CountResponses(), (*ri)->GetName().c_str());

    //     IncThreads();
    //     RenderThread sim(*ri);  // Render XML/CSV outputs
    //     boost::thread *thrd = new boost::thread(sim);

    //     // Delay until a thread is terminated, if we have reached the limit
    //     while (threads >= thread_limit)
    //         boost::thread::yield();

    //     // If a thread ended in error, abort the simulation
    //     if (error)
    //         throw std::runtime_error("Thread terminated with error. Aborting simulation.");

    //     // Add the thread pointers to a vector to be freed later
    //     running.push_back(thrd);
    // }

    // // Wait for all the render threads to finish
    // while (threads)
    // {
    //     boost::thread::yield();
    //     if (error)
    //         throw std::runtime_error("Thread terminated with error. Aborting simulation.");
    // }

    // // Clean all the thread pointers
    // for (std::vector<boost::thread*>::iterator i = running.begin(); i != running.end(); i++)
    //     delete *i;
    // running.clear();

    // Phase 3: Render receiver HDF5 responses serially (to avoid threading issues with C++ HDF5 library)
    for (ri = world->receivers.begin(); ri != world->receivers.end(); ri++)
    {
        (*ri)->RenderBin();
    }
}
