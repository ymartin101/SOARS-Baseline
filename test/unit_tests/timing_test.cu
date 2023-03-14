/// timing_test.cpp - Test the clock model code from SOARS
/// Marc Brooker
/// Edited by Yaaseen Martin, 28 August 2019

#include <config.h>
#include <iostream>
#include <boost/thread.hpp>
#include <sstream>
#include "rsnoise.cuh"
#include "rstiming.cuh"
#include "rsparameters.cuh"

using namespace rs;
using namespace std;

unsigned int processors = 4;

int main(int argc, char *argv[])
{
    // Get command line args
    if (argc != 4) {
        cout << "Usage: timing_test alpha weight count" << endl;
    }
    else {
        istringstream iss1(argv[1]);
        double alpha;
        iss1 >> alpha;
        double weight;
        istringstream iss2(argv[2]);
        iss2 >> weight;
        double count;
        istringstream iss3(argv[3]);
        iss3 >> count;
        cerr << "Alpha: " << alpha << " weight: " << weight << " count: " << count << endl;

        // Initialize noise generation
        rsNoise::InitializeNoise();
        rsParameters::modify_parms()->SetRate(1e3);
        PrototypeTiming proto("test");
        proto.AddAlpha(alpha, weight);
        proto.SetFrequency(1e6);
        ClockModelTiming timing("test");
        timing.InitializeModel(&proto);
        for (int i = 0; i < count; i++)
            cout << timing.NextNoiseSample() << endl;
    }
}
