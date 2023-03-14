/// noisegen_test.cpp - Test the noise generation code from SOARS
/// Marc Brooker
/// Edited by Yaaseen Martin, 28 August 2019

#include <config.h>
#include <iostream>
#include <boost/thread.hpp>
#include "rsnoise.cuh"

using namespace rs;
using namespace std;

unsigned int processors = 4;

int main()
{
    // Initialize noise generation
    rsNoise::InitializeNoise();

    // Create the generator (for alpha = 2)
    MultirateGenerator *gen = new MultirateGenerator(1, 4);
    for (int i = 0; i < 1e7; i++)
        gen->GetSample();
    delete gen;
}
