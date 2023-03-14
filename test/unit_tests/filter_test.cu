/// filter_test.cpp - Test the noise generation code from SOARS
/// Marc Brooker
/// Edited by Yaaseen Martin, 28 August 2019

#include <config.h>
#include <iostream>
#include <boost/thread.hpp>
#include "rsnoise.cuh"
#include "rsdsp.cuh"

using namespace rs;
using namespace std;

unsigned int processors = 4;

int main()
{
    rsNoise::InitializeNoise();

    DecadeUpsampler upsample;
    for (int i = 0; i < 1e4; i++) {
        rsFloat buf[10];
        upsample.Upsample(rsNoise::WGNSample(1), buf);
        for (int j = 0; j < 10; j++)
            cout << buf[j] << endl;
    }
}
