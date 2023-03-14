/// signal_test.cpp - Test program for the signal class
/// Marc Brooker, 25 May 2006
/// Edited by Yaaseen Martin, 28 August 2019

#include <iostream>
#include <math.h>
#include "rssignal.cuh"
#include "fftwcpp.h"

using namespace rsSignal;
using namespace std;

void dump(Signal &sig)
{
    cDbl *data = sig.DataPtr();
    for (int i = 0; i < sig.Size(); i++)
        cout << data[i] << "\n";
    cout << endl;
}

void fillsignal(Signal &sig, double freq, double sr, double time)
{
    int size = sr*time;
    double *data = new double[size];

    // Fill the data with a cos wave
    for (int i = 0; i < size; i++)
        data[i] = cos(i/sr*2*M_PI*freq);

    // Load the signal
    sig.Load(data, size, sr);
    delete[] data;
}


int main()
{
    Signal sig1, sig2;
    fillsignal(sig1, 1e3, 3e3, 1e-1);
    fillsignal(sig2, 1e3, 3e3, 1e-1);

    sig1.Decimate(3);
    sig1.Interpolate(3);
    sig2 -= sig1;

    FFTManager *manager = FFTManager::Instance();
    manager->Clean();
}
