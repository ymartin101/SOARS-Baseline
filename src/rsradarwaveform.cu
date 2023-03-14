/// rsradarwaveform.cpp - Classes for different types of radar waveforms
/// Marc Brooker, 24 May 2006
/// Edited by Yaaseen Martin, 27 August 2019

#include <math.h>
#include <stdexcept>
#include <fstream>
#include <boost/scoped_array.hpp>
#include <ctype.h>
#include "rsradarwaveform.cuh"
#include "rsparameters.cuh"
#include "rsdebug.cuh"
#include "rssignal.cuh"

using namespace rs;

/// RadarWaveform implementation

// Default constructor
RadarSignal::RadarSignal(std::string name, rsFloat power, rsFloat carrierfreq, rsFloat ext_temp, rsFloat length, Signal* signal):
    name(name),
    power(power),
    carrierfreq(carrierfreq),
    ext_temp(ext_temp),
    length(length),
    signal(signal),
    polar(1.0,0) // Default to horizontal polarization
{
    if (!signal)
        printf("RadarSignal cannot be constructed with NULL signal");
}

// Destructor
RadarSignal::~RadarSignal()
{
    delete signal;
}

// Return the power of the signal
rsFloat RadarSignal::GetPower() const
{
    return power;
}

// Get the carrier frequency
rsFloat RadarSignal::GetCarrier() const
{
    return carrierfreq;
}

// Get the external noise temperature (K)
rsFloat RadarSignal::GetTemp() const
{
    return ext_temp;
}

// Get the name of the pulse
std::string RadarSignal::GetName() const
{
  return name;
}

// Get the native sample rate of the pulse
rsFloat RadarSignal::GetRate() const
{
    return signal->Rate();
}

// Return the length of the pulse
rsFloat RadarSignal::GetLength() const
{
    return length;
}

// Render the waveform to the target buffer
boost::shared_array<rsComplex> RadarSignal::Render(const std::vector<InterpPoint> &points, unsigned int &size, rsFloat frac_win_delay) const
{
    //Render the return pulse
    boost::shared_array<rsComplex> data = signal->Render(points, power, size, frac_win_delay);

    //Scale the return pulse by the signal power
    rsFloat scale = std::sqrt(power);
    for (unsigned int i = 0; i < size; i++)
        data[i] *= scale;
    return data;
}

// Get the signal polarization
JonesVector RadarSignal::GetPolarization()
{
    return polar;
}

// Set the signal polarization
void RadarSignal::SetPolarization(const JonesVector &in)
{
    polar = in;
}

/// Functions to load signal data from files

// Create the pulse from parameters
rs::RadarSignal* rsPulseFactory::CreatePulse(const std::string& name, rsFloat power, rsFloat carrierfreq, rsFloat bw, rsFloat ext_temp, rsFloat length, rsFloat rate)
{
    // Create pulse time vector using an equivalent to "linspace" function (adapted from code by Damith Jinasena, 2014)
    std::vector<double> t_pulse;
    int iterate = 0;
    double pi = 3.14159265358979323846; // Define pi
    int n = int(length * rate);
    for (int i = 0; i <= (n-2); i++)
    {
        double temp = 0 + i*(length - 0)/(floor((double)n) - 1);
        t_pulse.insert(t_pulse.begin() + iterate, temp);
        iterate += 1;
    }
    t_pulse.insert(t_pulse.begin() + iterate, length);

    // Create the pulse's I and Q channels
    int tsize = t_pulse.size();
    std::vector<double> I(tsize);
    std::vector<double> Q(tsize);
    std::complex<double> j(0, 1);
    rsComplex *data = new std::complex<rsFloat>[tsize];

    for (int i = 0; i < tsize; i++){
        // Linear FM/chirp pulse
        // I[i] = cos(M_PI * (bw/length) * pow(t_pulse[i],2));
        // Q[i] = sin(M_PI * (bw/length) * pow(t_pulse[i],2));
        // data[i] = std::complex<double>(I[i], Q[i]);

        // Rectangular pulse
        data[i] = std::complex<double>(1, 0);
    }

    // Create the signal object
    Signal *signal = new Signal();

    // Load the pulse into the signal object
    signal->Load(data, tsize, rate);
    delete[] data;

    // Create the RadarSignal
    rs::RadarSignal *any = new rs::RadarSignal(name, power, carrierfreq, ext_temp, tsize/rate, signal);
    return any;
}
