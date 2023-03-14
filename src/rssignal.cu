/// rssignal.cpp - Class to contain an arbitrary frequency domain signal
/// Marc Brooker, 24 May 2006
/// Edited by Craig Tong
/// Edited by Yaaseen Martin, 27 August 2019

#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <stdio.h>
#include <boost/thread/mutex.hpp>
#include "rssignal.cuh"
#include "rsdebug.cuh"
#include "rsnoise.cuh"
#include "rsparameters.cuh"
#include "rsdsp.cuh"
#include "rsportable.cuh"

using namespace rsSignal;
using namespace rs;

// Global mutex to protect the InterpFilter filter calculation function
boost::mutex interp_mutex;

namespace
{

    // Compute the zeroth order modified bessel function of the first kind
    // Use the polynomial approximation from section 9.8 of "Handbook of Mathematical Functions" by Abramowitz and Stegun
    rsFloat BesselI0(rsFloat x)
    {
        rsFloat I0;
        rsFloat t = (x / 3.75);
        if (t < 0.0)
        {
            throw std::logic_error("Modified Bessel approximation only valid for x > 0");
        }
        else if (t <= 1.0)
        {
            t *= t;
            I0 = 1.0 + t * (3.5156229 + t * (3.0899424 + t * (1.2067492 + t * (0.2659732 + t * (0.0360768 + t * 0.0045813)))));
            // Error bounded to 1.6e-7
        }
        else
        {
            // t > 1;
            I0 = 0.39894228 + t * (0.01328592 + t * (0.00225319 + t * (-0.00157565 + t * (0.00916281 + t * (-0.02057706 + t * (0.02635537 + t * (-0.01647633 + t * 0.00392377)))))));
            I0 *= std::exp(x) / std::sqrt(x);
            // Error bounded to 1.9e-7
        }
        return I0;
    }

    class InterpFilter
    {
        public:
            // Compute the sinc function at the specified x
            inline rsFloat
            Sinc(rsFloat x) const;

            // Compute the value of a Kaiser Window at the given x
            inline rsFloat
            kaiser_win_compute(rsFloat x) const;

            // Calculate the value of the interpolation filter at time x
            inline rsFloat
            interp_filter(rsFloat x) const;

            // Get a pointer to the filter with approximately the specified delay
            const rsFloat*
            GetFilter(rsFloat delay) const;

            // Get a pointer to the class instance
            static InterpFilter*
            GetInstance()
            {

            // Protect this with a mutex --- all other operations are const
            boost::mutex::scoped_lock lock(interp_mutex);
            if (!instance)
                instance = new InterpFilter();
            return instance;
            }

        private:
            // Default constructor
            InterpFilter();

            // Copy Constructor
            InterpFilter(const InterpFilter& ifilt);

            // Assignment operator
            InterpFilter&
            operator=(const InterpFilter& ifilt);

            // Pointer to a single instance of the class
            static InterpFilter *instance;

            rsFloat alpha; // 'alpha' parameter
            rsFloat beta; // 'beta' parameter
            rsFloat bessel_beta; // I0(beta)
            int length;
            int table_filters; // Number of filters in the filter table
            rsFloat *filter_table; // Table of precalculated filters
    };

    // Interpfilter class constructor
    InterpFilter::InterpFilter()
    {
        length = rsParameters::render_filter_length();

        // Size of the table to use for interpolation
        table_filters = 1000;

        // Allocate memory for the table
        filter_table = new rsFloat[table_filters * length];

        // Alpha is half the filter length
        alpha = std::floor(rsParameters::render_filter_length() / 2.0);

        // Beta sets the window shape
        beta = 5;
        bessel_beta = BesselI0(beta);
        int hfilt = table_filters / 2;

        // Fill the table of filters; delay is the fraction of time ellapsed between samples
        for (int i = -hfilt; i < hfilt; i++)
        {
            rsFloat delay = i / rsFloat(hfilt);
            for (int j = -alpha; j < alpha; j++)
                filter_table[int((i + hfilt) * length + j + alpha)] = interp_filter(j - delay);
        }
        //rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "Built filter table of %d filters\n", table_filters);
    }

    // Get a pointer to the filter with approximately the specified delay
    const rsFloat*
    InterpFilter::GetFilter(rsFloat delay) const
    {
        int filt = (delay + 1) * (table_filters / 2);

        if ((delay <= -1) || (delay >= 1))
        {
            rsDebug::printf(rsDebug::RS_VERY_VERBOSE, "GetFilter %f %d\n", delay, filt);
            throw std::runtime_error("[BUG] Requested delay filter value out of range");
        }
        return &(filter_table[filt * length]);
    }

    // Look-up the value of the interpolation filter at time x
    rsFloat
    InterpFilter::interp_filter(rsFloat x) const
    {
        rsFloat w = kaiser_win_compute(x + alpha);
        rsFloat s = Sinc(x); // The filter value
        rsFloat filt = w * s;
        return filt;
    }

    // Compute the sinc function at the specified x
    rsFloat
    InterpFilter::Sinc(rsFloat x) const
    {
    if (x == 0)
        return 1.0;
    return std::sin(x * M_PI) / (x * M_PI);
    }

    // Compute the value of a Kaiser Window at the given x
    rsFloat
    InterpFilter::kaiser_win_compute(rsFloat x) const
    {
    if ((x < 0) || (x > (alpha * 2)))
        return 0;
    else
        return BesselI0(beta * std::sqrt(1 - std::pow((x - alpha) / alpha, 2))) / bessel_beta;
    }

    // Create an instance of InterpFilter
    InterpFilter* InterpFilter::instance = 0;

}

/// Simulate the effect of and ADC converter on the signal
void rsSignal::ADCSimulate(complex *data, unsigned int size, int bits, rsFloat fullscale)
{
    // Get the number of levels associated with the number of bits
    rsFloat levels = pow(2, bits - 1);
    for (unsigned int i = 0; i < size; i++)
    {
        // Simulate the ADC effect on the I and Q samples
        rsFloat I = std::floor(levels * data[i].real() / fullscale) / levels;
        rsFloat Q = std::floor(levels * data[i].imag() / fullscale) / levels;

        // Clamp I and Q to the range, simulating saturation in the adc
        if (I > 1)
            I = 1;
        else if (I < -1)
            I = -1;
        if (Q > 1)
            Q = 1;
        else if (Q < -1)
            Q = -1;

        // Overwrite data with the results
        data[i] = complex(I, Q);
    }
}

/// Signal Implementation

// Default constructor for brute signal
Signal::Signal() :
    data(0), size(0), rate(0)
{
}

// Default destructor for brutesignal
Signal::~Signal()
{
    delete[] data;
}

// Clear the data array, emptying the signal and freeing memory
void
Signal::Clear()
{
    delete[] data;
    size = 0;
    rate = 0;
    data = 0; // Set data to zero to prevent multiple frees
}

// Load data into the signal, with the given sample rate and size
void Signal::Load(const rsFloat *indata, unsigned int samples, rsFloat samplerate)
{
    // Remove the previous data
    Clear();

    // Set the size and samples attributes
    size = samples;
    rate = samplerate;

    // Allocate memory for the signal
    data = new complex[samples];
    // Copy the data

    for (unsigned int i = 0; i < samples; i++)
        data[i] = complex(indata[i], 0.0);
}

// Load data into the signal (time domain, complex)
void Signal::Load(const complex *indata, unsigned int samples, rsFloat samplerate)
{
    // Remove the previous data
    Clear();

    // Get the oversampling ratio
    unsigned int ratio = rsParameters::oversample_ratio();

    // Allocate memory for the signal
    data = new complex[samples * ratio];

    // Set the size and samples attributes
    size = samples * ratio;
    rate = samplerate * ratio;
    if (ratio == 1)
    {
        // Copy the data (using a loop for now, not sure memcpy() will always work on complex)
        for (unsigned int i = 0; i < samples; i++)
            data[i] = indata[i];
    }
    else
    {
        // Upsample the data into the target buffer
        Upsample(indata, samples, data, ratio);
    }
}

// Return the sample rate of the signal
rsFloat
Signal::Rate() const
{
    return rate;
}

// Return the size, in samples, of the signal
unsigned int
Signal::Size() const
{
    return size;
}

// Get a copy of the signal domain data
rsFloat*
Signal::CopyData() const
{
    rsFloat* result = new rsFloat[size];

    //Copy the data into result
    std::memcpy(result, data, sizeof(rsFloat) * size);
    return result;
}

// Get the number of samples of padding at the beginning of the pulse
unsigned int
Signal::Pad() const
{
    return pad;
}

// Render the pulse with the specified doppler, delay and amplitude
boost::shared_array<rsComplex>
Signal::Render(const std::vector<InterpPoint> &points, rsFloat trans_power, unsigned int &out_size, rsFloat frac_win_delay) const
{
    // Allocate memory for rendering
    rsComplex *out = new rsComplex[size];
    out_size = size;

    // Get the sample interval
    rsFloat timestep = 1.0 / rate;

    // Create the rendering window
    const int filt_length = rsParameters::render_filter_length();
    InterpFilter* interp = InterpFilter::GetInstance();

    // Loop through the interp points, rendering each in time
    std::vector<rs::InterpPoint>::const_iterator iter = points.begin();
    std::vector<rs::InterpPoint>::const_iterator next = iter + 1;
    if (next == points.end())
        next = iter;

    // Get the delay of the first point; idelay is measured in number of receiver samples (possibly with a fractional part)
    rsFloat idelay = rsPortable::rsRound(rate * (*iter).delay);

    // Memory to store the filter in
    const rsFloat *filt;

    // Loop over the pulse, performing the rendering
    rsFloat sample_time = (*iter).time;
    for (int i = 0; i < (int) size; i++, sample_time += timestep)
    {
        // Check if we should move on to the next set of interp points
        if ((sample_time > (*next).time))
        {
            iter = next;
            if ((next + 1) != points.end())
                next++;
        }

        // Get the weightings for the parameters
        rsFloat aw = 1, bw = 0;
        if (iter < next)
        {
          bw = (sample_time - (*iter).time) / ((*next).time - (*iter).time); // Fraction of time elapsed between 2 samples (interpolated)
          aw = 1 - bw; // Fraction of time remaining
        }

        // Now calculate the current sample parameters
        rsFloat amplitude = std::sqrt((*iter).power) * aw + std::sqrt((*next).power) * bw;
        rsFloat phase = (*iter).phase * aw + (*next).phase * bw;
        rsFloat fdelay = -(((*iter).delay * aw + (*next).delay * bw) * rate - idelay + frac_win_delay);
        int iSampleUnwrap = floor(fdelay); // Number of samples to unwrap by.
        fdelay -= iSampleUnwrap; // Re-calculate the delay filter for the given delay
        filt = interp->GetFilter(fdelay);

        // Get the start and end times of interpolation
        int start = -filt_length / 2;
        if ((i + start) < 0)
            start = -i;

        int end = filt_length / 2;
        if ((i + end) >= size)
            end = size - i;

        // Apply the filter
        complex accum = 0;
        for (int j = start; j < end; j++)
        {
              // Check that unwrapping doesn't put us out of bounds.
              if (i + j + iSampleUnwrap >= size || i + j + iSampleUnwrap < 0)
                    continue;
              accum += amplitude * data[i + j + iSampleUnwrap] * filt[j + filt_length/2]; // Apply unwrapping to Tx samples
        }

        // Perform IQ demodulation
        rs::rsComplex ph = exp(rs::rsComplex(0.0, 1.0) * phase);
        out[i] = ph * accum;
    }

    // Return the result
    return boost::shared_array<rs::rsComplex>(out);
}
