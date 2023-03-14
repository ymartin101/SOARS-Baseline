/// rsresponse.h - Classes for responses created by simulation
/// Marc Brooker, 3 August 2006
/// Edited by Yaaseen Martin, 27 August 2019

#ifndef __RSRESPONSE_H
#define __RSRESPONSE_H

#include <config.h>
#include <vector>
#include <string>
#include "rsradarwaveform.cuh"
#include <boost/utility.hpp>
#include <boost/shared_array.hpp>
#include <iosfwd>
#include "pugixml.hpp"

namespace rs {

    // Forward definition of Antenna (see rsantenna.h)
    class Antenna;
    // Forward definition of Transmitter (in this file)
    class Transmitter;

    // Class for both pulsed and CW responses
    class Response: boost::noncopyable {
        public:
            // Constructor
            Response(RadarSignal* wave, const Transmitter *transmitter);

            // Get the start time
            rsFloat StartTime() const;
            rsFloat EndTime() const;

            // Destructor
            ~Response();

            // Render the response to an XML file
            void RenderXML(pugi::xml_node &root);

            // Render the response to a CSV file
            void RenderCSV(std::ofstream &of);

            // Render the response to a binary file
            boost::shared_array<rsComplex> RenderBinary(rsFloat& rate, unsigned int &size, rsFloat frac_win_delay);

            // Get the length of the pulse
            rsFloat GetLength() const;

            // Get a pointer to the wave
            const rs::RadarSignal* GetWave() const;

            // Get the name of the transmitter that started this response
            std::string GetTransmitterName() const;

            // Add an interp point to the vector
            void AddInterpPoint(InterpPoint &point);

            // Return the number of Points
            int CountPoints() const;

        protected:
            const Transmitter *transmitter; // The transmitter that caused this response
            void RenderResponseXML(pugi::xml_node &root, const InterpPoint &point);
            void RenderResponseCSV(std::ofstream &of, const InterpPoint &point); // Render a InterpPoint as CSV
            const rs::RadarSignal *wave; // The waveform sent out by the transmitter
            std::vector<InterpPoint> points; // Waypoints from which the response parameters are interpolated
    };

}

#endif
