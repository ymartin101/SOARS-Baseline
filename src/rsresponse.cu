/// rsresponse.cpp - Implementation of ResponseBase and derived classes
/// Marc Brooker, 3 August 2006
/// Edited by Yaaseen Martin, 27 August 2019

#include <sstream>
#include <fstream>
#include <iomanip>
#include "rsresponse.cuh"
#include "rsradar.cuh"

/// Attach a text node to an XML element for creating structures like this:
// <node>
//   <name>text</name>
// </node>

using namespace rs;

namespace{

    // Attach a text node to an XML element
    void AttachTextNode(pugi::xml_node &root, std::string name, std::string text)
    {
        pugi::xml_node element = root.append_child(name.c_str());
        element.append_child(pugi::node_pcdata).set_value(text.c_str());
    }

    // Attach a text node to an XML element, getting the text by converting a rsFloat to a string
    void AttachRsFloatNode(pugi::xml_node &root, std::string name, rsFloat data, bool scientific = true, int precision = 10)
    {
        std::ostringstream oss;
        if (scientific)
            oss.setf(std::ios::scientific);
        oss << std::setprecision(precision) << data;
        AttachTextNode(root, name, oss.str());
    }

    // Attach a text node to an XML element, getting the text by converting an int to a string
    void AttachIntNode(pugi::xml_node &root, std::string name, int data)
    {
        std::ostringstream oss;
        oss << data;
        AttachTextNode(root, name, oss.str());
    }

}

/// Interppoint Implementation

// InterpPoint Constructor
InterpPoint::InterpPoint(rsFloat power, rsFloat start, rsFloat delay, rsFloat doppler, rsFloat phase, rsFloat noise_temperature):
    power(power),
    time(start),
    delay(delay),
    doppler(doppler),
    phase(phase),
    noise_temperature(noise_temperature)
{
}

/// ResponseBase Implementation

Response::Response(RadarSignal* wave, const Transmitter* transmitter):
    transmitter(transmitter),
    wave(wave)
{
}

Response::~Response()
{
}

// Return the time the pulse's energy starts
rsFloat Response::StartTime() const
{
    if (points.empty())
        return 0;
    return points.front().time;
}

// Return the time the pulse's energy ends
rsFloat Response::EndTime() const
{
    if (points.empty())
        return 0;
    return points.back().time;
}

// Return the length of the pulse
rsFloat Response::GetLength() const
{
    return EndTime()-StartTime();
}

// Get the name of the transmitter which caused this response
std::string Response::GetTransmitterName() const
{
    return transmitter->GetName();
}

// Return a pointer to the waveform
const rs::RadarSignal* Response::GetWave() const
{
    return wave;
}

// Render a single response point to XML
void Response::RenderResponseXML(pugi::xml_node &resp, const InterpPoint &point)
{
    // Create a node for the response
    pugi::xml_node element = resp.append_child("InterpolationPoint");

    // Attach nodes for properties of the response
    AttachRsFloatNode(element, "time", point.time, false);
    AttachRsFloatNode(element, "delay", point.delay, false);
    // AttachRsFloatNode(element, "amplitude", std::sqrt(point.power*wave->GetPower()), false);
    AttachRsFloatNode(element, "phase", point.phase, false);
    AttachRsFloatNode(element, "doppler", point.doppler, false); // Fd = (Fr - Ft), i.e. +V == +Fd target moves towards radar
    AttachRsFloatNode(element, "power", point.power*wave->GetPower());
    // AttachRsFloatNode(element, "Iamplitude", std::cos(point.phase)*std::sqrt(point.power*wave->GetPower()));
    // AttachRsFloatNode(element, "Qamplitude", std::sin(point.phase)*std::sqrt(point.power*wave->GetPower()));
    // AttachRsFloatNode(element, "noise_temperature", point.noise_temperature);
    // AttachRsFloatNode(element, "phasedeg", point.phase/M_PI*180);
}

// Render the response to an XML file
void Response::RenderXML(pugi::xml_node &rec)
{
    // Create a node for the response
    pugi::xml_node element = rec.append_child("Response");
    element.append_attribute("transmitter") = GetTransmitterName().c_str();

    // Attach nodes for properties of the response
    ::AttachRsFloatNode(element, "start", StartTime(), false);
    AttachTextNode(element, "name", wave->GetName());

    // Render each interpolation point in turn
    std::vector<InterpPoint>::iterator i;
    for (i = points.begin(); i != points.end(); i++)
        RenderResponseXML(element, *i);
}

// Render a InterpPoint as CSV
void Response::RenderResponseCSV(std::ofstream &of, const InterpPoint &point)
{
    of << point.time << ", " << point.power << ", " << point.phase << ", " << point.doppler << "\n";
}

// Render the response to a CSV file
void Response::RenderCSV(std::ofstream &of)
{
    //Render each interpolation point
    std::vector<InterpPoint>::const_iterator i;
    for (i = points.begin(); i != points.end(); i++)
        RenderResponseCSV(of, *i);
}

// Add an interp point to the vector
void Response::AddInterpPoint(InterpPoint &point)
{
    // Check that points are being added in order
    if ((!points.empty()) && (point.time < points.back().time))
        throw std::logic_error("BUG: Interpolation points not being added in order");

    // This method does not need a mutex as only one thread owns any non-const Response object
    points.push_back(point);
}

// Return the number of Points
int Response::CountPoints() const
{
    return points.size();
}

// Render the response to an array
boost::shared_array<rsComplex> Response::RenderBinary(rsFloat& rate, unsigned int &size, rsFloat frac_win_delay)
{
    rate = wave->GetRate();
    return wave->Render(points, size, frac_win_delay);
}
