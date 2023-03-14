/// Implementation of rsdebug.h - Debug messaging
/// Marc Brooker, 20 March 2006
/// Edited by Yaaseen Martin, 27 August 2019

#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstdarg>
#include "rsdebug.cuh"

// For boost::mutex
#include <boost/thread/mutex.hpp>

// The current debug level
rsDebug::Level debug_level = rsDebug::RS_VERY_VERBOSE;

// Use a mutex for log information printing (to stop messages from getting mangled)
boost::mutex debugMutex;

// Print out a debug message including the line and file name
void rsDebug::print(const rsDebug::Level level, const std::string &str, const std::string file, const int line) {
    if (level >= debug_level) {
        boost::mutex::scoped_lock lock(debugMutex); // Lock the mutex
        std::ostringstream oss;
        oss << "[" << file << " " << line << "] ";
        oss << str << "\n";
        std::cerr << oss.str();
        // Mutex will automatically be unlocked here by scoped_lock
    }
}

// Formatted print of the current debug level - does not include filename and line
// Uses the cstdarg variable arguments system and the vfprintf function to handle the arguments
void rsDebug::printf(const rsDebug::Level level, const char *format, ...)
{
    if (level >= debug_level) {
        boost::mutex::scoped_lock lock(debugMutex); // Lock the mutex
        va_list ap;
        va_start(ap, format);
        vfprintf(stderr, format, ap);
        va_end(ap);
        // Mutex will automatically be unlocked here by scoped_lock
    }
}

// See comments for printf(Level, char *)
void rsDebug::printf(const rsDebug::Level level, const std::string &format, ...){
    if (level >= debug_level) {
        boost::mutex::scoped_lock lock(debugMutex); // Lock the mutex
        va_list ap;
        va_start(ap, format);
        vfprintf(stderr, format.c_str(), ap);
        va_end(ap);
        // Mutex will automatically be unlocked here by scoped_lock
    }
}

// Set the current debug level
void rsDebug::setDebugLevel(rsDebug::Level level) {
    if (level <= rsDebug::RS_EXTREMELY_CRITICAL)
        debug_level = level;
}
