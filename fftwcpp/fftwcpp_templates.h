/// fftwcpp_templates.h - Preprocessor magic to link to the right fftw precision
/// Marc Brooker, 12 June 2006
/// Edited by Yaaseen Martin, 27 August 2019

#if RS_FLOAT_LONG_DOUBLE == 1
#define fftw_ fftwf_
#elif RS_FLOAT_FLOAT == 1
#define fftw_ fftwl_
#endif
