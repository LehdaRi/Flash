#ifndef FOURIER_HPP
#define FOURIER_HPP

#include <complex>
#define PI2 6.28318530718

typedef std::complex<float> compf;

namespace Fourier {

    void DFT(const float* src, compf* output, unsigned size, bool printTime = false);
    void getSpectrum(const compf* fourier, float* amp, float* phase, unsigned size);
    void FFT(const float* src, compf* output, unsigned size, bool printTime = false, unsigned stride = 1u);
    void invDFT(const float* amp, const float* phase, float* dest, unsigned size, bool printTime);

    //void dft_2d(const float* src, float* amp, float* phase, unsigned xsize, unsigned ysize);

} // namespace Fourier


#endif // FOURIER_HPP
