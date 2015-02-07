#ifndef FOURIER_HPP
#define FOURIER_HPP


#define PI2 6.28318530718


namespace Fourier {

    void DFT(const float* src, float* reArr, float* imgArr, unsigned size, bool printTime = false);
    void getSpectrum(const float* re, const float* img, float* amp, float* phase, unsigned size);
    void FFT(const float* src, float* reArr, float* imgArr, unsigned size, bool printTime = false, unsigned stride = 1u);
    void invDFT(const float* amp, const float* phase, float* dest, unsigned size, bool printTime);

    //void dft_2d(const float* src, float* amp, float* phase, unsigned xsize, unsigned ysize);

} // namespace Fourier


#endif // FOURIER_HPP
