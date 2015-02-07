#ifndef FOURIER_HPP
#define FOURIER_HPP

#include <complex>
#include <iostream>
#include <vector>
#define PI2 6.28318530718

typedef std::complex<float> compf;

namespace Fourier {

    void DFT(const float* src, compf* output, unsigned size, bool printTime = false);
    void getSpectrum(const compf* fourier, float* amp, float* phase, unsigned size);
    void invDFT(const float* amp, const float* phase, float* dest, unsigned size, bool printTime);

    template <typename T>
	void FFT(const T* src, std::complex<float>* out, size_t size, bool printTime, size_t stride) {
        std::chrono::steady_clock::time_point start, end;

        if (printTime)
            start = std::chrono::steady_clock::now();

        if (size == 1) {
            out[0].real(src[0]);
            out[0].imag(0.0f);
            return;
        }

        FFT<T>(src, out, size/2, false, 2*stride);
        FFT<T>(src+stride, out+size/2, size/2, false, 2*stride);

        for (auto k=0u; k<size/2; ++k) {
            std::complex<float> t = out[k];
            std::complex<float> exponent;
            exponent.imag(-PI2*k/size);

            out[k] = t + std::exp(exponent)*out[k+size/2];
            out[k+size/2] = t - std::exp(exponent)*out[k+size/2];
        }
     
        if (printTime) {
            end = std::chrono::steady_clock::now();
            std::cout << "Calculating FFT took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;
        }
    }

} // namespace Fourier


#endif // FOURIER_HPP
