#include "Fourier.hpp"

#include <cmath>
#include <complex>
#include <chrono>
#include <iostream>

void Fourier::DFT(const float* src, compf* output, unsigned size, bool printTime) {
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

    for (auto k=0u; k<size; ++k) {
        float re(0.0f), img(0.0f);

        for (auto n=0u; n<size; ++n) {
            float a = (float)(PI2*k*n)/size;
            re += src[n]*cosf(a);
            img += src[n]*sinf(-a);
        }

        output[k].real(re);
        output[k].imag(img);
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Calculating DFT took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;
}

void Fourier::getSpectrum(const compf* fourier, float* amp, float* phase, unsigned size) {
    for (auto k=0u; k<size; ++k) {
        amp[k] = sqrtf(powf(fourier[k].real(), 2.0f) + powf(fourier[k].imag(), 2.0f)) / size;
        phase[k] = atan2f(fourier[k].imag(), fourier[k].real());
    }
}

void Fourier::FFT(const float* src, compf* out, unsigned size, bool printTime, unsigned stride) {
    std::chrono::steady_clock::time_point start, end;

    if (printTime)
        start = std::chrono::steady_clock::now();

    if (size == 1) {
        out[0].real(src[0]);
        out[0].imag(0.0f);
        return;
    }

    FFT(src, out, size/2, false, 2*stride);
    FFT(src+stride, out+size/2, size/2, false, 2*stride);

    for (auto k=0u; k<size/2; ++k) {
        compf t = out[k];
        compf exponent;
        exponent.imag(-PI2*k/size);

        out[k] = t + std::exp(exponent)*out[k+size/2];
        out[k+size/2] = t - std::exp(exponent)*out[k+size/2];
    }

    if (printTime) {
        end = std::chrono::steady_clock::now();
        std::cout << "Calculating FFT took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;
    }
}

void Fourier::invDFT(const float* amp, const float* phase, float* dest, unsigned size, bool printTime = false) {
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

    for(auto k = 0u; k < size; k++) {
        for(auto n = 0u; n < size; n++) {
            dest[k] += amp[n]*cosf(PI2*k*(float)n/size+phase[n]);
        }
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Calculating IDFT took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;
}

/*void Fourier::dft2d(const float* src, float* amp, float* phase, unsigned size) {
    for (auto k=0u; k<size; ++k) {
        float re(0.0f), img(0.0f);

        for (auto n=0u; n<size; ++n) {
            float a = (float)(PI2*k*n)/size;
            re += src[n]*cosf(a);
            img += src[n]*sinf(-a);
        }

        amp[k] = sqrtf(powf(re, 2.0f) + powf(img, 2.0f)) / size;
        phase[k] = atan2f(img, re);
    }
}*/
