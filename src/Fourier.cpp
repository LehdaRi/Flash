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

void Fourier::getSpectrum(const compf* fourier, float* amp, float* phase, unsigned size, bool logarithmic) {
    for (auto k=0u; k<size; ++k) {
        amp[k] = logarithmic ? 10*log10(pow((std::norm(fourier[k]) / size), 2) / pow(1, 2)) : std::norm(fourier[k]) / float(size);
        phase[k] = std::arg(fourier[k]);

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
