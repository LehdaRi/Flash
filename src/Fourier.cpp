#include "Fourier.hpp"

#include <cmath>
#include <chrono>
#include <iostream>


void Fourier::DFT(const float* src, float* reArr, float* imgArr, unsigned size, bool printTime) {
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

    for (auto k=0u; k<size; ++k) {
        float re(0.0f), img(0.0f);

        for (auto n=0u; n<size; ++n) {
            float a = (float)(PI2*k*n)/size;
            re += src[n]*cosf(a);
            img += src[n]*sinf(-a);
        }

        reArr[k] = re;
        imgArr[k] = img;
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Calculating DFT took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;
}

void Fourier::getSpectrum(const float* re, const float* img, float* amp, float* phase, unsigned size) {
    for (auto k=0u; k<size; ++k) {
        amp[k] = sqrtf(powf(re[k], 2.0f) + powf(img[k], 2.0f)) / size;
        phase[k] = atan2f(img[k], re[k]);
    }
}

void Fourier::FFT(const float* src, float* reArr, float* imgArr, unsigned size, bool printTime, unsigned stride) {
    std::chrono::steady_clock::time_point start, end;

    if (printTime)
        start = std::chrono::steady_clock::now();

    if (size == 1) {
        reArr[0] = src[0];
        imgArr[0] = 0.0f;
        return;
    }

    FFT(src, reArr, imgArr, size/2, false, 2*stride);
    FFT(src+stride, reArr+size/2, imgArr+size/2, size/2, false, 2*stride);

    for (auto k=0u; k<size/2; ++k) {
        float tre = reArr[k];
        float timg = imgArr[k];
        float arg = (float)(PI2*k/size);
        float a = cosf(arg);
        float b = -sinf(arg);
        float c = reArr[k+size/2];
        float d = imgArr[k+size/2];

        reArr[k] = tre + a*c-b*d;
        imgArr[k] = timg + b*c+a*d;

        reArr[k+size/2] = tre - a*c-b*d;
        imgArr[k+size/2] = timg - b*c+a*d;

        /*
        float a = (float)(PI2/size*k);
        float tre = reArr[k];
        float timg = imgArr[k];

        reArr[k] = tre + cosf(a)*reArr[k+size/2]-sinf(-a)*imgArr[k+size/2];
        imgArr[k] = timg + sinf(-a)*reArr[k+size/2]+cosf(a)*imgArr[k+size/2];

        reArr[k+size/2] = tre - cosf(a)*reArr[k+size/2]-sinf(-a)*imgArr[k+size/2];
        imgArr[k+size/2] = timg - sinf(-a)*reArr[k+size/2]+cosf(a)*imgArr[k+size/2];
        */
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
