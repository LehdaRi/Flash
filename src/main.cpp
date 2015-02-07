#include "Fourier.hpp"

#include <SFML/Graphics.hpp>
#include <cmath>
#include <complex>


sf::VertexArray getGraph(const float* signal,
                         unsigned size,
                         unsigned xstart,
                         unsigned ystart,
                         float xscale,
                         float yscale) {
    sf::VertexArray signalGraph(sf::Lines, (size-1)*2);

    for (auto i=0u; i<size-1; ++i) {
        signalGraph[2*i].position = sf::Vector2f(xstart+xscale*i, ystart-yscale*signal[i]);
        signalGraph[2*i].color = sf::Color(255, 255, 255);
        signalGraph[2*i+1].position = sf::Vector2f(xstart+xscale*(i+1), ystart-yscale*signal[i+1]);
        signalGraph[2*i+1].color = sf::Color(255, 255, 255);
    }

    return signalGraph;
}


int main(void) {
    const unsigned signalSize = 1024;

    sf::RenderWindow window(sf::VideoMode(1024, 600), "SFML window");
    window.setFramerateLimit(60);

    //  Signal
    float* signal = new float[signalSize];
    for (auto i=0u; i<signalSize; ++i) {
        //signal[i] = powf(cosf(PI2*i*(1.0f/signalSize)) * cosf(PI2*i*(3.0f/signalSize) * sinf(PI2*i*(7.0f/signalSize))), 2.0f);
        signal[i] = (i/64)%2;
    }
    auto signalGraph = getGraph(signal, signalSize, 0, 100, 1.0f, 25.0f);

    //  Signal spectri
    compf signalFourier[signalSize];
    float* signalAmpSpectrum = new float[signalSize];
    float* signalPhaseSpectrum = new float[signalSize];
    Fourier::FFT(signal, signalFourier, signalSize, true, 1);
    //Fourier::DFT(signal, signalFourier, signalSize, true);
    Fourier::getSpectrum(signalFourier, signalAmpSpectrum, signalPhaseSpectrum, signalSize);
    auto signalAmpSpectrumGraph = getGraph(signalAmpSpectrum, signalSize, 0, 300, 1.0f, 100.0f);
    auto signalPhaseSpectrumGraph = getGraph(signalPhaseSpectrum, signalSize, 0, 400, 1.0f, 5.0f);

    //  Fourier series (inverse DFT)
    float* invSignal = new float[signalSize];
    Fourier::invDFT(signalAmpSpectrum, signalPhaseSpectrum, invSignal, signalSize, true);
    auto invSignalGraph = getGraph(invSignal, signalSize, 0, 500, 1.0f, 25.0f);

    while (window.isOpen())
    {
        // Event processing
        sf::Event event;
        while (window.pollEvent(event))
        {
           // Request for closing the window
           if (event.type == sf::Event::Closed)
               window.close();
        }

        window.clear();

        window.draw(signalGraph);
        window.draw(signalAmpSpectrumGraph);
        window.draw(signalPhaseSpectrumGraph);
        window.draw(invSignalGraph);

        window.display();
    }

    delete[] signal;
    delete[] signalAmpSpectrum;
    delete[] signalPhaseSpectrum;
    delete[] invSignal;
}
