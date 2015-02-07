#include "Fourier.hpp"

#include <SFML/Graphics.hpp>
#include <SFML/Audio.hpp>
#include <SFML/System.hpp>
#include <cmath>
#include <complex>
#include <vector>

/*
TODO

- Paddaa lähimpään kakkosen potenssiin
- Graafifunktio inteille
*/

sf::VertexArray getGraph(const float* signal, unsigned size, unsigned xstart, unsigned ystart, float xscale, float yscale) {
    sf::VertexArray signalGraph(sf::Lines, (size-1)*2);

    for (auto i=0u; i<size-1; ++i) {
        signalGraph[2*i].position = sf::Vector2f(xstart+xscale*i, ystart-yscale*signal[i]);
        signalGraph[2*i].color = sf::Color(255, 255, 255);
        signalGraph[2*i+1].position = sf::Vector2f(xstart+xscale*(i+1), ystart-yscale*signal[i+1]);
        signalGraph[2*i+1].color = sf::Color(255, 255, 255);
    }

    return signalGraph;
}

class ChunkRecorder: public sf::SoundRecorder
{
    std::vector<sf::Int16> samples;
    sf::Mutex mutex;

    virtual bool onStart() // optional
    {
        return true;
    }

    virtual bool onProcessSamples(const sf::Int16* samples, std::size_t sampleCount)
    {
        sf::Lock lock(this->mutex);
        this->samples.resize(sampleCount);
        std::copy(samples, samples+sampleCount, this->samples.begin());
        return true;
    }

    virtual void onStop() // optional
    {
        // clean up whatever has to be done after the capture is finished
    }

public:

    std::vector<sf::Int16> getChunk() {
        sf::Lock lock(this->mutex);
        return this->samples;
    }
};

int main(void) {

    sf::RenderWindow window(sf::VideoMode(1024, 800), "SFML window");
    window.setFramerateLimit(60);    

    // Vectors
    std::vector<std::complex<float>> chunkFourierVect(0);
    std::vector<float> chunkAmpSpectVect(0);
    std::vector<float> chunkPhaseSpectVect(0);

    // Pointers to beginning of vectors
    std::complex<float>* chunkFourier;
    sf::Int16* samples;
    float* chunkAmpSpect;
    float* chunkPhaseSpect;
    float* inverse;

    // Other variables
    std::size_t chunkSize;
    unsigned paddedChunkSize;
    std::vector<float> inverseVect;

    // Start recording
    ChunkRecorder recorder;
    recorder.start(16385);
    //recorder.start(1024);

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

        auto samplesVect = recorder.getChunk();
        chunkSize = samplesVect.size();
        paddedChunkSize = pow(2, ceil(log2(chunkSize))); // Round to next highest power of 2

        //std::cout << chunkSize << " " << paddedChunkSize << std::endl;
        if(chunkSize > 0) {
            samplesVect.resize(paddedChunkSize, 0);
            samples = &samplesVect[0];

            chunkFourierVect.resize(paddedChunkSize);
            chunkFourier = &chunkFourierVect[0];

            chunkAmpSpectVect.resize(paddedChunkSize);
            chunkPhaseSpectVect.resize(paddedChunkSize);
            chunkAmpSpect = &chunkAmpSpectVect[0];
            chunkPhaseSpect = &chunkPhaseSpectVect[0];

            Fourier::FFT(samples, chunkFourier, paddedChunkSize, false, 1);
            Fourier::getSpectrum(chunkFourier, chunkAmpSpect, chunkPhaseSpect, paddedChunkSize/4);

            /*inverseVect.resize(1024);
            inverse = &inverseVect[0];
            Fourier::invDFT(chunkAmpSpect, chunkPhaseSpect, inverse, paddedChunkSize, false);*/

            auto chunkAmpSpectGraph = getGraph(chunkAmpSpect, paddedChunkSize, 0, 750, 2.0f, 0.0001f);
            //auto inverseGraph = getGraph(inverse, 1024, 0, 500, 1.0f, 1.0f);
            //auto phasePhaseSpectGraph = getGraph(chunkPhaseSpect, paddedChunkSize, 0, 300, 1.0f, 100.0f);

            for(auto val : chunkFourierVect) {
                //std::cout << val << std::endl;
            }
            window.clear();

            window.draw(chunkAmpSpectGraph);
            //window.draw(inverseGraph);
            //makewindow.draw(invSignalGraph);
        }

        window.display();
    }
    recorder.stop();
}
