#include "Fourier.hpp"

#include <SFML/Graphics.hpp>
#include <SFML/Audio.hpp>
#include <SFML/System.hpp>
#include <cmath>
#include <complex>
#include <vector>

/*
TODO
-
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

    virtual bool onProcessSamples(const sf::Int16* samples, std::size_t sampleCount)
    {
        /*
        Copies received samples to a vector. Thread safe.
        */

        sf::Lock lock(this->mutex);
        this->samples.resize(sampleCount);
        std::copy(samples, samples+sampleCount, this->samples.begin());
        return true;
    }

public:

    std::vector<sf::Int16> getChunk() {
        /*
        Returns a vector of audio samples. Thread safe.
        */

        sf::Lock lock(this->mutex);
        return this->samples;
    }
};


int main(void) {
    // Create an SFML window
    sf::RenderWindow window(sf::VideoMode(1024, 200), "SFML window");
    window.setFramerateLimit(60);    

    // Vectors for samples, F-transform, and spectrums
    std::vector<sf::Int16> samplesVect;
    std::vector<std::complex<float>> chunkFourierVect(0);
    std::vector<float> chunkAmpSpectVect(0);
    std::vector<float> chunkPhaseSpectVect(0);
    sf::VertexArray chunkAmpSpectGraph;

    // Pointers to the vectors for C-style array indexing
    sf::Int16* samples;
    std::complex<float>* chunkFourier;
    float* chunkAmpSpect;
    float* chunkPhaseSpect;

    // Other variables
    std::size_t chunkSize;
    std::size_t paddedChunkSize;

    // Start recording
    ChunkRecorder recorder;
    recorder.start(32770); // 16385, 32770

    while (window.isOpen())
    {
        // Event processing
        sf::Event event;
        while (window.pollEvent(event))
        {
           // Window close request
           if (event.type == sf::Event::Closed)
               window.close();
        }

        samplesVect = recorder.getChunk(); // Get a chunk of samples from the mic
        chunkSize = samplesVect.size();
        paddedChunkSize = pow(2, ceil(log2(chunkSize))); // Round to next highest power of 2
        //std::cout << chunkSize << " " << paddedChunkSize << std::endl;

        if(chunkSize > 0) {

            // Resize vectors
            samplesVect.resize(paddedChunkSize, 0); // This MUST be padded with zeroes!
            chunkFourierVect.resize(paddedChunkSize);
            chunkAmpSpectVect.resize(paddedChunkSize);
            //chunkPhaseSpectVect.resize(paddedChunkSize);

            // Assign pointers
            samples = &samplesVect[0];
            chunkFourier = &chunkFourierVect[0];
            chunkAmpSpect = &chunkAmpSpectVect[0];
            //chunkPhaseSpect = &chunkPhaseSpectVect[0];

            // Do the F-transform
            Fourier::FFT(samples, chunkFourier, paddedChunkSize, false, 1);
            //Get amplitudes and phases for each sinusoid. Because of symmetry, we can discard the second half.
            Fourier::getAmpSpectrum(chunkFourier, chunkAmpSpect, paddedChunkSize/2, true);

            /*for(auto val : chunkAmpSpectVect) {
                std::cout << val << std::endl;
            }*/

            // Get the graph
            chunkAmpSpectGraph = getGraph(chunkAmpSpect, paddedChunkSize, 0, 200, 0.5f, 1.0f);

            // Draw
            window.clear();
            //window.draw(borders);
            window.draw(chunkAmpSpectGraph);
        }

        window.display();
    }
    recorder.stop();
}
