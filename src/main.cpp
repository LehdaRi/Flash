#include <SFML/Graphics.hpp>
#include <SFML/Audio.hpp>
#include <SFML/System.hpp>
#include <unsupported/Eigen/FFT>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <complex>
#include <vector>
#include <deque>


sf::VertexArray getGraph(const float* signal, unsigned size, unsigned xstart, unsigned ystart, float xscale, float yscale) {
    if (size < 1)
        return sf::VertexArray(sf::Lines, 0);

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
public:

    //  fetch a chunk of audio data from sample deque
    std::vector<float>& drawSamples(unsigned chunkSize, unsigned newDataAmount) {
        sf::Lock lock(mutex_);

        for (auto i=0u; i<newDataAmount && samples_.size()>0; ++i) {
            timevec_.push_back(samples_.front());
            samples_.pop_front();
        }

        if (timevec_.size() > chunkSize)
            timevec_.erase(timevec_.begin(), timevec_.begin()+(timevec_.size()-chunkSize));

        return timevec_;
    }

    //  returns amplitude spectrum of last fetched chunk
    void getAmpSpectrum(unsigned size, std::vector<float>& spectrum, bool logarithmic) {
        if (timevec_.size() != size)
            return;

        fft_.fwd(freqvec_, timevec_);

        spectrum.clear();
        spectrum.reserve(freqvec_.size());

        for (auto& f : freqvec_)
            spectrum.push_back( logarithmic ? 10*log10(pow((std::norm(f)/freqvec_.size()), 2) / pow(512, 2)) : std::norm(f)/freqvec_.size() );
    }

private:
    std::deque<float> samples_;
    sf::Mutex mutex_;
    Eigen::FFT<float> fft_;
    std::vector<float> timevec_;
    std::vector<std::complex<float>> freqvec_;

    virtual bool onProcessSamples(const sf::Int16* samples, std::size_t sampleCount) {
        sf::Lock lock(mutex_);

        for (auto i=0u; i<sampleCount; ++i)
            samples_.push_back((float)samples[i]);

        if (samples_.size() > 441000)
            samples_.clear();

        return true;
    }
};


class BufferFilter {
public:
    BufferFilter(unsigned size) :
        min_(0.0f), max_(0.0f)
    {
        buffer_.resize(size, 0.0f);
    }

    void pushSample(float sample) {
        buffer_.push_back(sample);

        float front = buffer_.front();
        buffer_.pop_front();

        if (abs(min_ - front) < 0.00001f)
            updateMin();

        if (abs(max_ - front) < 0.00001f)
            updateMax();
    }

    float getMin(void) {
        return min_;
    }

    float getMax(void) {
        return max_;
    }

    float getAverage(void) {
        float avg = 0.0f;
        for (auto& sample : buffer_)
            avg += sample;
        return avg / buffer_.size();
    }

    float scale(float in) {
        return (in - min_) / (max_ - min_);
    }

private:
    std::deque<float> buffer_;
    float min_;
    float max_;

    void updateMin(void) {
        min_ = buffer_.front();
        for (auto& sample : buffer_)
            min_ = sample < min_ ? sample : min_;
    }

    void updateMax(void) {
        max_ = buffer_.front();
        for (auto& sample : buffer_)
            max_ = sample > max_ ? sample : max_;
    }
};


void readConfigFile(int& triggerFreqBegin, int& triggerFreqNum, int& triggerFreqStrife, float& flashThreshold, float& gamma) {
    FILE* f = fopen("flashconfig.txt", "r");

    if (f) {
        if (fscanf(f, "%i %i %i %f %f", &triggerFreqBegin, &triggerFreqNum, &triggerFreqStrife, &flashThreshold, &gamma) < 5)
            std::cout << "Unable to read config values from flashconfig.txt" << std::endl;
        fclose(f);
        printf("flashconfig.txt read succesfully:\ntriggerFreqBegin: %i\ntriggerFreqNum: %i\ntriggerFreqStrife: %i\nflashThreshold: %0.3f\ngamma: %0.3f\n",
               triggerFreqBegin, triggerFreqNum, triggerFreqStrife, flashThreshold, gamma);
    }
    else
        std::cout << "Unable to read flashconfig.txt" << std::endl;
}


void writeConfigFile(int triggerFreqBegin, int triggerFreqNum, int triggerFreqStrife, float flashThreshold, float gamma) {
    FILE* f = fopen("flashconfig.txt", "w+");

    if (f) {
        fprintf(f, "%i %i %i %f %f", triggerFreqBegin, triggerFreqNum, triggerFreqStrife, flashThreshold, gamma);
        fclose(f);
    }
    else
        std::cout << "Unable to write flashconfig.txt" << std::endl;
}


int main(void) {
    // Create an SFML window
    sf::RenderWindow window(sf::VideoMode(1024, 768), "SFML window",  sf::Style::Fullscreen);
    window.setFramerateLimit(60);
    window.setVerticalSyncEnabled(true);

    // Vectors for samples, F-transform, and spectrums
    std::vector<float> chunkAmpSpect(0);
    sf::VertexArray chunkAmpSpectGraph;

    // Start recording
    ChunkRecorder recorder;
    recorder.start(); // 16385, 32770

    BufferFilter ampFilter(60);
    float lastTriggerValue = 0.0f;

    unsigned rgbMasks[] = {
        0xff0000,
        0xffff00,
        0x00ff00,
        0x00ffff,
        0x0000ff,
        0xff00ff
    };
    int activeMask = 0;

    //  control variables
    int triggerFreqBegin = 5;
    int triggerFreqNum = 10;
    int triggerFreqStrife = 1;
    float flashThreshold = 1.0f;
    float gamma = 1.0f;
    readConfigFile(triggerFreqBegin, triggerFreqNum, triggerFreqStrife, flashThreshold, gamma);

    bool showSpectrum = false;

    while (window.isOpen())
    {
        // Event processing
        sf::Event event;
        while (window.pollEvent(event))
        {
            // Window close request
            switch (event.type) {
            case sf::Event::Closed:
                window.close();
            break;
            case sf::Event::KeyPressed:
                switch (event.key.code) {
                case sf::Keyboard::Escape:
                    window.close();
                break;
                case sf::Keyboard::Space:
                    showSpectrum = !showSpectrum;
                break;
                case sf::Keyboard::Q:
                    if (triggerFreqBegin > 0)
                        printf("triggerFreqBegin: %i\n", --triggerFreqBegin);
                break;
                case sf::Keyboard::W:
                    if (triggerFreqBegin < 1000)
                        printf("triggerFreqBegin: %i\n", ++triggerFreqBegin);
                break;
                case sf::Keyboard::A:
                    if (triggerFreqNum > 1)
                        printf("triggerFreqNum: %i\n", --triggerFreqNum);
                break;
                case sf::Keyboard::S:
                    if (triggerFreqNum < 100)
                        printf("triggerFreqNum: %i\n", ++triggerFreqNum);
                break;
                case sf::Keyboard::Z:
                    if (triggerFreqStrife > 1)
                        printf("triggerFreqStrife: %i\n", --triggerFreqStrife);
                break;
                case sf::Keyboard::X:
                    if (triggerFreqStrife < 100)
                        printf("triggerFreqStrife: %i\n", ++triggerFreqStrife);
                break;
                case sf::Keyboard::E:
                    if (flashThreshold > -100.0f) {
                        flashThreshold -= 1.0f;
                        printf("flashThreshold: %f\n", flashThreshold);
                    }
                break;
                case sf::Keyboard::R:
                    if (flashThreshold < 100.0f) {
                        flashThreshold += 1.0f;
                        printf("flashThreshold: %f\n", flashThreshold);
                    }
                break;
                case sf::Keyboard::D:
                    if (gamma > 0.1f) {
                        gamma -= 0.1f;
                        printf("gamma: %f\n", gamma);
                    }
                break;
                case sf::Keyboard::F:
                    if (gamma < 3.0f) {
                        gamma += 0.1f;
                        printf("gamma: %f\n", gamma);
                    }
                break;
                default:
                break;
                }
            break;
            default:
            break;
            }
        }

        recorder.drawSamples(2205, 800); //  draw 800 new samples every frame (44100Hz/60Fps = 735 new samples/frame)
        recorder.getAmpSpectrum(2205, chunkAmpSpect, true);

        if(chunkAmpSpect.size() > 0) {


            // Draw
            float a = 0.0f;
            for (auto i=0; i<triggerFreqNum; i++) {
                a += chunkAmpSpect.at(triggerFreqBegin + i*triggerFreqStrife);
            }
            a /= triggerFreqNum;

            ampFilter.pushSample(a);

            a = a*0.6+lastTriggerValue*0.4;
            lastTriggerValue = a;

            float b = (a-ampFilter.getAverage()) / (ampFilter.getMax()-ampFilter.getAverage()+flashThreshold);
            if (b < 0.0f) b = 0.0f;
            if (b > 1.0f) b = 1.0f;
            b = pow(b, 1.0f/gamma);
            unsigned char c = b*255;

            unsigned mask = rgbMasks[activeMask];

            static int maskTimer = 0;
            maskTimer++;
            if (maskTimer > 30) {
                activeMask = (activeMask+1)%6;
                maskTimer = 0;
            }

            window.clear(sf::Color(((c << 16) & mask) >> 16, ((c << 8) & mask) >> 8, c & mask));

            if (showSpectrum) {
                chunkAmpSpectGraph = getGraph(&chunkAmpSpect[0], chunkAmpSpect.size(), 0, 700, 1.0f, 2.0f);
                window.draw(chunkAmpSpectGraph);

                //  draw trigger visualization
                sf::VertexArray triggerLines(sf::Lines, triggerFreqNum*2);

                for (auto i=0; i<triggerFreqNum; ++i) {
                    triggerLines[2*i].position = sf::Vector2f(triggerFreqBegin + i*triggerFreqStrife, 450);
                    triggerLines[2*i].color = sf::Color(255, 255, 255);
                    triggerLines[2*i+1].position = sf::Vector2f(triggerFreqBegin + i*triggerFreqStrife, 768);
                    triggerLines[2*i+1].color = sf::Color(255, 255, 255);
                }

                window.draw(triggerLines);
                }
        }

        window.display();
    }

    recorder.stop();
    writeConfigFile(triggerFreqBegin, triggerFreqNum, triggerFreqStrife, flashThreshold, gamma);
}
