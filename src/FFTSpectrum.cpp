
class FFTSpectrum : public sf::SoundRecorder
{
    virtual bool onStart()
    {
        // INIT
        return true;
    }

    virtual bool onProcessSamples(const Int16* samples, std::size_t sampleCount)
    {
        
        return true;
    }

    virtual void onStop() // optional
    {
        // clean up whatever has to be done after the capture is finished
        ...
    }
}