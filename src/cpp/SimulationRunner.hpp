#ifndef SIMULATIONRUNNER_HPP
#define SIMULATIONRUNNER_HPP

#include "GuitarAmp.hpp"
#include "InputProvider.hpp"
#include "OutputCollector.hpp"

class SimulationRunner
{
    const char *outputPath;
    InputProvider<GuitarAmp> *input;
    OutputCollector<GuitarAmp> *output;
    GuitarAmp *amp;

    std::thread thread;
    bool running;

    void _run()
    {
        running = true;
        amp->run();
        input->start(amp);
        output->start(outputPath);
        output->join();
        amp->join();
        input->join();
        running = false;
    }

    public:
        void run(const char *inputPath, const char *impulseResponsePath, const char *_outputPath)
        {
            outputPath = _outputPath;
            input = new InputProvider<GuitarAmp>(inputPath);
            amp = new GuitarAmp(input->getFS(), impulseResponsePath);
            output = new OutputCollector<GuitarAmp>(*amp);
            thread = std::thread(&SimulationRunner::_run, this);
        }

        void join()
        {
            thread.join();
        }

        double getProgress()
        {
            return (double)output->getCollected() / input->getFrames();
        }
        bool isRunning()
        {
            return running;
        }
};

#endif