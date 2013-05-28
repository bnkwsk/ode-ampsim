#pragma once

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
        delete input;
        delete output;
        delete amp;
        running = false;
    }

    public:
        SimulationRunner() : running(false)
        {
            //
        }

        void run(const char *inputPath, const char *impulseResponsePath, const char *_outputPath)
        {
            outputPath = _outputPath;
            input = new InputProvider<GuitarAmp>(inputPath);
            amp = new GuitarAmp(input->getFS(), impulseResponsePath);
            output = new OutputCollector<GuitarAmp>(*amp);
            if(thread.joinable())
                thread.join();
            thread = std::thread(&SimulationRunner::_run, this);
        }

        void join()
        {
            //  if(running)
                thread.join();
        }

        void interrupt()
        {
            if(running)
                input->interrupt();
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
