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
    bool running, interrupted;

    void _run()
    {
        running = true;
        interrupted = false;
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
        if(thread.joinable())
            thread.join();
    }

    void interrupt()
    {
        interrupted = true;
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
    bool isInterrupted()
    {
        return interrupted;
    }
};
