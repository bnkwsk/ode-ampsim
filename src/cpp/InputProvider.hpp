#pragma once

#include <thread>
#include <atomic>

#include <sndfile.hh>
#include <sndfile.h>

#include "Sample.hpp"

template<typename P>
class InputProvider
{
    std::thread thread;
    P *processor;
    SndfileHandle in_file;
    bool interrupted;

    void run()
    {
        if(!in_file)
            return;

        int length = getFrames();

        float u_in;
        bool last = false;
        in_file.seek(0, SEEK_SET);
        for(int sample = 0; sample < length; ++sample)
        {
            in_file.readf(&u_in, 1);
            last = interrupted || (sample == length - 1);
            processor->setInput(Sample(u_in, last));
            if(last)
                break;
        }
        in_file = SndfileHandle();
    }

public:
    InputProvider(const char *inputPath) : interrupted(false)
    {
        in_file = SndfileHandle(inputPath);
    }
    double getFS()
    {
        return in_file.samplerate();
    }
    long long getFrames()
    {
        return in_file.frames();
    }
    void start(P *_processor)
    {
        processor = _processor;
        thread = std::thread(&InputProvider::run, this);
    }
    void join()
    {
        thread.join();
    }

    void interrupt()
    {
        interrupted = true;
    }
};
