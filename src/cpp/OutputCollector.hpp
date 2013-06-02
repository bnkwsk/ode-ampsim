#pragma once

#include <thread>

#include <sndfile.hh>
#include <sndfile.h>

#include "Sample.hpp"

template<typename P>
class OutputCollector
{
    std::thread thread;
    P &processor;
    const char *path;
    bool process;
    long long collected;
    
    void run()
    {
        double f_s = processor.getFS();
        SndfileHandle out_file(path, SFM_WRITE, SF_FORMAT_WAV | SF_FORMAT_PCM_16, 1, int(f_s));

        if(!out_file)
            return;
        Sample out;
        float u_out;
        bool last = false;
        collected = 0;

        while(!last) {
            out = processor.getOutput();
            u_out = out.getValue();
            if(u_out > 1.0f || u_out < -1.0f)
            {
                if(u_out > 0.0f)
                    u_out = 1.0f;
                else
                    u_out = -1.0f;
            }
            last = out.isLast();
            out_file.write(&u_out, 1);
            ++collected;
        }
    }
public:
    OutputCollector(P &_processor) : processor(_processor)
    {
        //
    }
    void start(const char *pth)
    {
        path = pth;
        thread = std::thread(&OutputCollector::run, this);
    }
    long long getCollected()
    {
        return collected;
    }
    void join()
    {
        thread.join();
    }
    void stop()
    {
        process = false;
    }
};