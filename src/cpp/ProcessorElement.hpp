#pragma once

#include<thread>

#include "LockingQueue.hpp"
#include "Sample.hpp"

class ProcessorElement
{
    protected:
        typedef LockingQueue<Sample> QueueType;
        QueueType *input;
        QueueType *output;
        std::thread thread;

        virtual void run() {}
        
    public:
        void setInput(QueueType *in)
        {
            input = in;
        }
        void setOutput(QueueType *out)
        {
            output = out;
        }
        void start()
        {
            thread = std::thread(&ProcessorElement::run, this);
        }
        void join()
        {
            thread.join();
        }
};