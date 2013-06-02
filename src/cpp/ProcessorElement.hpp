#pragma once

#include<thread>

#include "LockingQueue.hpp"
#include "Sample.hpp"

class ProcessorElement
{
private:
    std::thread thread;
protected:
    typedef LockingQueue<Sample> QueueType;
    QueueType *input;
    QueueType *output;
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