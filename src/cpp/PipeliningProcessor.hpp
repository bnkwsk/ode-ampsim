#ifndef PIPELININGPROCESSOR_H
#define PIPELININGPROCESSOR_H

#include <iostream>
#include <vector>

#include <thread>

#include "LockingQueue.hpp"


template<typename SampleType, typename ElementType>
class PipeliningProcessor
{
protected:
    std::vector<ElementType *> elements;
    std::vector<LockingQueue<SampleType> *> queues;
    LockingQueue<SampleType> *input;
    LockingQueue<SampleType> *output;
public:
    PipeliningProcessor()
    {
        //
    }

    ~PipeliningProcessor()
    {
        for(LockingQueue<SampleType> *queue : queues)
            delete queue;
    }

    void setup()
    {
        if(!elements.size())
            return;
        queues.push_back(new LockingQueue<SampleType>());
        elements[0]->setInput(queues[0]);
        for(int i = 1; i < elements.size(); ++i)
        {
            queues.push_back(new LockingQueue<SampleType>());
            elements[i - 1]->setOutput(queues[i]);
            elements[i]->setInput(queues[i]);
        }
        queues.push_back(new LockingQueue<SampleType>());
        elements[elements.size() - 1]->setOutput(queues[queues.size() - 1]);
        input = queues[0];
        output = queues[queues.size() - 1];
    }

    void run()
    {
        for(ElementType *element : elements)
            element->start();
    }

    void join()
    {
        for(ElementType *element : elements)
            element->join();
    }

    void setInput(SampleType in)
    {
        input->push(in);
    }

    SampleType getOutput()
    {
        return output->pop();
    }
};

#endif