#include <deque>
#include <vector>
#include <algorithm>

#include <sndfile.hh>
#include <sndfile.h>
#include <viennacl/vector.hpp>
#include <viennacl/linalg/inner_prod.hpp>

#include "ProcessorElement.hpp"

class CabinetSimulator : public ProcessorElement
{
    protected:
        typedef viennacl::vector<float> GPUVector;
        typedef viennacl::scalar<float> GPUScalar;
        std::deque<float> inputBuffer;
        std::vector<float> cpuInput;
        long long impulseLength;
        float *impulse;
        GPUVector gpuImpulse;
        GPUVector gpuInput;
        GPUScalar gpuResult;

        void run()
        {
            Sample in, out;
            double u_in;
            bool last = false;

            while(!last)
            {
                in = input->pop();
                u_in = in.getValue();
                last = in.isLast();
                inputBuffer.pop_back();
                inputBuffer.push_front(u_in);
                output->push(Sample(process(), last));
            }
        }

        float process()
        {
            for(long long i = 0; i < impulseLength; ++i)
                cpuInput[i] = inputBuffer[i];
            viennacl::fast_copy(cpuInput, gpuInput);
            gpuResult = viennacl::linalg::inner_prod(gpuInput, gpuImpulse);
            return float(gpuResult);
        }
    public:
        CabinetSimulator(const char *path) : ProcessorElement()
        {
            SndfileHandle impulseFile = SndfileHandle(path);
            impulseLength = impulseFile.frames();
            impulse = new float[impulseLength];
            impulseFile.readf(impulse, impulseLength);
            for(long long i = 0; i < impulseLength; ++i)
                inputBuffer.push_front(0.0f);
            gpuImpulse = GPUVector(impulseLength);
            viennacl::copy(std::vector<float>(impulse, impulse + impulseLength), gpuImpulse);
            gpuInput = GPUVector(impulseLength);
            cpuInput = std::vector<float>(impulseLength);
        }
        ~CabinetSimulator()
        {
            delete[] impulse;
        }
};