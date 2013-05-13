#pragma once

#include <deque>
#include <vector>
#include <algorithm>

#include <sndfile.hh>
#include <sndfile.h>
#include <viennacl/vector.hpp>
#include <viennacl/matrix.hpp>
#include <viennacl/fft.hpp>
#include <viennacl/linalg/inner_prod.hpp>
#include "viennacl/ocl/backend.hpp"

#include "ProcessorElement.hpp"
#include "Sample.hpp"

class CabinetSimulator : public ProcessorElement
{
protected:
    typedef float ScalarType;
    typedef std::vector<ScalarType> CPUVector;
    typedef viennacl::matrix<ScalarType> GPUMatrix;
    typedef viennacl::vector<ScalarType> GPUVector;
    typedef viennacl::scalar<ScalarType> GPUScalar;

    const int blockSize = 1024; // it has to be a power of 2
    int filterPartQuantity;

    CPUVector cpuInput;
    CPUVector cpuOutput;
    CPUVector zeroBlock;
    
    GPUMatrix filterSpectra;
    GPUVector gpuInput;
    GPUVector gpuOutput;

    void run()
    {
        Sample in, out;
        double u_in;
        bool last = false;

        viennacl::ocl::program &program = viennacl::ocl::current_context().get_program("ComplexPointwiseMultiplyAdd");
        viennacl::ocl::kernel &ComplexPointwiseMultiplyAdd = program.get_kernel("ComplexPointwiseMultiplyAdd");
        ComplexPointwiseMultiplyAdd.global_work_size(0, blockSize);
        ComplexPointwiseMultiplyAdd.local_work_size(0, 2);

        while(!last)
        {
            // pack the input data [previousInput:newInput]
            // -- move the right hand side of the input buffer to the left
            for(int i = 0; i < blockSize; i += 2)
            {
                cpuInput[i] = cpuInput[i + blockSize];
            }
            // -- pack the new input data
            for(int i = blockSize; i < 2 * blockSize; i += 2)
            {
                in = input->pop();
                u_in = in.getValue();
                last = in.isLast();
                cpuInput[i] = u_in;
                if(last)
                {
                    for(; i < 2 * blockSize; i += 2)
                        cpuInput[i] = 0.0f;
                    break;
                }
            }
            processBlock(ComplexPointwiseMultiplyAdd);
            for(int i = blockSize; i < 2 * blockSize; i += 2)
                output->push(Sample(cpuOutput[i] / 20.0, last ? (i == 2 * blockSize - 2) : false));
        }
    }

    void processBlock(viennacl::ocl::kernel &kernel)
    {
        viennacl::fast_copy(cpuInput, gpuInput);

        viennacl::inplace_fft(gpuInput);

        gpuOutput.clear();

        // compute the COMPLEX product
        viennacl::ocl::enqueue(kernel(gpuInput, filterSpectra, gpuOutput, static_cast<cl_uint>(filterSpectra.size1())));

        viennacl::inplace_ifft(gpuOutput);

        viennacl::fast_copy(gpuOutput, cpuOutput);
    }
public:
    CabinetSimulator(const char *path) : ProcessorElement()
    {
        std::vector<ScalarType> impulse;
        int filterLength;
        std::vector<GPUVector> filterParts;

        SndfileHandle impulseFile = SndfileHandle(path);
        filterLength = impulseFile.frames();
        ScalarType *_impulse = new ScalarType[filterLength];
        impulseFile.readf(_impulse, filterLength);
        impulse = CPUVector(_impulse, _impulse + filterLength);
        cpuOutput = CPUVector(2 * blockSize);
        gpuInput = GPUVector(2 * blockSize);
        gpuOutput = GPUVector(2 * blockSize);

        for(int i = 0; i < 2 * blockSize; ++i)
        {
            zeroBlock.push_back(0.0f);
            cpuInput.push_back(0.0f);
        }
        
        filterPartQuantity = 2 * ceil(filterLength / (double)blockSize);
        filterSpectra = GPUMatrix(filterPartQuantity, 2 * blockSize);

        // fill the filter parts in and perform the FFT
        CPUVector cpuPart = CPUVector(2 * blockSize);
        for(int i = blockSize; i < 2 * blockSize; ++i)
            cpuPart[i] = 0.0f;
        for(int p = 0; p < filterPartQuantity; ++p)
        {
            for(int i = 0; i < blockSize; ++i)
                if(i % 2 == 0)
                    if(p * blockSize / 2 + i / 2 < filterLength)
                        cpuPart[i] = impulse[p * blockSize / 2 + i / 2];
                    else
                        cpuPart[i] = 0.0f;
                else
                    cpuPart[i] = 0.0f;
            filterParts.push_back(GPUVector(2 * blockSize));
            GPUVector &lastPart = filterParts[filterParts.size() - 1];
            viennacl::fast_copy(cpuPart, lastPart);
            viennacl::inplace_fft(lastPart);
            for(int i = 0; i < 2 * blockSize; ++i)
                filterSpectra(p, i) = lastPart[i];
        }
    }
};