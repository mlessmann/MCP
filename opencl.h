#pragma once

#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>
#include <string>
#include <vector>

// OpenCL kernel code
static const std::string kernel_code = R"kernel_code(

__kernel void life(__global char *data, // TODO: Wie kann die Eingabefunktion übergeben werden?
                            uint  cols,
                            uint  rows)
{
    size_t pos = get_global_id(0);
    uint   col = (pos % (cols + 1));
    uint   row = (pos / (cols + 1)) + 1;
    
    if (col < cols)
        data[pos] = (col % (row + 1) == 0) ? 'X' : '_';
    else
        data[pos] = '\n'; // line feed in last column
}

)kernel_code";

typedef std::vector<std::vector<double>> vector_t;

template <typename Func>
vector_t jakobiOpenCL(vector_t     u,                // Eingabevector, mit Rand
                      Func         f,                // Eingabefunktion
                      const double h,                // Feinheit des Gitters
                      int          &iteration_count, // Out-Variable (profiling)
                      const double change_threshold, // Abbruch, wenn Änderung kleiner Wert
                      const int    max_iterations)   // Abbruch, wenn Anzahl der Iterationen erreicht
{
    try {
        // Get platform and device (old OpenCL 1.1 approach)
        std::vector<cl::Platform> platforms;
        cl::Platform::get(&platforms);
        std::vector<cl::Device> devices;
        try {
            platforms.at(0).getDevices(CL_DEVICE_TYPE_GPU, &devices);
        } catch (cl::Error& e) {
            if (e.err() == CL_DEVICE_NOT_FOUND) {
                std::cerr << "No GPU device found, falling back to CPU..." << std::endl;
                platforms.at(0).getDevices(CL_DEVICE_TYPE_CPU, &devices);
            } else throw;
        }
        
        // Build kernel
        cl::Context context(devices);
        cl::Program::Sources sources;
        sources.push_back({kernel_code.c_str(), kernel_code.length()});
        cl::Program program(context, sources);
        program.build(devices);
        cl::Kernel kernel(program, "gaussseidel");
        
        // Initialize command queue and buffers
        cl::CommandQueue queue(context, devices.at(0));
        cl::Buffer print_buffer(context, CL_MEM_READ_WRITE,
                                u.size() * u.size() * sizeof(double));
        // TODO: Copy to buffer
        
        // Set kernel parameters
        kernel.setArg(0, print_buffer);
        kernel.setArg(1, (cl_uint) cols);
        kernel.setArg(2, (cl_uint) rows);
        
        // Queue: Run kernel
        std::vector<cl::Event> kernel_finished(1);
        queue.enqueueNDRangeKernel(kernel, cl::NullRange,
            print_data.size(), cl::NullRange,
            nullptr, &kernel_finished.at(0));
        
        // Queue: Read results
        queue.enqueueReadBuffer(print_buffer, false,
            0, print_data.size() * sizeof(char), print_data.data(),
            &kernel_finished, nullptr);
        
        // Wait for queue to finish
        queue.finish();
    } catch (cl::Error& e) {
        std::cerr << "Failure: " << e.what() << " (Error code: " << e.err() << ")" << std::endl;
    } catch (std::exception& e) {
        std::cerr << "Failure: " << e.what() << std::endl;
    }
}
