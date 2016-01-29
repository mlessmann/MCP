#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>
#include <iostream>
#include <string>

int main(int argc, char **argv) {
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
        
        // Some debug output
        std::cout <<   "Platform:      " << platforms.at(0).getInfo<CL_PLATFORM_NAME>()
                  << "\nDevice 0:      " << devices.at(0).getInfo<CL_DEVICE_NAME>()
                  << "\nCompute Units: " << devices.at(0).getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>()
                      << " (x " << devices.at(0).getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>() << ")"
                  << "\nGlobal Memory: " << (devices.at(0).getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>() >> 20)
                      << " MByte" << std::endl;
    } catch (cl::Error& e) {
        std::cerr << "Failure: " << e.what() << " (Error code: " << e.err() << ")" << std::endl;
    } catch (std::exception& e) {
        std::cerr << "Failure: " << e.what() << std::endl;
    }

    return EXIT_SUCCESS;
}
