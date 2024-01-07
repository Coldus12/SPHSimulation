#ifndef SPHSIMULATION_SESPH_H
#define SPHSIMULATION_SESPH_H

#include <memory>
#include <vector>
#include <iterator>
#include "vulkan/vulkan.hpp"
#include "ComputeShader.hpp"
#include "SPH.h"
#include "VulkanWrapper.h"
#include "CPUSim.hpp"

namespace Vltava {
    class SESPH : public SPH {
    public:
        // other variables

        // variables for GPU sim

        // variables for CPU sim

        // other functions
        ~SESPH() {
            std::remove("sesph_data.txt");
            std::ofstream sesphData("sesph_data.txt");
            std::copy(times.begin(), times.end(), std::ostream_iterator<float>(sesphData, "\n"));
            sesphData.close();
            cleanup();
        }

        // functions for GPU sim
        void initGpuSim(std::vector<Buffer>* uBuffers, std::vector<Buffer>* sBuffers);
        void gpuTimeStep() override;

        // functions for CPU sim
        void cpuTimeStep() override;
        void setNeighbour(bool val) {
            neighbour = val;
        }

    private:
        // other variables
        bool data_has_been_set = false;
        bool simProps_has_been_set = false;
        std::vector<float> times;

        // variables for GPU sim
        std::unique_ptr<ComputeShader> densityComp;
        std::unique_ptr<ComputeShader> particleIterComp;
        std::unique_ptr<ComputeShader> gridPlacementComp;
        std::unique_ptr<ComputeShader> cleanGridComp;
        vk::Fence compFence;

        // variables for CPU sim

        // other functions
        void cleanup();

        // functions for GPU sim
        void initComputeShaders(std::vector<Buffer>* uBuffers, std::vector<Buffer>* sBuffers);
        void createSyncObjects();

        // functions for CPU sim
        void originalCalculateRhoAndP(float dt);
        void neighbourCalculateRhoAndP(float dt);
        void originalIter(float dt);
        void neighbourIter(float dt);
    };
}

#endif //SPHSIMULATION_SESPH_H
