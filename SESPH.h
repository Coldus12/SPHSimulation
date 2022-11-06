#ifndef SPHSIMULATION_SESPH_H
#define SPHSIMULATION_SESPH_H

#include <memory>
#include <vector>
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
        std::string log() override;

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
        void originalCalculateRhoAndP();
        void neighbourCalculateRhoAndP();
        void originalIter();
        void neighbourIter();
    };
}

#endif //SPHSIMULATION_SESPH_H
