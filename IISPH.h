#ifndef SPHSIMULATION_IISPH_H
#define SPHSIMULATION_IISPH_H

#include "SPH.h"
#include "ComputeShader.hpp"

namespace Vltava {

    class IISPH : public SPH {
    public:
        // other variables

        // variables for GPU sim

        // variables for CPU sim

        // other functions
        void setBuffers(std::vector<Buffer>* uBuffers, std::vector<Buffer>* sBuffers) override;

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

        // variables for GPU sim
        /*std::unique_ptr<ComputeShader> densityComp;
        std::unique_ptr<ComputeShader> particleIterComp;
        std::unique_ptr<ComputeShader> gridPlacementComp;
        std::unique_ptr<ComputeShader> cleanGridComp;
        vk::Fence compFence;*/

        // variables for CPU sim
        std::vector<glm::vec3> dii; // displacement
        std::vector<float> aii;
        std::vector<glm::vec3> sumDijPj;
        std::vector<float> rho_adv;
        std::vector<float> rho_pred;
        std::vector<glm::vec3> v_adv;

        // other functions
        void cleanup();

        // functions for GPU sim
        void initComputeShaders(std::vector<Buffer>* uBuffers, std::vector<Buffer>* sBuffers);
        void createSyncObjects();

        // functions for CPU sim
        void calculateRho();
        float calculateAverageError();
        void computeVadvAndDii(float dt);
        void computeRhoadvAndAii(float dt);
        void predictAdvection(float dt);
        void computeSumDijPj(float dt);
        void updatePressure(float dt);
        void pressureSolve(float dt);
        void integrate(float dt);
    };

} // Vltava

#endif //SPHSIMULATION_IISPH_H
