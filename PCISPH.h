#ifndef SPHSIMULATION_PCISPH_H
#define SPHSIMULATION_PCISPH_H

#include "SPH.h"
#include "ComputeShader.hpp"
#include "VulkanWrapper.h"

namespace Vltava {
    struct AdditionalPCISPHData {
        glm::vec3 v_adv = glm::vec3(0);
        float rho_pred = 0;
        glm::vec3 pred_pos = glm::vec3(0);
        float padding1 = 0;
        glm::vec3 prev_p_acc = glm::vec3(0);
        float padding2 = 0;
    };

    class PCISPH : public SPH {
    public:
        // other variables

        // variables for GPU sim

        // variables for CPU sim

        // other functions
        void setBuffers(std::vector<Buffer>* uBuffers, std::vector<Buffer>* sBuffers) override;
        ~PCISPH() {
            cleanup();
        }

        // functions for GPU sim
        void initGpuSim(std::vector<Buffer>* uBuffers, std::vector<Buffer>* sBuffers);
        void gpuTimeStep() override;

        // functions for CPU sim
        void cpuTimeStep() override;

    private:
        // other variables
        float kPCI = 0.0;

        // variables for GPU sim
        std::unique_ptr<ComputeShader> cleanGridComp;
        std::unique_ptr<ComputeShader> gridPlacementComp;
        std::unique_ptr<ComputeShader> advVelocityComp;
        std::unique_ptr<ComputeShader> predictedRhoComp;
        std::unique_ptr<ComputeShader> pressureAccelerationComp;
        std::unique_ptr<ComputeShader> densityAndPressureUpdateComp;
        std::unique_ptr<ComputeShader> particleIterComp;
        vk::Fence compFence;

        // variables for CPU sim
        std::vector<float> rho_pred;
        std::vector<glm::vec3> v_adv;
        std::vector<glm::vec3> pred_pos;
        std::vector<glm::vec3> prev_p_acc;
        //float error = 10.0;

        // other functions
        void compute_kPCI();
        void cleanup();

        // functions for GPU sim
        void initComputeShaders(std::vector<Buffer>* uBuffers, std::vector<Buffer>* sBuffers);
        void createSyncObjects();
        void gpuPredictAdvection();
        void gpuPressureSolve();
        void gpuIntegrate();

        // functions for CPU sim
        void pressureSolve(float dt);
        void computeAdvectedVelocity(float dt);
        void computePredictedRho(float dt);
        void computePressureAcceleration(float dt);
        void updateDensityAndPressure(float dt);
        float calculateError(float dt);
        void particleAdvection(float dt);
    };

} // Vltava

#endif //SPHSIMULATION_PCISPH_H
