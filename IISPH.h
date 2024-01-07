#ifndef SPHSIMULATION_IISPH_H
#define SPHSIMULATION_IISPH_H

#include <iterator>
#include "SPH.h"
#include "ComputeShader.hpp"
#include "VulkanWrapper.h"

namespace Vltava {
    struct AdditionalIISPHData {
        glm::vec3 dii = glm::vec3(0);
        float aii = 0;
        glm::vec3 sumDijPj = glm::vec3(0);
        float rho_adv = 0;
        glm::vec3 v_adv = glm::vec3(0);
        float rho_pred = 0;
    };

    class IISPH : public SPH {
    public:
        // other variables

        // variables for GPU sim

        // variables for CPU sim

        // other functions
        ~IISPH() {
            std::remove("iisph_data.txt");
            std::ofstream iisphData("iisph_data.txt");
            std::copy(times.begin(), times.end(), std::ostream_iterator<float>(iisphData, "\n"));
            iisphData.close();
            cleanup();
        }
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
        float time = 0.0F;
        std::vector<float> times;

        // variables for GPU sim
        std::unique_ptr<ComputeShader> gridPlacementComp;
        std::unique_ptr<ComputeShader> cleanGridComp;
        // predict advection
        std::unique_ptr<ComputeShader> densityComp;
        std::unique_ptr<ComputeShader> advVelocityAndDiiComp;
        std::unique_ptr<ComputeShader> advDensityAndAiiComp;

        // pressure solve
        std::unique_ptr<ComputeShader> sumDijPjComp;
        std::unique_ptr<ComputeShader> pressureUpdateComp;
        //std::unique_ptr<ComputeShader> pressureClampComp;

        // integrate
        std::unique_ptr<ComputeShader> particleIterComp;

        vk::Fence compFence;

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
        void gpuPredictAdvection();
        void gpuPressureSolveUpdate();
        void gpuIntegrate();

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
