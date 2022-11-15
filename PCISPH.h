#ifndef SPHSIMULATION_PCISPH_H
#define SPHSIMULATION_PCISPH_H

#include "SPH.h"
#include "ComputeShader.hpp"

namespace Vltava {

    class PCISPH : public SPH {
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

    private:
        // other variables
        float kPCI = 0.0;

        // variables for GPU sim

        // variables for CPU sim
        std::vector<float> rho_pred;
        std::vector<glm::vec3> v_adv;
        std::vector<glm::vec3> pred_pos;
        std::vector<glm::vec3> prev_p_acc;
        //float error = 10.0;

        // other functions
        void compute_kPCI();

        // functions for GPU sim

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
