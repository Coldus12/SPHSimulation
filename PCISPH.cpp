#include "PCISPH.h"

namespace Vltava {
    // GPU functions
    //------------------------------------------------------------------------------------------------------------------
    void PCISPH::gpuTimeStep() {
        logCpuData = false;
    }

    // CPU functions
    //------------------------------------------------------------------------------------------------------------------

    // PCISPH time step
    //----------------------------------------------------
    void PCISPH::cpuTimeStep() {
        logCpuData = true;
        float dt = 0.01;

        pressureSolve(dt);
        particleAdvection(dt);
    }

    // PCISPH iterative pressure solver
    //----------------------------------------------------
    void PCISPH::pressureSolve(float dt) {
        computeAdvectedVelocity(dt);
        computePredictedRho(dt);

        float error = 10;
        int nr = 0;

        while (error > 0.1 && nr < 50) {
            computePressureAcceleration(dt);
            updateDensityAndPressure(dt);

            nr++;
            error = abs(calculateError(dt));
            std::cout << "nr :" << nr << " error: " << error << std::endl;
        }
    }

    void PCISPH::computeAdvectedVelocity(float dt) {
        auto& particles = first ? particles1 : particles2;

        // Going through all particles
        for (int i = 0; i < particles.size(); i++) {
            // Forces
            glm::vec3 f(0);

            // dii
            //glm::vec3 dii(0);

            // Calculating viscosity
            glm::vec3 viscosity = glm::vec3(0);

            for (int j = 0; j < particles.size(); j++) {
                if (j == i)
                    continue;

                glm::vec3 xij = particles[i].x - particles[j].x;

                float pval = 0;
                if (particles[j].rho != 0)
                    pval = (particles[j].m / particles[j].rho) * (dot(xij, gradKernel(i, j)) / (dot(xij, xij) + 0.01 * props.kernelh));

                glm::vec3 vij = particles[i].v - particles[j].v;
                viscosity += pval * vij;
            }

            float nu = 0.01;
            viscosity *= 2 * nu * particles[i].m;

            if (particles[i].staticP == 0)
                f += viscosity;

            f += glm::vec3(0, 0, -9.81 * particles[i].m);
            //f += prev_p_acc[i] * particles[i].m;

            // Calculating v_adv
            v_adv[i] = particles[i].v + dt * f / particles[i].m;
            pred_pos[i] = particles[i].x + dt * (v_adv[i] + prev_p_acc[i] * dt);
        }
    }

    void PCISPH::computePredictedRho(float dt) {
        auto& particles = first ? particles1 : particles2;

        for (int i = 0; i < particles.size(); i++) {
            float density = 0;
            float real_density = 0;

            for (int j = 0; j < particles.size(); j++) {
                if (j == i) continue;

                density += particles[j].m * kernel(pred_pos[i], pred_pos[j]);
                real_density += particles[j].m * kernel(i, j);
            }

            rho_pred[i] = density;
            particles[i].rho = real_density;

            // Pressure init
            particles[i].p = kPCI * ((density - props.desired_density)/props.desired_density);
        }
    }

    void PCISPH::computePressureAcceleration(float dt) {
        auto& particles = first ? particles1 : particles2;

        for (int i = 0; i < particles.size(); i++) {
            glm::vec3 a_p(0);
            float dpi = particles[i].p / (particles[i].rho * particles[i].rho);

            for (int j = 0; j < particles.size(); j++) {
                if (j == i) continue;

                float dpj = particles[j].p / (particles[j].rho * particles[j].rho);
                a_p += -particles[j].m * (dpi + dpj) * gradKernel(i, j);
            }

            prev_p_acc[i] = a_p;
        }
    }

    void PCISPH::updateDensityAndPressure(float dt) {
        auto& particles = first ? particles1 : particles2;

        std::cout <<"\n\n";
        for (int i = 0; i < particles.size(); i++) {
            float rho_change = 0;

            for (int j = 0; j < particles.size(); j++) {
                if (j == i) continue;

                rho_change += particles[j].m * dot((dt * prev_p_acc[i] - dt * prev_p_acc[j]), gradKernel(i, j));
            }

            rho_change *= dt;
            rho_pred[i] += rho_change;
            particles[i].rho = rho_pred[i];
            particles[i].p += kPCI * ((rho_pred[i] - props.desired_density)/props.desired_density);
        }
    }

    float PCISPH::calculateError(float dt) {
        auto& particles = first ? particles1 : particles2;
        float error = 0.0f;

        for (int i = 0; i < particles.size(); i++) {
            error += particles[i].rho;
        }

        error /= particles.size();
        return (error - props.desired_density)/props.desired_density;
    }

    // PCISPH particle advection
    //----------------------------------------------------
    void PCISPH::particleAdvection(float dt) {
        const auto& p1 = first ? particles1 : particles2;
        auto& p2 = first ? particles2 : particles1;

        for (int i = 0; i < p1.size(); i++) {

            if (p1.at(i).staticP == 0) {
                glm::vec3 viNext = v_adv[i];
                glm::vec3 xiNext = p1[i].x;
                viNext += prev_p_acc[i] * dt;
                xiNext += viNext * dt;

                p2[i].x = xiNext;
                p2[i].v = viNext;
                p2[i].rho = p1[i].rho;
                p2[i].p = p1[i].p;

                container(i);
            } else {
                p2.at(i) = p1.at(i);
            }
        }

        first = !first;
    }

    // Other functions
    //------------------------------------------------------------------------------------------------------------------
    void PCISPH::compute_kPCI() {
        auto& particles = first ? particles1 : particles2;
        float dt = 0.01f;

        glm::vec3 pi(0);
        glm::vec3 pj(0);

        float particleR = props.kernelh/4;
        glm::vec3 sumGradWij(0);
        float sumDotGradWij = 0;

        glm::vec3 gk(0); //gradKernel

        // There should be at maximum 5 particles between
        // one "edge" of the sampling sphere and another "edge".
        // Therefore if I go through a 5x5x5 particle block,
        // I should have sampled more than enough particles
        // for a good kPCI.
        for (int i = -2; i < 3; i++) {
            for (int j = -2; j < 3; j++) {
                for (int k = -2; k < 3; k++) {
                    pj = glm::vec3(i * particleR, j * particleR, k * particleR);

                    gk = gradKernel(pi, pj);

                    sumGradWij += gk;
                    sumDotGradWij += dot(gk, gk);
                }
            }
        }

        kPCI = (0.5f * props.desired_density * props.desired_density) / (dt * dt * particles[0].m * particles[0].m);
        kPCI /= dot(sumGradWij, sumGradWij) + sumDotGradWij;

        std::cout << kPCI << std::endl;
    }

    void PCISPH::setBuffers(std::vector<Buffer> *uBuffers, std::vector<Buffer> *sBuffers) {
        SPH::setBuffers(uBuffers, sBuffers);

        rho_pred = std::vector<float>(particles1.size(), 0.0f);
        v_adv = std::vector<glm::vec3>(particles1.size(), glm::vec3(0));
        pred_pos = std::vector<glm::vec3>(particles1.size(), glm::vec3(0));
        prev_p_acc = std::vector<glm::vec3>(particles1.size(), glm::vec3(0));

        compute_kPCI();
    }
} // Vltava