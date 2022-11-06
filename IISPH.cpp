#include <chrono>
#include <thread>
#include "IISPH.h"

namespace Vltava {
    // GPU functions
    //------------------------------------------------------------------------------------------------------------------
    void IISPH::gpuTimeStep() {
        logCpuData = false;

        gpuPredictAdvection();

        int nr = 0;
        float error = 10;

        //while (nr < 1) {
        while (error > 0.01 && nr < 100) {
            gpuPressureSolveUpdate();

            // Density error calculation
            error = 0;

            auto data = sBuffers->at(3).getData<AdditionalIISPHData>();
            //std::this_thread::sleep_for(std::chrono::milliseconds(1000));
            for (auto& add: data) {
                //std::cout << add.rho_pred << " " << add.rho_adv << " " << add.aii << std::endl;
                error += add.rho_pred;
            }
            error /= data.size();
            error -= props.desired_density;

            nr++;

            std::cout << "[IISPH GPU] error: " << error << " nr: " << nr << std::endl;
        }

        gpuIntegrate();
    }

    void IISPH::gpuPredictAdvection() {
        vk::CommandBufferBeginInfo beginInfo({}, nullptr);
        auto& computeCmdBuffer = VulkanWrapper::getInstance().getCompCmdBuffer();
        computeCmdBuffer.begin(beginInfo);

        std::vector<vk::BufferMemoryBarrier> membarriers;
        for (auto& buffer: *sBuffers) {
            membarriers.emplace_back(
                    vk::AccessFlagBits::eShaderWrite,
                    vk::AccessFlagBits::eShaderRead,
                    VulkanWrapper::getInstance().getComputeQueueFamily(),
                    VulkanWrapper::getInstance().getComputeQueueFamily(),
                    buffer.getBufferHandle(),
                    0,
                    VK_WHOLE_SIZE
            );
        }

        // Predict advection
        //--------------------
        densityComp->bindPipelineAndDescriptors(computeCmdBuffer);
        computeCmdBuffer.dispatch(particles1.size(), 1, 1);
        computeCmdBuffer.pipelineBarrier(
                vk::PipelineStageFlagBits::eComputeShader,
                vk::PipelineStageFlagBits::eComputeShader,
                {},
                {},
                membarriers,
                {}
        );

        advVelocityAndDiiComp->bindPipelineAndDescriptors(computeCmdBuffer);
        computeCmdBuffer.dispatch(particles1.size(), 1, 1);
        computeCmdBuffer.pipelineBarrier(
                vk::PipelineStageFlagBits::eComputeShader,
                vk::PipelineStageFlagBits::eComputeShader,
                {},
                {},
                membarriers,
                {}
        );

        advDensityAndAiiComp->bindPipelineAndDescriptors(computeCmdBuffer);
        computeCmdBuffer.dispatch(particles1.size(), 1, 1);
        computeCmdBuffer.pipelineBarrier(
                vk::PipelineStageFlagBits::eComputeShader,
                vk::PipelineStageFlagBits::eComputeShader,
                {},
                {},
                membarriers,
                {}
        );

        computeCmdBuffer.end();

        vk::SubmitInfo submitInfo(
                {},
                {},
                computeCmdBuffer,
                {}
        );

        //VulkanResources::getInstance().computeQueue->submit(submitInfo, compFence);
        VulkanResources::getInstance().computeQueue->submit(submitInfo);
        VulkanResources::getInstance().logDev->getHandle().waitIdle();
    }

    void IISPH::gpuPressureSolveUpdate() {
        vk::CommandBufferBeginInfo beginInfo({}, nullptr);
        auto& computeCmdBuffer = VulkanWrapper::getInstance().getCompCmdBuffer();
        computeCmdBuffer.begin(beginInfo);

        std::vector<vk::BufferMemoryBarrier> membarriers;
        for (auto& buffer: *sBuffers) {
            membarriers.emplace_back(
                    vk::AccessFlagBits::eShaderWrite,
                    vk::AccessFlagBits::eShaderRead,
                    VulkanWrapper::getInstance().getComputeQueueFamily(),
                    VulkanWrapper::getInstance().getComputeQueueFamily(),
                    buffer.getBufferHandle(),
                    0,
                    VK_WHOLE_SIZE
            );
        }

        // Pressure solve
        //--------------------
        sumDijPjComp->bindPipelineAndDescriptors(computeCmdBuffer);
        computeCmdBuffer.dispatch(particles1.size(), 1, 1);
        computeCmdBuffer.pipelineBarrier(
                vk::PipelineStageFlagBits::eComputeShader,
                vk::PipelineStageFlagBits::eComputeShader,
                {},
                {},
                membarriers,
                {}
        );

        pressureUpdateComp->bindPipelineAndDescriptors(computeCmdBuffer);
        computeCmdBuffer.dispatch(particles1.size(), 1, 1);
        computeCmdBuffer.pipelineBarrier(
                vk::PipelineStageFlagBits::eComputeShader,
                vk::PipelineStageFlagBits::eComputeShader,
                {},
                {},
                membarriers,
                {}
        );

        computeCmdBuffer.end();

        vk::SubmitInfo submitInfo(
                {},
                {},
                computeCmdBuffer,
                {}
        );

        //VulkanResources::getInstance().computeQueue->submit(submitInfo, compFence);
        VulkanResources::getInstance().computeQueue->submit(submitInfo);
        VulkanResources::getInstance().logDev->getHandle().waitIdle();
    }

    void IISPH::gpuIntegrate() {
        vk::CommandBufferBeginInfo beginInfo({}, nullptr);
        auto& computeCmdBuffer = VulkanWrapper::getInstance().getCompCmdBuffer();
        computeCmdBuffer.begin(beginInfo);

        std::vector<vk::BufferMemoryBarrier> membarriers;
        for (auto& buffer: *sBuffers) {
            membarriers.emplace_back(
                    vk::AccessFlagBits::eShaderWrite,
                    vk::AccessFlagBits::eShaderRead,
                    VulkanWrapper::getInstance().getComputeQueueFamily(),
                    VulkanWrapper::getInstance().getComputeQueueFamily(),
                    buffer.getBufferHandle(),
                    0,
                    VK_WHOLE_SIZE
            );
        }

        // Integrate
        //--------------------
        particleIterComp->bindPipelineAndDescriptors(computeCmdBuffer);
        computeCmdBuffer.dispatch(particles1.size(), 1, 1);
        computeCmdBuffer.pipelineBarrier(
                vk::PipelineStageFlagBits::eComputeShader,
                vk::PipelineStageFlagBits::eComputeShader,
                {},
                {},
                membarriers,
                {}
        );

        computeCmdBuffer.end();

        vk::SubmitInfo submitInfo(
                {},
                {},
                computeCmdBuffer,
                {}
        );

        //VulkanResources::getInstance().computeQueue->submit(submitInfo, compFence);
        VulkanResources::getInstance().computeQueue->submit(submitInfo);
        VulkanResources::getInstance().logDev->getHandle().waitIdle();
    }

    void IISPH::initGpuSim(std::vector<Buffer> *uBuffers, std::vector<Buffer> *sBuffers) {
        densityComp.reset();
        advVelocityAndDiiComp.reset();
        advDensityAndAiiComp.reset();

        sumDijPjComp.reset();
        pressureUpdateComp.reset();

        particleIterComp.reset();

        createSyncObjects();
        initComputeShaders(uBuffers, sBuffers);
    }

    void IISPH::createSyncObjects() {
        vk::FenceCreateInfo fenceInfo(vk::FenceCreateFlagBits::eSignaled);
        compFence =  VulkanResources::getInstance().logDev->getHandle().createFence(fenceInfo);
    }

    void IISPH::initComputeShaders(std::vector<Buffer> *uBuffers, std::vector<Buffer> *sBuffers) {
        setBuffers(uBuffers, sBuffers);

        Buffer additionalDataBuffer(
                particles1.size() * sizeof(AdditionalIISPHData),
                vk::BufferUsageFlagBits::eStorageBuffer,
                vk::MemoryPropertyFlagBits::eHostVisible | vk::MemoryPropertyFlagBits::eHostCoherent
        );
        additionalDataBuffer.bind(0);
        std::vector<AdditionalIISPHData> add_data = std::vector<AdditionalIISPHData>(particles1.size(), AdditionalIISPHData());
        additionalDataBuffer.writeToBuffer(add_data.data(), add_data.size() * sizeof(AdditionalIISPHData));
        sBuffers->push_back(std::move(additionalDataBuffer));

        densityComp = std::make_unique<ComputeShader>("shaders/IISPH_calculateRho.spv");
        densityComp->setBuffers(uBuffers, sBuffers);
        densityComp->createPipeline();

        advVelocityAndDiiComp = std::make_unique<ComputeShader>("shaders/IISPH_computeVadvAndDii.spv");
        advVelocityAndDiiComp->setBuffers(uBuffers, sBuffers);
        advVelocityAndDiiComp->createPipeline();

        advDensityAndAiiComp = std::make_unique<ComputeShader>("shaders/IISPH_computeRhoadvAndAii.spv");
        advDensityAndAiiComp->setBuffers(uBuffers, sBuffers);
        advDensityAndAiiComp->createPipeline();

        sumDijPjComp = std::make_unique<ComputeShader>("shaders/IISPH_computeSumDijPj.spv");
        sumDijPjComp->setBuffers(uBuffers, sBuffers);
        sumDijPjComp->createPipeline();

        pressureUpdateComp = std::make_unique<ComputeShader>("shaders/IISPH_updatePressure.spv");
        pressureUpdateComp->setBuffers(uBuffers, sBuffers);
        pressureUpdateComp->createPipeline();

        particleIterComp = std::make_unique<ComputeShader>("shaders/IISPH_iterateParticles.spv");
        particleIterComp->setBuffers(uBuffers, sBuffers);
        particleIterComp->createPipeline();
    }

    // CPU functions
    //------------------------------------------------------------------------------------------------------------------

    // IISPH time step
    //----------------------------------------------------
    void IISPH::cpuTimeStep() {
        logCpuData = true;

        float dt = 0.01;

        // TODO iter
        for (int iter = 0; iter < 1; iter++) {
            predictAdvection(dt);
            pressureSolve(dt);
            integrate(dt);
        }
    }

    // Predict advection
    //----------------------------------------------------
    void IISPH::predictAdvection(float dt) {
        calculateRho();
        computeVadvAndDii(dt);
        computeRhoadvAndAii(dt);
    }

    void IISPH::calculateRho() {
        auto& particles = first ? particles1 : particles2;

        for (int i = 0; i < particles.size(); i++) {

            // Calculating rho
            float rho = 0;

            for (int j = 0; j < particles.size(); j++) {
                if (i == j)
                    continue;

                rho += particles[j].m * kernel(i, j);
            }

            particles[i].rho = rho;
        }
    }

    void IISPH::computeVadvAndDii(float dt) {
        auto& particles = first ? particles1 : particles2;

        // Going through all particles
        for (int i = 0; i < particles.size(); i++) {
            // Forces
            glm::vec3 f(0);

            // dii
            glm::vec3 dii(0);

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

                dii += -particles[j].m * gradKernel(i,j);
            }

            float nu = 0.01;
            viscosity *= 2 * nu * particles[i].m;

            if (particles[i].staticP == 0)
                f += viscosity;

            f += glm::vec3(0, 0, -9.81 * particles[i].m);

            // Calculating v_adv
            v_adv[i] = particles[i].v + dt * f / particles[i].m;

            // Calculating dii
            dii *= dt * dt / (particles[i].rho * particles[i].rho); // TODO, what happens if rho == 0?
            this->dii[i] = dii;
        }
    }

    void IISPH::computeRhoadvAndAii(float dt) {
        auto& particles = first ? particles1 : particles2;

        for (int i = 0; i < particles.size(); i++) {
            float rho_adv = 0;
            float aii = 0;

            for (int j = 0; j < particles.size(); j++) {
                if (i == j)
                    continue;

                // dji = ??????
                // dji = mi / rho_i^2 * grad(j,i)?
                //TODO rho == 0?
                glm::vec3 dji = -dt * dt * particles[i].m / (particles[i].rho * particles[i].rho) * gradKernel(j, i);
                aii += particles[j].m * glm::dot(dii[i] - dji, gradKernel(i,j));

                glm::vec3 vij_adv = v_adv[i] - v_adv[j];
                rho_adv += particles[j].m * glm::dot(vij_adv, gradKernel(i,j));
            }

            rho_adv *= dt;
            rho_adv += particles[i].rho;
            this->rho_adv[i] = rho_adv;
            this->aii[i] = aii;

            // pressure init
            particles[i].p *= 0.5;
        }
    }

    // Pressure solve
    //----------------------------------------------------
    void IISPH::pressureSolve(float dt) {
        int nr = 0;
        float error = 10;

        while (error > 0.01 && nr < 100) {
            computeSumDijPj(dt);
            updatePressure(dt);

            error = abs(calculateAverageError());
            nr++;

            std::cout << "[IISPH CPU] error: " << error << " nr " << std::endl;
        }
        auto& particles = first ? particles1 : particles2;
    }

    void IISPH::computeSumDijPj(float dt) {
        auto& particles = first ? particles1 : particles2;

        for (int i = 0; i < particles.size(); i++) {
            glm::vec3 dijpj(0);

            for (int j = 0; j < particles.size(); j++) {
                if (i == j)
                    continue;

                // TODO rho == 0?
                dijpj += particles[j].m * particles[j].p / (particles[j].rho * particles[j].rho) * gradKernel(i,j);
            }

            dijpj *= -dt*dt;

            sumDijPj[i] = dijpj;
        }
    }

    void IISPH::updatePressure(float dt) {
        auto& particles = first ? particles1 : particles2;

        float omega = 0.5f;

        for (int i = 0; i < particles.size(); i++) {
            float sum = 0;

            for (int j = 0; j < particles.size(); j++) {
                if (i == j)
                    continue;

                glm::vec3 dji = -dt*dt*particles[i].m / (particles[i].rho * particles[i].rho) * gradKernel(j, i);
                sum += particles[j].m * glm::dot(sumDijPj[i] - dii[j] * particles[j].p - (sumDijPj[j] - dji * particles[i].p),
                                                 gradKernel(i,j));
            }

            rho_pred[i] = rho_adv[i] + particles[i].p * aii[i] + sum;
            particles[i].p = (1 - omega) * particles[i].p + (omega/aii[i]) * (props.desired_density - rho_adv[i] - sum);
        }
    }

    float IISPH::calculateAverageError() {
        float rho_avg = 0;

        for (auto& rho : rho_pred) {
            rho_avg += rho;
        }

        rho_avg /= rho_pred.size();
        return rho_avg - props.desired_density;
    }

    // Integrate
    //----------------------------------------------------
    void IISPH::integrate(float dt) {
        const auto& p1 = first ? particles1 : particles2;
        auto& p2 = first ? particles2 : particles1;

        for (int i = 0; i < p1.size(); i++) {

            if (p1.at(i).staticP == 0) {

                // Calculate pressure
                glm::vec3 pressure = glm::vec3(0);

                // Original
                for (int j = 0; j < p1.size(); j++) {
                    if (j == i)
                        continue;

                    // Note to self: as the particles get further and further from each other the density decreases which means rho --> 0
                    // whihc leads to something/0^2, which is either inf or -inf ----> nan or -nan
                    float val = 0;
                    if (p1[j].rho != 0 && p1[i].rho != 0)
                        val = p1[j].m * ((p1[i].p / pow(p1[i].rho, 2)) + (p1[j].p / pow(p1[j].rho, 2)));

                    pressure += -val * gradKernel(i, j);
                }

                /*pressure *= -p1[i].m;
                glm::vec3 acc = (pressure) / p1[i].m;*/

                glm::vec3 viNext = v_adv[i];
                glm::vec3 xiNext = p1[i].x;
                viNext += pressure * dt;
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
    void IISPH::setBuffers(std::vector<Buffer> *uBuffers, std::vector<Buffer> *sBuffers) {
        SPH::setBuffers(uBuffers, sBuffers);

        dii = std::vector<glm::vec3>(particles1.size(), glm::vec3(0));
        sumDijPj = std::vector<glm::vec3>(particles1.size(), glm::vec3(0));
        v_adv = std::vector<glm::vec3>(particles1.size(), glm::vec3(0));

        aii = std::vector<float>(particles1.size(), 0.0f);
        rho_adv = std::vector<float>(particles1.size(), 0.0f);
        rho_pred = std::vector<float>(particles1.size(), 0.0f);
    }

    void IISPH::cleanup() {
        VulkanResources::getInstance().logDev->getHandle().destroyFence(compFence);
    }
} // Vltava