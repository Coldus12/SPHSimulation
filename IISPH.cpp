#include <chrono>
#include <thread>
#include "IISPH.h"

#define show_error false

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

            auto data = sBuffers->at(4).getData<AdditionalIISPHData>();
            for (auto& add: data) {
                error += add.rho_pred;
            }
            error /= data.size();
            error -= props.desired_density;
            error = std::abs(error);

            nr++;
#if show_error
            std::cout << "[IISPH GPU] error: " << error << " nr: " << nr << std::endl;
#endif
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

        // Grid stuff
        //--------------------
        cleanGridComp->bindPipelineAndDescriptors(computeCmdBuffer);
        computeCmdBuffer.dispatch(list_size * cellx * celly * cellz / 64 + 1, 1, 1);
        computeCmdBuffer.pipelineBarrier(
                vk::PipelineStageFlagBits::eComputeShader,
                vk::PipelineStageFlagBits::eComputeShader,
                {},
                {},
                membarriers,
                {}
        );

        gridPlacementComp->bindPipelineAndDescriptors(computeCmdBuffer);
        computeCmdBuffer.dispatch(particles1.size() / 64 + 1, 1, 1);
        computeCmdBuffer.pipelineBarrier(
                vk::PipelineStageFlagBits::eComputeShader,
                vk::PipelineStageFlagBits::eComputeShader,
                {},
                {},
                membarriers,
                {}
        );

        // Predict advection
        //--------------------
        densityComp->bindPipelineAndDescriptors(computeCmdBuffer);
        computeCmdBuffer.dispatch(particles1.size() / 64 + 1, 1, 1);
        computeCmdBuffer.pipelineBarrier(
                vk::PipelineStageFlagBits::eComputeShader,
                vk::PipelineStageFlagBits::eComputeShader,
                {},
                {},
                membarriers,
                {}
        );

        advVelocityAndDiiComp->bindPipelineAndDescriptors(computeCmdBuffer);
        computeCmdBuffer.dispatch(particles1.size() / 64 + 1, 1, 1);
        computeCmdBuffer.pipelineBarrier(
                vk::PipelineStageFlagBits::eComputeShader,
                vk::PipelineStageFlagBits::eComputeShader,
                {},
                {},
                membarriers,
                {}
        );

        advDensityAndAiiComp->bindPipelineAndDescriptors(computeCmdBuffer);
        computeCmdBuffer.dispatch(particles1.size() / 64 + 1, 1, 1);
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
        computeCmdBuffer.dispatch(particles1.size() / 64 + 1, 1, 1);
        computeCmdBuffer.pipelineBarrier(
                vk::PipelineStageFlagBits::eComputeShader,
                vk::PipelineStageFlagBits::eComputeShader,
                {},
                {},
                membarriers,
                {}
        );

        pressureUpdateComp->bindPipelineAndDescriptors(computeCmdBuffer);
        computeCmdBuffer.dispatch(particles1.size() / 64 + 1, 1, 1);
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

        // Clamp pressure
        //--------------------
        /*pressureClampComp->bindPipelineAndDescriptors(computeCmdBuffer);
        computeCmdBuffer.dispatch(particles1.size(), 1, 1);
        computeCmdBuffer.pipelineBarrier(
                vk::PipelineStageFlagBits::eComputeShader,
                vk::PipelineStageFlagBits::eComputeShader,
                {},
                {},
                membarriers,
                {}
        );*/

        // Integrate
        //--------------------
        particleIterComp->bindPipelineAndDescriptors(computeCmdBuffer);
        computeCmdBuffer.dispatch(particles1.size() / 64 + 1, 1, 1);
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

        /*pressureClampComp = std::make_unique<ComputeShader>("shaders/IISPH_clampPressure.spv");
        pressureClampComp->setBuffers(uBuffers, sBuffers);
        pressureClampComp->createPipeline();*/

        particleIterComp = std::make_unique<ComputeShader>("shaders/IISPH_iterateParticles.spv");
        particleIterComp->setBuffers(uBuffers, sBuffers);
        particleIterComp->createPipeline();

        gridPlacementComp = std::make_unique<ComputeShader>("shaders/gridPlacement.spv");
        gridPlacementComp->setBuffers(uBuffers, sBuffers);
        gridPlacementComp->createPipeline();

        cleanGridComp = std::make_unique<ComputeShader>("shaders/resetgrid.spv");
        cleanGridComp->setBuffers(uBuffers, sBuffers);
        cleanGridComp->createPipeline();
    }

    // CPU functions
    //------------------------------------------------------------------------------------------------------------------

    // IISPH time step
    //----------------------------------------------------
    void IISPH::cpuTimeStep() {
        logCpuData = true;

        // TODO iter
        for (int iter = 0; iter < 1; iter++) {
            place();

            predictAdvection(props.dt);
            pressureSolve(props.dt);
            integrate(props.dt);
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
            float og_rho = 0;

            if (!props.neighbour) {
                for (int j = 0; j < particles.size(); j++) {
                    if (i == j)
                        continue;

                    rho += particles[j].m * kernel(i, j);
                }
            } else {
                glm::vec3 tuple = determineGridTuple(i);
                Neighbourhood n = getNeighbouringCells(tuple);
                for (int nr = 0; nr < 27; nr++) {
                    glm::vec3 current = n.neighbour[nr];
                    if (!checkBounds(current)) continue;

                    int idx = getStartIdxOfCell(current);
                    if (idx >= 0) {
                        int size = grid_data.at(idx);

                        int iterIdx = 0;
                        for (int j = 1; j < size + 1; j++) {
                            iterIdx = grid_data.at(idx + j);
                            if (i == iterIdx) continue;

                            //density += in_data.p[gID].m * kernel(in_data.p[gID].x, in_data.p[iterIdx].x);
                            rho += particles[i].m * kernel(i, iterIdx);
                        }
                    }
                }
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

            if (!props.neighbour) {
                for (int j = 0; j < particles.size(); j++) {
                    if (j == i)
                        continue;

                    float pval = 0;
                    if (particles[j].rho != 0)
                        pval = (particles[j].m / particles[j].rho) * 1.0 / props.dt * 0.05 *
                               kernel(particles.at(i).x, particles.at(j).x);

                    glm::vec3 vij = particles[i].v - particles[j].v;
                    viscosity -= pval * vij;

                    dii -= particles[j].m * gradKernel(i, j);
                }
            } else {
                // Neighbour
                glm::vec3 tuple = determineGridTuple(i);
                Neighbourhood n = getNeighbouringCells(tuple);
                for (auto current : n.neighbour) {
                    if (!checkBounds(current)) continue;

                    int idx = getStartIdxOfCell(current);
                    int size = grid_data.at(idx);

                    int iterIdx = 0;
                    for (int j = 1; j < size+1; j++) {
                        iterIdx = grid_data.at(idx + j);
                        if (i == iterIdx) continue;

                        float pval = 0;
                        if (particles.at(iterIdx).rho != 0)
                            pval = (particles.at(iterIdx).m / particles.at(iterIdx).rho) * 1.0/props.dt * 0.05 * kernel(particles.at(i).x, particles.at(iterIdx).x);

                        glm::vec3 vij = particles[i].v - particles.at(iterIdx).v;
                        viscosity -= pval * vij;

                        dii -= particles[j].m * gradKernel(i, iterIdx);
                    }
                }
            }

            if (particles[i].staticP == 0)
                f += viscosity;

            f += glm::vec3(0, 0, -9.81 * particles[i].m);

            // Calculating v_adv
            v_adv[i] = particles[i].v + dt * f / particles[i].m;

            // Calculating dii
            if (abs(particles[i].rho) > 0.001)
                dii *= dt * dt / (particles[i].rho * particles[i].rho);
            else dii = glm::vec3(0);

            this->dii[i] = dii;
        }
    }

    void IISPH::computeRhoadvAndAii(float dt) {
        auto& particles = first ? particles1 : particles2;

        for (int i = 0; i < particles.size(); i++) {
            float rho_adv = 0;
            float aii = 0;

            if (!props.neighbour) {
                for (int j = 0; j < particles.size(); j++) {
                    if (i == j)
                        continue;

                    // dji = ??????
                    // dji = mi / rho_i^2 * grad(j,i)?
                    glm::vec3 dji(0);
                    if (abs(particles[i].rho) > 0.001)
                        dji = -dt * dt * particles[i].m / (particles[i].rho * particles[i].rho) * gradKernel(j, i);

                    aii += particles[j].m * glm::dot(dii[i] - dji, gradKernel(i, j));

                    glm::vec3 vij_adv = v_adv[i] - v_adv[j];
                    rho_adv += particles[j].m * glm::dot(vij_adv, gradKernel(i, j));
                }
            } else {
                glm::vec3 tuple = determineGridTuple(i);
                Neighbourhood n = getNeighbouringCells(tuple);
                for (auto current : n.neighbour) {
                    if (!checkBounds(current)) continue;

                    int idx = getStartIdxOfCell(current);
                    int size = grid_data.at(idx);

                    int iterIdx = 0;
                    for (int j = 1; j < size+1; j++) {
                        iterIdx = grid_data.at(idx + j);
                        if (i == iterIdx) continue;

                        glm::vec3 dji(0);
                        if (abs(particles[i].rho) > 0.001)
                            dji = -dt * dt * particles[i].m / (particles[i].rho * particles[i].rho) * gradKernel(iterIdx, i);

                        aii += particles[iterIdx].m * glm::dot(dii[i] - dji, gradKernel(i, iterIdx));

                        glm::vec3 vij_adv = v_adv[i] - v_adv[iterIdx];
                        rho_adv += particles[iterIdx].m * glm::dot(vij_adv, gradKernel(i, iterIdx));
                    }
                }
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
#if show_error
            std::cout << "[IISPH CPU] error: " << error << " nr " << nr << std::endl;
#endif
        }
    }

    void IISPH::computeSumDijPj(float dt) {
        auto& particles = first ? particles1 : particles2;

        for (int i = 0; i < particles.size(); i++) {
            glm::vec3 dijpj(0);

            if (!props.neighbour) {
                for (int j = 0; j < particles.size(); j++) {
                    if (i == j)
                        continue;

                    if (abs(particles[j].rho) > 0.0001)
                        dijpj += particles[j].m * particles[j].p / (particles[j].rho * particles[j].rho) *
                                 gradKernel(i, j);
                }
            } else {
                glm::vec3 tuple = determineGridTuple(i);
                Neighbourhood n = getNeighbouringCells(tuple);
                for (auto current : n.neighbour) {
                    if (!checkBounds(current)) continue;

                    int idx = getStartIdxOfCell(current);
                    int size = grid_data.at(idx);

                    int iterIdx = 0;
                    for (int j = 1; j < size+1; j++) {
                        iterIdx = grid_data.at(idx + j);
                        if (i == iterIdx) continue;

                        if (abs(particles[iterIdx].rho) > 0.0001)
                            dijpj += particles[iterIdx].m * particles[iterIdx].p / (particles[iterIdx].rho * particles[iterIdx].rho) *
                                     gradKernel(i, iterIdx);
                    }
                }
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

            if (!props.neighbour) {
                for (int j = 0; j < particles.size(); j++) {
                    if (i == j)
                        continue;

                    glm::vec3 dji(0);
                    if (abs(particles[i].rho) > 0.001)
                        dji = -dt * dt * particles[i].m / (particles[i].rho * particles[i].rho) * gradKernel(j, i);

                    sum += particles[j].m *
                           glm::dot(sumDijPj[i] - dii[j] * particles[j].p - (sumDijPj[j] - dji * particles[i].p),
                                    gradKernel(i, j));
                }
            } else {
                glm::vec3 tuple = determineGridTuple(i);
                Neighbourhood n = getNeighbouringCells(tuple);
                for (auto current : n.neighbour) {
                    if (!checkBounds(current)) continue;

                    int idx = getStartIdxOfCell(current);
                    int size = grid_data.at(idx);

                    int iterIdx = 0;
                    for (int j = 1; j < size+1; j++) {
                        iterIdx = grid_data.at(idx + j);
                        if (i == iterIdx) continue;

                        glm::vec3 dji(0);
                        if (abs(particles[i].rho) > 0.001)
                            dji = -dt * dt * particles[i].m / (particles[i].rho * particles[i].rho) * gradKernel(iterIdx, i);

                        sum += particles[iterIdx].m *
                               glm::dot(sumDijPj[i] - dii[iterIdx] * particles[iterIdx].p - (sumDijPj[iterIdx] - dji * particles[i].p),
                                        gradKernel(i, iterIdx));
                    }
                }
            }

            rho_pred[i] = rho_adv[i] + particles[i].p * aii[i] + sum;
            rho_pred[i] = rho_pred[i] < props.desired_density ? props.desired_density : rho_pred[i];
            particles[i].p = (1 - omega) * particles[i].p + (omega/aii[i]) * (props.desired_density - rho_adv[i] - sum);
            particles[i].p = particles[i].p < 0 ? 0 : particles[i].p;
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

                float dpi = 0;
                if (p1.at(i).rho != 0)
                    dpi = p1.at(i).p / (p1.at(i).rho * p1.at(i).rho);

                if (!props.neighbour) {
                    // Original
                    for (int j = 0; j < p1.size(); j++) {
                        if (j == i)
                            continue;

                        float dpj = 0;
                        if (p1.at(j).rho != 0)
                            dpj = p1.at(j).p / (p1.at(j).rho * p1.at(j).rho);

                        pressure -= p1.at(j).m * (dpi + dpj) * gradKernel(p1.at(i).x, p1.at(j).x);
                    }
                } else {
                    glm::vec3 tuple = determineGridTuple(i);
                    Neighbourhood n = getNeighbouringCells(tuple);
                    for (auto current : n.neighbour) {
                        if (!checkBounds(current)) continue;

                        int idx = getStartIdxOfCell(current);
                        int size = grid_data.at(idx);

                        int iterIdx = 0;
                        for (int j = 1; j < size+1; j++) {
                            iterIdx = grid_data.at(idx + j);
                            if (i == iterIdx) continue;

                            float dpj = 0;
                            if (p1.at(iterIdx).rho != 0)
                                dpj = p1.at(iterIdx).p / (p1.at(iterIdx).rho * p1.at(iterIdx).rho);

                            pressure -= p1.at(iterIdx).m * (dpi + dpj) * gradKernel(p1.at(i).x, p1.at(iterIdx).x);
                        }
                    }
                }

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