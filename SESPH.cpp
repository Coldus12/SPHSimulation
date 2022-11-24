#include "SESPH.h"

namespace Vltava {
    // Other functions
    //------------------------------------------------------------------------------------------------------------------
    void SESPH::cleanup() {
        VulkanResources::getInstance().logDev->getHandle().destroyFence(compFence);
    }

    // GPU functions
    //------------------------------------------------------------------------------------------------------------------
    void SESPH::initGpuSim(std::vector<Buffer>* uBuffers, std::vector<Buffer>* sBuffers) {
        densityComp.reset();
        particleIterComp.reset();
        gridPlacementComp.reset();
        cleanGridComp.reset();

        createSyncObjects();
        initComputeShaders(uBuffers, sBuffers);
    }

    void SESPH::createSyncObjects() {
        vk::FenceCreateInfo fenceInfo(vk::FenceCreateFlagBits::eSignaled);
        compFence =  VulkanResources::getInstance().logDev->getHandle().createFence(fenceInfo);
    }

    void SESPH::initComputeShaders(std::vector<Buffer>* uBuffers, std::vector<Buffer>* sBuffers) {
        setBuffers(uBuffers, sBuffers);

        densityComp = std::make_unique<ComputeShader>("shaders/SESPH_calculateRhoAndP.spv");
        densityComp->setBuffers(uBuffers, sBuffers);
        densityComp->createPipeline();

        particleIterComp = std::make_unique<ComputeShader>("shaders/SESPH_iterateParticles.spv");
        particleIterComp->setBuffers(uBuffers, sBuffers);
        particleIterComp->createPipeline();

        gridPlacementComp = std::make_unique<ComputeShader>("shaders/gridPlacement.spv");
        gridPlacementComp->setBuffers(uBuffers, sBuffers);
        gridPlacementComp->createPipeline();

        cleanGridComp = std::make_unique<ComputeShader>("shaders/resetgrid.spv");
        cleanGridComp->setBuffers(uBuffers, sBuffers);
        cleanGridComp->createPipeline();
    }

    void SESPH::gpuTimeStep() {
        logCpuData = false;

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

        // TODO something with iter
        for (int i = 0; i < 1; i++) {
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
        }

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

    // CPU functions
    //------------------------------------------------------------------------------------------------------------------
    void SESPH::cpuTimeStep() {
        logCpuData = true;

        // TODO something with iter
        for (int i = 0; i < 1; i++) {
            place();
            if (neighbour) {
                neighbourCalculateRhoAndP(props.dt);
                neighbourIter(props.dt);
            } else {
                originalCalculateRhoAndP(props.dt);
                originalIter(props.dt);
            }
        }
    }

    void SESPH::neighbourCalculateRhoAndP(float dt) {
        auto& particles = first ? particles1 : particles2;

        // Density calculation
        for (int i = 0; i < particles.size(); i++) {
            float density = 0.0f;
            float ogd = 0.0f;

            // Neighbour
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
                        density += particles[i].m * kernel(i, iterIdx);
                    }
                }
            }

            particles[i].rho = density;
            particles[i].padding = ogd;

            // Pressure calculation
            float ratio = density / props.desired_density;
            float p = props.k * (pow(ratio, 7) - 1.0);
            particles[i].p = p < 0 ? 0 : p;
        }
    }

    void SESPH::originalCalculateRhoAndP(float dt) {
        auto& particles = first ? particles1 : particles2;

        // Density calculation
        for (int i = 0; i < particles.size(); i++) {
            float density = 0.0f;
            float ogd = 0.0f;

            // Original
            for (int j = 0; j < particles.size(); j++) {
                if (i == j)
                    continue;

                density += particles[j].m * kernel(i, j);
            }

            particles[i].rho = density;
            particles[i].padding = ogd;

            // Pressure calculation
            float ratio = density / props.desired_density;
            float p = props.k * (pow(ratio, 7) - 1.0);
            particles[i].p = p < 0 ? 0 : p;
        }
    }

    void SESPH::neighbourIter(float dt) {
        const auto& p1 = first ? particles1 : particles2;
        auto& p2 = first ? particles2 : particles1;
        for (int i = 0; i < p1.size(); i++) {

            if (p1.at(i).staticP == 0) {

                // Calculate pressure
                glm::vec3 pressure = glm::vec3(0);
                glm::vec3 viscosity = glm::vec3(0);

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

                        // Note to self: as the particles get further and further from each other the density decreases which means rho --> 0
                        // whihc leads to something/0^2, which is either inf or -inf ----> nan or -nan
                        float val = 0;
                        if (p1[iterIdx].rho != 0 && p1[i].rho != 0)
                            val = p1[iterIdx].m * ((p1[i].p / pow(p1[i].rho, 2)) + (p1[iterIdx].p / pow(p1[iterIdx].rho, 2)));

                        pressure += val * gradKernel(i, iterIdx);

                        glm::vec3 xij = p1[i].x - p1[iterIdx].x;

                        float pval = 0;
                        if (p1[iterIdx].rho != 0)
                            pval = (p1[iterIdx].m / p1[iterIdx].rho) * (dot(xij, gradKernel(i, iterIdx)) /  (dot(xij, xij) + 0.01 * props.kernelh));

                        glm::vec3 vij = p1[i].v - p1[iterIdx].v;
                        viscosity += pval * vij;
                    }
                }

                pressure *= -p1[i].m;

                float nu = 0.01;
                viscosity *= 2 * nu * p1[i].m;

                glm::vec3 gravity(0, 0, -9.81 * p1[i].m);
                glm::vec3 acc = (pressure + viscosity + gravity) / p1[i].m;

                glm::vec3 viNext = p1[i].v;
                glm::vec3 xiNext = p1[i].x;
                viNext += acc * (float) dt;

                viNext = speedBound(viNext);

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

    void SESPH::originalIter(float dt) {
        const auto& p1 = first ? particles1 : particles2;
        auto& p2 = first ? particles2 : particles1;
        for (int i = 0; i < p1.size(); i++) {

            if (p1.at(i).staticP == 0) {

                // Calculate pressure
                glm::vec3 pressure = glm::vec3(0);
                glm::vec3 viscosity = glm::vec3(0);

                // Original
                for (int j = 0; j < p1.size(); j++) {
                    if (j == i)
                        continue;

                    // Note to self: as the particles get further and further from each other the density decreases which means rho --> 0
                    // whihc leads to something/0^2, which is either inf or -inf ----> nan or -nan
                    float val = 0;
                    if (p1[j].rho != 0 && p1[i].rho != 0)
                        val = p1[j].m * ((p1[i].p / pow(p1[i].rho, 2)) + (p1[j].p / pow(p1[j].rho, 2)));

                    pressure += val * gradKernel(i, j);

                    glm::vec3 xij = p1[i].x - p1[j].x;

                    float pval = 0;
                    if (p1[j].rho != 0)
                        pval = (p1[j].m / p1[j].rho) *
                               (dot(xij, gradKernel(i, j)) / (dot(xij, xij) + 0.01 * props.kernelh));

                    glm::vec3 vij = p1[i].v - p1[j].v;
                    viscosity += pval * vij;
                }

                pressure *= -p1[i].m;

                float nu = 0.01;
                viscosity *= 2 * nu * p1[i].m;

                glm::vec3 gravity(0, 0, -9.81 * p1[i].m);
                glm::vec3 acc = (pressure + viscosity + gravity) / p1[i].m;

                glm::vec3 viNext = p1[i].v;
                glm::vec3 xiNext = p1[i].x;
                viNext += acc * (float) dt;
                viNext = speedBound(viNext);
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
}