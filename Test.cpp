#include <thread>
#include "Test.h"
#include "CPUSim.hpp"
#include "VulkanWrapper.h"
#include "ComputeShader.hpp"

#define log false

namespace Vltava {
#if 0
    void Test::runTests() {
        initVulkanWrapper();

        basicPlaceTest();
        tupleBoundsCheck();
        basicNeighbourhoodTest();
        cpuGpuPlaceCompare();
    }

    // Grid tests
    //------------------------------------------------------------------------------------------------------------------
    // Testing placement
    void Test::basicPlaceTest() {
        std::cout << "doing basicPlaceTest" << std::endl;

        CPUSim sim1;

        Particle p1 = {
                glm::vec3(0.45f,0.78f,0.1f),
                0.1f,
                glm::vec3(0,0,0),
                0.005236f,
                0.0f,
                0.0f,
                0.0f
        };

        // Set simProps
        SimProps props{
                1.0f,
                0.001f,
                1.0f,
                0.2f,

                glm::vec4(-1,-1,-1, 0),
                glm::vec4(1, 1, 1, 0)
        };

        std::vector<Particle> particles;
        particles.push_back(p1);

        sim1.setSimProps(props);
        sim1.setData(particles);

        sim1.place();

        // number of cells on a given axis = 10; list_size = 12
        /* --> cellwidth = 0.2
         * --> cellTuple = abs(floor((a.x - x)/0.2), floor((a.y - y)/0.2), floor((a.z - z)/0.2))
         * --> cellTuple = (-1 - 0.45)/0.2 = 7 etc.
         * --> cellTuple = (7, 8, 9)
         */
        glm::vec3 cellTuple = sim1.determineGridTuple(0);
        assert((cellTuple.x == 7) && "CellTuple x isn't correct");
        assert((cellTuple.y == 8) && "CellTuple y isn't correct");
        assert((cellTuple.z == 5) && "CellTuple z isn't correct");

        // Index should be: x * cellz * celly + y * cellz + z
        // Therefore idx = (7 * 10 * 10 + 8 * 10 + 5) * 12 = 785 * 25 = 19 625
        int idx = sim1.getStartIdxOfCell(cellTuple);
        assert((idx == 19625) && "getStartIdxOfCell returned a wrong number.");

        // Since we only put in 1 particle, size = 1
        int size = sim1.grid_data[idx];
        assert((size == 1) && "size isn't correct");
    }

    void Test::tupleBoundsCheck() {
        std::cout << "doing tupleBoundsCheck" << std::endl;
        CPUSim sim1;

        Particle p1 = {
                glm::vec3(-2.45f,0.78f,0.1f),
                0.1f,
                glm::vec3(0,0,0),
                0.005236f,
                0.0f,
                0.0f,
                0.0f
        };

        Particle p2 = {
                glm::vec3(2.45f,0.78f,0.1f),
                0.1f,
                glm::vec3(0,0,0),
                0.005236f,
                0.0f,
                0.0f,
                0.0f
        };

        // Set simProps
        SimProps props{
                1.0f,
                0.001f,
                2.0f,
                0.2f,

                glm::vec4(-1,-1,-1, 0),
                glm::vec4(1, 1, 1, 0)
        };

        std::vector<Particle> particles;
        particles.push_back(p1);
        particles.push_back(p2);

        sim1.setSimProps(props);
        sim1.setData(particles);

        sim1.place();

        // number of cells on a given axis = 10; list_size = 12
        /* --> cellwidth = 0.2
         * --> cellTuple = abs(floor((a.x - x)/0.2), floor((a.y - y)/0.2), floor((a.z - z)/0.2))
         * --> cellTuple = abs( (-1 - (-2.45))/0.2 ) = 7 etc. --> problem
         * --> cellTuple = should be (7, 8, 9), but -2.45 is out of the bounds so
         * --> cellTuple = (-1, -1, -1)
         */
        glm::vec3 cellTuple = sim1.determineGridTuple(0);
        assert((cellTuple.x == -1) && "CellTuple x isn't correct");
        assert((cellTuple.y == -1) && "CellTuple y isn't correct");
        assert((cellTuple.z == -1) && "CellTuple z isn't correct");

        // There shouldn't be a calculated index.
        int idx = sim1.getStartIdxOfCell(cellTuple);
        assert((idx == -1) && "getStartIdxOfCell returned a wrong number.");

        cellTuple = sim1.determineGridTuple(1);
        assert((cellTuple.x == -1) && "2.CellTuple x isn't correct");
        assert((cellTuple.y == -1) && "2.CellTuple y isn't correct");
        assert((cellTuple.z == -1) && "2.CellTuple z isn't correct");

        idx = sim1.getStartIdxOfCell(cellTuple);
        assert((idx == -1) && "2.getStartIdxOfCell returned a wrong number.");
    }

    void Test::basicNeighbourhoodTest() {
        std::cout << "doing basicNeighbourhoodTest" << std::endl;

        CPUSim sim1;
        std::vector<Particle> particles;

        float s = 0.2;
        int r = -1;
        int z = -1;
        for (int i = 0; i < 1000; i++) {
            Particle data;

            if (i % 10 == 0)
                r++;

            if (i % (100) == 0)
                z++;

            //data.x = glm::vec3((i%10) * s,(r%10) * s,z * s);
            data.x = glm::vec3((i%10) * s,(r%10) * s,z * s);
            data.x -= glm::vec3(1);
            data.h = 0.1; // Only sets size for visualization atm
            data.v = glm::vec3(0,0,0);
            data.m = 0.005236f;
            //data[i].m = 1.0f;

            data.rho = 1;
            data.p = 1;

            data.staticP = 0;
            data.padding = 0;

            particles.push_back(data);
        }

        // Set simProps
        SimProps props{
                1.0f,
                0.001f,
                (float) particles.size(),
                0.2f,

                glm::vec4(-1,-1,-1, 0),
                glm::vec4(1, 1, 1, 0)
        };

        sim1.setSimProps(props);
        sim1.setData(particles);

        sim1.place();

        // s = 0.2; gridA = (-1, -1, -1)
        // x = (555%10 = 5) * 0.2 - 1, y = (555%100 = 5) * 0.2 - 1, z = (floor(555/100)) * 0.2 - 1
        // x = y = z = 0 -> middle
        glm::vec3 tuple = sim1.determineGridTuple(555); // middle -> should be 5,5,5
        assert((tuple.x == 5) && "[Neighbourhood] tuple x isnt correct");
        assert((tuple.y == 5) && "[Neighbourhood] tuple y isnt correct");
        assert((tuple.z == 5) && "[Neighbourhood] tuple z isnt correct");

        int neighbournr = 0;
        Neighbourhood nh = sim1.getNeighbouringCells(tuple);
        for (auto& current: nh.neighbour) {
            int idx = sim1.getStartIdxOfCell(current);
            int size = sim1.grid_data[idx];
            neighbournr += size;
        }
        assert((neighbournr == 27) && "[Neighbourhood] the number of neighbours isn't correct");
    }

    bool errorMargin(float val1, float val2, float margin) {
        if ((val1 + margin > val2) && (val1 - margin < val2))
            return true;

        return false;
    }

    // Seems to work with the old comp/comp_it shaders and CPUSim
    // Also works with the new ones. NOTE: ALWAYS START MORE GRID-CLEANING "THREADS"
    void Test::cpuGpuPlaceCompare() {
        CPUSim sim;

        std::vector<Buffer> sBuffers;
        std::vector<Buffer> uBuffers;

        int list_size = 25;

        // Data
        //-----------------------------
        Buffer bprops(
                sizeof(SimProps),
                vk::BufferUsageFlagBits::eUniformBuffer,
                vk::MemoryPropertyFlagBits::eHostCoherent | vk::MemoryPropertyFlagBits::eHostVisible
        );
        uBuffers.push_back(std::move(bprops));

        std::vector<Particle> particles;
        float s=0.10;
        int r = -1;
        int z = -1;
        float mass = 0.005236f;
        int pnrAlongAxis = round(pow(64, 1.0/3.0));

        // Actual simulated particles
        //-----------------------------
        for (int i = 0; i < 64; i++) {
            Particle data;

            if (i % pnrAlongAxis == 0)
                r++;

            if (i % (pnrAlongAxis * pnrAlongAxis) == 0)
                z++;

            data.x = glm::vec3((i%pnrAlongAxis) * s,(r%pnrAlongAxis) * s, z * s);
            data.h = 0.1; // Only sets size for visualization atm
            data.v = glm::vec3(0,0,0);
            data.m = mass;

            data.rho = 1;
            data.p = 1;

            data.staticP = 0;
            data.padding = 0;

            particles.push_back(data);
        }
        vk::DeviceSize size = sizeof(Particle) * particles.size();

        SimProps props = {
                1.0f,
                0.001f,
                (float) particles.size(),
                0.2f,

                glm::vec4(-2, -2, -2, 0),
                glm::vec4(2, 2, 2, 0)
        };
        uBuffers[0].writeToBuffer(&props, sizeof(props));
        sim.setSimProps(props);
        sim.setData(particles);

        Buffer inBuffer(
                size,
                vk::BufferUsageFlagBits::eStorageBuffer,
                vk::MemoryPropertyFlagBits::eHostVisible | vk::MemoryPropertyFlagBits::eHostCoherent
        );
        inBuffer.bind(0);

        Buffer outBuffer(
                size,
                vk::BufferUsageFlagBits::eStorageBuffer,
                vk::MemoryPropertyFlagBits::eHostVisible | vk::MemoryPropertyFlagBits::eHostCoherent
        );
        outBuffer.setSize(size);
        outBuffer.bind(0);

        int cellx = int(ceil(abs((props.gridB.x - props.gridA.x)/props.kernelh))); // Number of cells in x direction
        int celly = int(ceil(abs((props.gridB.y - props.gridA.y)/props.kernelh))); // Number of cells in y direction
        int cellz = int(ceil(abs((props.gridB.z - props.gridA.z)/props.kernelh))); // Number of cells in z direction

        // (0.2/0.1)^3 * 1.5 = 12
        Buffer gridBuffer(
        list_size * cellx * celly * cellz * sizeof(int),
                vk::BufferUsageFlagBits::eStorageBuffer,
                vk::MemoryPropertyFlagBits::eHostVisible | vk::MemoryPropertyFlagBits::eHostCoherent
        );
        gridBuffer.bind(0);

        std::vector<int> val;
        for (int i = 0; i < list_size * cellx * celly * cellz; i++)
            val.push_back(0);

        sBuffers.push_back(std::move(inBuffer));
        sBuffers.push_back(std::move(outBuffer));

        for (auto& b: sBuffers)
            b.writeToBuffer(particles.data(), size);

        sBuffers.push_back(std::move(gridBuffer));
        sBuffers[2].writeToBuffer(val.data(), val.size() * sizeof(int));

        // ComputeShaders
        //-----------------------------
        auto densityComp = ComputeShader("shaders/comp.spv");
        densityComp.setBuffers(&uBuffers, &sBuffers);
        densityComp.createPipeline();

        auto particleIterComp = ComputeShader("shaders/comp_it.spv");
        particleIterComp.setBuffers(&uBuffers, &sBuffers);
        particleIterComp.createPipeline();

        auto gridPlacementComp = ComputeShader("shaders/atomic.spv");
        gridPlacementComp.setBuffers(&uBuffers, &sBuffers);
        gridPlacementComp.createPipeline();

        auto cleanGridComp = ComputeShader("shaders/resetgrid.spv");
        cleanGridComp.setBuffers(&uBuffers, &sBuffers);
        cleanGridComp.createPipeline();

        // Sync objects
        //-----------------------------
        std::vector<vk::Semaphore> imageAvailableSemaphores;
        std::vector<vk::Semaphore> renderFinishedSemaphores;
        std::vector<vk::Fence> inFlightFences;
        vk::Fence compFence;

        imageAvailableSemaphores.reserve(VulkanResources::getInstance().FRAMES_IN_FLIGHT);
        renderFinishedSemaphores.reserve(VulkanResources::getInstance().FRAMES_IN_FLIGHT);
        inFlightFences.reserve(VulkanResources::getInstance().FRAMES_IN_FLIGHT);

        vk::FenceCreateInfo fenceInfo(vk::FenceCreateFlagBits::eSignaled);
        vk::SemaphoreCreateInfo semaphoreCreateInfo;

        compFence =  VulkanResources::getInstance().logDev->getHandle().createFence(fenceInfo);
        for (int i = 0; i < VulkanResources::getInstance().FRAMES_IN_FLIGHT; i++) {
            imageAvailableSemaphores.push_back(VulkanResources::getInstance().logDev->getHandle().createSemaphore(semaphoreCreateInfo));
            renderFinishedSemaphores.push_back(VulkanResources::getInstance().logDev->getHandle().createSemaphore(semaphoreCreateInfo));
            inFlightFences.push_back(VulkanResources::getInstance().logDev->getHandle().createFence(fenceInfo));
        }

        // Dispatch compute
        //-----------------------------
        vk::CommandBufferBeginInfo beginInfo({}, nullptr);
        auto& computeCmdBuffer = VulkanWrapper::getInstance().getCompCmdBuffer();
        computeCmdBuffer.begin(beginInfo);

        std::vector<vk::BufferMemoryBarrier> membarriers;
        for (auto& buffer: sBuffers) {
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

        for (int i = 0; i < 300; i++) { /* Default: 150 */
            sim.runSESPH(1);

            cleanGridComp.bindPipelineAndDescriptors(computeCmdBuffer);
            computeCmdBuffer.dispatch(1000, 1, 1);
            computeCmdBuffer.pipelineBarrier(
                    vk::PipelineStageFlagBits::eComputeShader,
                    vk::PipelineStageFlagBits::eComputeShader,
                    {},
                    {},
                    membarriers,
                    {}
            );

            gridPlacementComp.bindPipelineAndDescriptors(computeCmdBuffer);
            computeCmdBuffer.dispatch(1, 1, 1);
            computeCmdBuffer.pipelineBarrier(
                    vk::PipelineStageFlagBits::eComputeShader,
                    vk::PipelineStageFlagBits::eComputeShader,
                    {},
                    {},
                    membarriers,
                    {}
            );

            densityComp.bindPipelineAndDescriptors(computeCmdBuffer);
            computeCmdBuffer.dispatch(1, 1, 1);
            computeCmdBuffer.pipelineBarrier(
                    vk::PipelineStageFlagBits::eComputeShader,
                    vk::PipelineStageFlagBits::eComputeShader,
                    {},
                    {},
                    membarriers,
                    {}
            );

            particleIterComp.bindPipelineAndDescriptors(computeCmdBuffer);
            computeCmdBuffer.dispatch(1, 1, 1);
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

        // Actual submit
        //-----------------------------
        VulkanResources::getInstance().computeQueue->submit(submitInfo);

        std::this_thread::sleep_for(std::chrono::milliseconds(500));

        // Checking results
        //-----------------------------
        // grid data
        auto gpuGrid = sBuffers[2].getData<int>();
        auto cpuGrid = sim.grid_data;

        // Checking grid sizes
        assert((gpuGrid.size() == cpuGrid.size()) && "Grid sizes are different!");

        int cpuSize = 0;
        int gpuSize = 0;
        int startIdx = 0;

        for (int i = 0; i < sim.cellx; i++) {
            for (int j = 0; j < sim.celly; j++) {
                for (int k = 0; k < sim.cellz; k++) {
                    startIdx = sim.getStartIdxOfCell(glm::vec3(i,j,k));
                    cpuSize = cpuGrid[startIdx];
                    gpuSize = cpuGrid[startIdx];

                    // Checking cell sizes
                    assert((cpuSize == gpuSize) && "Number of particles in cell differs.");

                    int sum = 0;
                    for (int s = 0; s < cpuSize; s++) {
                        sum += cpuGrid[startIdx + s + 1] - gpuGrid[startIdx + s + 1];
                    }

#if log
                    if (sum != 0) {
                        std::cout << i << " " << j << " " << k << "; sum: " << sum << std::endl;
                        for (int s = 0; s < cpuSize; s++) {
                            std::cout << cpuGrid[startIdx + s + 1] << " " << gpuGrid[startIdx + s + 1] << std::endl;
                        }
                        std::cout << "-----------" << std::endl;
                    }
#endif
                    assert((sum == 0) && "Sum isn't zero. -> There are different particles contained in the cell in CPU vs GPU version.");
                }
            }
        }

        auto gpuParticles = sBuffers[1].getData<Particle>();
        auto cpuParticles = sim.particles2;

        assert((gpuParticles.size() == cpuParticles.size()) && "Particle list sizes are different!");
        //std::cout << "rhos" << std::endl;
        for (int i = 0; i < gpuParticles.size(); i++) {
#if log
            std::cout << i << " gpu: " << gpuParticles[i].rho << "; cpu: " << cpuParticles[i].rho << ";" << std::endl;
            std::cout << "CPU pos: " << cpuParticles[i].x.x << " " << cpuParticles[i].x.y << " " << cpuParticles[i].x.z <<
                  "; CPU velocity: " << cpuParticles[i].v.x << " " << cpuParticles[i].v.y << " " << cpuParticles[i].v.z << std::endl;
            std::cout << "GPU pos: " << gpuParticles[i].x.x << " " << gpuParticles[i].x.y << " " << gpuParticles[i].x.z <<
                  "; GPU velocity: " << gpuParticles[i].v.x << " " << gpuParticles[i].v.y << " " << gpuParticles[i].v.z <<std::endl << std::endl;
#endif
            if (!errorMargin(gpuParticles[i].rho, cpuParticles[i].rho, 0.1)) std::cout << i << " gpu: " << gpuParticles[i].rho << "; parallelCpuSim: " << cpuParticles[i].rho << ";" << std::endl;
            //assert((errorMargin(gpuParticles[i].rho, cpuParticles[i].rho, 0.1)) && "Rhos are not within error margin");
        }
    }
#endif
}