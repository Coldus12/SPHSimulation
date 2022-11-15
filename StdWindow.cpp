#include "StdWindow.hpp"
#include "VltavaFunctions.hpp"
#include "ComputeShader.hpp"
#include "ParticleModel.hpp"
#include "CPUSim.hpp"
#include "ImLog.hpp"

#define IMGUI_ENABLED true
#define list_size 25

static Vltava::ImLog my_log;

namespace Vltava {
    bool StdWindow::run = false;
    bool StdWindow::rot = false;
    bool StdWindow::logB = false;
    bool StdWindow::parallelCpuSim = false;

    // Constructor
    //------------------------------------------------------------------------------------------------------------------
    StdWindow::StdWindow(int width, int height) : width(width), height(height) {
        VulkanResources::getInstance().FRAMES_IN_FLIGHT = 2;

        glfwInit();
        glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API);
        //glfwWindowHint(GLFW_RESIZABLE, GLFW_FALSE);

        window = glfwCreateWindow(width, height, "Window", nullptr, nullptr);
        glfwSetWindowUserPointer(window, this);
        glfwSetFramebufferSizeCallback(window, frameBufferResizeCallback);

        glfwSetKeyCallback(window, key_callback);

        initVulkan();
#if IMGUI_ENABLED
        initImgui();
#endif
        mainloop();
    }

    // Destructor
    //------------------------------------------------------------------------------------------------------------------
    StdWindow::~StdWindow() {
        cleanup();
    }

    // Cleaning up
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::cleanup() {
#if IMGUI_ENABLED
        vkDestroyDescriptorPool(
                VulkanResources::getInstance().logDev->getHandle(),
                imguiPool,
                nullptr
        );
        ImGui_ImplVulkan_Shutdown();
#endif
        for (int i = 0; i < VulkanResources::getInstance().FRAMES_IN_FLIGHT; i++) {
            VulkanResources::getInstance().logDev->getHandle().destroySemaphore(renderFinishedSemaphores[i]);
            VulkanResources::getInstance().logDev->getHandle().destroySemaphore(imageAvailableSemaphores[i]);
            VulkanResources::getInstance().logDev->getHandle().destroyFence(inFlightFences[i]);
        }
    }

    // Initializing Vulkan environment
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::initVulkan() {
        gridA = {-6,-1.5,-0.2};
        gridB = {2, 2, 2};

        VulkanWrapper::getInstance().createVulkanWindow(window);
        createSyncObjects();
        sesph = std::make_unique<SESPH>();
        iisph = std::make_unique<IISPH>();
        pcisph = std::make_unique<PCISPH>();

        setComputeData();

        loadModel();
    }

    // Sync objects
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::createSyncObjects() {
        imageAvailableSemaphores.reserve(VulkanResources::getInstance().FRAMES_IN_FLIGHT);
        renderFinishedSemaphores.reserve(VulkanResources::getInstance().FRAMES_IN_FLIGHT);
        inFlightFences.reserve(VulkanResources::getInstance().FRAMES_IN_FLIGHT);

        vk::FenceCreateInfo fenceInfo(vk::FenceCreateFlagBits::eSignaled);
        vk::SemaphoreCreateInfo semaphoreCreateInfo;

        for (int i = 0; i < VulkanResources::getInstance().FRAMES_IN_FLIGHT; i++) {
            imageAvailableSemaphores.push_back(VulkanResources::getInstance().logDev->getHandle().createSemaphore(semaphoreCreateInfo));
            renderFinishedSemaphores.push_back(VulkanResources::getInstance().logDev->getHandle().createSemaphore(semaphoreCreateInfo));
            inFlightFences.push_back(VulkanResources::getInstance().logDev->getHandle().createFence(fenceInfo));
        }
    }

    // Framebuffer resize callback function
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::frameBufferResizeCallback(GLFWwindow *window, int width, int height) {
        auto app = reinterpret_cast<StdWindow *>(glfwGetWindowUserPointer(window));
        app->framebufferResized = true;
    }

    // Try changing data to std::vec, and only allocate gpu memory after you are done filling said vector.
    void StdWindow::setComputeData() {
        //cpusim = std::make_unique<CPUSim>();

        // Water rest density = 1000 kg/m^3
        //
        // Particle r = 0.05 meter
        // r^3 = 0.000125
        // 4/3 * PI ~ 4.1888
        // Particle volume = 5.2359 * 10^-4 m^3 = 0.0005236 m^3 (4/3 * 0.05^3 * PI)
        // Particle mass = volume * density = 0.5236 kg
        //
        // Smoothing length := 3 * 2 * r = 0.1 * 3 meter (for now)
        int pnrAlongAxis = round(pow(particleNr, 1.0/3.0));
        std::vector<Particle> particles;

        nrOfP = pnrAlongAxis * pnrAlongAxis * pnrAlongAxis;
        particleNr = nrOfP;

        // Solid box
        //int bsize = 16;

        float s=0.1;
        //vk::DeviceSize size = sizeof(Particle) * nrOfP * bsize * bsize;
        //auto* data = new Particle[nrOfP + bsize * bsize];

        int r = -1;
        int z = -1;

        float mass = 0.5236f;

        float distAxis = s * pnrAlongAxis / 2.0f;

        // Actual simulated particles
        //--------------------------------------------------------------------------------------------------------------
        for (int i = 0; i < nrOfP; i++) {
            Particle data;

            if (i % pnrAlongAxis == 0)
                r++;

            if (i % (pnrAlongAxis * pnrAlongAxis) == 0)
                z++;

            data.x = glm::vec3((i%pnrAlongAxis) * s - distAxis,(r%pnrAlongAxis) * s - distAxis, z * s);
            data.h = 0.1; // Only sets size for visualization atm
            data.v = glm::vec3(0,0,0);
            data.m = mass;
            //data[i].m = 1.0f;

            data.rho = 1;
            data.p = 1;

            data.staticP = 0;
            data.padding = 0;

            particles.push_back(data);
        }

        float dist = 0.1;
        r = -1;
        z = -1;

        // Container
        //--------------------------------------------------------------------------------------------------------------
        if (makeContainer) {
            auto container = createContainerAt(containerPos, dist, containerSize);
            particles.insert(particles.end(), std::begin(container), std::end(container));
        }

        all_particle_nr = particles.size();

        vk::DeviceSize size = sizeof(Particle) * all_particle_nr;

        Buffer propsBuffer(
                sizeof(SimProps),
                vk::BufferUsageFlagBits::eUniformBuffer,
                vk::MemoryPropertyFlagBits::eHostCoherent | vk::MemoryPropertyFlagBits::eHostVisible
        );
        uBuffers.push_back(std::move(propsBuffer));

        props = {
                1000.0f,
                50.0f,
                (float) all_particle_nr,
                0.2f,

                glm::vec4(gridA, 0),
                glm::vec4(gridB, 0)
        };

        //cpusim->setSimProps(props);
        //cpusim->setData(particles);

        uBuffers[0].writeToBuffer(&props, sizeof(props));

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
        outBuffer.setSize(all_particle_nr * sizeof(Particle));
        outBuffer.bind(0);

        cellx = int(ceil(abs((props.gridB.x - props.gridA.x)/props.kernelh))); // Number of cells in x direction
        celly = int(ceil(abs((props.gridB.y - props.gridA.y)/props.kernelh))); // Number of cells in y direction
        cellz = int(ceil(abs((props.gridB.z - props.gridA.z)/props.kernelh))); // Number of cells in z direction
        
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

        sesph->initGpuSim(&uBuffers, &sBuffers);
        sesph->setCellSizes(cellx, celly, cellz, list_size);
        //iisph->setBuffers(&uBuffers, &sBuffers);

        iisph->initGpuSim(&uBuffers, &sBuffers);
        pcisph->setBuffers(&uBuffers, &sBuffers);
    }

    void StdWindow::resetData() {
        sBuffers.clear();
        uBuffers.clear();
        setComputeData();
        model->changeModel(all_particle_nr, &uBuffers, &sBuffers);
    }

    std::vector<Particle> StdWindow::createContainerAt(glm::vec3 pos, float dist, int sideLength) {
        std::vector<Particle> container;
        container.reserve(sideLength * sideLength * 5);

        float mass = 0.5236f;

        float m = (sideLength * dist) / 2.0f;
        glm::vec3 middle = glm::vec3(pos.x - m, pos.y - m, pos.z);

        int r = -1;
        for (int i = 0; i < sideLength * sideLength; i++) {
            if (i%sideLength == 0)
                r++;

            Particle data;

            // Bottom first
            data.x = glm::vec3((i % sideLength) * dist + middle.x,(r % sideLength) * dist + middle.y, middle.z);
            data.h = 0.1f;
            data.v = glm::vec3(0,0,0);
            data.m = mass;
            //data[i].m = 1.0f;

            data.rho = 1;
            data.p = 1;

            data.staticP = 1;
            data.padding = 0;

            container.push_back(data);

            // Second bottom
            data.x = glm::vec3((i % sideLength) * dist + middle.x,(r % sideLength) * dist + middle.y, middle.z - 1 * dist);
            container.push_back(data);

            // Wall1
            data.x = glm::vec3((i % sideLength) * dist + middle.x,dist + middle.y, (r % sideLength) * dist + middle.z);
            container.push_back(data);

            data.x = glm::vec3((i % sideLength) * dist + middle.x,-dist*0 + middle.y, (r % sideLength) * dist + middle.z);
            container.push_back(data);

            // Wall2
            data.x = glm::vec3(dist*1 + middle.x,(i % sideLength) * dist + middle.y, (r % sideLength) * dist + middle.z);
            container.push_back(data);

            data.x = glm::vec3(-dist*0 + middle.x,(i % sideLength) * dist + middle.y, (r % sideLength) * dist + middle.z);
            container.push_back(data);

            // Wall3
            data.x = glm::vec3((i % sideLength) * dist + middle.x, (sideLength - 1) * dist + middle.y, (r % sideLength) * dist + middle.z);
            container.push_back(data);
            data.x = glm::vec3((i % sideLength) * dist + middle.x, (sideLength - 2)* 1 * dist + middle.y, (r % sideLength) * dist + middle.z);
            container.push_back(data);

            // Wall4
            data.x = glm::vec3((sideLength + 1*0) * dist + middle.x, (i % sideLength) * dist + middle.y, (r % sideLength) * dist + middle.z);
            container.push_back(data);
            data.x = glm::vec3((sideLength - 1) * 1 * dist + middle.x, (i % sideLength) * dist + middle.y, (r % sideLength) * dist + middle.z);
            container.push_back(data);
        }

        return container;
    }

    void StdWindow::runComp() {
        if (run) {
            if (!cpuSim) {
                if (sph)
                    sesph->gpuTimeStep();
                else
                    iisph->gpuTimeStep();
                //dispatchCompute(nrOfP, 1, 1);
            } else
                runCpuSim(nrOfIter);

            log();
            run = false;
        }
    }

    void StdWindow::runCpuSim(int iterNr) {
        //cpusim->runIISPH(iterNr);
        //cpusim->runSESPH(iterNr);
        std::vector<Particle> particles;
        if (sph) {
            sesph->cpuTimeStep();
            particles = sesph->getFirst() ? sesph->getData1() : sesph->getData2();
        } else {
            //iisph->cpuTimeStep();
            //particles = iisph->getFirst() ? iisph->getData1() : iisph->getData2();
            //uBuffers.clear();

            pcisph->cpuTimeStep();
            particles = pcisph->getFirst() ? pcisph->getData1() : pcisph->getData2();
        }
        sBuffers.clear();

        uint64_t size = particles.size() * sizeof(Particle);

        //uBuffers[0].writeToBuffer(&props, sizeof(props));

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
        outBuffer.setSize(all_particle_nr * sizeof(Particle));
        outBuffer.bind(0);

        sBuffers.push_back(std::move(inBuffer));
        sBuffers.push_back(std::move(outBuffer));

        for (auto& b: sBuffers)
            b.writeToBuffer(particles.data(), size);

        model->changeModel(particles.size(), &uBuffers, &sBuffers);
    }

    void StdWindow::log() {
        if (write_log) {
            std::string str = "";

            if (sph)
                str = sesph->log();
            else
                str = iisph->log();

            my_log.addLog(str.c_str());

            if (console_log)
                std::cout << str << std::endl;
        }
    }

    // Main loop
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::mainloop() {
        resetData();
        start = std::chrono::high_resolution_clock::now();

        while (!glfwWindowShouldClose(window)) {
            glfwPollEvents();
            if (logB) {
                logB = false;
                log();
            }

            if (realTime) {
                auto now = std::chrono::high_resolution_clock::now();

                // "If timeout is zero, then vkWaitForFences does not wait, but simply returns the current state of the fences."
                //if (VulkanResources::getInstance().logDev->getHandle().waitForFences(compFence, false, 0) != vk::Result::eTimeout) {
                    // Sadly the fence does not seem to work the way i hoped it would
                    float time = std::chrono::duration<float, std::chrono::milliseconds::period>(now - start).count();

                    time /= 10;
                    //if (time < 1) time = 1;
                    if (time > 1) {
                        //VulkanResources::getInstance().logDev->getHandle().resetFences(compFence);
                        std::cout << "time: " << time << std::endl;
                        if (time > 1)
                            time = 1;

                        //std::cout << "dispatch nr: " << time << std::endl;

                        nrOfIter = time;

                        if (!cpuSim)
                            //dispatchCompute(nrOfP, 1, 1);
                            if (sph)
                                sesph->gpuTimeStep();
                            else
                                iisph->gpuTimeStep();
                        else
                            runCpuSim(nrOfIter);
                        //sesph->gpuTimeStep();


                        log();
                        start = std::chrono::high_resolution_clock::now();
                    }
                //}
            }

#if IMGUI_ENABLED
            ImGui_ImplVulkan_NewFrame();
            ImGui_ImplGlfw_NewFrame();
            ImGui::NewFrame();
            //ImGui::ShowDemoWindow();
            ImGui::Begin("Test");

            ImGui::Text("Reset data??");
            if (ImGui::Button("Yes"))
                resetData();

            ImGui::Checkbox("SESPH or IISPH, checked means SE", &sph);

            if (ImGui::Button("Toggle rotation"))
                rot = !rot;

            if (ImGui::Button("Iter"))
                run = true;

            if (ImGui::Button("CPU")) {
                parallelCpuSim = !parallelCpuSim;
            }

            ImGui::Text("Number of iterations: "); ImGui::SameLine(); ImGui::InputInt("##", &nrOfIter);
            ImGui::Text("Number of particles: "); ImGui::SameLine(); ImGui::InputInt("###", &particleNr);

            ImGui::Checkbox("Show log?", &show_log);
            ImGui::Checkbox("Write log?", &write_log);
            ImGui::Checkbox("Write log into cli?", &console_log);
            //ImGui::Checkbox("Real time?", &realTime);

            if (ImGui::Button("Real time flip")) {
                realTime = !realTime;
                start = std::chrono::high_resolution_clock::now();
            }

            ImGui::Checkbox("Try cpuSim?", &cpuSim);

            ImGui::Checkbox("Set boundaries?", &setBoundaries);
            if (setBoundaries) {
                ImGui::Text("gridA: "); ImGui::SameLine(); ImGui::InputFloat3("##", (float*)&gridA);
                ImGui::Text("gridB: "); ImGui::SameLine(); ImGui::InputFloat3("s", (float*)&gridB);
            }

            ImGui::Checkbox("Create container?", &makeContainer);
            if (makeContainer) {
                ImGui::Text("Container position: "); ImGui::SameLine(); ImGui::InputFloat3("##", (float*)&containerPos);
                ImGui::Text("Container size: "); ImGui::SameLine(); ImGui::InputInt("s", &containerSize);

                props = {
                        1000.0f,
                        50.0f,
                        (float) all_particle_nr,
                        0.2f,

                        glm::vec4(gridA, 0),
                        glm::vec4(gridB, 0)
                };

                //cpusim->setSimProps(props);
                uBuffers[0].writeToBuffer(&props, sizeof(props));
            }

            if (show_log)
                my_log.draw("Log");

            ImGui::End();
#endif
            runComp();
            drawFrame();
        }

        VulkanResources::getInstance().logDev->getHandle().waitIdle();
    }

    // GLFW window keyboard callback function
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
        if (key == GLFW_KEY_W && action == GLFW_PRESS) {
            run = true;
        }

        if (key == GLFW_KEY_Q && action == GLFW_PRESS) {
            rot = !rot;
        }

        if (key == GLFW_KEY_E && action == GLFW_PRESS) {
            logB = true;
        }

        if (key == GLFW_KEY_C && action == GLFW_PRESS) {
            parallelCpuSim = !parallelCpuSim;
        }
    }

    // DearImgui initialization
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::initImgui() {
        auto& commandBuffers = VulkanWrapper::getInstance().getCmdBuffers();

        //1: create descriptor pool for IMGUI
        // the size of the pool is very oversize, but it's copied from imgui demo itself.
        VkDescriptorPoolSize pool_sizes[] =
                {
                        { VK_DESCRIPTOR_TYPE_SAMPLER, 1000 },
                        { VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1000 },
                        { VK_DESCRIPTOR_TYPE_SAMPLED_IMAGE, 1000 },
                        { VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, 1000 },
                        { VK_DESCRIPTOR_TYPE_UNIFORM_TEXEL_BUFFER, 1000 },
                        { VK_DESCRIPTOR_TYPE_STORAGE_TEXEL_BUFFER, 1000 },
                        { VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1000 },
                        { VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1000 },
                        { VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER_DYNAMIC, 1000 },
                        { VK_DESCRIPTOR_TYPE_STORAGE_BUFFER_DYNAMIC, 1000 },
                        { VK_DESCRIPTOR_TYPE_INPUT_ATTACHMENT, 1000 }
                };

        VkDescriptorPoolCreateInfo pool_info = {};
        pool_info.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
        pool_info.flags = VK_DESCRIPTOR_POOL_CREATE_FREE_DESCRIPTOR_SET_BIT;
        pool_info.maxSets = 1000;
        pool_info.poolSizeCount = std::size(pool_sizes);
        pool_info.pPoolSizes = pool_sizes;

        vkCreateDescriptorPool(
                VulkanResources::getInstance().logDev->getHandle(),
                &pool_info,
                nullptr,
                &imguiPool
        );

        // 2: initialize imgui library

        //this initializes the core structures of imgui
        ImGui::CreateContext();

        //this initializes imgui for GLFW
        ImGui_ImplGlfw_InitForVulkan(window, true);

        //this initializes imgui for Vulkan
        ImGui_ImplVulkan_InitInfo init_info = {};
        init_info.Instance = VulkanResources::getInstance().instance->getHandle();
        init_info.PhysicalDevice = VulkanResources::getInstance().physDev->getHandle();
        init_info.Device = VulkanResources::getInstance().logDev->getHandle();
        init_info.Queue = *VulkanResources::getInstance().graphicsQueue;
        init_info.DescriptorPool = imguiPool;
        init_info.MinImageCount = 3;
        init_info.ImageCount = 3;
        init_info.MSAASamples = VK_SAMPLE_COUNT_1_BIT;

        ImGui_ImplVulkan_Init(&init_info, *VulkanResources::getInstance().renderPass);

        //execute a gpu command to upload imgui font textures
        vk::CommandBufferBeginInfo beginInfo({}, nullptr);
        commandBuffers[currentFrame].begin(beginInfo);
        ImGui_ImplVulkan_CreateFontsTexture(commandBuffers[currentFrame]);
        commandBuffers[currentFrame].end();

        vk::SubmitInfo submitInfo({},{}, commandBuffers[currentFrame]);
        std::array<vk::SubmitInfo, 1> infos = {submitInfo};

        VulkanResources::getInstance().graphicsQueue->submit(infos, inFlightFences[currentFrame]);
        VulkanResources::getInstance().logDev->getHandle().waitIdle();

        //clear font textures from parallelCpuSim data
        ImGui_ImplVulkan_DestroyFontUploadObjects();
    }

    // Load model
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::loadModel() {
        model = std::make_unique<ParticleModel>(&uBuffers, &sBuffers);
    }

    // Draw frame
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::drawFrame() {
        auto& commandBuffers = VulkanWrapper::getInstance().getCmdBuffers();

#if IMGUI_ENABLED
        ImGui::Render();
#endif
        VulkanResources::getInstance().logDev->getHandle().waitForFences(inFlightFences[currentFrame], true, UINT64_MAX);
        //VulkanResources::getInstance().logDev->getHandle().waitForFences(inFlightFences[currentFrame], true, UINT64_MAX);
        //VulkanResources::getInstance().logDev->getHandle().waitForFences(inFlightFences[currentFrame], true, 10);

        vk::Result result;
        uint32_t imageIndex;

        std::tie(result, imageIndex) = VulkanResources::getInstance().logDev->getHandle().acquireNextImageKHR(
                VulkanResources::getInstance().swapChain->getHandle(),
                UINT64_MAX,
                imageAvailableSemaphores[currentFrame]
        );

        if (result == vk::Result::eErrorOutOfDateKHR) {
            VulkanWrapper::getInstance().recreateSwapChain(window);
            model->recreatePipeline();
            return;
        } else if (result != vk::Result::eSuccess && result != vk::Result::eSuboptimalKHR) {
            throw std::runtime_error("Failed to acquire swap chain image!");
        }

        VulkanResources::getInstance().logDev->getHandle().resetFences(inFlightFences[currentFrame]);

        commandBuffers[currentFrame].reset();
        recordCommandBuffer(imageIndex);

        vk::PipelineStageFlags flags(vk::PipelineStageFlagBits::eColorAttachmentOutput);
        vk::SubmitInfo submitInfo(
                (uint32_t) 1,
                &imageAvailableSemaphores[currentFrame],
                &flags,
                (uint32_t) 1,
                &commandBuffers[currentFrame],
                (uint32_t) 1,
                &renderFinishedSemaphores[currentFrame]
        );
        std::array<vk::SubmitInfo, 1> infos = {submitInfo};

        VulkanResources::getInstance().graphicsQueue->submit(infos, inFlightFences[currentFrame]);

        std::array<vk::SwapchainKHR, 1> swapChains = {VulkanResources::getInstance().swapChain->getHandle()};
        vk::PresentInfoKHR presentInfo(renderFinishedSemaphores[currentFrame], swapChains, imageIndex);

        /** For some reason instead of returning the vk::eErrorOutOfDateKHR result,
         * this just throws an exception with the this error. But this try/catch block
         * seems to fix the problem.
         * */
        try {
            result = VulkanResources::getInstance().presentQueue->presentKHR(presentInfo);
        } catch (vk::SystemError &err) {
            framebufferResized = false;
            VulkanWrapper::getInstance().recreateSwapChain(window);
            model->recreatePipeline();
        }

        if (result == vk::Result::eErrorOutOfDateKHR || result == vk::Result::eSuboptimalKHR || framebufferResized) {
            framebufferResized = false;
            VulkanWrapper::getInstance().recreateSwapChain(window);
            model->recreatePipeline();
        } else if (result != vk::Result::eSuccess) {
            throw std::runtime_error("Failed to present swap chain image!");
        }

        currentFrame = (currentFrame + 1) % VulkanResources::getInstance().FRAMES_IN_FLIGHT;
    }

    void StdWindow::recordCommandBuffer(uint32_t imageIndex) {
        auto& commandBuffers = VulkanWrapper::getInstance().getCmdBuffers();
        auto& swapChainFramebuffers = VulkanWrapper::getInstance().getSCFrameBuffers();

        vk::CommandBufferBeginInfo beginInfo({}, nullptr);
        commandBuffers[currentFrame].begin(beginInfo);

        std::array<vk::ClearValue, 2> clearValues{};
        std::array<float, 4> color = {0.0f, 0.0f, 0.0f, 1.0f};
        clearValues[0].color = vk::ClearColorValue(color);
        clearValues[1].depthStencil = vk::ClearDepthStencilValue(1.0f, 0);

        vk::RenderPassBeginInfo renderPassInfo(
                *VulkanResources::getInstance().renderPass,
                swapChainFramebuffers[imageIndex],
                {{0, 0}, VulkanResources::getInstance().extent},
                static_cast<uint32_t>(clearValues.size()),
                clearValues.data()
        );

        commandBuffers[currentFrame].beginRenderPass(renderPassInfo, vk::SubpassContents::eInline);

        // Drawing the models
        model->draw(commandBuffers[currentFrame],currentFrame);
        //model->draw(commandBuffers->at(currentFrame));

#if IMGUI_ENABLED
        // Recording imgui draw cmds into the commandBuffer
        ImGui_ImplVulkan_RenderDrawData(ImGui::GetDrawData(), commandBuffers[currentFrame]);
#endif

        commandBuffers[currentFrame].endRenderPass();
        commandBuffers[currentFrame].end();
    }
}