#include "StdWindow.hpp"
#include "VltavaFunctions.hpp"
#include "ComputeShader.hpp"
#include "ParticleModel.hpp"

namespace Vltava {
    bool StdWindow::run1 = false;
    bool StdWindow::run2 = false;
    bool StdWindow::rot = true;

    // Main loop
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::mainloop() {
        while (!glfwWindowShouldClose(window)) {
            glfwPollEvents();
            runComp();
            drawFrame();
        }

        logicalDevice->getHandle().waitIdle();
    }

    void StdWindow::runComp() {
        if (run1) {
            for (int i = 0; i < 100; i++) {
                comp1->dispatch(computeCmdBuffer, 64, 1, 1);
                comp2->dispatch(computeCmdBuffer, 64, 1, 1);
            }

            auto spheres = sBuffers[1].getData<Particle>();
            std::cout << spheres.size() << std::endl;

            for (int i = 0; i < 64; i++) {
                std::cout << "Density: " << spheres[i].rho << "; Pressure: " << spheres[i].p << " ; Position: " << spheres[i].x.x << " " << spheres[i].x.y << " " << spheres[i].x.z << " ; Mass: " << spheres[i].m  << " ; Padding "  << spheres[i].padding1 <<";\n";
            }
            std::cout << std::endl;
            std::cout << "Done" << std::endl;

            run1 = false;
        }
    }

    void StdWindow::computeStuff() {
        VulkanResources updated{
                &renderPass,
                &*physicalDevice,
                &*logicalDevice,
                &*instance,
                &commandPool,
                &*graphicsQueue,
                &*computeQueue,
                swapChainExtent,
                MAX_FRAMES_IN_FLIGHT
        };

        // density / pressure calculation
        SimProps props{
            1.0f,
            1.0f,
            64.0f,
            3.0f
        };

        Buffer UBO(
                updated,
                sizeof(SimProps),
                vk::BufferUsageFlagBits::eUniformBuffer,
                vk::MemoryPropertyFlagBits::eHostCoherent | vk::MemoryPropertyFlagBits::eHostVisible
         );
        UBO.writeToBuffer(&props, sizeof(SimProps));

        //std::vector<Buffer> uBuffers;
        uBuffers.push_back(std::move(UBO));

        int nrOfP = 64;
        float s=0.1;
        vk::DeviceSize size = sizeof(Particle) * nrOfP;
        auto* data = new Particle[nrOfP];

        int r = -1;
        int z = -1;

        float mass = (2.0f/3.0f * 1) * (2.0f/3.0f * 1) * (2.0f/3.0f * 1);

        for (int i = 0; i < nrOfP; i++) {
            if (i%4 == 0)
                r++;

            if (i%16 == 0)
                z++;

            data[i].x = glm::vec3((i%4) * s,(r%4) * s, z * s);

            data[i].h = 1;
            data[i].v = glm::vec3(0,0,0);
            data[i].m = mass;
            //data[i].m = 1.0f;

            data[i].rho = 0;
            data[i].p = 0;

            data[i].padding1 = 0;
            data[i].padding2 = 0;
        }

        Buffer inBuffer(
                updated,
                size,
                vk::BufferUsageFlagBits::eStorageBuffer,
                vk::MemoryPropertyFlagBits::eHostVisible | vk::MemoryPropertyFlagBits::eHostCoherent
        );
        inBuffer.bind(0);

        Buffer outBuffer(
                updated,
                size,
                vk::BufferUsageFlagBits::eStorageBuffer,
                vk::MemoryPropertyFlagBits::eHostVisible | vk::MemoryPropertyFlagBits::eHostCoherent
        );
        outBuffer.setSize(nrOfP * sizeof(Particle));
        outBuffer.bind(0);

        inBuffer.writeToBuffer(data, size);
        delete[] data;

        //std::vector<Buffer> sBuffers;
        sBuffers.push_back(std::move(inBuffer));
        sBuffers.push_back(std::move(outBuffer));

        //ComputeShader comp(updated, "shaders/comp.spv");
        comp1 = std::make_unique<ComputeShader>(updated, "shaders/comp.spv");
        comp1->setBuffers(&uBuffers, &sBuffers);
        comp1->createPipeline();
        //comp.createCommandBuffer(computeQueueFamily);
        comp1->dispatch(computeCmdBuffer, size / sizeof(Particle) + 1, 1, 1);

        auto spheres = sBuffers[1].getData<Particle>();
        std::cout << spheres.size() << std::endl;

        for (int i = 0; i < nrOfP; i++) {
            std::cout << "Density: " << spheres[i].rho << "; Pressure: " << spheres[i].p << ";\n";
        }
        std::cout << std::endl;
        std::cout << "Done" << std::endl;

        //ComputeShader comp2(updated, "shaders/comp_it.spv");
        comp2 = std::make_unique<ComputeShader>(updated, "shaders/comp_it.spv");
        comp2->setBuffers(&uBuffers, &sBuffers);
        comp2->createPipeline();
        //comp.createCommandBuffer(computeQueueFamily);
        comp2->dispatch(computeCmdBuffer, size / sizeof(Particle) + 1, 1, 1);
    }

    void StdWindow::key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
        if (key == GLFW_KEY_W && action == GLFW_PRESS) {
            run1 = true;
        }

        /*if (key == GLFW_KEY_S && action == GLFW_PRESS) {
            run2 = true;
        }*/

        if (key == GLFW_KEY_Q && action == GLFW_PRESS) {
            rot = !rot;
        }
    }

    // Constructor
    //------------------------------------------------------------------------------------------------------------------
    StdWindow::StdWindow(int width, int height) : width(width), height(height) {
        glfwInit();
        glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API);
        //glfwWindowHint(GLFW_RESIZABLE, GLFW_FALSE);

        window = glfwCreateWindow(width, height, "Window", nullptr, nullptr);
        glfwSetWindowUserPointer(window, this);
        glfwSetFramebufferSizeCallback(window, frameBufferResizeCallback);

        glfwSetKeyCallback(window, key_callback);

        initVulkan();
        mainloop();
    }

    void StdWindow::initVulkan() {
        createInstance("Test");
        createSurface();
        selectPhysicalDevice();
        selectQueues();
        createLogicalDevice();
        createSwapChain();
        createImageViews();
        createRenderPass();
        createFramebuffers();
        createCommandPool();
        createCommandBuffers();
        createSyncObjects();

        computeStuff();

        loadModel();
    }

    // Cleaning up
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::cleanup() {
        cleanupSwapChain();

        for (int i = 0; i < MAX_FRAMES_IN_FLIGHT; i++) {
            logicalDevice->getHandle().destroySemaphore(renderFinishedSemaphores[i]);
            logicalDevice->getHandle().destroySemaphore(imageAvailableSemaphores[i]);
            logicalDevice->getHandle().destroyFence(inFlightFences[i]);
        }

        glfwDestroyWindow(window);
        //glfwTerminate();
    }

    // Destructor
    //------------------------------------------------------------------------------------------------------------------
    StdWindow::~StdWindow() {
        logicalDevice->getHandle().waitIdle();
        cleanup();
        logicalDevice->getHandle().waitIdle();
    }

    // Instance creation
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::createInstance(std::string app_name) {
        instance = std::make_unique<MInstance>(app_name, enableValidationLayers);
    }

    // Creating a surface
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::createSurface() {
        surface = std::make_unique<MSurface>(instance->getHandle(), window);
    }

    // Selecting the best gpu for the task
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::selectPhysicalDevice() {
        physicalDevice = std::make_unique<MPhysDev>(instance->getHandle(), deviceExtensions);
    }

    // Selecting queue(s)
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::selectQueues() {
        std::vector<vk::QueueFamilyProperties> queueFamilies = physicalDevice->getHandle().getQueueFamilyProperties();

        uint32_t i = 0;
        for (const auto &queueFamily: queueFamilies) {
            if (physicalDevice->getHandle().getSurfaceSupportKHR(i, surface->getHandle())) {
                presentQueueFamily = i;
            }

            if (queueFamily.queueFlags & vk::QueueFlagBits::eGraphics) {
                graphicsQueueFamily = i;
            }

            if (queueFamily.queueFlags & vk::QueueFlagBits::eCompute) {
                computeQueueFamily = i;
            }
            i++;
        }
    }

    // Creating a logical device
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::createLogicalDevice() {
        std::set<uint32_t> uniqueQueueFamilies({graphicsQueueFamily, presentQueueFamily, computeQueueFamily});

        logicalDevice = std::make_unique<MLogDev>(
                uniqueQueueFamilies,
                physicalDevice->getHandle(),
                deviceExtensions,
                &validationLayers
        );

        // Queue creation
        //--------------------------------------------------------------------------------------------------------------
        graphicsQueue = std::make_unique<vk::Queue>(logicalDevice->getHandle().getQueue(graphicsQueueFamily, 0));
        presentQueue = std::make_unique<vk::Queue>(logicalDevice->getHandle().getQueue(presentQueueFamily, 0));
        computeQueue = std::make_unique<vk::Queue>(logicalDevice->getHandle().getQueue(computeQueueFamily, 0));
    }

    // SwapChain
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::createSwapChain() {
        swapChain = std::make_unique<MSwapChain>(
                window,
                physicalDevice->getHandle(),
                logicalDevice->getHandle(),
                surface->getHandle(),
                vk::Format::eB8G8R8A8Srgb,
                vk::ColorSpaceKHR::eSrgbNonlinear,
                vk::PresentModeKHR::eMailbox,
                graphicsQueueFamily,
                presentQueueFamily
        );

        swapChainImages = logicalDevice->getHandle().getSwapchainImagesKHR(swapChain->getHandle());
        imageViews.reserve(swapChainImages.size());

        swapChainImageFormat = swapChain->getFormat();
        swapChainExtent = swapChain->getExtent();
    }

    void StdWindow::recreateSwapChain() {
        int w = 0, h = 0;
        glfwGetFramebufferSize(window, &w, &h);
        while (w == 0 || h == 0) {
            glfwGetFramebufferSize(window, &w, &h);
            glfwWaitEvents();
        }

        logicalDevice->getHandle().waitIdle();
        cleanupSwapChain();
        createSwapChain();
        createImageViews();
        createRenderPass();

        // Updating / recreating the pipelines of the models
        VulkanResources updated{
            &renderPass,
            &*physicalDevice,
            &*logicalDevice,
            &*instance,
            &commandPool,
            &*graphicsQueue,
            &*computeQueue,
            swapChainExtent,
            MAX_FRAMES_IN_FLIGHT
        };

        model->updateResources(updated);

        createFramebuffers();
    }

    void StdWindow::cleanupSwapChain() {
        for (auto framebuffer : swapChainFramebuffers) {
            logicalDevice->getHandle().destroyFramebuffer(framebuffer);
        }
        swapChainFramebuffers.clear();

        logicalDevice->getHandle().destroyRenderPass(renderPass);

        for (auto imgView: imageViews) {
            logicalDevice->getHandle().destroyImageView(imgView);
        }
        imageViews.clear();

        swapChain.reset();
    }

    // Framebuffer resize callback function
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::frameBufferResizeCallback(GLFWwindow *window, int width, int height) {
        auto app = reinterpret_cast<StdWindow *>(glfwGetWindowUserPointer(window));
        app->framebufferResized = true;
    }

    // ImageView
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::createImageViews() {
        imageViews.reserve(swapChainImages.size());

        vk::ImageViewCreateInfo createInfo(
                {},
                {},
                vk::ImageViewType::e2D,
                swapChainImageFormat,
                {vk::ComponentSwizzle::eIdentity, vk::ComponentSwizzle::eIdentity, vk::ComponentSwizzle::eIdentity,
                 vk::ComponentSwizzle::eIdentity},
                {vk::ImageAspectFlagBits::eColor, 0, 1, 0, 1}
        );

        for (auto image: swapChainImages) {
            createInfo.image = image;
            imageViews.push_back(logicalDevice->getHandle().createImageView(createInfo));
        }
    }

    // RenderPass
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::createRenderPass() {
        vk::SubpassDependency dependency(
                VK_SUBPASS_EXTERNAL,
                0,
                vk::PipelineStageFlagBits::eColorAttachmentOutput,
                vk::PipelineStageFlagBits::eColorAttachmentOutput,
                {},
                vk::AccessFlagBits::eColorAttachmentWrite
        );

        vk::AttachmentDescription colorAttachment(
                {},
                swapChainImageFormat,
                vk::SampleCountFlagBits::e1,
                vk::AttachmentLoadOp::eClear,
                vk::AttachmentStoreOp::eStore,
                vk::AttachmentLoadOp::eDontCare,
                vk::AttachmentStoreOp::eDontCare,
                vk::ImageLayout::eUndefined,
                vk::ImageLayout::ePresentSrcKHR
        );

        vk::AttachmentReference colorReference(0, vk::ImageLayout::eColorAttachmentOptimal);

        vk::SubpassDescription subpass({}, vk::PipelineBindPoint::eGraphics, {}, colorReference);
        vk::RenderPassCreateInfo createInfo(
                {},
                1,
                &colorAttachment,
                1,
                &subpass,
                1,
                &dependency
        );

        renderPass = logicalDevice->getHandle().createRenderPass(createInfo);
    }

    // Framebuffers
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::createFramebuffers() {
        swapChainFramebuffers.reserve(imageViews.size());

        for (int i = 0; i < imageViews.size(); i++) {
            vk::FramebufferCreateInfo framebufferInfo(
                    {},
                    renderPass,
                    1,
                    &imageViews[i],
                    swapChainExtent.width,
                    swapChainExtent.height, 1
            );

            swapChainFramebuffers.push_back(logicalDevice->getHandle().createFramebuffer(framebufferInfo));
        }
    }

    // Commandpool
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::createCommandPool() {
        vk::CommandPoolCreateInfo poolInfo(vk::CommandPoolCreateFlagBits::eResetCommandBuffer, graphicsQueueFamily);
        commandPool = logicalDevice->getHandle().createCommandPool(poolInfo);

        vk::CommandPoolCreateInfo cmdPoolInfo(vk::CommandPoolCreateFlagBits::eResetCommandBuffer, computeQueueFamily);
        computeCmdPool = logicalDevice->getHandle().createCommandPool(cmdPoolInfo);
    }

    // Command buffers
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::createCommandBuffers() {
        vk::CommandBufferAllocateInfo allocInfo(commandPool, vk::CommandBufferLevel::ePrimary, MAX_FRAMES_IN_FLIGHT);
        commandBuffers = logicalDevice->getHandle().allocateCommandBuffers(allocInfo);

        vk::CommandBufferAllocateInfo cmdBufferInfo(computeCmdPool, vk::CommandBufferLevel::ePrimary, 1);
        computeCmdBuffer = logicalDevice->getHandle().allocateCommandBuffers(cmdBufferInfo).front();
    }

    // Sync objects
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::createSyncObjects() {
        imageAvailableSemaphores.reserve(MAX_FRAMES_IN_FLIGHT);
        renderFinishedSemaphores.reserve(MAX_FRAMES_IN_FLIGHT);
        inFlightFences.reserve(MAX_FRAMES_IN_FLIGHT);

        vk::FenceCreateInfo fenceInfo(vk::FenceCreateFlagBits::eSignaled);
        vk::SemaphoreCreateInfo semaphoreCreateInfo;

        for (int i = 0; i < MAX_FRAMES_IN_FLIGHT; i++) {
            imageAvailableSemaphores.push_back(logicalDevice->getHandle().createSemaphore(semaphoreCreateInfo));
            renderFinishedSemaphores.push_back(logicalDevice->getHandle().createSemaphore(semaphoreCreateInfo));
            inFlightFences.push_back(logicalDevice->getHandle().createFence(fenceInfo));
        }
    }

    // Load model
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::loadModel() {
        VulkanResources res{
            &renderPass,
            &*physicalDevice,
            &*logicalDevice,
            &*instance,
            &commandPool,
            &*graphicsQueue,
            &*computeQueue,
            swapChainExtent,
            MAX_FRAMES_IN_FLIGHT
        };

        //model = std::make_unique<Model>(res);
        //model->loadModel(""); // Actual load

        model = std::make_unique<ParticleModel>(res, &uBuffers, &sBuffers);
    }

    // Draw frame
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::drawFrame() {
        logicalDevice->getHandle().waitForFences(inFlightFences[currentFrame], true, UINT64_MAX);

        vk::Result result;
        uint32_t imageIndex;

        std::tie(result, imageIndex) = logicalDevice->getHandle().acquireNextImageKHR(swapChain->getHandle(), UINT64_MAX, imageAvailableSemaphores[currentFrame]);

        if (result == vk::Result::eErrorOutOfDateKHR) {
            recreateSwapChain();
            return;
        } else if (result != vk::Result::eSuccess && result != vk::Result::eSuboptimalKHR) {
            throw std::runtime_error("Failed to acquire swap chain image!");
        }

        logicalDevice->getHandle().resetFences(inFlightFences[currentFrame]);

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

        graphicsQueue->submit(infos, inFlightFences[currentFrame]);

        std::array<vk::SwapchainKHR, 1> swapChains = {swapChain->getHandle()};
        vk::PresentInfoKHR presentInfo(renderFinishedSemaphores[currentFrame], swapChains, imageIndex);

        /** For some reason instead of returning the vk::eErrorOutOfDateKHR result,
         * this just throws an exception with the this error. But this try/catch block
         * seems to fix the problem.
         * */
        try {
            result = presentQueue->presentKHR(presentInfo);
        } catch (vk::SystemError &err) {
            framebufferResized = false;
            recreateSwapChain();
        }

        if (result == vk::Result::eErrorOutOfDateKHR || result == vk::Result::eSuboptimalKHR || framebufferResized) {
            framebufferResized = false;
            recreateSwapChain();
        } else if (result != vk::Result::eSuccess) {
            throw std::runtime_error("Failed to present swap chain image!");
        }

        currentFrame = (currentFrame + 1) % MAX_FRAMES_IN_FLIGHT;
    }

    void StdWindow::recordCommandBuffer(uint32_t imageIndex) {
        vk::CommandBufferBeginInfo beginInfo({}, nullptr);
        commandBuffers[currentFrame].begin(beginInfo);

        std::array<float, 4> color = {0.0f, 0.0f, 0.0f, 1.0f};
        vk::ClearValue clrVal((vk::ClearColorValue(color)));

        vk::RenderPassBeginInfo renderPassInfo(
                renderPass,
                swapChainFramebuffers[imageIndex],
                {{0, 0}, swapChainExtent},
                1, &clrVal
        );

        commandBuffers[currentFrame].beginRenderPass(renderPassInfo, vk::SubpassContents::eInline);

        // Drawing the models
        model->draw(commandBuffers[currentFrame],currentFrame);
        //model->draw(commandBuffers->at(currentFrame));

        commandBuffers[currentFrame].endRenderPass();
        commandBuffers[currentFrame].end();
    }
}