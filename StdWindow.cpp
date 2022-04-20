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
                comp1->dispatch(computeCmdBuffer->getBuffers()[0], 64, 1, 1);
                comp2->dispatch(computeCmdBuffer->getBuffers()[0], 64, 1, 1);
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

        /*if (run1) {
            if (comp1 != nullptr) {
                comp1->dispatch(*computeCmdBuffer, 64 / sizeof(Particle) + 1, 1, 1);
                run1 = false;
            }
        }

        if (run2) {
            if (comp2 != nullptr) {
                comp2->dispatch(*computeCmdBuffer, 64 / sizeof(Particle) + 1, 1, 1);
                auto spheres = sBuffers[1].getData<Particle>();
                std::cout << spheres.size() << std::endl;

                for (int i = 0; i < 64; i++) {
                    std::cout << "Density: " << spheres[i].rho << "; Pressure: " << spheres[i].p << " ; Position: " << spheres[i].x.x << " " << spheres[i].x.y << " " << spheres[i].x.z <<";\n";
                }
                std::cout << std::endl;
                std::cout << "Done" << std::endl;
                run2 = false;
            }
        }*/
    }

    void StdWindow::computeStuff() {
        VulkanResources updated{
                &*renderPass,
                &*physicalDevice,
                &*logicalDevice,
                &*instance,
                &*commandPool,
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
        comp1->dispatch(computeCmdBuffer->getBuffers()[0], size / sizeof(Particle) + 1, 1, 1);

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
        comp2->dispatch(computeCmdBuffer->getBuffers()[0], size / sizeof(Particle) + 1, 1, 1);
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
        glfwDestroyWindow(window);
    }

    // Destructor
    //------------------------------------------------------------------------------------------------------------------
    StdWindow::~StdWindow() {
        cleanup();
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
            &*renderPass,
            &*physicalDevice,
            &*logicalDevice,
            &*instance,
            &*commandPool,
            &*graphicsQueue,
            &*computeQueue,
            swapChainExtent,
            MAX_FRAMES_IN_FLIGHT
        };

        model->updateResources(updated);

        createFramebuffers();
    }

    void StdWindow::cleanupSwapChain() {
        swapChainFramebuffers.clear();
        renderPass.reset();
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
            imageViews.emplace_back(logicalDevice->getHandle(), createInfo);
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

        renderPass = std::make_unique<MRenderPass>(logicalDevice->getHandle(), subpass, colorAttachment);
    }

    // Framebuffers
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::createFramebuffers() {
        swapChainFramebuffers.reserve(imageViews.size());

        for (int i = 0; i < imageViews.size(); i++) {
            vk::FramebufferCreateInfo framebufferInfo(
                    {},
                    renderPass->getHandle(),
                    1,
                    imageViews[i].getAddress(),
                    swapChainExtent.width,
                    swapChainExtent.height, 1
            );

            swapChainFramebuffers.emplace_back(logicalDevice->getHandle(), framebufferInfo);
        }
    }

    // Commandpool
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::createCommandPool() {
        vk::CommandPoolCreateInfo poolInfo(vk::CommandPoolCreateFlagBits::eResetCommandBuffer, graphicsQueueFamily);
        commandPool = std::make_unique<MCommandPool>(logicalDevice->getHandle(), poolInfo);

        vk::CommandPoolCreateInfo cmdPoolInfo(vk::CommandPoolCreateFlagBits::eResetCommandBuffer, computeQueueFamily);
        computeCmdPool = std::make_unique<MCommandPool>(logicalDevice->getHandle(), cmdPoolInfo);
    }

    // Command buffers
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::createCommandBuffers() {
        vk::CommandBufferAllocateInfo allocInfo(commandPool->getHandle(), vk::CommandBufferLevel::ePrimary, MAX_FRAMES_IN_FLIGHT);
        commandBuffers = std::make_unique<MCommandBuffers>(logicalDevice->getHandle(), allocInfo);

        vk::CommandBufferAllocateInfo cmdBufferInfo(computeCmdPool->getHandle(), vk::CommandBufferLevel::ePrimary, 1);
        computeCmdBuffer = std::make_unique<MCommandBuffers>(logicalDevice->getHandle(), allocInfo);
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
            imageAvailableSemaphores.emplace_back(logicalDevice->getHandle(), semaphoreCreateInfo);
            renderFinishedSemaphores.emplace_back(logicalDevice->getHandle(), semaphoreCreateInfo);
            inFlightFences.emplace_back(logicalDevice->getHandle(), fenceInfo);
        }
    }

    // Load model
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::loadModel() {
        VulkanResources res{
            &*renderPass,
            &*physicalDevice,
            &*logicalDevice,
            &*instance,
            &*commandPool,
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
        logicalDevice->getHandle().waitForFences(inFlightFences[currentFrame].getHandle(), true, UINT64_MAX);

        vk::Result result;
        uint32_t imageIndex;

        std::tie(result, imageIndex) = logicalDevice->getHandle().acquireNextImageKHR(swapChain->getHandle(), UINT64_MAX, imageAvailableSemaphores[currentFrame].getHandle());

        if (result == vk::Result::eErrorOutOfDateKHR) {
            recreateSwapChain();
            return;
        } else if (result != vk::Result::eSuccess && result != vk::Result::eSuboptimalKHR) {
            throw std::runtime_error("Failed to acquire swap chain image!");
        }

        logicalDevice->getHandle().resetFences(inFlightFences[currentFrame].getHandle());

        commandBuffers->getBuffers().at(currentFrame).reset();
        recordCommandBuffer(imageIndex);

        vk::PipelineStageFlags flags(vk::PipelineStageFlagBits::eColorAttachmentOutput);
        vk::SubmitInfo submitInfo(
                (uint32_t) 1,
                imageAvailableSemaphores[currentFrame].getAddress(),
                &flags,
                (uint32_t) 1,
                &commandBuffers->getBuffers().at(currentFrame),
                (uint32_t) 1,
                renderFinishedSemaphores[currentFrame].getAddress()
        );
        std::array<vk::SubmitInfo, 1> infos = {submitInfo};

        graphicsQueue->submit(infos, inFlightFences[currentFrame].getHandle());

        std::array<vk::SwapchainKHR, 1> swapChains = {swapChain->getHandle()};
        vk::PresentInfoKHR presentInfo(*renderFinishedSemaphores[currentFrame].getAddress(), swapChains, imageIndex);

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
        commandBuffers->getBuffers().at(currentFrame).begin(beginInfo);

        std::array<float, 4> color = {0.0f, 0.0f, 0.0f, 1.0f};
        vk::ClearValue clrVal((vk::ClearColorValue(color)));

        vk::RenderPassBeginInfo renderPassInfo(
                renderPass->getHandle(),
                swapChainFramebuffers[imageIndex].getHandle(),
                {{0, 0}, swapChainExtent},
                1, &clrVal
        );

        commandBuffers->getBuffers().at(currentFrame).beginRenderPass(renderPassInfo, vk::SubpassContents::eInline);

        // Drawing the models
        model->draw(commandBuffers->getBuffers().at(currentFrame),currentFrame);
        //model->draw(commandBuffers->at(currentFrame));

        commandBuffers->getBuffers().at(currentFrame).endRenderPass();
        commandBuffers->getBuffers().at(currentFrame).end();
    }
}