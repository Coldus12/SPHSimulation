#include "StdWindow.hpp"
#include "VltavaFunctions.hpp"
#include "ComputeShader.hpp"

namespace Vltava {
    // Main loop
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::mainloop() {
        while (!glfwWindowShouldClose(window)) {
            glfwPollEvents();
            drawFrame();
        }

        device->waitIdle();
    }

    struct Stuff {
        glm::vec3 pos;
        float rad;
    };

    void StdWindow::computeStuff() {
        VulkanResources updated{
                &*renderPass,
                &*physicalDevice,
                &*device,
                &*instance,
                &*commandPool,
                &*graphicsQueue,
                &*computeQueue,
                swapChainExtent,
                MAX_FRAMES_IN_FLIGHT
        };

        vk::DeviceSize size = sizeof(Stuff) * 512;
        auto* data = new Stuff[512];
        for (int i = 0; i < 512; i++) {
            data[i].pos = glm::vec3(-i,i,-i);
            data[i].rad = sin(i);
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
        outBuffer.setSize(512 * sizeof(Stuff));
        outBuffer.bind(0);

        inBuffer.writeToBuffer(data, size);
        delete[] data;

        std::vector<Buffer> buffers;
        buffers.push_back(std::move(inBuffer));
        buffers.push_back(std::move(outBuffer));

        ComputeShader comp(updated, "shaders/comp.spv");
        comp.setStorageBuffers(buffers);
        comp.createPipeline();
        comp.createCommandBuffer(computeQueueFamily);
        comp.dispatch(size / sizeof(Stuff) + 1, 1, 1);

        auto spheres = buffers[1].getData<Stuff>();
        std::cout << spheres.size() << std::endl;

        for (int i = 0; i < 512; i++) {
            std::cout << spheres[i].rad << " ";
        }
        std::cout << std::endl;
        std::cout << "Done" << std::endl;
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

        initVulkan();
        computeStuff();
        mainloop();
    }

    void StdWindow::initVulkan() {
        createInstance("Test");
        setupDebugMessenger();
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

        loadModel();
    }

    // Destructor
    //------------------------------------------------------------------------------------------------------------------
    StdWindow::~StdWindow() {
        cleanup();
    }

    // Instance creation
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::createInstance(std::string app_name) {
        vk::raii::Context context;
        vk::ApplicationInfo appInfo(app_name.c_str(), VK_MAKE_VERSION(1, 0, 0), "No Engine", VK_MAKE_VERSION(1, 0, 0),
                                    VK_API_VERSION_1_0);
        vk::InstanceCreateInfo createInfo({}, &appInfo);

        auto glfwExtensions = getRequiredExtensions();

        if (enableValidationLayers) {
            if (!checkValidationLayerSupport())
                throw std::runtime_error("Validation layers requested, but not available!");

            createInfo.enabledLayerCount = static_cast<uint32_t>(validationLayers.size());
            createInfo.ppEnabledLayerNames = validationLayers.data();

            vk::DebugUtilsMessengerCreateInfoEXT debugCreateInfo(
                    {},
                    vk::DebugUtilsMessageSeverityFlagBitsEXT::eVerbose |
                    vk::DebugUtilsMessageSeverityFlagBitsEXT::eWarning |
                    vk::DebugUtilsMessageSeverityFlagBitsEXT::eError,
                    vk::DebugUtilsMessageTypeFlagBitsEXT::eGeneral |
                    vk::DebugUtilsMessageTypeFlagBitsEXT::eValidation |
                    vk::DebugUtilsMessageTypeFlagBitsEXT::ePerformance,
                    debugCallback
            );

            createInfo.pNext = (vk::DebugUtilsMessengerCreateInfoEXT *) &debugCreateInfo;
        } else {
            createInfo.enabledLayerCount = 0;
        }

        createInfo.enabledExtensionCount = static_cast<uint32_t>(glfwExtensions.size());
        createInfo.ppEnabledExtensionNames = glfwExtensions.data();

        createInfo.pNext = nullptr;

        instance = std::make_unique<vk::raii::Instance>(std::move(vk::raii::Instance(context, createInfo)));
    }

    // Validation layer
    //------------------------------------------------------------------------------------------------------------------
    bool StdWindow::checkValidationLayerSupport() {
        uint32_t layerCount;
        vk::enumerateInstanceLayerProperties(&layerCount, nullptr);

        std::vector<vk::LayerProperties> availableLayers(layerCount);
        vk::enumerateInstanceLayerProperties(&layerCount, availableLayers.data());

        for (auto layerName: validationLayers) {
            bool layerFound = false;

            for (const auto &layerProperties: validationLayers) {
                if (strcmp(layerProperties, layerName) == 0) {
                    layerFound = true;
                    break;
                }
            }

            if (!layerFound)
                return false;
        }

        return true;
    }

    void StdWindow::setupDebugMessenger() {
        if (!enableValidationLayers) return;

        vk::DebugUtilsMessengerCreateInfoEXT debugCreateInfo(
                {},
                vk::DebugUtilsMessageSeverityFlagBitsEXT::eVerbose |
                vk::DebugUtilsMessageSeverityFlagBitsEXT::eWarning |
                vk::DebugUtilsMessageSeverityFlagBitsEXT::eError,
                vk::DebugUtilsMessageTypeFlagBitsEXT::eGeneral |
                vk::DebugUtilsMessageTypeFlagBitsEXT::eValidation |
                vk::DebugUtilsMessageTypeFlagBitsEXT::ePerformance,
                debugCallback
        );

        pfnVkCreateDebugUtilsMessengerEXT = reinterpret_cast<PFN_vkCreateDebugUtilsMessengerEXT>( instance->getProcAddr(
                "vkCreateDebugUtilsMessengerEXT"));
        if (!pfnVkCreateDebugUtilsMessengerEXT) {
            std::cout << "GetInstanceProcAddr: Unable to find pfnVkCreateDebugUtilsMessengerEXT function." << std::endl;
            exit(1);
        }

        pfnVkDestroyDebugUtilsMessengerEXT = reinterpret_cast<PFN_vkDestroyDebugUtilsMessengerEXT>( instance->getProcAddr(
                "vkDestroyDebugUtilsMessengerEXT"));
        if (!pfnVkDestroyDebugUtilsMessengerEXT) {
            std::cout << "GetInstanceProcAddr: Unable to find pfnVkDestroyDebugUtilsMessengerEXT function."
                      << std::endl;
            exit(1);
        }

        debugUtilsMessenger = std::make_unique<vk::raii::DebugUtilsMessengerEXT>(*instance, debugCreateInfo);
    }

    VKAPI_ATTR VkBool32 VKAPI_CALL StdWindow::debugCallback(VkDebugUtilsMessageSeverityFlagBitsEXT severity,
                                                            VkDebugUtilsMessageTypeFlagsEXT type,
                                                            const VkDebugUtilsMessengerCallbackDataEXT *pCallbackData,
                                                            void *pUserData) {

        std::cerr << "[Validation layer] " << pCallbackData->pMessage << std::endl;
        return VK_FALSE;
    }

    VKAPI_ATTR VkResult VKAPI_CALL StdWindow::vkCreateDebugUtilsMessengerEXT(
            VkInstance instance,
            const VkDebugUtilsMessengerCreateInfoEXT *pCreateInfo,
            const VkAllocationCallbacks *pAllocator,
            VkDebugUtilsMessengerEXT *pDebugMessenger) {

        return pfnVkCreateDebugUtilsMessengerEXT(instance, pCreateInfo, pAllocator, pDebugMessenger);
    }

    VKAPI_ATTR void VKAPI_CALL StdWindow::vkDestroyDebugUtilsMessengerEXT(
            VkInstance instance,
            VkDebugUtilsMessengerEXT debugMessenger,
            const VkAllocationCallbacks *pAllocator) {

        return pfnVkDestroyDebugUtilsMessengerEXT(instance, debugMessenger, pAllocator);
    }

    std::vector<const char *> StdWindow::getRequiredExtensions() {
        uint32_t extensionCount = 0;
        const char **p_ext = glfwGetRequiredInstanceExtensions(&extensionCount);
        std::vector<const char *> extensions(p_ext, p_ext + extensionCount);

        if (enableValidationLayers)
            extensions.push_back(VK_EXT_DEBUG_UTILS_EXTENSION_NAME);

        return extensions;
    }

    // Selecting the best gpu for the task
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::selectPhysicalDevice() {
        vk::raii::PhysicalDevices devs(*instance);

        if (devs.empty())
            throw std::runtime_error("Failed to find GPUs with Vulkan support!");

        int bestNr = getBestDevice(devs);
        std::cout << "BEST NUMBER IS: " << bestNr << std::endl;
        physicalDevice = std::make_unique<vk::raii::PhysicalDevice>(std::move(devs.at(bestNr)));
    }

    int StdWindow::getBestDevice(const vk::raii::PhysicalDevices &devs) {
        int maxScore = 0;
        int currentBestNr = 0;

        for (int i = 0; i < devs.size(); i++) {
            int score = rateDeviceSuitability(devs.at(i));

            if (score >= maxScore) {
                maxScore = score;
                currentBestNr = i;
            }
        }

        std::cout << "[GPU Score] " << maxScore << std::endl;

        return currentBestNr;
    }

    uint32_t StdWindow::rateDeviceSuitability(const vk::raii::PhysicalDevice &dev) {
        uint32_t score = 0;

        vk::PhysicalDeviceProperties properties = dev.getProperties();
        vk::PhysicalDeviceFeatures features = dev.getFeatures();

        if (properties.deviceType == vk::PhysicalDeviceType::eDiscreteGpu)
            score += 1000;

        score += properties.limits.maxImageDimension2D;

        if (!features.geometryShader)
            return 0;

        if (!checkDeviceExtensionSupport(dev))
            score = 0;

        return score;
    }

    bool StdWindow::checkDeviceExtensionSupport(const vk::raii::PhysicalDevice &dev) {
        uint32_t extensionCount = 0;
        std::vector<vk::ExtensionProperties> availableExtensions = dev.enumerateDeviceExtensionProperties();

        std::set<std::string> requiredExtensions(deviceExtensions.begin(), deviceExtensions.end());

        for (const auto &extension: availableExtensions) {
            requiredExtensions.erase(extension.extensionName);
        }

        return requiredExtensions.empty();
    }

    // Selecting queue(s)
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::selectQueues() {
        std::vector<vk::QueueFamilyProperties> queueFamilies = physicalDevice->getQueueFamilyProperties();

        uint32_t i = 0;
        for (const auto &queueFamily: queueFamilies) {
            if (physicalDevice->getSurfaceSupportKHR(i, **surface)) {
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
        std::vector<vk::DeviceQueueCreateInfo> queueCreateInfos;
        float queuePriority = 1.0f;

        for (auto &queueFamily: uniqueQueueFamilies) {
            vk::DeviceQueueCreateInfo createInfo({}, queueFamily, 1, &queuePriority);
            queueCreateInfos.push_back(createInfo);
        }

        vk::PhysicalDeviceFeatures features;

        vk::DeviceCreateInfo createInfo(
                {},
                queueCreateInfos.size(),
                queueCreateInfos.data(),
                0,
                nullptr,
                static_cast<uint32_t>(deviceExtensions.size()),
                deviceExtensions.data(),
                &features
        );

        if (enableValidationLayers) {
            // Add validation layer stuff
            createInfo.enabledLayerCount = static_cast<uint32_t>(validationLayers.size());
            createInfo.ppEnabledLayerNames = validationLayers.data();
        }

        device = std::make_unique<vk::raii::Device>(*physicalDevice, createInfo);
        graphicsQueue = std::make_unique<vk::raii::Queue>(device->getQueue(graphicsQueueFamily, 0));
        presentQueue = std::make_unique<vk::raii::Queue>(device->getQueue(presentQueueFamily, 0));
        computeQueue = std::make_unique<vk::raii::Queue>(device->getQueue(computeQueueFamily, 0));
    }

    // Cleaning up
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::cleanup() {
        cleanupSwapChain();
        glfwDestroyWindow(window);
    }

    // Creating a surface
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::createSurface() {
        vk::SurfaceKHR surf(nullptr);

        surface = std::make_unique<vk::raii::SurfaceKHR>(*instance, surf);
        glfwCreateWindowSurface(**instance, window, nullptr, (VkSurfaceKHR *) &**surface);
    }

    // SwapChain
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::createSwapChain() {

        // Surface format
        std::vector<vk::SurfaceFormatKHR> formats = physicalDevice->getSurfaceFormatsKHR(**surface);
        vk::SurfaceFormatKHR &chosenFormat = formats[0];

        for (const auto &format: formats) {
            if (format.format == vk::Format::eB8G8R8A8Srgb &&
                format.colorSpace == vk::ColorSpaceKHR::eSrgbNonlinear) {
                chosenFormat = format;
                break;
            }
        }

        // Surface present mode
        std::vector<vk::PresentModeKHR> presentModes = physicalDevice->getSurfacePresentModesKHR(**surface);
        vk::PresentModeKHR &chosenMode = presentModes[0];

        for (const auto &mode: presentModes) {
            if (mode == vk::PresentModeKHR::eMailbox) {
                chosenMode = mode;
                break;
            }
        }

        // Extent
        vk::SurfaceCapabilitiesKHR capabilities = physicalDevice->getSurfaceCapabilitiesKHR(**surface);
        vk::Extent2D extent = capabilities.currentExtent;
        if (capabilities.currentExtent.width == std::numeric_limits<uint32_t>::max()) {
            int width, height;
            glfwGetFramebufferSize(window, &width, &height);

            extent.width = static_cast<uint32_t>(width);
            extent.height = static_cast<uint32_t>(height);

            extent.width = std::clamp(extent.width, capabilities.minImageExtent.width,
                                      capabilities.maxImageExtent.width);
            extent.height = std::clamp(extent.height, capabilities.minImageExtent.height,
                                       capabilities.maxImageExtent.height);
        }

        // Actual Swap Chain creation
        uint32_t imageCount = capabilities.minImageCount + 1;
        if (capabilities.maxImageCount > 0 && imageCount > capabilities.maxImageCount)
            imageCount = capabilities.maxImageCount;

        uint32_t queueFamilyIndices[] = {graphicsQueueFamily, presentQueueFamily};
        vk::SwapchainCreateInfoKHR createInfo(
                {},
                **surface,
                imageCount,
                chosenFormat.format,
                chosenFormat.colorSpace,
                extent,
                1,
                vk::ImageUsageFlagBits::eColorAttachment,
                vk::SharingMode::eExclusive,
                0,
                nullptr,
                capabilities.currentTransform,
                vk::CompositeAlphaFlagBitsKHR::eOpaque,
                chosenMode,
                true,
                nullptr
        );

        if (graphicsQueueFamily != presentQueueFamily) {
            createInfo.imageSharingMode = vk::SharingMode::eConcurrent;
            createInfo.queueFamilyIndexCount = 2;
            createInfo.pQueueFamilyIndices = queueFamilyIndices;
        }

        swapChain = std::make_unique<vk::raii::SwapchainKHR>(*device, createInfo);
        swapChainImages = swapChain->getImages();
        imageViews.reserve(swapChainImages.size());

        swapChainImageFormat = chosenFormat.format;
        swapChainExtent = extent;
    }

    void StdWindow::recreateSwapChain() {
        int w = 0, h = 0;
        glfwGetFramebufferSize(window, &w, &h);
        while (w == 0 || h == 0) {
            glfwGetFramebufferSize(window, &w, &h);
            glfwWaitEvents();
        }

        device->waitIdle();
        cleanupSwapChain();
        createSwapChain();
        createImageViews();
        createRenderPass();

        // Updating / recreating the pipelines of the models
        VulkanResources updated{
            &*renderPass,
            &*physicalDevice,
            &*device,
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
            imageViews.push_back(vk::raii::ImageView(*device, createInfo));
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
        vk::RenderPassCreateInfo createInfo({}, 1, &colorAttachment, 1, &subpass, 1, &dependency);

        renderPass = std::make_unique<vk::raii::RenderPass>(*device, createInfo);
    }

    // Framebuffers
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::createFramebuffers() {
        swapChainFramebuffers.reserve(imageViews.size());

        for (auto const &view: imageViews) {
            vk::FramebufferCreateInfo framebufferInfo(
                    {},
                    **renderPass,
                    1,
                    &*view,
                    swapChainExtent.width,
                    swapChainExtent.height, 1
            );

            swapChainFramebuffers.push_back(vk::raii::Framebuffer(*device, framebufferInfo));
        }
    }

    // Commandpool
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::createCommandPool() {
        vk::CommandPoolCreateInfo poolInfo(vk::CommandPoolCreateFlagBits::eResetCommandBuffer, graphicsQueueFamily);
        commandPool = std::make_unique<vk::raii::CommandPool>(*device, poolInfo);
    }

    // Command buffers
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::createCommandBuffers() {
        vk::CommandBufferAllocateInfo allocInfo(**commandPool, vk::CommandBufferLevel::ePrimary, MAX_FRAMES_IN_FLIGHT);
        commandBuffers = std::make_unique<vk::raii::CommandBuffers>(*device, allocInfo);
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
            imageAvailableSemaphores.push_back(vk::raii::Semaphore(*device, semaphoreCreateInfo));
            renderFinishedSemaphores.push_back(vk::raii::Semaphore(*device, semaphoreCreateInfo));
            inFlightFences.push_back(vk::raii::Fence(*device, fenceInfo));
        }
    }

    // Load model
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::loadModel() {
        VulkanResources res{
            &*renderPass,
            &*physicalDevice,
            &*device,
            &*instance,
            &*commandPool,
            &*graphicsQueue,
            &*computeQueue,
            swapChainExtent,
            MAX_FRAMES_IN_FLIGHT
        };

        model = std::make_unique<Model>(res);
    }

    // Draw frame
    //------------------------------------------------------------------------------------------------------------------
    void StdWindow::drawFrame() {
        device->waitForFences(*inFlightFences[currentFrame], true, UINT64_MAX);

        vk::Result result;
        uint32_t imageIndex;

        std::tie(result, imageIndex) = swapChain->acquireNextImage(UINT64_MAX, *imageAvailableSemaphores[currentFrame]);
        if (result == vk::Result::eErrorOutOfDateKHR) {
            recreateSwapChain();
            return;
        } else if (result != vk::Result::eSuccess && result != vk::Result::eSuboptimalKHR) {
            throw std::runtime_error("Failed to acquire swap chain image!");
        }

        device->resetFences(*inFlightFences[currentFrame]);

        commandBuffers->at(currentFrame).reset();
        recordCommandBuffer(imageIndex);

        vk::PipelineStageFlags flags(vk::PipelineStageFlagBits::eColorAttachmentOutput);
        vk::SubmitInfo submitInfo(
                (uint32_t) 1,
                &*imageAvailableSemaphores[currentFrame],
                &flags,
                (uint32_t) 1,
                &*commandBuffers->at(currentFrame),
                (uint32_t) 1,
                &*renderFinishedSemaphores[currentFrame]
        );
        std::array<vk::SubmitInfo, 1> infos = {submitInfo};

        graphicsQueue->submit(infos, *inFlightFences[currentFrame]);

        std::array<vk::SwapchainKHR, 1> swapChains = {**swapChain};
        vk::PresentInfoKHR presentInfo(*renderFinishedSemaphores[currentFrame], swapChains, imageIndex);

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
        commandBuffers->at(currentFrame).begin(beginInfo);

        std::array<float, 4> color = {0.0f, 0.0f, 0.0f, 1.0f};
        vk::ClearValue clrVal((vk::ClearColorValue(color)));

        vk::RenderPassBeginInfo renderPassInfo(
                **renderPass,
                *swapChainFramebuffers[imageIndex],
                {{0, 0}, swapChainExtent},
                1, &clrVal
        );

        commandBuffers->at(currentFrame).beginRenderPass(renderPassInfo, vk::SubpassContents::eInline);

        // Drawing the models
        model->draw(commandBuffers->at(currentFrame),currentFrame);
        //model->draw(commandBuffers->at(currentFrame));

        commandBuffers->at(currentFrame).endRenderPass();
        commandBuffers->at(currentFrame).end();
    }
}