#include "VulkanWrapper.h"

namespace Vltava {
    void VulkanWrapper::createVulkanWindow(GLFWwindow* window) {
        createInstance("Test");
        createSurface(window);
        selectPhysicalDevice();
        selectQueues();
        createLogicalDevice();
        createSwapChain(window);
        createImageViews();
        createRenderPass();
        createCommandPool();
        createDepthResources();
        createFramebuffers();
        createCommandBuffers();
    }

    // Instance creation
    //------------------------------------------------------------------------------------------------------------------
    void VulkanWrapper::createInstance(std::string app_name) {
        VulkanResources::getInstance().instance = std::make_unique<MInstance>(
                app_name,
                enableValidationLayers
        );
    }

    // Creating a surface
    //------------------------------------------------------------------------------------------------------------------
    void VulkanWrapper::createSurface(GLFWwindow* window) {
        VulkanResources::getInstance().surface = std::make_unique<MSurface>(
                VulkanResources::getInstance().instance->getHandle(),
                window
        );
    }

    // Selecting the best gpu for the task
    //------------------------------------------------------------------------------------------------------------------
    void VulkanWrapper::selectPhysicalDevice() {
        VulkanResources::getInstance().physDev = std::make_unique<MPhysDev>(
                VulkanResources::getInstance().instance->getHandle(),
                deviceExtensions
        );
    }

    // Selecting queue(s)
    //------------------------------------------------------------------------------------------------------------------
    void VulkanWrapper::selectQueues() {
        std::vector<vk::QueueFamilyProperties> queueFamilies = VulkanResources::getInstance().physDev->getHandle().getQueueFamilyProperties();

        uint32_t i = 0;
        for (const auto &queueFamily: queueFamilies) {
            if (VulkanResources::getInstance().physDev->getHandle().getSurfaceSupportKHR(i, VulkanResources::getInstance().surface->getHandle())) {
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
    void VulkanWrapper::createLogicalDevice() {
        std::set<uint32_t> uniqueQueueFamilies({graphicsQueueFamily, presentQueueFamily, computeQueueFamily});

        VulkanResources::getInstance().logDev = std::make_unique<MLogDev>(
                uniqueQueueFamilies,
                VulkanResources::getInstance().physDev->getHandle(),
                deviceExtensions,
                &validationLayers
        );

        // Queue creation
        //--------------------------------------------------------------------------------------------------------------
        //Clearing previous queues
        VulkanResources::getInstance().graphicsQueue.reset();
        VulkanResources::getInstance().presentQueue.reset();
        VulkanResources::getInstance().computeQueue.reset();

        VulkanResources::getInstance().graphicsQueue = std::make_unique<vk::Queue>(VulkanResources::getInstance().logDev->getHandle().getQueue(graphicsQueueFamily, 0));
        VulkanResources::getInstance().presentQueue = std::make_unique<vk::Queue>(VulkanResources::getInstance().logDev->getHandle().getQueue(presentQueueFamily, 0));
        VulkanResources::getInstance().computeQueue = std::make_unique<vk::Queue>(VulkanResources::getInstance().logDev->getHandle().getQueue(computeQueueFamily, 0));
    }

    // SwapChain
    //------------------------------------------------------------------------------------------------------------------
    void VulkanWrapper::createSwapChain(GLFWwindow* window) {
        VulkanResources::getInstance().swapChain = std::make_unique<MSwapChain>(
                window,
                VulkanResources::getInstance().physDev->getHandle(),
                VulkanResources::getInstance().logDev->getHandle(),
                VulkanResources::getInstance().surface->getHandle(),
                vk::Format::eB8G8R8A8Srgb,
                vk::ColorSpaceKHR::eSrgbNonlinear,
                vk::PresentModeKHR::eMailbox,
                graphicsQueueFamily,
                presentQueueFamily
        );

        swapChainImages = VulkanResources::getInstance().logDev->getHandle().getSwapchainImagesKHR(VulkanResources::getInstance().swapChain->getHandle());
        imageViews.reserve(swapChainImages.size());

        swapChainImageFormat = VulkanResources::getInstance().swapChain->getFormat();
        VulkanResources::getInstance().extent = VulkanResources::getInstance().swapChain->getExtent();
    }

    void VulkanWrapper::recreateSwapChain(GLFWwindow* window) {
        int w = 0, h = 0;
        glfwGetFramebufferSize(window, &w, &h);
        while (w == 0 || h == 0) {
            glfwGetFramebufferSize(window, &w, &h);
            glfwWaitEvents();
        }

        VulkanResources::getInstance().logDev->getHandle().waitIdle();
        cleanupSwapChain();
        createSwapChain(window);
        createImageViews();
        createRenderPass();
        createDepthResources();
        createFramebuffers();
    }

    void VulkanWrapper::cleanupSwapChain() {
        VulkanResources::getInstance().logDev->getHandle().destroyImageView(depthImageView);
        VulkanResources::getInstance().logDev->getHandle().destroyImage(depthImage);
        VulkanResources::getInstance().logDev->getHandle().freeMemory(depthImageMemory);

        for (auto framebuffer : swapChainFramebuffers) {
            VulkanResources::getInstance().logDev->getHandle().destroyFramebuffer(framebuffer);
        }
        swapChainFramebuffers.clear();

        VulkanResources::getInstance().logDev->getHandle().destroyRenderPass(*VulkanResources::getInstance().renderPass);

        for (auto imgView: imageViews) {
            VulkanResources::getInstance().logDev->getHandle().destroyImageView(imgView);
        }
        imageViews.clear();

        VulkanResources::getInstance().swapChain.reset();
    }

    // ImageView
    //------------------------------------------------------------------------------------------------------------------
    void VulkanWrapper::createImageViews() {
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
            imageViews.push_back(VulkanResources::getInstance().logDev->getHandle().createImageView(createInfo));
        }
    }

    // RenderPass
    //------------------------------------------------------------------------------------------------------------------
    void VulkanWrapper::createRenderPass() {
        vk::SubpassDependency dependency(
                VK_SUBPASS_EXTERNAL,
                0,
                vk::PipelineStageFlagBits::eColorAttachmentOutput | vk::PipelineStageFlagBits::eEarlyFragmentTests,
                vk::PipelineStageFlagBits::eColorAttachmentOutput | vk::PipelineStageFlagBits::eEarlyFragmentTests,
                {},
                vk::AccessFlagBits::eColorAttachmentWrite | vk::AccessFlagBits::eDepthStencilAttachmentWrite
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

        vk::AttachmentDescription depthAttachment(
                {},
                vk::Format::eD32Sfloat,
                vk::SampleCountFlagBits::e1,
                vk::AttachmentLoadOp::eClear,
                vk::AttachmentStoreOp::eDontCare,
                vk::AttachmentLoadOp::eDontCare,
                vk::AttachmentStoreOp::eDontCare,
                vk::ImageLayout::eUndefined,
                vk::ImageLayout::eDepthStencilAttachmentOptimal
        );

        vk::AttachmentReference depthReference(1, vk::ImageLayout::eDepthStencilAttachmentOptimal);

        vk::SubpassDescription subpass({}, vk::PipelineBindPoint::eGraphics, {}, colorReference, {}, &depthReference);
        std::array<vk::AttachmentDescription, 2> attachments = {colorAttachment, depthAttachment};
        vk::RenderPassCreateInfo createInfo(
                {},
                static_cast<uint32_t>(attachments.size()),
                attachments.data(),
                1,
                &subpass,
                1,
                &dependency
        );

        VulkanResources::getInstance().renderPass = std::make_unique<vk::RenderPass>(
                VulkanResources::getInstance().logDev->getHandle().createRenderPass(createInfo)
        );
    }

    // Framebuffers
    //------------------------------------------------------------------------------------------------------------------
    void VulkanWrapper::createFramebuffers() {
        swapChainFramebuffers.reserve(imageViews.size());

        for (int i = 0; i < imageViews.size(); i++) {
            std::array<vk::ImageView, 2> attachemnts = {
                    imageViews[i],
                    depthImageView
            };

            vk::FramebufferCreateInfo framebufferInfo(
                    {},
                    *VulkanResources::getInstance().renderPass,
                    static_cast<uint32_t>(attachemnts.size()),
                    attachemnts.data(),
                    VulkanResources::getInstance().extent.width,
                    VulkanResources::getInstance().extent.height,
                    1
            );

            swapChainFramebuffers.push_back(VulkanResources::getInstance().logDev->getHandle().createFramebuffer(framebufferInfo));
        }
    }

    // Commandpool
    //------------------------------------------------------------------------------------------------------------------
    void VulkanWrapper::createCommandPool() {
        vk::CommandPoolCreateInfo poolInfo(vk::CommandPoolCreateFlagBits::eResetCommandBuffer, graphicsQueueFamily);
        VulkanResources::getInstance().graphicalCmdPool = std::make_unique<vk::CommandPool>(
                VulkanResources::getInstance().logDev->getHandle().createCommandPool(poolInfo)
        );

        vk::CommandPoolCreateInfo cmdPoolInfo(vk::CommandPoolCreateFlagBits::eResetCommandBuffer, computeQueueFamily);
        VulkanResources::getInstance().computeCmdPool = std::make_unique<vk::CommandPool>(
                VulkanResources::getInstance().logDev->getHandle().createCommandPool(cmdPoolInfo)
        );
    }

    // Command buffers
    //------------------------------------------------------------------------------------------------------------------
    void VulkanWrapper::createCommandBuffers() {
        vk::CommandBufferAllocateInfo allocInfo(*VulkanResources::getInstance().graphicalCmdPool, vk::CommandBufferLevel::ePrimary, VulkanResources::getInstance().FRAMES_IN_FLIGHT);
        commandBuffers = VulkanResources::getInstance().logDev->getHandle().allocateCommandBuffers(allocInfo);

        vk::CommandBufferAllocateInfo cmdBufferInfo(*VulkanResources::getInstance().computeCmdPool, vk::CommandBufferLevel::ePrimary, 1);
        computeCmdBuffer = VulkanResources::getInstance().logDev->getHandle().allocateCommandBuffers(cmdBufferInfo).front();
    }

    // DepthBuffering and related functions
    //------------------------------------------------------------------------------------------------------------------
    void VulkanWrapper::createImage(uint32_t width,
                                uint32_t height,
                                vk::Format format,
                                vk::ImageTiling tiling,
                                vk::ImageUsageFlags usage,
                                vk::MemoryPropertyFlags properties,
                                vk::Image& image,
                                vk::DeviceMemory& imageMemory) {

        vk::ImageCreateInfo imageInfo(
                {},
                vk::ImageType::e2D,
                format,
                vk::Extent3D(width, height, 1),
                1,
                1,
                vk::SampleCountFlagBits::e1,
                tiling,
                usage,
                vk::SharingMode::eExclusive,
                {},
                vk::ImageLayout::eUndefined
        );

        image = VulkanResources::getInstance().logDev->getHandle().createImage(imageInfo);

        vk::MemoryRequirements memReq = VulkanResources::getInstance().logDev->getHandle().getImageMemoryRequirements(image);
        vk::MemoryAllocateInfo allocInfo(memReq.size, findMemoryType(memReq.memoryTypeBits, properties));
        imageMemory = VulkanResources::getInstance().logDev->getHandle().allocateMemory(allocInfo);
        VulkanResources::getInstance().logDev->getHandle().bindImageMemory(image, imageMemory, 0);
    }

    vk::ImageView VulkanWrapper::createImageView(vk::Image image, vk::Format format, vk::ImageAspectFlags aspectFlags) {
        vk::ImageViewCreateInfo viewInfo(
                {},
                image,
                vk::ImageViewType::e2D,
                format,
                {},
                vk::ImageSubresourceRange(aspectFlags,0,1,0,1)
        );

        vk::ImageView ret = VulkanResources::getInstance().logDev->getHandle().createImageView(viewInfo);
        return ret;
    }

    void VulkanWrapper::createDepthResources() {
        vk::Format depthFormat = vk::Format::eD32Sfloat;

        createImage(
                VulkanResources::getInstance().extent.width,
                VulkanResources::getInstance().extent.height,
                depthFormat,
                vk::ImageTiling::eOptimal,
                vk::ImageUsageFlagBits::eDepthStencilAttachment,
                vk::MemoryPropertyFlagBits::eDeviceLocal,
                depthImage,
                depthImageMemory
        );

        depthImageView = createImageView(depthImage, depthFormat, vk::ImageAspectFlagBits::eDepth);
    }

    // Cleaning up
    //------------------------------------------------------------------------------------------------------------------
    void VulkanWrapper::cleanup() {
        cleanupSwapChain();
        VulkanResources::getInstance().logDev->getHandle().destroyCommandPool(*VulkanResources::getInstance().graphicalCmdPool);
        VulkanResources::getInstance().logDev->getHandle().destroyCommandPool(*VulkanResources::getInstance().computeCmdPool);
    }
}