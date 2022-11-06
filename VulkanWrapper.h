#ifndef SPHSIMULATION_VULKANWRAPPER_H
#define SPHSIMULATION_VULKANWRAPPER_H

#define GLFW_INCLUDE_NONE
#define GLFW_INCLUDE_VULKAN
#include "GLFW/glfw3.h"

#include <cstdint>
#include <string>
#include "vulkan/vulkan.hpp"

#include "Managed/Managed.hpp"
#include "VltavaFunctions.hpp"

namespace Vltava {

    class VulkanWrapper {
    public:
        // Functions
        //--------------------------------------------------------------------------------------------------------------
        static VulkanWrapper& getInstance() {
            static VulkanWrapper instance;
            return instance;
        }

        ~VulkanWrapper() { cleanup(); }

        virtual void createVulkanWindow(GLFWwindow* window);
        virtual void createWindowless();

        void cleanupSwapChain();
        virtual void recreateSwapChain(GLFWwindow* window);
        void cleanup();

        uint32_t getGraphicsQueueFamily() const { return graphicsQueueFamily; }
        uint32_t getPresentQueueFamily() const { return presentQueueFamily; }
        uint32_t getComputeQueueFamily() const { return computeQueueFamily; }

        vk::CommandBuffer& getCompCmdBuffer() { return computeCmdBuffer; }
        std::vector<vk::CommandBuffer>& getCmdBuffers() { return commandBuffers; }
        std::vector<vk::Framebuffer>& getSCFrameBuffers() { return swapChainFramebuffers; }

        // Variables
        //--------------------------------------------------------------------------------------------------------------
    private:
        enum State {
            uninitialized,
            windowless,
            vulkanWindow
        };

        State currentState = uninitialized;

        VulkanWrapper() = default;
        // Functions
        //--------------------------------------------------------------------------------------------------------------
        void createDepthResources();
        vk::ImageView createImageView(vk::Image image, vk::Format format, vk::ImageAspectFlags aspectFlags);
        void createImage(
                uint32_t width,
                uint32_t height,
                vk::Format format,
                vk::ImageTiling tiling,
                vk::ImageUsageFlags usage,
                vk::MemoryPropertyFlags properties,
                vk::Image& image,
                vk::DeviceMemory& imagememory
        );

        virtual void createInstance(std::string app_name);
        virtual void createSurface(GLFWwindow* window);
        virtual void selectPhysicalDevice();
        virtual void selectQueues(bool graphical = true);
        virtual void createLogicalDevice();
        virtual void createSwapChain(GLFWwindow* window);
        virtual void createImageViews();
        virtual void createRenderPass();
        virtual void createFramebuffers();
        virtual void createCommandPool();
        virtual void createCommandBuffers();

        // Variables
        //--------------------------------------------------------------------------------------------------------------
        bool enableValidationLayers = false;

        const std::vector<const char *> deviceExtensions = {
                VK_KHR_SWAPCHAIN_EXTENSION_NAME
        };

        const std::vector<const char *> validationLayers = {
                "VK_LAYER_KHRONOS_validation"
        };

        uint32_t graphicsQueueFamily = (uint32_t) -1;
        uint32_t presentQueueFamily = (uint32_t) -1;
        uint32_t computeQueueFamily = (uint32_t) -1;

        std::vector<vk::Image> swapChainImages;
        std::vector<vk::ImageView> imageViews;
        std::vector<vk::Framebuffer> swapChainFramebuffers;
        vk::Format swapChainImageFormat;
        vk::Image depthImage;
        vk::DeviceMemory depthImageMemory;
        vk::ImageView depthImageView;
        std::vector<vk::CommandBuffer> commandBuffers;
        vk::CommandBuffer computeCmdBuffer;
    };

}

#endif //SPHSIMULATION_VULKANWRAPPER_H
