#ifndef SPHSIMULATION_STDWINDOW_HPP
#define SPHSIMULATION_STDWINDOW_HPP

#define GLFW_INCLUDE_VULKAN
#include "GLFW/glfw3.h"

#include "vulkan/vulkan.hpp"
#include "vulkan/vulkan_raii.hpp"

#include "Model.hpp"
#include "ParticleModel.hpp"
#include "ComputeShader.hpp"

#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <fstream>

namespace Vltava {
    class StdWindow {

    public:
        // Functions
        //--------------------------------------------------------------------------------------------------------------
        explicit StdWindow(int width = 1280, int height = 720);
        ~StdWindow();
        static bool rot;

    private:
        // Variables
        //--------------------------------------------------------------------------------------------------------------
        int width;
        int height;
        int MAX_FRAMES_IN_FLIGHT = 2;
        bool framebufferResized = false;
        GLFWwindow *window;

        static bool run1;
        static bool run2;
        void runComp();

        std::unique_ptr<vk::raii::Instance> instance;
        std::unique_ptr<vk::raii::PhysicalDevice> physicalDevice;
        std::unique_ptr<vk::raii::Device> device;
        std::unique_ptr<vk::raii::DebugUtilsMessengerEXT> debugUtilsMessenger;
        std::unique_ptr<vk::raii::SwapchainKHR> swapChain;
        std::unique_ptr<vk::raii::SurfaceKHR> surface;
        std::unique_ptr<vk::raii::RenderPass> renderPass;
        std::unique_ptr<vk::raii::CommandPool> commandPool;
        std::unique_ptr<vk::raii::CommandBuffers> commandBuffers;
        std::unique_ptr<vk::raii::CommandPool> computeCmdPool;
        std::unique_ptr<vk::raii::CommandBuffer> computeCmdBuffer;

        std::unique_ptr<vk::raii::Queue> graphicsQueue;
        std::unique_ptr<vk::raii::Queue> presentQueue;
        std::unique_ptr<vk::raii::Queue> computeQueue;

        std::unique_ptr<Model> model;
        //std::unique_ptr<ParticleModel> model;
        std::vector<Buffer> sBuffers;
        std::vector<Buffer> uBuffers;

        std::unique_ptr<ComputeShader> comp1;
        std::unique_ptr<ComputeShader> comp2;

        uint32_t currentFrame = 0;
        uint32_t graphicsQueueFamily = (uint32_t) -1;
        uint32_t presentQueueFamily = (uint32_t) -1;
        uint32_t computeQueueFamily = (uint32_t) -1;

        vk::Format swapChainImageFormat;
        vk::Extent2D swapChainExtent;

        std::vector<VkImage> swapChainImages; // even the sample uses VkImage and not vk::Image (i guess because swapchain.getImages() returns a vector of VkImage and not vk::Image)
        std::vector<vk::raii::ImageView> imageViews;
        std::vector<vk::raii::Framebuffer> swapChainFramebuffers;
        std::vector<vk::raii::Semaphore> imageAvailableSemaphores;
        std::vector<vk::raii::Semaphore> renderFinishedSemaphores;
        std::vector<vk::raii::Fence> inFlightFences;

        PFN_vkCreateDebugUtilsMessengerEXT pfnVkCreateDebugUtilsMessengerEXT = VK_NULL_HANDLE;
        PFN_vkDestroyDebugUtilsMessengerEXT pfnVkDestroyDebugUtilsMessengerEXT = VK_NULL_HANDLE;

        bool enableValidationLayers = true;

        const std::vector<const char *> deviceExtensions = {
                VK_KHR_SWAPCHAIN_EXTENSION_NAME
        };

        const std::vector<const char *> validationLayers = {
                "VK_LAYER_KHRONOS_validation"
        };

        // Functions
        //--------------------------------------------------------------------------------------------------------------
        virtual void initVulkan();
        virtual void createInstance(std::string app_name);
        virtual void setupDebugMessenger();
        virtual void createSurface();
        virtual void selectPhysicalDevice();
        virtual void selectQueues();
        virtual void createLogicalDevice();
        virtual void recreateSwapChain();
        virtual void createSwapChain();
        virtual void createImageViews();
        virtual void createRenderPass();
        virtual void createFramebuffers();
        virtual void createCommandPool();
        virtual void createCommandBuffers();
        virtual void createSyncObjects();
        virtual void drawFrame();
        virtual void recordCommandBuffer(uint32_t imageIndex);
        void loadModel();
        void cleanup();

        void computeStuff();

        void mainloop();

        // Helper functions
        //--------------------------------------------------------------------------------------------------------------
        // GPU selection helper functions
        int getBestDevice(const vk::raii::PhysicalDevices &devs);
        uint32_t rateDeviceSuitability(const vk::raii::PhysicalDevice &dev);
        bool checkDeviceExtensionSupport(const vk::raii::PhysicalDevice &dev);

        // Validation layer helper function
        bool checkValidationLayerSupport();
        std::vector<const char *> getRequiredExtensions();

        VKAPI_ATTR VkResult VKAPI_CALL vkCreateDebugUtilsMessengerEXT(
                VkInstance instance,
                const VkDebugUtilsMessengerCreateInfoEXT *pCreateInfo,
                const VkAllocationCallbacks *pAllocator,
                VkDebugUtilsMessengerEXT *pDebugMessenger);

        VKAPI_ATTR void VKAPI_CALL vkDestroyDebugUtilsMessengerEXT(
                VkInstance instance,
                VkDebugUtilsMessengerEXT debugMessenger,
                const VkAllocationCallbacks *pAllocator);

        static VKAPI_ATTR VkBool32 VKAPI_CALL debugCallback(VkDebugUtilsMessageSeverityFlagBitsEXT severity,
                                                            VkDebugUtilsMessageTypeFlagsEXT type,
                                                            const VkDebugUtilsMessengerCallbackDataEXT *pCallbackData,
                                                            void *pUserData);

        // Swapchain helper options
        void cleanupSwapChain();
        static void frameBufferResizeCallback(GLFWwindow *window, int width, int height);

        static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
    };
}

#endif //SPHSIMULATION_STDWINDOW_HPP
