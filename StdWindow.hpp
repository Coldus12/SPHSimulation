#ifndef SPHSIMULATION_STDWINDOW_HPP
#define SPHSIMULATION_STDWINDOW_HPP

#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_vulkan.h"

#define GLFW_INCLUDE_NONE
#define GLFW_INCLUDE_VULKAN
#include "GLFW/glfw3.h"

#include "vulkan/vulkan.hpp"

#include "Model.hpp"
#include "ParticleModel.hpp"
#include "ComputeShader.hpp"
#include "Managed/Managed.hpp"

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
        //int MAX_FRAMES_IN_FLIGHT = 2;

        bool framebufferResized = false;
        GLFWwindow *window;

        static bool run;
        int nrOfIter = 100;
        void runComp();

        uint32_t currentFrame = 0;
        uint32_t graphicsQueueFamily = (uint32_t) -1;
        uint32_t presentQueueFamily = (uint32_t) -1;
        uint32_t computeQueueFamily = (uint32_t) -1;

        std::vector<vk::CommandBuffer> commandBuffers;
        vk::CommandBuffer computeCmdBuffer;

        //vk::Extent2D swapChainExtent;
        std::vector<vk::Image> swapChainImages;
        std::vector<vk::ImageView> imageViews;
        std::vector<vk::Framebuffer> swapChainFramebuffers;

        std::unique_ptr<Model> model;
        // For the compute shaders.
        // These contain the particle data.
        std::vector<Buffer> sBuffers;
        std::vector<Buffer> uBuffers;

        std::unique_ptr<ComputeShader> comp1;
        std::unique_ptr<ComputeShader> comp2;

        vk::Format swapChainImageFormat;

        // even the sample uses VkImage and not vk::Image (i guess because swapchain.getImages() returns a vector of VkImage and not vk::Image)
        std::vector<vk::Semaphore> imageAvailableSemaphores;
        std::vector<vk::Semaphore> renderFinishedSemaphores;
        std::vector<vk::Fence> inFlightFences;

        bool enableValidationLayers = false;

        const std::vector<const char *> deviceExtensions = {
                VK_KHR_SWAPCHAIN_EXTENSION_NAME
        };

        const std::vector<const char *> validationLayers = {
                "VK_LAYER_KHRONOS_validation"
        };

        //IMGUI
        VkDescriptorPool imguiPool;

        // Functions
        //--------------------------------------------------------------------------------------------------------------
        void initImgui();

        virtual void initVulkan();
        virtual void createInstance(std::string app_name);
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
        void dispatchCompute(int groupCountX, int groupCountY, int groupCountZ);
        void resetData();

        void mainloop();

        // Swapchain helper options
        void cleanupSwapChain();
        static void frameBufferResizeCallback(GLFWwindow *window, int width, int height);

        static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
    };
}

#endif //SPHSIMULATION_STDWINDOW_HPP
