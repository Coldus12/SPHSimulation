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
#include <chrono>

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
        int particleNr = 512;
        int nrOfP = particleNr;
        int all_particle_nr = particleNr; // TODO rename all variations of "particle nr" to reflect what the difference is between them
        bool show_log = false;
        bool write_log = false;
        bool console_log = false;
        bool realTime = false;
        bool makeContainer = false;
        int containerSize = 16;
        glm::vec3 containerPos = glm::vec3(0.0f,0.0f,-0.2f);
        std::chrono::time_point<std::chrono::high_resolution_clock> start;

        void runComp();
        std::vector<Particle> createContainerAt(glm::vec3 pos, float dist, int sideLength);
        void log();

        vk::Image depthImage;
        vk::DeviceMemory depthImageMemory;
        vk::ImageView depthImageView;
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

        std::unique_ptr<ParticleModel> model;
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

        vk::Fence compFence;

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

        void initCompute();
        void dispatchCompute(int groupCountX, int groupCountY, int groupCountZ);
        void setComputeData();
        void resetData();

        void mainloop();

        // Swapchain helper options
        void cleanupSwapChain();
        static void frameBufferResizeCallback(GLFWwindow *window, int width, int height);

        static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
    };
}

#endif //SPHSIMULATION_STDWINDOW_HPP
