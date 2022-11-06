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
#include "CPUSim.hpp"
#include "VulkanWrapper.h"
#include "SESPH.h"
#include "IISPH.h"

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
        SimProps props;
        glm::vec3 gridA;
        glm::vec3 gridB;

        int width;
        int height;
        //int MAX_FRAMES_IN_FLIGHT = 2;

        bool framebufferResized = false;
        GLFWwindow *window;

        static bool run;
        static bool parallelCpuSim;
        static bool logB;
        bool cpuSim = false;

        int nrOfIter = 100;
        int particleNr = 512;
        int nrOfP = particleNr;
        int all_particle_nr = particleNr; // TODO rename all variations of "particle nr" to reflect what the difference is between them
        bool show_log = false;
        bool write_log = false;
        bool console_log = false;
        bool realTime = false;
        bool makeContainer = false;
        bool setBoundaries = false;
        int containerSize = 16;
        glm::vec3 containerPos = glm::vec3(0.0f,0.0f,-0.2f);
        std::chrono::time_point<std::chrono::high_resolution_clock> start;

        int cellx, celly, cellz;

        void runComp();
        std::vector<Particle> createContainerAt(glm::vec3 pos, float dist, int sideLength);
        void log();

        uint32_t currentFrame = 0;

        //std::unique_ptr<CPUSim> cpusim;
        std::unique_ptr<ParticleModel> model;

        // For the compute shaders.
        // These contain the particle data.
        std::vector<Buffer> sBuffers;
        std::vector<Buffer> uBuffers;

        std::unique_ptr<SESPH> sesph;
        std::unique_ptr<IISPH> iisph;

        //IMGUI
        VkDescriptorPool imguiPool;

        std::vector<vk::Semaphore> imageAvailableSemaphores;
        std::vector<vk::Semaphore> renderFinishedSemaphores;
        std::vector<vk::Fence> inFlightFences;

        // Functions
        //--------------------------------------------------------------------------------------------------------------
        void initImgui();
        virtual void initVulkan();
        virtual void createSyncObjects();
        virtual void drawFrame();
        virtual void recordCommandBuffer(uint32_t imageIndex);
        void loadModel();
        void cleanup();
        
        void setComputeData();
        void resetData();

        void mainloop();

        void runCpuSim(int iterNr);

        static void frameBufferResizeCallback(GLFWwindow *window, int width, int height);
        static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
    };
}

#endif //SPHSIMULATION_STDWINDOW_HPP
