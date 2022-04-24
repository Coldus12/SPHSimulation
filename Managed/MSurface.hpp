#ifndef SPHSIMULATION_MSURFACE_HPP
#define SPHSIMULATION_MSURFACE_HPP

#define GLFW_INCLUDE_VULKAN
#include "GLFW/glfw3.h"
#include "vulkan/vulkan.hpp"

namespace Vltava {
    class MSurface {
    public:
        MSurface(VkInstance instance, GLFWwindow* w);
        ~MSurface();

        vk::SurfaceKHR getHandle();
    private:
        VkInstance instance;
        std::unique_ptr<vk::SurfaceKHR> surfaceHandle;
    };
}

#endif //SPHSIMULATION_MSURFACE_HPP
