#include "MSurface.hpp"

namespace Vltava {
    MSurface::MSurface(VkInstance instance, GLFWwindow *w) : instance(instance) {
        VkSurfaceKHR surf = {};

        if (glfwCreateWindowSurface(instance, w, nullptr, &surf) != VK_SUCCESS) {
            throw std::runtime_error("[ERROR] Failed to create window surface!\n[ERROR] Thrown from MSurface.hpp;");
        }

        surfaceHandle = std::make_unique<vk::SurfaceKHR>(surf);
    }

    MSurface::~MSurface() {
        vkDestroySurfaceKHR(instance, (VkSurfaceKHR) *surfaceHandle, nullptr);
    }

    vk::SurfaceKHR MSurface::getHandle() {
        return *surfaceHandle;
    }
}
