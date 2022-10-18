#ifndef SPHSIMULATION_MSWAPCHAIN_HPP
#define SPHSIMULATION_MSWAPCHAIN_HPP

#include <GLFW/glfw3.h>
#include <vulkan/vulkan.hpp>
#include <memory>

namespace Vltava {
    class MSwapChain {
    public:
        MSwapChain(GLFWwindow* w,
                   vk::PhysicalDevice physDev,
                   vk::Device dev,
                   vk::SurfaceKHR surface,
                   vk::Format format,
                   vk::ColorSpaceKHR colorSpace,
                   vk::PresentModeKHR presentMode,
                   uint32_t graphicsQFamily,
                   uint32_t presentQFamily);
        MSwapChain();
        ~MSwapChain();

        void cleanup();
        vk::Extent2D getExtent();
        vk::Format getFormat();
        vk::SwapchainKHR getHandle();
    private:
        std::unique_ptr<vk::SwapchainKHR> swapChainHandle;
        vk::Format format;
        vk::Format imageFormat;
        vk::ColorSpaceKHR colorSpace;
        vk::PresentModeKHR presentMode;

        vk::Device logDev;

        vk::Extent2D extent;

        void createSwapchain(GLFWwindow* window,
                             vk::PhysicalDevice physDev,
                             vk::Device logDev,
                             vk::SurfaceKHR surface,
                             uint32_t graphicsQFamily, uint32_t presentQFamily);
    };
}

#endif //SPHSIMULATION_MSWAPCHAIN_HPP
