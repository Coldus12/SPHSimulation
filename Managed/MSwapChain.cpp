#include <limits>
#include "MSwapChain.hpp"

namespace Vltava {
    MSwapChain::MSwapChain(GLFWwindow* w, vk::PhysicalDevice physDev, vk::Device logDev, vk::SurfaceKHR surface,
                           vk::Format format,vk::ColorSpaceKHR colorSpace, vk::PresentModeKHR presentMode,
                           uint32_t graphicsQFamily, uint32_t presentQFamily)
                           : logDev(logDev), format(format), colorSpace(colorSpace), presentMode(presentMode) {
        createSwapchain(w, physDev, logDev, surface, graphicsQFamily, presentQFamily);
    }

    MSwapChain::~MSwapChain() {
        cleanup();
    }

    vk::SwapchainKHR MSwapChain::getHandle() {
        return *swapChainHandle;
    }

    void MSwapChain::createSwapchain(GLFWwindow* window,
                                     vk::PhysicalDevice physDev,
                                     vk::Device logDev,
                                     vk::SurfaceKHR surface,
                                     uint32_t graphicsQFamily, uint32_t presentQFamily) {

        // Surface format
        std::vector<vk::SurfaceFormatKHR> formats = physDev.getSurfaceFormatsKHR(surface);
        vk::SurfaceFormatKHR &chosenFormat = formats[0];

        for (const auto &format: formats) {
            if (format.format == this->format &&
                format.colorSpace == this->colorSpace) {
                chosenFormat = format;
                break;
            }
        }

        // Surface present mode
        std::vector<vk::PresentModeKHR> presentModes = physDev.getSurfacePresentModesKHR(surface);
        vk::PresentModeKHR &chosenMode = presentModes[0];

        for (const auto &mode: presentModes) {
            if (mode == this->presentMode) {
                chosenMode = mode;
                break;
            }
        }

        // Extent
        vk::SurfaceCapabilitiesKHR capabilities = physDev.getSurfaceCapabilitiesKHR(surface);
        extent = capabilities.currentExtent;
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

        uint32_t queueFamilyIndices[] = {graphicsQFamily, presentQFamily};
        vk::SwapchainCreateInfoKHR createInfo(
                {},
                surface,
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

        if (graphicsQFamily != presentQFamily) {
            createInfo.imageSharingMode = vk::SharingMode::eConcurrent;
            createInfo.queueFamilyIndexCount = 2;
            createInfo.pQueueFamilyIndices = queueFamilyIndices;
        }

        imageFormat = chosenFormat.format;
        swapChainHandle = std::make_unique<vk::SwapchainKHR>(logDev.createSwapchainKHR(createInfo));
    }

    vk::Extent2D MSwapChain::getExtent() {
        return extent;
    }

    void MSwapChain::cleanup() {
        logDev.destroySwapchainKHR(*swapChainHandle);
        swapChainHandle.reset();
    }

    vk::Format MSwapChain::getFormat() {
        return imageFormat;
    }
}