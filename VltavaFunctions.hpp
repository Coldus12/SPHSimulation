#ifndef SPHSIMULATION_VLTAVAFUNCTIONS_HPP
#define SPHSIMULATION_VLTAVAFUNCTIONS_HPP

#include <iostream>
#include <fstream>
#include "vulkan/vulkan.hpp"
#include "Managed/Managed.hpp"

namespace Vltava {
    /*
     * VulkanResources::getInstance().renderPass = renderPass.get();
        VulkanResources::getInstance().physDev = physicalDevice.get();
        VulkanResources::getInstance().dev = logicalDevice.get();
        VulkanResources::getInstance().instance = instance.get();
        VulkanResources::getInstance().commandPool = commandPool.get();
        VulkanResources::getInstance().graphicsQueue = graphicsQueue.get();
        VulkanResources::getInstance().computeQueue = computeQueue.get();
        VulkanResources::getInstance().extent = swapChainExtent;
        VulkanResources::getInstance().FRAMES_IN_FLIGHT = MAX_FRAMES_IN_FLIGHT;
     *
     * */

    /*class VulkanResources {
    public:
        static VulkanResources& getInstance() {
            static VulkanResources instance;
            return instance;
        }

        VulkanResources(VulkanResources const&) = delete;
        void operator=(VulkanResources const&) = delete;

        MRenderPass* renderPass = nullptr;
        MPhysDev* physDev = nullptr;
        MLogDev* dev = nullptr;
        MInstance* instance = nullptr;
        MCommandPool* commandPool = nullptr;
        vk::Queue* graphicsQueue = nullptr;
        vk::Queue* computeQueue = nullptr;

        vk::Extent2D extent;
        int FRAMES_IN_FLIGHT = 0;

    private:
        VulkanResources() = default;
    };*/

    struct VulkanResources {
        vk::RenderPass* renderPass = nullptr;
        MPhysDev* physDev = nullptr;
        MLogDev* dev = nullptr;
        MInstance* instance = nullptr;
        vk::CommandPool* commandPool = nullptr;
        vk::Queue* graphicsQueue = nullptr;
        vk::Queue* computeQueue = nullptr;

        vk::Extent2D extent;
        int FRAMES_IN_FLIGHT = 0;
    };

    uint32_t findMemoryType(uint32_t typeFilter, const vk::PhysicalDeviceMemoryProperties &properties, vk::MemoryPropertyFlags flags);
    std::vector<char> readFile(const std::string &filename);
    void sanityCheck(const std::string &str);
}

#endif //SPHSIMULATION_VLTAVAFUNCTIONS_HPP
