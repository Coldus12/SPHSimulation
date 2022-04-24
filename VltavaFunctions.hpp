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
        VulkanResources::getInstance().graphicalCmdPool = graphicalCmdPool.get();
        VulkanResources::getInstance().graphicsQueue = graphicsQueue.get();
        VulkanResources::getInstance().computeQueue = computeQueue.get();
        VulkanResources::getInstance().extent = swapChainExtent;
        VulkanResources::getInstance().FRAMES_IN_FLIGHT = MAX_FRAMES_IN_FLIGHT;
     *
     * */

    class VulkanResources {
    public:
        static VulkanResources& getInstance() {
            static VulkanResources instance;
            return instance;
        }

        VulkanResources(VulkanResources const&) = delete;
        void operator=(VulkanResources const&) = delete;

        std::unique_ptr<MInstance> instance;
        std::unique_ptr<MSurface> surface;
        std::unique_ptr<MPhysDev> physDev;
        std::unique_ptr<MLogDev> logDev;
        std::unique_ptr<MSwapChain> swapChain;
        std::unique_ptr<vk::RenderPass> renderPass;
        std::unique_ptr<vk::CommandPool> graphicalCmdPool;
        std::unique_ptr<vk::CommandPool> computeCmdPool;

        std::unique_ptr<vk::Queue> graphicsQueue;
        std::unique_ptr<vk::Queue> presentQueue;
        std::unique_ptr<vk::Queue> computeQueue;

        vk::Extent2D extent;
        int FRAMES_IN_FLIGHT = 0;
    private:
        VulkanResources() = default;
    };

    uint32_t findMemoryType(uint32_t typeFilter, const vk::PhysicalDeviceMemoryProperties &properties, vk::MemoryPropertyFlags flags);
    std::vector<char> readFile(const std::string &filename);
    void sanityCheck(const std::string &str);
}

#endif //SPHSIMULATION_VLTAVAFUNCTIONS_HPP
