#ifndef SPHSIMULATION_VLTAVAFUNCTIONS_HPP
#define SPHSIMULATION_VLTAVAFUNCTIONS_HPP

#include <iostream>
#include <fstream>
#include "vulkan/vulkan.hpp"
#include "vulkan/vulkan_raii.hpp"

namespace Vltava {
    struct VulkanResources {
        vk::raii::RenderPass* renderPass;
        vk::raii::PhysicalDevice* physDev;
        vk::raii::Device* dev;
        vk::raii::Instance* instance;
        vk::raii::CommandPool* commandPool;
        vk::raii::Queue* graphicsQueue;

        vk::Extent2D extent;
        int FRAMES_IN_FLIGHT = 0;
    };

    uint32_t findMemoryType(uint32_t typeFilter, const vk::PhysicalDeviceMemoryProperties &properties, vk::MemoryPropertyFlags flags);
    std::vector<char> readFile(const std::string &filename);
    void sanityCheck(const std::string &str);
}

#endif //SPHSIMULATION_VLTAVAFUNCTIONS_HPP
