//
// Created by PCF112021 on 3/8/2022.
//

#ifndef SPHSIMULATION_VLTAVAFUNCTIONS_H
#define SPHSIMULATION_VLTAVAFUNCTIONS_H

#include <iostream>
#include <fstream>
#include "vulkan/vulkan.hpp"
#include "vulkan/vulkan_raii.hpp"

namespace Vltava {
    uint32_t findMemoryType(uint32_t typeFilter, const vk::PhysicalDeviceMemoryProperties &properties, vk::MemoryPropertyFlags flags);
    std::vector<char> readFile(const std::string &filename);
    void sanityCheck(const std::string &str);
}

#endif //SPHSIMULATION_VLTAVAFUNCTIONS_H
