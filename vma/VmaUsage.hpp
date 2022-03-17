//
// Created by PCF112021 on 3/17/2022.
//

#ifndef SPHSIMULATION_VMAUSAGE_HPP
#define SPHSIMULATION_VMAUSAGE_HPP

#ifdef _WIN32

#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#define VK_USE_PLATFORM_WIN32_KHR

#else  // #ifdef _WIN32

#include <vulkan/vulkan.h>

#endif  // #ifdef _WIN32

#include "vk_mem_alloc.h"

#endif //SPHSIMULATION_VMAUSAGE_HPP
