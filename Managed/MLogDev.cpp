#include "MLogDev.hpp"

namespace Vltava {
    MLogDev::MLogDev(std::set<uint32_t> uniqueQueueFamilies,
                     vk::PhysicalDevice physDev,
                     const std::vector<const char *> deviceExtensions,
                     const std::vector<const char*>* validationLayers) {

        std::vector<vk::DeviceQueueCreateInfo> queueCreateInfos;
        float queuePriority = 1.0f;

        for (auto &queueFamily: uniqueQueueFamilies) {
            vk::DeviceQueueCreateInfo createInfo({}, queueFamily, 1, &queuePriority);
            queueCreateInfos.push_back(createInfo);
        }

        vk::PhysicalDeviceFeatures features;

        vk::DeviceCreateInfo createInfo(
                {},
                queueCreateInfos.size(),
                queueCreateInfos.data(),
                0,
                nullptr,
                static_cast<uint32_t>(deviceExtensions.size()),
                deviceExtensions.data(),
                &features
        );

        if (validationLayers != nullptr) {
            // Add validation layer stuff
            createInfo.enabledLayerCount = static_cast<uint32_t>(validationLayers->size());
            createInfo.ppEnabledLayerNames = validationLayers->data();
        }

        device = std::make_unique<vk::Device>(physDev.createDevice(createInfo));
    }

    MLogDev::~MLogDev() {
        device->destroy();
    }

    vk::Device MLogDev::getHandle() {
        return *device;
    }
}