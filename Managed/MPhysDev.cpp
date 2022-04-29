#include <set>
#include "MPhysDev.hpp"
#include "MInstance.hpp"
#include <iostream>

namespace Vltava {
    MPhysDev::MPhysDev(vk::Instance instance, const std::vector<const char *> deviceExtensions)
    : instance(instance), deviceExtensions(deviceExtensions) {
        selectPhysicalDevice();
    }

    MPhysDev::~MPhysDev() {
        // No need to manually destroy / free anything.
    }

    void MPhysDev::selectPhysicalDevice() {
        auto devices = instance.enumeratePhysicalDevices();

        if (devices.empty())
            throw std::runtime_error("[ERROR] Failed to find GPUs with Vulkan support!\n [ERROR] Thrown from MPhysDev.");

        physDev = std::make_unique<vk::PhysicalDevice>(getBestDevice(devices));

        std::cout << "[Device name] " << physDev->getProperties().deviceName << std::endl;
    }

    vk::PhysicalDevice MPhysDev::getBestDevice(const std::vector<vk::PhysicalDevice> &devs) {
        int maxScore = 0;
        vk::PhysicalDevice currentBest;

        for (const auto& device: devs) {
            int score = rateDeviceSuitability(device);

            if (score >= maxScore) {
                maxScore = score;
                currentBest = device;
            }
        }

        if (maxScore == 0)
            throw std::runtime_error("[ERROR] No suitable GPU found!\n[ERROR] Thrown from MPhysDev.");

        return currentBest;
    }

    uint32_t MPhysDev::rateDeviceSuitability(const vk::PhysicalDevice &dev) {
        uint32_t score = 0;

        vk::PhysicalDeviceProperties properties = dev.getProperties();
        vk::PhysicalDeviceFeatures features = dev.getFeatures();

        if (properties.deviceType == vk::PhysicalDeviceType::eDiscreteGpu)
            score += 1000;

        score += properties.limits.maxImageDimension2D;

        if (!features.geometryShader)
            return 0;

        if (!checkDeviceExtensionSupport(dev))
            score = 0;

        return score;
    }

    bool MPhysDev::checkDeviceExtensionSupport(const vk::PhysicalDevice &dev) {
        uint32_t extensionCount = 0;
        std::vector<vk::ExtensionProperties> availableExtensions = dev.enumerateDeviceExtensionProperties();

        std::set<std::string> requiredExtensions(deviceExtensions.begin(), deviceExtensions.end());

        for (const auto &extension: availableExtensions) {
            requiredExtensions.erase(extension.extensionName);
        }

        return requiredExtensions.empty();
    }

    vk::PhysicalDevice MPhysDev::getHandle() {
        return *physDev;
    }
}