#ifndef SPHSIMULATION_MPHYSDEV_HPP
#define SPHSIMULATION_MPHYSDEV_HPP


#include <vulkan/vulkan.hpp>

namespace Vltava {
    class MPhysDev {
    public:
        MPhysDev(vk::Instance instance, const std::vector<const char *> deviceExtensions);
        ~MPhysDev();

        vk::PhysicalDevice getHandle();
    private:
        vk::Instance instance;
        std::unique_ptr<vk::PhysicalDevice> physDev;
        const std::vector<const char *> deviceExtensions;

        void selectPhysicalDevice();
        vk::PhysicalDevice getBestDevice(const std::vector<vk::PhysicalDevice> &devs);
        uint32_t rateDeviceSuitability(const vk::PhysicalDevice &dev);
        bool checkDeviceExtensionSupport(const vk::PhysicalDevice &dev);
    };
}


#endif //SPHSIMULATION_MPHYSDEV_HPP
