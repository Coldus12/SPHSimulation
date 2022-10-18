#ifndef SPHSIMULATION_MLOGDEV_HPP
#define SPHSIMULATION_MLOGDEV_HPP

#include <vulkan/vulkan.hpp>
#include <set>
#include <memory>

namespace Vltava {
    class MLogDev {
    public:
        MLogDev(std::set<uint32_t> uniqueQueueFamilies,
                vk::PhysicalDevice physDev,
                const std::vector<const char *> deviceExtensions,
                const std::vector<const char *>* validationLayers = nullptr);
        ~MLogDev();

        vk::Device getHandle();
        vk::Device* getAddress() {
            return device.get();
        }
    private:
        std::unique_ptr<vk::Device> device;
    };
}

#endif //SPHSIMULATION_MLOGDEV_HPP
