#ifndef SPHSIMULATION_MDESCRIPTORPOOL_HPP
#define SPHSIMULATION_MDESCRIPTORPOOL_HPP

#include <vulkan/vulkan.hpp>

namespace Vltava {
    class MDescriptorPool {
    public:
        MDescriptorPool(vk::Device dev, vk::DescriptorPoolCreateInfo info) : dev(dev) {
            handle = std::make_unique<vk::DescriptorPool>(dev.createDescriptorPool(info));
        }

        ~MDescriptorPool() {
            dev.destroyDescriptorPool(*handle);
        }

        vk::DescriptorPool getHandle() {
            return *handle;
        }

        vk::DescriptorPool* getAddress() {
            return handle.get();
        }

    private:
        std::unique_ptr<vk::DescriptorPool> handle;
        vk::Device dev;
    };
}

#endif //SPHSIMULATION_MDESCRIPTORPOOL_HPP
