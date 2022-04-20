#ifndef SPHSIMULATION_MDESCRIPTORSETLAYOUT_HPP
#define SPHSIMULATION_MDESCRIPTORSETLAYOUT_HPP

#include <vulkan/vulkan.hpp>

namespace Vltava {
    class MDescriptorSetLayout {
    public:
        MDescriptorSetLayout(vk::Device dev, vk::DescriptorSetLayoutCreateInfo info) : dev(dev) {
            handle = std::make_unique<vk::DescriptorSetLayout>(dev.createDescriptorSetLayout(info));
        }

        ~MDescriptorSetLayout() {
            dev.destroyDescriptorSetLayout(*handle);
        }

        vk::DescriptorSetLayout getHandle() {
            return *handle;
        }

        vk::DescriptorSetLayout* getAddress() {
            return handle.get();
        }

    private:
        std::unique_ptr<vk::DescriptorSetLayout> handle;
        vk::Device dev;
    };
}

#endif //SPHSIMULATION_MDESCRIPTORSETLAYOUT_HPP
