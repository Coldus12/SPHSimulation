#ifndef SPHSIMULATION_MDESCRIPTORSETS_HPP
#define SPHSIMULATION_MDESCRIPTORSETS_HPP

#include <vulkan/vulkan.hpp>

namespace Vltava {
    class MDescriptorSets {
    public:
        MDescriptorSets(vk::Device dev, vk::DescriptorSetAllocateInfo allocInfo) : dev(dev), info(allocInfo) {
            descriptorSets = dev.allocateDescriptorSets(allocInfo);
        }

        ~MDescriptorSets() {
            dev.freeDescriptorSets(info.descriptorPool, descriptorSets);
        }

        std::vector<vk::DescriptorSet> getSets() {
            return descriptorSets;
        }

        std::vector<vk::DescriptorSet>* getAddress() {
            return &descriptorSets;
        }

    private:
        vk::Device dev;
        vk::DescriptorSetAllocateInfo info;
        std::vector<vk::DescriptorSet> descriptorSets;
    };
}

#endif //SPHSIMULATION_MDESCRIPTORSETS_HPP
