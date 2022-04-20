#ifndef SPHSIMULATION_MDESCRIPTORSET_HPP
#define SPHSIMULATION_MDESCRIPTORSET_HPP

#include <vulkan/vulkan.hpp>

namespace Vltava {
    class MDescriptorSet {
    public:
        MDescriptorSet(vk::Device dev, vk::DescriptorSet set, vk::DescriptorPool p) : dev(dev), poolie(p) {
            handle = std::make_unique<vk::DescriptorSet>(set);
        }

        MDescriptorSet(MDescriptorSet&& set)  noexcept {
            handle = std::move(set.handle);
            this->dev = set.dev;
            this->poolie = set.poolie;
        }

        ~MDescriptorSet() {
            dev.freeDescriptorSets(poolie, *handle);
        }

        vk::DescriptorSet getHandle() {
            return *handle;
        }

        vk::DescriptorSet* getAddress() {
            return handle.get();
        }

    private:
        std::unique_ptr<vk::DescriptorSet> handle;
        vk::Device dev;
        vk::DescriptorPool poolie;
    };
}

#endif //SPHSIMULATION_MDESCRIPTORSET_HPP
