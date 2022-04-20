#ifndef SPHSIMULATION_MPIPELINELAYOUT_HPP
#define SPHSIMULATION_MPIPELINELAYOUT_HPP

#include <vulkan/vulkan.hpp>

namespace Vltava {
    class MPipelineLayout {
    public:
        MPipelineLayout(vk::Device dev, vk::PipelineLayoutCreateInfo info) : dev(dev) {
            handle = std::make_unique<vk::PipelineLayout>(dev.createPipelineLayout(info));
        }

        ~MPipelineLayout() {
            dev.destroyPipelineLayout(*handle);
        }

        vk::PipelineLayout getHandle() {
            return *handle;
        }

        vk::PipelineLayout* getAddress() {
            return handle.get();
        }

    private:
        std::unique_ptr<vk::PipelineLayout> handle;
        vk::Device dev;
    };
}

#endif //SPHSIMULATION_MPIPELINELAYOUT_HPP
