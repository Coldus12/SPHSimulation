#ifndef SPHSIMULATION_MPIPELINE_HPP
#define SPHSIMULATION_MPIPELINE_HPP

#include <vulkan/vulkan.hpp>

namespace Vltava {
    class MPipeline {
    public:
        MPipeline(vk::Device dev, vk::ComputePipelineCreateInfo info) : dev(dev) {
            auto tmp = dev.createComputePipeline({}, info);
            handle = std::make_unique<vk::Pipeline>(tmp);
        }

        MPipeline(vk::Device dev, vk::GraphicsPipelineCreateInfo info) : dev(dev) {
            auto tmp = dev.createGraphicsPipeline({}, info);
            handle = std::make_unique<vk::Pipeline>(tmp);
        }

        ~MPipeline() {
            dev.destroyPipeline(*handle);
        }

        vk::Pipeline getHandle() {
            return *handle;
        }

        vk::Pipeline* getAddress() {
            return handle.get();
        }

    private:
        std::unique_ptr<vk::Pipeline> handle;
        vk::Device dev;
    };
}

#endif //SPHSIMULATION_MPIPELINE_HPP
