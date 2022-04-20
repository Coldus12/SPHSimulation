#ifndef SPHSIMULATION_MCOMMANDBUFFERS_HPP
#define SPHSIMULATION_MCOMMANDBUFFERS_HPP

#include <vulkan/vulkan.hpp>

namespace Vltava {
    class MCommandBuffers {
    public:
        MCommandBuffers(vk::Device dev, vk::CommandBufferAllocateInfo allocInfo) : dev(dev), info(allocInfo) {
            commandBuffers = dev.allocateCommandBuffers(allocInfo);
        }

        ~MCommandBuffers() {
            dev.freeCommandBuffers(info.commandPool, commandBuffers);
        }

        std::vector<vk::CommandBuffer> getBuffers() {
            return commandBuffers;
        }

        std::vector<vk::CommandBuffer>* getAddress() {
            return &commandBuffers;
        }
    private:
        vk::Device dev;
        vk::CommandBufferAllocateInfo info;
        std::vector<vk::CommandBuffer> commandBuffers;
    };
}

#endif //SPHSIMULATION_MCOMMANDBUFFERS_HPP
