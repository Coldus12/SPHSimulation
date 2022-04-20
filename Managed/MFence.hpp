#ifndef SPHSIMULATION_MFENCE_HPP
#define SPHSIMULATION_MFENCE_HPP

#include <vulkan/vulkan.hpp>

namespace Vltava {
    class MFence {
    public:
        MFence(vk::Device dev, vk::FenceCreateInfo info) : dev(dev) {
            handle = dev.createFence(info);
        }

        ~MFence() {
            dev.destroyFence(handle);
        }

        vk::Fence getHandle() {
            return handle;
        }
    private:
        vk::Fence handle;
        vk::Device dev;
    };
}

#endif //SPHSIMULATION_MFENCE_HPP
