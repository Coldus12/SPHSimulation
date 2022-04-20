#ifndef SPHSIMULATION_MFRAMEBUFFER_HPP
#define SPHSIMULATION_MFRAMEBUFFER_HPP

#include <vulkan/vulkan.hpp>

namespace Vltava {
    class MFramebuffer {
    public:
        MFramebuffer(vk::Device dev, vk::FramebufferCreateInfo info) : dev(dev) {
            handle = dev.createFramebuffer(info);
        }

        ~MFramebuffer() {
            dev.destroyFramebuffer(handle);
        }

        vk::Framebuffer getHandle() {
            return handle;
        }
    private:
        vk::Framebuffer handle;
        vk::Device dev;
    };
}

#endif //SPHSIMULATION_MFRAMEBUFFER_HPP
