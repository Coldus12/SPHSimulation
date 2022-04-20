#ifndef SPHSIMULATION_MRENDERPASS_HPP
#define SPHSIMULATION_MRENDERPASS_HPP

#include <vulkan/vulkan.hpp>

namespace Vltava {
    class MRenderPass {
    public:
        MRenderPass(vk::Device dev, vk::SubpassDescription subpass, vk::AttachmentDescription colorAttachment) : dev(dev) {
            handle = dev.createRenderPass(vk::RenderPassCreateInfo(vk::RenderPassCreateFlags(), colorAttachment, subpass));
        }

        ~MRenderPass() {
            dev.destroyRenderPass(handle);
        }

        vk::RenderPass getHandle() {
            return handle;
        }
    private:
        vk::RenderPass handle;
        vk::Device dev;
    };
}

#endif //SPHSIMULATION_MRENDERPASS_HPP
