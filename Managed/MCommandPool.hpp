#ifndef SPHSIMULATION_MCOMMANDPOOL_HPP
#define SPHSIMULATION_MCOMMANDPOOL_HPP

#include <vulkan/vulkan.hpp>

namespace Vltava {
    class MCommandPool {
    public:
        MCommandPool(vk::Device dev, vk::CommandPoolCreateInfo info) : dev(dev) {
            handle = dev.createCommandPool(info);
        }

        ~MCommandPool() {
            dev.destroyCommandPool(handle);
        }

        vk::CommandPool getHandle() {
            return handle;
        }
    private:
        vk::CommandPool handle;
        vk::Device dev;
    };
}


#endif //SPHSIMULATION_MCOMMANDPOOL_HPP
