#ifndef SPHSIMULATION_MBUFFER_HPP
#define SPHSIMULATION_MBUFFER_HPP

#include <vulkan/vulkan.hpp>

namespace Vltava {
    class MBuffer {
    public:
        MBuffer(vk::Device dev, vk::BufferCreateInfo info) : dev(dev) {
            handle = std::make_unique<vk::Buffer>(dev.createBuffer(info));
        }

        MBuffer(vk::Device dev, vk::Buffer buffer) : dev(dev) {
            handle = std::make_unique<vk::Buffer>(buffer);
        }

        ~MBuffer() {
            dev.destroyBuffer(*handle);
        }

        vk::Buffer getHandle() {
            return *handle;
        }

        vk::Buffer* getAddress() {
            return handle.get();
        }

    private:
        std::unique_ptr<vk::Buffer> handle;
        vk::Device dev;
    };
}

#endif //SPHSIMULATION_MBUFFER_HPP
