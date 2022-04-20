#ifndef SPHSIMULATION_MDEVICEMEMORY_HPP
#define SPHSIMULATION_MDEVICEMEMORY_HPP

#include <vulkan/vulkan.hpp>

namespace Vltava {
    class MDeviceMemory {
    public:
        MDeviceMemory(vk::Device dev, vk::MemoryAllocateInfo info) : dev(dev) {
            handle = std::make_unique<vk::DeviceMemory>(dev.allocateMemory(info));
        }

        MDeviceMemory(vk::Device dev, vk::DeviceMemory mem) : dev(dev) {
            handle = std::make_unique<vk::DeviceMemory>(mem);
        }

        ~MDeviceMemory() {
            dev.freeMemory(*handle);
        }

        vk::DeviceMemory getHandle() {
            return *handle;
        }

        vk::DeviceMemory* getAddress() {
            return handle.get();
        }

    private:
        std::unique_ptr<vk::DeviceMemory> handle;
        vk::Device dev;
    };
}


#endif //SPHSIMULATION_MDEVICEMEMORY_HPP
