#ifndef SPHSIMULATION_BUFFER_HPP
#define SPHSIMULATION_BUFFER_HPP

#include <vulkan/vulkan.hpp>
#include <vulkan/vulkan_raii.hpp>
#include "VltavaFunctions.hpp"

namespace Vltava {
    class Buffer {
    public:
        Buffer() = delete;
        Buffer(const VulkanResources& vkResources,
               vk::DeviceSize bufferSize,
               vk::BufferUsageFlags usage,
               vk::MemoryPropertyFlags memFlags
        );

        void bind(int offset);
        void writeToBuffer(void* bufferData, size_t size);
        void writeData(void* data, size_t size);
        vk::Buffer getBufferHandle();

        static void copyBuffer(vk::Buffer src, vk::Buffer dst, vk::DeviceSize size);
        static void updateResources(const VulkanResources& vkResources);

    private:
        bool bound = false;
        inline static VulkanResources res;
        std::unique_ptr<vk::raii::Buffer> vkBuffer;
        std::unique_ptr<vk::raii::DeviceMemory> vkBufferMemory;

        std::pair<vk::raii::Buffer, vk::raii::DeviceMemory> createBuffer(vk::DeviceSize bufferSize,
                                                                         vk::BufferUsageFlags usage,
                                                                         vk::MemoryPropertyFlags memFlags);
    };
}


#endif //SPHSIMULATION_BUFFER_HPP
