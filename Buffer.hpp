#ifndef SPHSIMULATION_BUFFER_HPP
#define SPHSIMULATION_BUFFER_HPP

#include <vulkan/vulkan.hpp>
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

        template<typename T>
        std::vector<T> getData() {
            std::vector<T> retVector;
            size_t objectSize = sizeof(T);
            size_t nrOfObjects = dataSize / objectSize;
            size_t wholeSize = sizeof(T) * nrOfObjects;
            retVector.reserve(nrOfObjects);

            bind(0);
            //T* bufferData = (T*) res.dev->getHandle().mapMemory<T*>(vkBufferMemory->getHandle(), 0, wholeSize);
            T* bufferData = map<T>();

            for (int i = 0; i < nrOfObjects; i++) {
                retVector.push_back(bufferData[i]);
            }

            //res.dev->getHandle().unmapMemory<T*>();
            unmap();

            return retVector;
        }

        vk::Buffer getBufferHandle() const;
        size_t getSize() const;
        void setSize(size_t size); // Need this for output storage buffers which would have a size of 0

        static void copyBuffer(vk::Buffer src, vk::Buffer dst, vk::DeviceSize size);
        static void updateResources(const VulkanResources& vkResources);

        // Taken from: https://github.com/jherico/Vulkan/blob/cpp/base/vks/allocation.hpp
        void* mapped{ nullptr };
        template <typename T = void>
        inline T* map(size_t offset = 0, VkDeviceSize size = VK_WHOLE_SIZE) {
            mapped = res.dev->getHandle().mapMemory(vkBufferMemory->getHandle(), offset, size, vk::MemoryMapFlags());
            return (T*)mapped;
        }

        inline void unmap() {
            res.dev->getHandle().unmapMemory(vkBufferMemory->getHandle());
            mapped = nullptr;
        }

    private:
        size_t dataSize = 0;
        bool bound = false;
        inline static VulkanResources res;
        std::unique_ptr<MBuffer> vkBuffer;
        std::unique_ptr<MDeviceMemory> vkBufferMemory;

        std::pair<vk::Buffer, vk::DeviceMemory> createBuffer(vk::DeviceSize bufferSize,
                                                             vk::BufferUsageFlags usage,
                                                             vk::MemoryPropertyFlags memFlags);
    };
}


#endif //SPHSIMULATION_BUFFER_HPP
