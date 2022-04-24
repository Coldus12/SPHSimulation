#ifndef SPHSIMULATION_BUFFER_HPP
#define SPHSIMULATION_BUFFER_HPP

#include <vulkan/vulkan.hpp>
#include "VltavaFunctions.hpp"

namespace Vltava {
    class Buffer {
    public:
        Buffer() = delete;
        Buffer(vk::DeviceSize bufferSize,
               vk::BufferUsageFlags usage,
               vk::MemoryPropertyFlags memFlags
        );
        Buffer(Buffer&& buff) noexcept;
        ~Buffer();

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
            T* bufferData = (T*) VulkanResources::getInstance().logDev->getHandle().mapMemory(vkBufferMemory, 0, wholeSize);

            for (int i = 0; i < nrOfObjects; i++) {
                retVector.push_back(bufferData[i]);
            }

            VulkanResources::getInstance().logDev->getHandle().unmapMemory(vkBufferMemory);

            return retVector;
        }

        vk::Buffer getBufferHandle() const;
        size_t getSize() const;
        void setSize(size_t size); // Need this for output storage buffers which would have a size of 0
        static void copyBuffer(vk::Buffer src, vk::Buffer dst, vk::DeviceSize size);

    private:
        size_t dataSize = 0;
        bool bound = false;
        bool moved = false;
        vk::Buffer vkBuffer;
        vk::DeviceMemory vkBufferMemory;

        std::pair<vk::Buffer, vk::DeviceMemory> createBuffer(vk::DeviceSize bufferSize,
                                                             vk::BufferUsageFlags usage,
                                                             vk::MemoryPropertyFlags memFlags);
    };
}


#endif //SPHSIMULATION_BUFFER_HPP
