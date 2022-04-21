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
            //T* bufferData = (T*) res.dev->getHandle().mapMemory<T*>(vkBufferMemory->getHandle(), 0, wholeSize);
            //T* bufferData = map<T>();

            T* bufferData = (T*) res.dev->getHandle().mapMemory(vkBufferMemory, 0, wholeSize);

            for (int i = 0; i < nrOfObjects; i++) {
                retVector.push_back(bufferData[i]);
            }

            //res.dev->getHandle().unmapMemory<T*>();
            res.dev->getHandle().unmapMemory(vkBufferMemory);

            return retVector;
        }

        vk::Buffer getBufferHandle() const;
        size_t getSize() const;
        void setSize(size_t size); // Need this for output storage buffers which would have a size of 0

        static void copyBuffer(vk::Buffer src, vk::Buffer dst, vk::DeviceSize size);
        static void updateResources(const VulkanResources& vkResources);

    private:
        size_t dataSize = 0;
        bool bound = false;
        bool moved = false;
        inline static VulkanResources res;
        vk::Buffer vkBuffer;
        vk::DeviceMemory vkBufferMemory;

        std::pair<vk::Buffer, vk::DeviceMemory> createBuffer(vk::DeviceSize bufferSize,
                                                             vk::BufferUsageFlags usage,
                                                             vk::MemoryPropertyFlags memFlags);
    };
}


#endif //SPHSIMULATION_BUFFER_HPP
