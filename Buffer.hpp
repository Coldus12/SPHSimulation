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

        template<typename T>
        //std::vector<T> getData(size_t nrOfObjects) {
        std::vector<T> getData() {
            std::vector<T> retVector;
            size_t objectSize = sizeof(T);
            size_t nrOfObjects = dataSize / objectSize;
            size_t wholeSize = sizeof(T) * nrOfObjects;
            retVector.reserve(nrOfObjects);

            bind(0);
            T* bufferData = (T*) vkBufferMemory->mapMemory(0, wholeSize);

            for (int i = 0; i < nrOfObjects; i++) {
                retVector.push_back(bufferData[i]);
            }

            vkBufferMemory->unmapMemory();

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
        inline static VulkanResources res;
        std::unique_ptr<vk::raii::Buffer> vkBuffer;
        std::unique_ptr<vk::raii::DeviceMemory> vkBufferMemory;

        std::pair<vk::raii::Buffer, vk::raii::DeviceMemory> createBuffer(vk::DeviceSize bufferSize,
                                                                         vk::BufferUsageFlags usage,
                                                                         vk::MemoryPropertyFlags memFlags);
    };
}


#endif //SPHSIMULATION_BUFFER_HPP
