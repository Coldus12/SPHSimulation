//
// Created by PCF112021 on 3/19/2022.
//

#include "Buffer.hpp"

namespace Vltava {
    Buffer::Buffer(const VulkanResources &vkResources,
                   vk::DeviceSize bufferSize,
                   vk::BufferUsageFlags usage,
                   vk::MemoryPropertyFlags memFlags) {

        updateResources(vkResources);

        auto pair = createBuffer(bufferSize, usage, memFlags);
        vkBuffer = std::make_unique<vk::raii::Buffer>(std::move(pair.first));
        vkBufferMemory = std::make_unique<vk::raii::DeviceMemory>(std::move(pair.second));
    }

    std::pair<vk::raii::Buffer, vk::raii::DeviceMemory> Buffer::createBuffer(vk::DeviceSize bufferSize,
                                                                             vk::BufferUsageFlags usage,
                                                                             vk::MemoryPropertyFlags memFlags) {
        vk::BufferCreateInfo bufferInfo(
                {},
                bufferSize,
                usage,
                vk::SharingMode::eExclusive
        );

        vk::raii::Buffer localBuffer(*res.dev, bufferInfo);

        vk::MemoryRequirements memReq = localBuffer.getMemoryRequirements();
        vk::PhysicalDeviceMemoryProperties memProperties = res.physDev->getMemoryProperties();

        vk::MemoryAllocateInfo memAllocInfo(memReq.size, {});
        memAllocInfo.memoryTypeIndex = findMemoryType(
                memReq.memoryTypeBits,
                memProperties,
                memFlags
        );

        vk::raii::DeviceMemory devMem(*res.dev, memAllocInfo);

        return std::pair<vk::raii::Buffer, vk::raii::DeviceMemory>(std::move(localBuffer), std::move(devMem));
    }

    void Buffer::writeToBuffer(void *bufferData, size_t size) {
        bind(0);
        void* data = vkBufferMemory->mapMemory(0, VK_WHOLE_SIZE);
        memcpy(data, bufferData, size);
        vkBufferMemory->unmapMemory();

        dataSize = size;
    }

    void Buffer::copyBuffer(vk::Buffer src, vk::Buffer dst, vk::DeviceSize size) {
        if (res.dev == nullptr)
            throw std::runtime_error("Buffer's static VulkanResources is null!");

        vk::CommandBufferAllocateInfo allocateInfo(
                **res.commandPool,
                vk::CommandBufferLevel::ePrimary,
                1
        );

        std::vector<vk::raii::CommandBuffer> commandBuffers = res.dev->allocateCommandBuffers(allocateInfo);

        vk::CommandBufferBeginInfo beginInfo(vk::CommandBufferUsageFlagBits::eOneTimeSubmit);
        commandBuffers[0].begin(beginInfo);

        vk::BufferCopy copyRegion(0,0, size);
        commandBuffers[0].copyBuffer(src, dst, copyRegion);
        commandBuffers[0].end();

        vk::SubmitInfo submitInfo;
        submitInfo.commandBufferCount = 1;
        submitInfo.pCommandBuffers = &*commandBuffers[0];

        res.graphicsQueue->submit(submitInfo);
        res.graphicsQueue->waitIdle();
    }

    vk::Buffer Buffer::getBufferHandle() const {
        return **vkBuffer;
    }

    size_t Buffer::getSize() const {
        return dataSize;
    }

    void Buffer::setSize(size_t size) {
        dataSize = size;
    }

    void Buffer::bind(int offset) {
        if (!bound)
            vkBuffer->bindMemory(**vkBufferMemory, offset);

        bound = true;
    }

    void Buffer::updateResources(const VulkanResources &vkResources) {
        res.dev = vkResources.dev;
        res.extent = vkResources.extent;
        res.renderPass = vkResources.renderPass;
        res.physDev = vkResources.physDev;
        res.instance = vkResources.instance;
        res.commandPool = vkResources.commandPool;
        res.graphicsQueue = vkResources.graphicsQueue;
        res.FRAMES_IN_FLIGHT = vkResources.FRAMES_IN_FLIGHT;
    }
}