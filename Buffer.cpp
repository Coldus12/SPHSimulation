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
        vkBuffer = std::make_unique<MBuffer>(res.dev->getHandle(), pair.first);
        vkBufferMemory = std::make_unique<MDeviceMemory>(res.dev->getHandle(), pair.second);
    }

    std::pair<vk::Buffer, vk::DeviceMemory> Buffer::createBuffer(vk::DeviceSize bufferSize,
                                                                             vk::BufferUsageFlags usage,
                                                                             vk::MemoryPropertyFlags memFlags) {
        vk::BufferCreateInfo bufferInfo(
                {},
                bufferSize,
                usage,
                vk::SharingMode::eExclusive
        );

        MBuffer localBuffer(res.dev->getHandle(), bufferInfo);

        vk::MemoryRequirements memReq = res.dev->getHandle().getBufferMemoryRequirements(localBuffer.getHandle());
        vk::PhysicalDeviceMemoryProperties memProperties = res.physDev->getHandle().getMemoryProperties();

        vk::MemoryAllocateInfo memAllocInfo(memReq.size, {});
        memAllocInfo.memoryTypeIndex = findMemoryType(
                memReq.memoryTypeBits,
                memProperties,
                memFlags
        );

        MDeviceMemory devMem(res.dev->getHandle(), memAllocInfo);

        return std::pair<vk::Buffer, vk::DeviceMemory>(localBuffer.getHandle(), devMem.getHandle());
    }

    void Buffer::writeToBuffer(void *bufferData, size_t size) {
        bind(0);

        void* data = map<void*>(0, VK_WHOLE_SIZE);
        memcpy(data, bufferData, size);
        unmap();

        dataSize = size;
    }

    void Buffer::copyBuffer(vk::Buffer src, vk::Buffer dst, vk::DeviceSize size) {
        if (res.dev == nullptr)
            throw std::runtime_error("Buffer's static VulkanResources is null!");

        vk::CommandBufferAllocateInfo allocateInfo(
                res.commandPool->getHandle(),
                vk::CommandBufferLevel::ePrimary,
                1
        );

        MCommandBuffers commandBuffers(res.dev->getHandle(), allocateInfo);

        vk::CommandBufferBeginInfo beginInfo(vk::CommandBufferUsageFlagBits::eOneTimeSubmit);
        commandBuffers.getBuffers()[0].begin(beginInfo);

        vk::BufferCopy copyRegion(0,0, size);
        commandBuffers.getBuffers()[0].copyBuffer(src, dst, copyRegion);
        commandBuffers.getBuffers()[0].end();

        vk::SubmitInfo submitInfo;
        submitInfo.commandBufferCount = 1;
        submitInfo.pCommandBuffers = &commandBuffers.getBuffers()[0];

        res.graphicsQueue->submit(submitInfo);
        res.graphicsQueue->waitIdle();
    }

    vk::Buffer Buffer::getBufferHandle() const {
        return vkBuffer->getHandle();
    }

    size_t Buffer::getSize() const {
        return dataSize;
    }

    void Buffer::setSize(size_t size) {
        dataSize = size;
    }

    void Buffer::bind(int offset) {
        if (!bound)
            res.dev->getHandle().bindBufferMemory(vkBuffer->getHandle(), vkBufferMemory->getHandle(), offset);

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