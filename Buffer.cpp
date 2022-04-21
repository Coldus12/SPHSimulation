//
// Created by PCF112021 on 3/19/2022.
//

#include <glm/vec2.hpp>
#include "Buffer.hpp"

namespace Vltava {
    Buffer::Buffer(const VulkanResources &vkResources,
                   vk::DeviceSize bufferSize,
                   vk::BufferUsageFlags usage,
                   vk::MemoryPropertyFlags memFlags) {

        updateResources(vkResources);

        auto pair = createBuffer(bufferSize, usage, memFlags);
        vkBuffer = pair.first;
        vkBufferMemory = pair.second;
    }

    Buffer::Buffer(Buffer &&buff) noexcept {
        buff.moved = true;
        this->vkBuffer = buff.vkBuffer;
        this->vkBufferMemory = buff.vkBufferMemory;
        this->dataSize = buff.dataSize;
    }

    Buffer::~Buffer() {
        if (!moved) {
            res.dev->getHandle().destroyBuffer(vkBuffer);
            res.dev->getHandle().freeMemory(vkBufferMemory);
        }
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

        vk::Buffer localBuffer = res.dev->getHandle().createBuffer(bufferInfo);

        vk::MemoryRequirements memReq = res.dev->getHandle().getBufferMemoryRequirements(localBuffer);
        vk::PhysicalDeviceMemoryProperties memProperties = res.physDev->getHandle().getMemoryProperties();

        vk::MemoryAllocateInfo memAllocInfo(memReq.size, {});
        memAllocInfo.memoryTypeIndex = findMemoryType(
                memReq.memoryTypeBits,
                memProperties,
                memFlags
        );

        vk::DeviceMemory devMem = res.dev->getHandle().allocateMemory(memAllocInfo);

        return std::pair<vk::Buffer, vk::DeviceMemory>(localBuffer, devMem);
    }

    void Buffer::writeToBuffer(void *bufferData, size_t size) {
        bind(0);
        void* data = res.dev->getHandle().mapMemory(vkBufferMemory, 0, size);
        memcpy(data, bufferData, size);
        res.dev->getHandle().unmapMemory(vkBufferMemory);

        dataSize = size;
    }

    void Buffer::copyBuffer(vk::Buffer src, vk::Buffer dst, vk::DeviceSize size) {
        if (res.graphicsQueue == nullptr)
            throw std::runtime_error("Buffer's static VulkanResources is null!");

        vk::CommandBufferAllocateInfo allocateInfo(
                *res.commandPool,
                vk::CommandBufferLevel::ePrimary,
                1
        );

        vk::CommandBuffer commandBuffer = res.dev->getHandle().allocateCommandBuffers(allocateInfo).front();

        vk::CommandBufferBeginInfo beginInfo(vk::CommandBufferUsageFlagBits::eOneTimeSubmit);
        commandBuffer.begin(beginInfo);

        vk::BufferCopy copyRegion(0,0, size);
        commandBuffer.copyBuffer(src, dst, copyRegion);
        commandBuffer.end();

        vk::SubmitInfo submitInfo;
        submitInfo.commandBufferCount = 1;
        submitInfo.pCommandBuffers = &commandBuffer;

        //vk::SubmitInfo submitInfo({}, {}, commandBuffer);

        res.graphicsQueue->submit(submitInfo);
        res.graphicsQueue->waitIdle();

        res.dev->getHandle().freeCommandBuffers(*res.commandPool, commandBuffer);
    }

    vk::Buffer Buffer::getBufferHandle() const {
        return vkBuffer;
    }

    size_t Buffer::getSize() const {
        return dataSize;
    }

    void Buffer::setSize(size_t size) {
        dataSize = size;
    }

    void Buffer::bind(int offset) {
        if (!bound)
            res.dev->getHandle().bindBufferMemory(vkBuffer, vkBufferMemory, offset);

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