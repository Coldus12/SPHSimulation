//
// Created by PCF112021 on 3/19/2022.
//

#include <glm/vec2.hpp>
#include "Buffer.hpp"

namespace Vltava {
    Buffer::Buffer(vk::DeviceSize bufferSize,
                   vk::BufferUsageFlags usage,
                   vk::MemoryPropertyFlags memFlags) {

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
            VulkanResources::getInstance().logDev->getHandle().destroyBuffer(vkBuffer);
            VulkanResources::getInstance().logDev->getHandle().freeMemory(vkBufferMemory);
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

        vk::Buffer localBuffer = VulkanResources::getInstance().logDev->getHandle().createBuffer(bufferInfo);

        vk::MemoryRequirements memReq = VulkanResources::getInstance().logDev->getHandle().getBufferMemoryRequirements(localBuffer);
        vk::PhysicalDeviceMemoryProperties memProperties = VulkanResources::getInstance().physDev->getHandle().getMemoryProperties();

        vk::MemoryAllocateInfo memAllocInfo(memReq.size, {});
        memAllocInfo.memoryTypeIndex = findMemoryType(
                memReq.memoryTypeBits,
                memProperties,
                memFlags
        );

        vk::DeviceMemory devMem = VulkanResources::getInstance().logDev->getHandle().allocateMemory(memAllocInfo);

        return std::pair<vk::Buffer, vk::DeviceMemory>(localBuffer, devMem);
    }

    void Buffer::writeToBuffer(void *bufferData, size_t size) {
        bind(0);
        void* data = VulkanResources::getInstance().logDev->getHandle().mapMemory(vkBufferMemory, 0, size);
        memcpy(data, bufferData, size);
        VulkanResources::getInstance().logDev->getHandle().unmapMemory(vkBufferMemory);

        dataSize = size;
    }

    void Buffer::copyBuffer(vk::Buffer src, vk::Buffer dst, vk::DeviceSize size) {
        if (VulkanResources::getInstance().graphicsQueue == nullptr)
            throw std::runtime_error("[ERROR] VulkanResources' graphicsQueue is null!\n[ERROR] Thrown from Buffer.copyBuffer().");

        vk::CommandBufferAllocateInfo allocateInfo(
                *VulkanResources::getInstance().graphicalCmdPool,
                vk::CommandBufferLevel::ePrimary,
                1
        );

        vk::CommandBuffer commandBuffer = VulkanResources::getInstance().logDev->getHandle().allocateCommandBuffers(allocateInfo).front();

        vk::CommandBufferBeginInfo beginInfo(vk::CommandBufferUsageFlagBits::eOneTimeSubmit);
        commandBuffer.begin(beginInfo);

        vk::BufferCopy copyRegion(0,0, size);
        commandBuffer.copyBuffer(src, dst, copyRegion);
        commandBuffer.end();

        vk::SubmitInfo submitInfo;
        submitInfo.commandBufferCount = 1;
        submitInfo.pCommandBuffers = &commandBuffer;

        //vk::SubmitInfo submitInfo({}, {}, commandBuffer);

        VulkanResources::getInstance().graphicsQueue->submit(submitInfo);
        VulkanResources::getInstance().graphicsQueue->waitIdle();

        VulkanResources::getInstance().logDev->getHandle().freeCommandBuffers(*VulkanResources::getInstance().graphicalCmdPool, commandBuffer);
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
            VulkanResources::getInstance().logDev->getHandle().bindBufferMemory(vkBuffer, vkBufferMemory, offset);

        bound = true;
    }
}