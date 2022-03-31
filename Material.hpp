//
// Created by PCF112021 on 3/28/2022.
//

#ifndef SPHSIMULATION_MATERIAL_HPP
#define SPHSIMULATION_MATERIAL_HPP

#include "VltavaFunctions.hpp"
#include "Buffer.hpp"

namespace Vltava {
    class Material {
    public:
        Material(VulkanResources &resources, std::string pathToVert, std::string pathToFrag);
        void updateResources(const VulkanResources& resources);

        template<typename T>
        void uploadVertexData(std::vector<T> vertices) {
            vk::DeviceSize bufferSize = sizeof(T) * vertices.size();
            nrOfVertices = vertices.size();

            // Staging buffer creation
            //---------------------------------
            Buffer stagingBuffer(
                    res,
                    bufferSize,
                    vk::BufferUsageFlagBits::eVertexBuffer | vk::BufferUsageFlagBits::eTransferSrc,
                    vk::MemoryPropertyFlagBits::eHostVisible | vk::MemoryPropertyFlagBits::eHostCoherent
            );

            stagingBuffer.writeToBuffer(vertices.data(), bufferSize);

            // Device local buffer creation
            //---------------------------------
            vertexBuffer = std::make_unique<Buffer>(
                    res,
                    bufferSize,
                    vk::BufferUsageFlagBits::eTransferDst | vk::BufferUsageFlagBits::eVertexBuffer,
                    vk::MemoryPropertyFlagBits::eDeviceLocal
            );

            vertexBuffer->bind(0);

            // Copying data from staging buffer to local buffer
            //---------------------------------
            Buffer::copyBuffer(stagingBuffer.getBufferHandle(),vertexBuffer->getBufferHandle(),bufferSize);
        }

        void uploadIndexData(std::vector<uint16_t> indices);
        void setBuffers(const std::vector<Buffer>* uniformBuffers, const std::vector<Buffer>* storageBuffers);
        void setBuffers(const std::vector<Buffer*>* uniformBuffers, const std::vector<Buffer*>* storageBuffers);
        void createPipeline(std::vector<vk::VertexInputBindingDescription> bindingDescriptions,
                            std::vector<vk::VertexInputAttributeDescription> attributeDescriptions);
        void draw(const vk::raii::CommandBuffer& buffer, uint32_t currentFrame);

        void recreatePipeline();

    private:
        VulkanResources& res;
        float aspect;
        int nrOfIndices = 0;
        int nrOfVertices = 0;
        std::string vertPath;
        std::string fragPath;

        std::unique_ptr<vk::raii::PipelineLayout> pipelineLayout;
        std::unique_ptr<vk::raii::Pipeline> pipeline;
        std::unique_ptr<vk::raii::DescriptorSetLayout> setLayout;
        std::unique_ptr<vk::raii::DescriptorPool> setPool;
        std::vector<vk::raii::DescriptorSet> sets;

        std::unique_ptr<Buffer> vertexBuffer;
        std::unique_ptr<Buffer> indexBuffer;

        std::vector<vk::VertexInputBindingDescription> bindings;
        std::vector<vk::VertexInputAttributeDescription> attribs;
    };
}

#endif //SPHSIMULATION_MATERIAL_HPP