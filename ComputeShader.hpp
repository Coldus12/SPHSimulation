#ifndef SPHSIMULATION_COMPUTESHADER_HPP
#define SPHSIMULATION_COMPUTESHADER_HPP

#include <string>
#include "VltavaFunctions.hpp"
#include "Buffer.hpp"

namespace Vltava {
    class ComputeShader {
    public:
        ComputeShader(VulkanResources &resources, std::string pathToShader);
        void updateResources(const VulkanResources& resources);
        void createPipeline();
        vk::Pipeline getPipelineHandle();
        vk::PipelineLayout getPipelineLayoutHandle();

        void setStorageBuffers(const std::vector<Buffer>& buffers);
        void createCommandBuffer(uint32_t computeQueueFamily);
        void dispatch(uint32_t groupCountX, uint32_t groupCountY, uint32_t groupCountZ);

    private:
        std::string computeShaderPath;
        VulkanResources& resources;

        std::unique_ptr<vk::raii::PipelineLayout> pipelineLayout;
        std::unique_ptr<vk::raii::Pipeline> computePipeline;
        std::unique_ptr<vk::raii::DescriptorSetLayout> setLayout;
        std::unique_ptr<vk::raii::DescriptorPool> setPool;
        std::unique_ptr<vk::raii::DescriptorSet> set;
        std::unique_ptr<vk::raii::CommandPool> cmdPool;
        std::unique_ptr<vk::raii::CommandBuffer> cmdBuffer;
    };
}


#endif //SPHSIMULATION_COMPUTESHADER_HPP
