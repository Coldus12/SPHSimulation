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

        void setBuffers(const std::vector<Buffer>* uniformBuffers, const std::vector<Buffer>* storageBuffers);
        //void createCommandBuffer(uint32_t computeQueueFamily);
        void dispatch(const vk::CommandBuffer &cmdBuffer, uint32_t groupCountX, uint32_t groupCountY, uint32_t groupCountZ);

    private:
        std::string computeShaderPath;
        VulkanResources resources;

        std::unique_ptr<MPipelineLayout> pipelineLayout;
        std::unique_ptr<MPipeline> computePipeline;
        std::unique_ptr<MDescriptorSetLayout> setLayout;
        std::unique_ptr<MDescriptorPool> setPool;
        std::unique_ptr<MDescriptorSets> set;
    };
}


#endif //SPHSIMULATION_COMPUTESHADER_HPP
