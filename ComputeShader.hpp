#ifndef SPHSIMULATION_COMPUTESHADER_HPP
#define SPHSIMULATION_COMPUTESHADER_HPP

#include <string>
#include "VltavaFunctions.hpp"
#include "Buffer.hpp"

namespace Vltava {
    class ComputeShader {
    public:
        ComputeShader(VulkanResources &resources, std::string pathToShader);
        ~ComputeShader();
        void updateResources(const VulkanResources& resources);
        void createPipeline();
        vk::Pipeline getPipelineHandle();
        vk::PipelineLayout getPipelineLayoutHandle();

        void setBuffers(const std::vector<Buffer>* uniformBuffers, const std::vector<Buffer>* storageBuffers);
        void dispatch(const vk::CommandBuffer &cmdBuffer, uint32_t groupCountX, uint32_t groupCountY, uint32_t groupCountZ);

    private:
        std::string computeShaderPath;
        VulkanResources resources;

        vk::PipelineLayout pipelineLayout;
        vk::Pipeline computePipeline;
        vk::DescriptorSetLayout setLayout;
        vk::DescriptorPool setPool;
        vk::DescriptorSet set;
    };
}


#endif //SPHSIMULATION_COMPUTESHADER_HPP
