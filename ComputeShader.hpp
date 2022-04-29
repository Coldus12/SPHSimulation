#ifndef SPHSIMULATION_COMPUTESHADER_HPP
#define SPHSIMULATION_COMPUTESHADER_HPP

#include <string>
#include "VltavaFunctions.hpp"
#include "Buffer.hpp"

namespace Vltava {
    class ComputeShader {
    public:
        ComputeShader(std::string pathToShader);
        ~ComputeShader();
        void recreatePipeline();
        void createPipeline();
        vk::Pipeline getPipelineHandle();
        vk::PipelineLayout getPipelineLayoutHandle();

        void setBuffers(const std::vector<Buffer>* uniformBuffers, const std::vector<Buffer>* storageBuffers);
        //void dispatch(const vk::CommandBuffer &cmdBuffer, uint32_t groupCountX, uint32_t groupCountY, uint32_t groupCountZ);

        void bindPipelineAndDescriptors(vk::CommandBuffer& cmdBuffer);

    private:
        std::string computeShaderPath;

        vk::PipelineLayout pipelineLayout;
        vk::Pipeline computePipeline;
        vk::DescriptorSetLayout setLayout;
        vk::DescriptorPool setPool;
        vk::DescriptorSet set;
    };
}


#endif //SPHSIMULATION_COMPUTESHADER_HPP
