//
// Created by PCF112021 on 3/22/2022.
//

#ifndef SPHSIMULATION_COMPUTESHADER_H
#define SPHSIMULATION_COMPUTESHADER_H

#include <string>
#include "VltavaFunctions.hpp"

namespace Vltava {
    class ComputeShader {
    public:
        ComputeShader(VulkanResources &resources, std::string pathToShader);
        void updateResources(const VulkanResources& resources);
        void createPipeline(uint32_t layoutCount, vk::DescriptorSetLayout* pLayouts);
        vk::Pipeline getPipelineHandle();
        vk::PipelineLayout getPipelineLayoutHandle();

    private:
        uint32_t layoutCount;
        vk::DescriptorSetLayout* pLayouts;

        std::string computeShaderPath;
        VulkanResources& resources;

        std::unique_ptr<vk::raii::PipelineLayout> pipelineLayout;
        std::unique_ptr<vk::raii::Pipeline> computePipeline;
        std::unique_ptr<vk::raii::Queue> computeQueue;
    };
}


#endif //SPHSIMULATION_COMPUTESHADER_H
