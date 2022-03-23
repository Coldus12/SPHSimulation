//
// Created by PCF112021 on 3/22/2022.
//

#include "ComputeShader.h"

namespace Vltava {

    ComputeShader::ComputeShader(VulkanResources &resources, std::string pathToShader) : resources(resources), computeShaderPath(pathToShader) {
        //createPipeline();
    }

    void ComputeShader::updateResources(const VulkanResources &res) {
        // Updating the resources
        resources.dev = res.dev;
        resources.extent = res.extent;
        resources.renderPass = res.renderPass;
        resources.physDev = res.physDev;
        resources.instance = res.instance;
        resources.commandPool = res.commandPool;
        resources.graphicsQueue = res.graphicsQueue;
        resources.computeQueue = res.computeQueue;
        resources.FRAMES_IN_FLIGHT = res.FRAMES_IN_FLIGHT;

        // Deleting the out of date pipeline
        computePipeline.reset();

        // Creating new pipeline with the updated resources
        createPipeline(layoutCount, pLayouts);
    }

    vk::Pipeline ComputeShader::getPipelineHandle() {
        return **computePipeline;
    }

    vk::PipelineLayout ComputeShader::getPipelineLayoutHandle() {
        return **pipelineLayout;
    }

    void ComputeShader::createPipeline(uint32_t layoutCount, vk::DescriptorSetLayout* pLayouts) {
        this->layoutCount = layoutCount;
        this->pLayouts = pLayouts;

        auto code = readFile(computeShaderPath);
        vk::ShaderModuleCreateInfo compInfo({}, code.size(), reinterpret_cast<uint32_t*>(code.data()));
        vk::raii::ShaderModule compModule(*resources.dev, compInfo);
        vk::PipelineShaderStageCreateInfo compStageInfo({}, vk::ShaderStageFlagBits::eCompute, *compModule, "main");
        vk::PipelineShaderStageCreateInfo stages[] = {compStageInfo};

        vk::PipelineLayoutCreateInfo layoutInfo(
                {},
                layoutCount,
                pLayouts,
                0,
                nullptr
        );

        pipelineLayout = std::make_unique<vk::raii::PipelineLayout>(*resources.dev, layoutInfo);

        vk::Pipeline plHandle;
        vk::ComputePipelineCreateInfo computeInfo(
                {},
                compStageInfo,
                **pipelineLayout,
                plHandle,
                0
        );

        computePipeline = std::make_unique<vk::raii::Pipeline>(*resources.dev, nullptr, computeInfo);
    }
};