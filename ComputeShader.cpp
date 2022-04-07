#include "ComputeShader.hpp"

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
        createPipeline();
    }

    vk::Pipeline ComputeShader::getPipelineHandle() {
        return **computePipeline;
    }

    vk::PipelineLayout ComputeShader::getPipelineLayoutHandle() {
        return **pipelineLayout;
    }

    void ComputeShader::createPipeline() {
        auto code = readFile(computeShaderPath);
        vk::ShaderModuleCreateInfo compInfo({}, code.size(), reinterpret_cast<uint32_t*>(code.data()));
        vk::raii::ShaderModule compModule(*resources.dev, compInfo);
        vk::PipelineShaderStageCreateInfo compStageInfo({}, vk::ShaderStageFlagBits::eCompute, *compModule, "main");
        vk::PipelineShaderStageCreateInfo stages[] = {compStageInfo};

        vk::PipelineLayoutCreateInfo layoutInfo(
                {},
                1,
                &**setLayout,
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

    void ComputeShader::setBuffers(const std::vector<Buffer>* uniformBuffers, const std::vector<Buffer>* storageBuffers) {
        // Descriptor set layout creation
        //---------------------------------
        std::vector<vk::DescriptorSetLayoutBinding> bindings;
        bindings.reserve(storageBuffers->size() + uniformBuffers->size());

        for (int i = 0; i < uniformBuffers->size(); i++) {
            bindings.emplace_back(i, vk::DescriptorType::eUniformBuffer, 1, vk::ShaderStageFlagBits::eCompute);
        }

        for (int i = uniformBuffers->size(); i < uniformBuffers->size() + storageBuffers->size(); i++) {
            bindings.emplace_back(i, vk::DescriptorType::eStorageBuffer, 1, vk::ShaderStageFlagBits::eCompute);
        }

        vk::DescriptorSetLayoutCreateInfo layoutInfo({}, bindings);
        setLayout = std::make_unique<vk::raii::DescriptorSetLayout>(*resources.dev, layoutInfo);

        // Descriptor pool creation
        //---------------------------------
        std::vector<vk::DescriptorPoolSize> sizes = {
                {vk::DescriptorType::eUniformBuffer, static_cast<uint32_t>(uniformBuffers->size())},
                {vk::DescriptorType::eStorageBuffer, static_cast<uint32_t>(storageBuffers->size())}
        };

        vk::DescriptorPoolCreateInfo poolInfo(
                vk::DescriptorPoolCreateFlagBits::eFreeDescriptorSet,
                20,
                sizes
        );

        setPool = std::make_unique<vk::raii::DescriptorPool>(*resources.dev, poolInfo);

        // Descriptor set creation
        //---------------------------------
        vk::DescriptorSetAllocateInfo allocInfo(**setPool, 1, &**setLayout);
        set = std::make_unique<vk::raii::DescriptorSet>(
                std::move(resources.dev->allocateDescriptorSets(allocInfo).front())
        );

        std::vector<vk::DescriptorBufferInfo> bufferInfos;
        std::vector<vk::WriteDescriptorSet> writeSets;
        bufferInfos.reserve(storageBuffers->size() + uniformBuffers->size());
        writeSets.reserve(storageBuffers->size() + uniformBuffers->size());

        for (int i = 0; i < uniformBuffers->size(); i++) {
            bufferInfos.emplace_back(
                    uniformBuffers->at(i).getBufferHandle(),
                    0,
                    uniformBuffers->at(i).getSize()
            );

            writeSets.emplace_back(
                    **set,
                    i,
                    0,
                    1,
                    vk::DescriptorType::eUniformBuffer,
                    nullptr,
                    &bufferInfos[i],
                    nullptr
            );
        }

        for (int i = 0; i < storageBuffers->size(); i++) {
            bufferInfos.emplace_back(
                    storageBuffers->at(i).getBufferHandle(),
                    0,
                    storageBuffers->at(i).getSize()
            );

            writeSets.emplace_back(
                    **set,
                    i + uniformBuffers->size(),
                    0,
                    1,
                    vk::DescriptorType::eStorageBuffer,
                    nullptr,
                    &bufferInfos[i + uniformBuffers->size()],
                    nullptr
            );
        }

        resources.dev->updateDescriptorSets(writeSets, nullptr);
    }

    /*void ComputeShader::createCommandBuffer(uint32_t computeQueueFamily) {
        vk::CommandPoolCreateInfo cmdPoolInfo(vk::CommandPoolCreateFlagBits::eResetCommandBuffer, computeQueueFamily);
        cmdPool = std::make_unique<vk::raii::CommandPool>(*resources.dev, cmdPoolInfo);

        vk::CommandBufferAllocateInfo cmdBufferInfo(**cmdPool, vk::CommandBufferLevel::ePrimary, 1);
        vk::raii::CommandBuffers cmdBuffers(*resources.dev, cmdBufferInfo);

        cmdBuffer = std::make_unique<vk::raii::CommandBuffer>(std::move(cmdBuffers.front()));
    }*/

    void ComputeShader::dispatch(const vk::raii::CommandBuffer &cmdBuffer, uint32_t groupCountX, uint32_t groupCountY, uint32_t groupCountZ) {
        vk::CommandBufferBeginInfo beginInfo({}, nullptr);
        cmdBuffer.begin(beginInfo);
        cmdBuffer.bindPipeline(vk::PipelineBindPoint::eCompute, **computePipeline);
        cmdBuffer.bindDescriptorSets(
                vk::PipelineBindPoint::eCompute,
                **pipelineLayout,
                0,
                **set,
                nullptr
        );
        cmdBuffer.dispatch(groupCountX, groupCountY, groupCountZ);
        cmdBuffer.end();

        vk::SubmitInfo submitInfo(
                {},
                {},
                *cmdBuffer,
                {}
        );

        /*vk::SubmitInfo submitInfo(
                (uint32_t) 0,
                {},
                {},
                (uint32_t) 1,
                &*cmdBuffer,
                (uint32_t) 0,
                {}
        );
        std::array<vk::SubmitInfo, 1> infos = {submitInfo};*/
        resources.computeQueue->submit(submitInfo);
        //resources.computeQueue->submit(infos);
        resources.computeQueue->waitIdle();
    }
}