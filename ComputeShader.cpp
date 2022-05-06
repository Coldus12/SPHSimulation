#include "ComputeShader.hpp"

namespace Vltava {

    ComputeShader::ComputeShader(std::string pathToShader) : computeShaderPath(pathToShader) {}

    ComputeShader::~ComputeShader() {
        VulkanResources::getInstance().logDev->getHandle().destroyDescriptorPool(setPool);
        VulkanResources::getInstance().logDev->getHandle().destroyDescriptorSetLayout(setLayout);
        VulkanResources::getInstance().logDev->getHandle().destroyPipelineLayout(pipelineLayout);
        VulkanResources::getInstance().logDev->getHandle().destroyPipeline(computePipeline);
    }

    void ComputeShader::recreatePipeline() {
        VulkanResources::getInstance().logDev->getHandle().destroyPipelineLayout(pipelineLayout);
        VulkanResources::getInstance().logDev->getHandle().destroyPipeline(computePipeline);
        createPipeline();
    }

    vk::Pipeline ComputeShader::getPipelineHandle() {
        return computePipeline;
    }

    vk::PipelineLayout ComputeShader::getPipelineLayoutHandle() {
        return pipelineLayout;
    }

    void ComputeShader::createPipeline() {
        auto code = readFile(computeShaderPath);
        vk::ShaderModuleCreateInfo compInfo({}, code.size(), reinterpret_cast<uint32_t*>(code.data()));
        vk::ShaderModule compModule = VulkanResources::getInstance().logDev->getHandle().createShaderModule(compInfo);
        vk::PipelineShaderStageCreateInfo compStageInfo({}, vk::ShaderStageFlagBits::eCompute, compModule, "main");
        vk::PipelineShaderStageCreateInfo stages[] = {compStageInfo};

        vk::PipelineLayoutCreateInfo layoutInfo(
                {},
                1,
                &setLayout,
                0,
                nullptr
        );

        pipelineLayout = VulkanResources::getInstance().logDev->getHandle().createPipelineLayout(layoutInfo);

        vk::Pipeline plHandle;
        vk::ComputePipelineCreateInfo computeInfo(
                {},
                compStageInfo,
                pipelineLayout,
                plHandle,
                0
        );

        computePipeline = VulkanResources::getInstance().logDev->getHandle().createComputePipeline({}, computeInfo).value;
        VulkanResources::getInstance().logDev->getHandle().destroyShaderModule(compModule);
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
        setLayout = VulkanResources::getInstance().logDev->getHandle().createDescriptorSetLayout(layoutInfo);

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

        setPool = VulkanResources::getInstance().logDev->getHandle().createDescriptorPool(poolInfo);

        // Descriptor set creation
        //---------------------------------
        vk::DescriptorSetAllocateInfo allocInfo(setPool, 1, &setLayout);
        set = VulkanResources::getInstance().logDev->getHandle().allocateDescriptorSets(allocInfo).front();

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
                    set,
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
                    set,
                    i + uniformBuffers->size(),
                    0,
                    1,
                    vk::DescriptorType::eStorageBuffer,
                    nullptr,
                    &bufferInfos[i + uniformBuffers->size()],
                    nullptr
            );
        }

        VulkanResources::getInstance().logDev->getHandle().updateDescriptorSets(writeSets, nullptr);
    }

    void ComputeShader::bindPipelineAndDescriptors(vk::CommandBuffer& cmdBuffer) {
        cmdBuffer.bindPipeline(vk::PipelineBindPoint::eCompute, computePipeline);
        cmdBuffer.bindDescriptorSets(
                vk::PipelineBindPoint::eCompute,
                pipelineLayout,
                0,
                set,
                nullptr
        );
    }
}