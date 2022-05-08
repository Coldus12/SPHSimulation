//
// Created by PCF112021 on 3/28/2022.
//

#include "Material.hpp"

#include <utility>

namespace Vltava {
    Material::Material(std::string pathToVert, std::string pathToFrag) :
    vertPath(std::move(pathToVert)), fragPath(std::move(pathToFrag)) {

        aspect = VulkanResources::getInstance().extent.width / (float) VulkanResources::getInstance().extent.height;
    }

    Material::~Material() {
        cleanup();
    }

    void Material::uploadIndexData(std::vector<uint16_t> indices) {
        indexBuffer.reset();

        vk::DeviceSize bufferSize = sizeof(uint16_t) * indices.size();
        nrOfIndices = indices.size();

        // Staging buffer creation
        //---------------------------------
        Buffer stagingBuffer(
                bufferSize,
                vk::BufferUsageFlagBits::eTransferSrc,
                vk::MemoryPropertyFlagBits::eHostVisible | vk::MemoryPropertyFlagBits::eHostCoherent
        );

        stagingBuffer.writeToBuffer(indices.data(), bufferSize);

        // Device local buffer creation
        //---------------------------------
        indexBuffer = std::make_unique<Buffer>(
                bufferSize,
                vk::BufferUsageFlagBits::eTransferDst | vk::BufferUsageFlagBits::eIndexBuffer,
                vk::MemoryPropertyFlagBits::eDeviceLocal
        );

        indexBuffer->bind(0);

        // Copying data from staging buffer to local buffer
        //---------------------------------
        Buffer::copyBuffer(stagingBuffer.getBufferHandle(),indexBuffer->getBufferHandle(),bufferSize);
    }

    /**
     * Vector of pointers version.
     *
     * This version stores the pointers to the buffers and as a result it is now easier "scissor" the needed /
     * required vectors together.
     *
     * @uniformBuffers This vector should contain FRAMES_IN_FLIGHT * (UBO(s) in the shader) nr of buffers.
     * For example if the shader has 3 UBOs and 2 FRAMES_IN_FLIGHT, then this vector should contain
     * 6 (UBO * FRAMES_IN_FLIGHT) buffers. Each bufferNr%FRAMES_IN_FLIGHT has the same binding.
     * In the case of our example buffer nr 0 and buffer nr 3 has the same binding, and so do nr 1 - nr 4,
     * and nr 2 - nr 5.
     * */
    void Material::setBuffers(const std::vector<Buffer*>* uniformBuffers,
                              const std::vector<Buffer*>* storageBuffers,
                              vk::ShaderStageFlags shaderStages) {

        std::vector<Buffer*> placeHolder;
        if (uniformBuffers == nullptr)
            uniformBuffers = &placeHolder;

        if (storageBuffers == nullptr)
            storageBuffers = &placeHolder;

        // Descriptor set layout creation
        //---------------------------------
        std::vector<vk::DescriptorSetLayoutBinding> bindings;
        bindings.reserve(storageBuffers->size() + uniformBuffers->size());

        int bufferNrPerFrame = (uniformBuffers->size() / VulkanResources::getInstance().FRAMES_IN_FLIGHT);

        for (int i = 0; i < bufferNrPerFrame; i++) {
            bindings.emplace_back(i, vk::DescriptorType::eUniformBuffer, 1, shaderStages);
        }

        for (int i = bufferNrPerFrame; i < bufferNrPerFrame + storageBuffers->size(); i++) {
            bindings.emplace_back(i, vk::DescriptorType::eStorageBuffer, 1, shaderStages);
        }

        vk::DescriptorSetLayoutCreateInfo layoutInfo({}, bindings);
        setLayout = VulkanResources::getInstance().logDev->getHandle().createDescriptorSetLayout(layoutInfo);

        // Descriptor pool creation
        //---------------------------------
        std::vector<vk::DescriptorPoolSize> sizes;
        if (!uniformBuffers->empty())
            sizes.emplace_back(
                    vk::DescriptorType::eUniformBuffer,
                    static_cast<uint32_t>(uniformBuffers->size()) * VulkanResources::getInstance().FRAMES_IN_FLIGHT
            );

        if (!storageBuffers->empty())
            sizes.emplace_back(
                    vk::DescriptorType::eStorageBuffer,
                    static_cast<uint32_t>(storageBuffers->size()) * VulkanResources::getInstance().FRAMES_IN_FLIGHT
            );

        vk::DescriptorPoolCreateInfo poolInfo(
                vk::DescriptorPoolCreateFlagBits::eFreeDescriptorSet,
                20,
                sizes
        );

        setPool = VulkanResources::getInstance().logDev->getHandle().createDescriptorPool(poolInfo);

        // Descriptor set creation
        //---------------------------------
        vk::DescriptorSetAllocateInfo allocInfo(setPool, 1, &setLayout);
        for (int j = 0; j < VulkanResources::getInstance().FRAMES_IN_FLIGHT; j++) {
            sets.push_back(VulkanResources::getInstance().logDev->getHandle().allocateDescriptorSets(allocInfo).front());

            std::vector<vk::DescriptorBufferInfo> bufferInfos;
            std::vector<vk::WriteDescriptorSet> writeSets;
            bufferInfos.reserve(storageBuffers->size() + bufferNrPerFrame);
            writeSets.reserve(storageBuffers->size() + bufferNrPerFrame);

            for (int i = j * bufferNrPerFrame; i < (j+1) * bufferNrPerFrame; i++) {
                bufferInfos.emplace_back(
                        uniformBuffers->at(i)->getBufferHandle(),
                        0,
                        uniformBuffers->at(i)->getSize()
                );

                writeSets.emplace_back(
                        sets[j],
                        i % bufferNrPerFrame,
                        0,
                        1,
                        vk::DescriptorType::eUniformBuffer,
                        nullptr,
                        &bufferInfos[i % bufferNrPerFrame],
                        nullptr
                );
            }

            for (int i = 0; i < storageBuffers->size(); i++) {
                bufferInfos.emplace_back(
                        storageBuffers->at(i)->getBufferHandle(),
                        0,
                        storageBuffers->at(i)->getSize()
                );

                writeSets.emplace_back(
                        sets[j],
                        i + bufferNrPerFrame,
                        0,
                        1,
                        vk::DescriptorType::eStorageBuffer,
                        nullptr,
                        &bufferInfos[i + bufferNrPerFrame],
                        nullptr
                );
            }

            VulkanResources::getInstance().logDev->getHandle().updateDescriptorSets(writeSets, nullptr);
        }
    }

    /** @uniformBuffers This vector should contain FRAMES_IN_FLIGHT * (UBO(s) in the shader) nr of buffers.
     * For example if the shader has 3 UBOs and 2 FRAMES_IN_FLIGHT, then this vector should contain
     * 6 (UBO * FRAMES_IN_FLIGHT) buffers. Each bufferNr%FRAMES_IN_FLIGHT has the same binding.
     * In the case of our example buffer nr 0 and buffer nr 3 has the same binding, and so do nr 1 - nr 4,
     * and nr 2 - nr 5.
     * */
    void Material::setBuffers(const std::vector<Buffer>* uniformBuffers,
                              const std::vector<Buffer>* storageBuffers,
                              vk::ShaderStageFlags shaderStages) {

        std::vector<Buffer> placeHolder;
        if (uniformBuffers == nullptr)
            uniformBuffers = &placeHolder;

        if (storageBuffers == nullptr)
            storageBuffers = &placeHolder;

        // Descriptor set layout creation
        //---------------------------------
        std::vector<vk::DescriptorSetLayoutBinding> bindings;
        bindings.reserve(storageBuffers->size() + uniformBuffers->size());

        int bufferNrPerFrame = (uniformBuffers->size() / VulkanResources::getInstance().FRAMES_IN_FLIGHT);

        for (int i = 0; i < bufferNrPerFrame; i++) {
            bindings.emplace_back(i, vk::DescriptorType::eUniformBuffer, 1, shaderStages);
        }

        for (int i = bufferNrPerFrame; i < bufferNrPerFrame + storageBuffers->size(); i++) {
            bindings.emplace_back(i, vk::DescriptorType::eStorageBuffer, 1, shaderStages);
        }

        vk::DescriptorSetLayoutCreateInfo layoutInfo({}, bindings);
        setLayout = VulkanResources::getInstance().logDev->getHandle().createDescriptorSetLayout(layoutInfo);

        // Descriptor pool creation
        //---------------------------------
        std::vector<vk::DescriptorPoolSize> sizes;
        if (!uniformBuffers->empty())
            sizes.emplace_back(
                    vk::DescriptorType::eUniformBuffer,
                    static_cast<uint32_t>(uniformBuffers->size()) * VulkanResources::getInstance().FRAMES_IN_FLIGHT
            );

        if (!storageBuffers->empty())
            sizes.emplace_back(
                    vk::DescriptorType::eStorageBuffer,
                    static_cast<uint32_t>(storageBuffers->size()) * VulkanResources::getInstance().FRAMES_IN_FLIGHT
            );

        vk::DescriptorPoolCreateInfo poolInfo(
                vk::DescriptorPoolCreateFlagBits::eFreeDescriptorSet,
                20,
                sizes
        );

        setPool = VulkanResources::getInstance().logDev->getHandle().createDescriptorPool(poolInfo);

        // Descriptor set creation
        //---------------------------------
        vk::DescriptorSetAllocateInfo allocInfo(setPool, 1, &setLayout);
        for (int j = 0; j < VulkanResources::getInstance().FRAMES_IN_FLIGHT; j++) {
            sets.push_back(VulkanResources::getInstance().logDev->getHandle().allocateDescriptorSets(allocInfo).front());

            std::vector<vk::DescriptorBufferInfo> bufferInfos;
            std::vector<vk::WriteDescriptorSet> writeSets;
            bufferInfos.reserve(storageBuffers->size() + bufferNrPerFrame);
            writeSets.reserve(storageBuffers->size() + bufferNrPerFrame);

            for (int i = j * bufferNrPerFrame; i < (j+1) * bufferNrPerFrame; i++) {
                bufferInfos.emplace_back(
                        uniformBuffers->at(i).getBufferHandle(),
                        0,
                        uniformBuffers->at(i).getSize()
                );

                writeSets.emplace_back(
                        sets[j],
                        i % bufferNrPerFrame,
                        0,
                        1,
                        vk::DescriptorType::eUniformBuffer,
                        nullptr,
                        &bufferInfos[i % bufferNrPerFrame],
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
                        sets[j],
                        i + bufferNrPerFrame,
                        0,
                        1,
                        vk::DescriptorType::eStorageBuffer,
                        nullptr,
                        &bufferInfos[i + bufferNrPerFrame],
                        nullptr
                );
            }

            VulkanResources::getInstance().logDev->getHandle().updateDescriptorSets(writeSets, nullptr);
        }
    }

    void Material::createPipeline(std::vector<vk::VertexInputBindingDescription> bindingDescriptions,
                                  std::vector<vk::VertexInputAttributeDescription> attributeDescriptions,
                                  vk::PrimitiveTopology topology) {
        bindings = bindingDescriptions;
        attribs = attributeDescriptions;

        auto vertCode = readFile(vertPath);
        auto fragCode = readFile(fragPath);

        vk::ShaderModuleCreateInfo vertInfo({}, vertCode.size(), reinterpret_cast<uint32_t*>(vertCode.data()));
        vk::ShaderModuleCreateInfo fragInfo({}, fragCode.size(), reinterpret_cast<uint32_t*>(fragCode.data()));

        vk::ShaderModule vertModule = VulkanResources::getInstance().logDev->getHandle().createShaderModule(vertInfo);
        vk::ShaderModule fragModule = VulkanResources::getInstance().logDev->getHandle().createShaderModule(fragInfo);

        vk::PipelineShaderStageCreateInfo vertStageInfo({}, vk::ShaderStageFlagBits::eVertex, vertModule, "main");
        vk::PipelineShaderStageCreateInfo fragStageInfo({}, vk::ShaderStageFlagBits::eFragment, fragModule, "main");

        vk::PipelineShaderStageCreateInfo shaderStages[] = {vertStageInfo, fragStageInfo};

        vk::PipelineVertexInputStateCreateInfo vertexInputInfo(
                {},
                bindingDescriptions.size(),
                bindingDescriptions.data(),
                attributeDescriptions.size(),
                attributeDescriptions.data()
        );

        vk::PipelineInputAssemblyStateCreateInfo inputAssembly({}, topology, false);

        vk::Viewport viewport(
                0.0f,
                0.0f,
                (float) VulkanResources::getInstance().extent.width,
                (float) VulkanResources::getInstance().extent.height,
                0.0f,
                1.0f
        );
        vk::Rect2D scissor({0, 0}, VulkanResources::getInstance().extent);
        vk::PipelineViewportStateCreateInfo viewportState({}, 1, &viewport, 1, &scissor);

        vk::PipelineRasterizationStateCreateInfo rasterizer(
                {},
                false,
                false,
                vk::PolygonMode::eFill,
                vk::CullModeFlagBits::eBack,
                vk::FrontFace::eCounterClockwise,
                false,
                0.0f,
                0.0f,
                0.0f,
                1.0f
        );

        vk::PipelineMultisampleStateCreateInfo multisampling(
                {},
                vk::SampleCountFlagBits::e1,
                false, 1.0f,
                nullptr,
                false,
                false
        );

        vk::PipelineColorBlendAttachmentState colorBlendAttachment(
                true,
                vk::BlendFactor::eSrcAlpha,
                vk::BlendFactor::eOneMinusSrcAlpha,
                vk::BlendOp::eAdd,
                vk::BlendFactor::eOne,
                vk::BlendFactor::eZero,
                vk::BlendOp::eAdd
        );

        colorBlendAttachment.colorWriteMask = {
                vk::ColorComponentFlagBits::eR |
                vk::ColorComponentFlagBits::eG |
                vk::ColorComponentFlagBits::eB |
                vk::ColorComponentFlagBits::eA
        };

        vk::PipelineColorBlendStateCreateInfo colorBlending(
                {},
                false,
                vk::LogicOp::eCopy,
                1,
                &colorBlendAttachment,
                {0.0f, 0.0f, 0.0f, 0.0f}
        );

        vk::PipelineLayoutCreateInfo pipelineLayoutInfo(
                {},
                1,
                &setLayout,
                0,
                nullptr
        );

        pipelineLayout = VulkanResources::getInstance().logDev->getHandle().createPipelineLayout(pipelineLayoutInfo);

        vk::GraphicsPipelineCreateInfo pipelineInfo(
                {},
                2,
                shaderStages,
                &vertexInputInfo,
                &inputAssembly,
                nullptr,
                &viewportState,
                &rasterizer,
                &multisampling,
                nullptr,
                &colorBlending,
                nullptr,
                pipelineLayout,
                *VulkanResources::getInstance().renderPass,
                0
        );

        pipeline = VulkanResources::getInstance().logDev->getHandle().createGraphicsPipeline({}, pipelineInfo).value;

        VulkanResources::getInstance().logDev->getHandle().destroyShaderModule(vertModule);
        VulkanResources::getInstance().logDev->getHandle().destroyShaderModule(fragModule);
    }

    void Material::draw(const vk::CommandBuffer &cmdBuffer, uint32_t currentFrame) {
        cmdBuffer.bindPipeline(vk::PipelineBindPoint::eGraphics, pipeline);

        cmdBuffer.bindDescriptorSets(
                vk::PipelineBindPoint::eGraphics,
                pipelineLayout,
                0,
                {sets[currentFrame]},
                nullptr
        );

        cmdBuffer.bindVertexBuffers(0, vertexBuffer->getBufferHandle(), {0});
        if (nrOfIndices != 0)
            cmdBuffer.bindIndexBuffer(indexBuffer->getBufferHandle(), 0, vk::IndexType::eUint16); // Extra

        if (nrOfIndices == 0)
            cmdBuffer.draw(nrOfVertices, 1, 0, 0); // --> changed to drawIndexed for indexed use-cases
        else
            cmdBuffer.drawIndexed(nrOfIndices, 1, 0, 0, 0);

    }

    void Material::recreatePipeline() {
        VulkanResources::getInstance().logDev->getHandle().destroyPipeline(pipeline);
        aspect = VulkanResources::getInstance().extent.width / (float) VulkanResources::getInstance().extent.height;
        createPipeline(bindings, attribs);
    }

    void Material::cleanup() {
        VulkanResources::getInstance().logDev->getHandle().freeDescriptorSets(setPool, sets);
        VulkanResources::getInstance().logDev->getHandle().destroyDescriptorPool(setPool);
        VulkanResources::getInstance().logDev->getHandle().destroyDescriptorSetLayout(setLayout);
        VulkanResources::getInstance().logDev->getHandle().destroyPipelineLayout(pipelineLayout);
        VulkanResources::getInstance().logDev->getHandle().destroyPipeline(pipeline);
    }
}