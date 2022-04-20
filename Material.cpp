//
// Created by PCF112021 on 3/28/2022.
//

#include "Material.hpp"

#include <utility>

namespace Vltava {
    Material::Material(VulkanResources &resources, std::string pathToVert, std::string pathToFrag) :
    res(resources), vertPath(std::move(pathToVert)), fragPath(std::move(pathToFrag)) {

        aspect = res.extent.width / (float) res.extent.height;
    }

    void Material::updateResources(const VulkanResources &resources) {
        // Updating the resources
        res.dev = resources.dev;
        res.extent = resources.extent;
        res.renderPass = resources.renderPass;
        res.physDev = resources.physDev;
        res.instance = resources.instance;
        res.commandPool = resources.commandPool;
        res.graphicsQueue = resources.graphicsQueue;
        res.computeQueue = resources.computeQueue;
        res.FRAMES_IN_FLIGHT = resources.FRAMES_IN_FLIGHT;

        aspect = res.extent.width / (float) resources.extent.height;
    }

    void Material::uploadIndexData(std::vector<uint16_t> indices) {
        vk::DeviceSize bufferSize = sizeof(uint16_t) * indices.size();
        nrOfIndices = indices.size();

        // Staging buffer creation
        //---------------------------------
        Buffer stagingBuffer(
                res,
                bufferSize,
                vk::BufferUsageFlagBits::eTransferSrc,
                vk::MemoryPropertyFlagBits::eHostVisible | vk::MemoryPropertyFlagBits::eHostCoherent
        );

        stagingBuffer.writeToBuffer(indices.data(), bufferSize);

        // Device local buffer creation
        //---------------------------------
        indexBuffer = std::make_unique<Buffer>(
                res,
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

        int bufferNrPerFrame = (uniformBuffers->size() / res.FRAMES_IN_FLIGHT);

        for (int i = 0; i < bufferNrPerFrame; i++) {
            bindings.emplace_back(i, vk::DescriptorType::eUniformBuffer, 1, shaderStages);
        }

        for (int i = bufferNrPerFrame; i < bufferNrPerFrame + storageBuffers->size(); i++) {
            bindings.emplace_back(i, vk::DescriptorType::eStorageBuffer, 1, shaderStages);
        }

        vk::DescriptorSetLayoutCreateInfo layoutInfo({}, bindings);
        setLayout = std::make_unique<MDescriptorSetLayout>(res.dev->getHandle(), layoutInfo);

        // Descriptor pool creation
        //---------------------------------
        std::vector<vk::DescriptorPoolSize> sizes;
        if (!uniformBuffers->empty())
            sizes.emplace_back(vk::DescriptorType::eUniformBuffer, static_cast<uint32_t>(uniformBuffers->size()) * res.FRAMES_IN_FLIGHT);

        if (!storageBuffers->empty())
            sizes.emplace_back(vk::DescriptorType::eStorageBuffer, static_cast<uint32_t>(storageBuffers->size()) * res.FRAMES_IN_FLIGHT);

        vk::DescriptorPoolCreateInfo poolInfo(
                vk::DescriptorPoolCreateFlagBits::eFreeDescriptorSet,
                20,
                sizes
        );

        setPool = std::make_unique<MDescriptorPool>(res.dev->getHandle(), poolInfo);

        // Descriptor set creation
        //---------------------------------
        vk::DescriptorSetAllocateInfo allocInfo(setPool->getHandle(), 1, setLayout->getAddress());
        for (int j = 0; j < res.FRAMES_IN_FLIGHT; j++) {
            sets.emplace_back(res.dev->getHandle(), res.dev->getHandle().allocateDescriptorSets(allocInfo).front(), setPool->getHandle());

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
                        sets[j].getHandle(),
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
                        sets[j].getHandle(),
                        i + bufferNrPerFrame,
                        0,
                        1,
                        vk::DescriptorType::eStorageBuffer,
                        nullptr,
                        &bufferInfos[i + bufferNrPerFrame],
                        nullptr
                );
            }

            res.dev->getHandle().updateDescriptorSets(writeSets, nullptr);
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

        int bufferNrPerFrame = (uniformBuffers->size() / res.FRAMES_IN_FLIGHT);

        for (int i = 0; i < bufferNrPerFrame; i++) {
            bindings.emplace_back(i, vk::DescriptorType::eUniformBuffer, 1, shaderStages);
        }

        for (int i = bufferNrPerFrame; i < bufferNrPerFrame + storageBuffers->size(); i++) {
            bindings.emplace_back(i, vk::DescriptorType::eStorageBuffer, 1, shaderStages);
        }

        vk::DescriptorSetLayoutCreateInfo layoutInfo({}, bindings);
        setLayout = std::make_unique<MDescriptorSetLayout>(res.dev->getHandle(), layoutInfo);

        // Descriptor pool creation
        //---------------------------------
        std::vector<vk::DescriptorPoolSize> sizes;
        if (!uniformBuffers->empty())
            sizes.emplace_back(vk::DescriptorType::eUniformBuffer, static_cast<uint32_t>(uniformBuffers->size()) * res.FRAMES_IN_FLIGHT);

        if (!storageBuffers->empty())
            sizes.emplace_back(vk::DescriptorType::eStorageBuffer, static_cast<uint32_t>(storageBuffers->size()) * res.FRAMES_IN_FLIGHT);

        vk::DescriptorPoolCreateInfo poolInfo(
                vk::DescriptorPoolCreateFlagBits::eFreeDescriptorSet,
                20,
                sizes
        );

        setPool = std::make_unique<MDescriptorPool>(res.dev->getHandle(), poolInfo);

        // Descriptor set creation
        //---------------------------------
        vk::DescriptorSetAllocateInfo allocInfo(setPool->getHandle(), 1, setLayout->getAddress());
        for (int j = 0; j < res.FRAMES_IN_FLIGHT; j++) {
            sets.emplace_back(res.dev->getHandle(), res.dev->getHandle().allocateDescriptorSets(allocInfo).front(), setPool->getHandle());

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
                        sets[j].getHandle(),
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
                        sets[j].getHandle(),
                        i + bufferNrPerFrame,
                        0,
                        1,
                        vk::DescriptorType::eStorageBuffer,
                        nullptr,
                        &bufferInfos[i + bufferNrPerFrame],
                        nullptr
                );
            }

            res.dev->getHandle().updateDescriptorSets(writeSets, nullptr);
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

        MShaderModule vertModule(res.dev->getHandle(), vertInfo);
        MShaderModule fragModule(res.dev->getHandle(), fragInfo);

        vk::PipelineShaderStageCreateInfo vertStageInfo({}, vk::ShaderStageFlagBits::eVertex, vertModule.getHandle(), "main");
        vk::PipelineShaderStageCreateInfo fragStageInfo({}, vk::ShaderStageFlagBits::eFragment, fragModule.getHandle(), "main");

        vk::PipelineShaderStageCreateInfo shaderStages[] = {vertStageInfo, fragStageInfo};

        vk::PipelineVertexInputStateCreateInfo vertexInputInfo(
                {},
                bindingDescriptions.size(),
                bindingDescriptions.data(),
                attributeDescriptions.size(),
                attributeDescriptions.data()
        );

        vk::PipelineInputAssemblyStateCreateInfo inputAssembly({}, topology, false);

        vk::Viewport viewport(0.0f, 0.0f, (float) res.extent.width, (float) res.extent.height, 0.0f, 1.0f);
        vk::Rect2D scissor({0, 0}, res.extent);
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
                false,
                vk::BlendFactor::eOne,
                vk::BlendFactor::eZero,
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
                setLayout->getAddress(),
                0,
                nullptr
        );

        pipelineLayout = std::make_unique<MPipelineLayout>(res.dev->getHandle(), pipelineLayoutInfo);
        //vk::raii::PipelineLayout pipelineLayout(*resources.dev, pipelineLayoutInfo);

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
                pipelineLayout->getHandle(),
                res.renderPass->getHandle(),
                0
        );

        pipeline = std::make_unique<MPipeline>(res.dev->getHandle(), pipelineInfo);
    }

    void Material::draw(const vk::CommandBuffer &cmdBuffer, uint32_t currentFrame) {
        cmdBuffer.bindPipeline(vk::PipelineBindPoint::eGraphics, pipeline->getHandle());

        cmdBuffer.bindDescriptorSets(
                vk::PipelineBindPoint::eGraphics,
                pipelineLayout->getHandle(),
                0,
                {sets[currentFrame].getHandle()},
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
        pipeline.reset();
        aspect = res.extent.width / (float) res.extent.height;
        createPipeline(bindings, attribs);
    }
}