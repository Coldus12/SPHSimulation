//
// Created by coldus on 3/8/22.
//

#include "Model.h"
#include "VltavaFunctions.h"

#define VMA_IMPLEMENTATION
#include "vma/vk_mem_alloc.hpp"

Vltava::Model::Model(PipelineCreationResources &resources) : resources(resources) {
    loadModel("");
    loadShaders("shaders/vert.spv", "shaders/frag.spv");
    createPipeline();
}

void Vltava::Model::updateResources(const PipelineCreationResources &res) {
    // Updating the resources
    resources.dev = res.dev;
    resources.extent = res.extent;
    resources.renderPass = res.renderPass;
    resources.physDev = res.physDev;

    // Deleting the out of date pipeline
    graphicsPipeline.reset();

    // Creating new pipeline with the updated resources
    createPipeline();
}

std::pair<vk::raii::Buffer, vk::raii::DeviceMemory> Vltava::Model::createBuffer(vk::DeviceSize bufferSize,
                                                                                vk::BufferUsageFlags usage,
                                                                                vk::MemoryPropertyFlags memFlags) {
    vk::BufferCreateInfo bufferInfo(
            {},
            bufferSize,
            usage,
            vk::SharingMode::eExclusive
    );

    vk::raii::Buffer localBuffer(*resources.dev, bufferInfo);

    vk::MemoryRequirements memReq = localBuffer.getMemoryRequirements();
    vk::PhysicalDeviceMemoryProperties memProperties = resources.physDev->getMemoryProperties();

    vk::MemoryAllocateInfo memAllocInfo(memReq.size, {});
    memAllocInfo.memoryTypeIndex = findMemoryType(
            memReq.memoryTypeBits,
            memProperties,
            memFlags
    );

    vk::raii::DeviceMemory devMem(*resources.dev, memAllocInfo);

    return std::pair<vk::raii::Buffer, vk::raii::DeviceMemory>(std::move(localBuffer), std::move(devMem));
}

void Vltava::Model::loadModel(const std::string &path) {
    vertices.reserve(3);
    vertices.push_back({{0.0f, -0.5f},{1.0f, 0.0f, 0.0f}});
    vertices.push_back({{0.5f, 0.5f},{0.0f, 1.0f, 0.0f}});
    vertices.push_back({{-0.5f, 0.5f},{0.0f, 0.0f, 1.0f}});

    vk::DeviceSize bufferSize = sizeof(vertices[0]) * vertices.size();

    auto hostlocal = createBuffer(
            bufferSize,
            vk::BufferUsageFlagBits::eVertexBuffer,
            vk::MemoryPropertyFlagBits::eHostVisible | vk::MemoryPropertyFlagBits::eHostCoherent
    );

    //vk::raii::Buffer stagingBuffer = std::move(hostlocal.first);
    //vk::raii::DeviceMemory stagingMemory = std::move(hostlocal.second);

    buffer = std::make_unique<vk::raii::Buffer>(std::move(hostlocal.first));
    vertexBuffer = std::make_unique<vk::raii::DeviceMemory>(std::move(hostlocal.second));

    buffer->bindMemory(**vertexBuffer, 0);
    void *data = vertexBuffer->mapMemory(0, VK_WHOLE_SIZE);
    memcpy(data, vertices.data(), bufferSize);
    vertexBuffer->unmapMemory();

    /*stagingBuffer.bindMemory(*stagingMemory, 0);
    void *data = stagingMemory.mapMemory(0, VK_WHOLE_SIZE);
    memcpy(data, vertices.data(), bufferSize);
    stagingMemory.unmapMemory();

    auto devicelocal = createBuffer(
            bufferSize,
            vk::BufferUsageFlagBits::eTransferDst | vk::BufferUsageFlagBits::eVertexBuffer,
            vk::MemoryPropertyFlagBits::eDeviceLocal
    );*/
}

void Vltava::Model::loadShaders(const std::string &vertShaderPath, const std::string &fragShaderPath) {
    this->vertShaderPath = vertShaderPath;
    this->fragShaderPath = fragShaderPath;
}

void Vltava::Model::createPipeline() {
    auto vertCode = readFile(vertShaderPath);
    auto fragCode = readFile(fragShaderPath);

    vk::ShaderModuleCreateInfo vertInfo({}, vertCode.size(), reinterpret_cast<uint32_t*>(vertCode.data()));
    vk::ShaderModuleCreateInfo fragInfo({}, fragCode.size(), reinterpret_cast<uint32_t*>(fragCode.data()));

    vk::raii::ShaderModule vertModule(*resources.dev, vertInfo);
    vk::raii::ShaderModule fragModule(*resources.dev, fragInfo);

    vk::PipelineShaderStageCreateInfo vertStageInfo({}, vk::ShaderStageFlagBits::eVertex, *vertModule, "main");
    vk::PipelineShaderStageCreateInfo fragStageInfo({}, vk::ShaderStageFlagBits::eFragment, *fragModule, "main");

    vk::PipelineShaderStageCreateInfo shaderStages[] = {vertStageInfo, fragStageInfo};

    auto bindingDescription = Vertex::getBindingDescription();
    auto attributeDescription = Vertex::getAttributeDescriptins();

    vk::PipelineVertexInputStateCreateInfo vertexInputInfo(
            {},
            1,
            &bindingDescription,
            static_cast<uint32_t>(attributeDescription.size()),
            attributeDescription.data()
    );

    vk::PipelineInputAssemblyStateCreateInfo inputAssembly({}, vk::PrimitiveTopology::eTriangleList, false);

    vk::Viewport viewport(0.0f, 0.0f, (float) resources.extent->width, (float) resources.extent->height, 0.0f, 1.0f);
    vk::Rect2D scissor({0, 0}, *resources.extent);
    vk::PipelineViewportStateCreateInfo viewportState({}, 1, &viewport, 1, &scissor);

    vk::PipelineRasterizationStateCreateInfo rasterizer(
            {},
            false,
            false,
            vk::PolygonMode::eFill,
            vk::CullModeFlagBits::eBack, vk::FrontFace::eClockwise,
            false, 0.0f, 0.0f, 0.0f, 1.0f
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
            0,
            nullptr,
            0,
            nullptr
    );

    //pipelineLayout = std::make_unique<vk::raii::PipelineLayout>(resources.dev, pipelineLayoutInfo);
    vk::raii::PipelineLayout pipelineLayout(*resources.dev, pipelineLayoutInfo);

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
            *pipelineLayout,
            **resources.renderPass,
            0
    );

    graphicsPipeline = std::make_unique<vk::raii::Pipeline>(*resources.dev, nullptr, pipelineInfo);
}

void Vltava::Model::draw(const vk::raii::CommandBuffer& cmdBuffer) {
    cmdBuffer.bindPipeline(vk::PipelineBindPoint::eGraphics, **graphicsPipeline);
    cmdBuffer.bindVertexBuffers(0, **buffer, {0});
    cmdBuffer.draw(static_cast<uint32_t>(vertices.size()), 1, 0, 0);
}