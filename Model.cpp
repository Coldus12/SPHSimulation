#include "Model.hpp"
#include "VltavaFunctions.hpp"

#include <glm/gtc/matrix_transform.hpp>
#include <chrono>

Vltava::Model::Model(VulkanResources &resources) : resources(resources) {
    aspect = resources.extent.width / (float) resources.extent.height;

    loadModel("");
    loadShaders("shaders/vert.spv", "shaders/frag.spv");
    createPipeline();
}

void Vltava::Model::updateResources(const VulkanResources &res) {
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
    graphicsPipeline.reset();

    aspect = resources.extent.width / (float) resources.extent.height;

    // Creating new pipeline with the updated resources
    createPipeline();
}

void Vltava::Model::loadModel(const std::string &path) {
    vertices.reserve(3);
    vertices.push_back({{-0.5f, -0.5f},{0.0f, 0.0f, 1.0f}});
    vertices.push_back({{0.5f, -0.5f},{0.0f, 1.0f, 0.0f}});
    vertices.push_back({{0.5f, 0.5f},{0.0f, 0.0f, 1.0f}});
    // For index buffer showcase
    vertices.push_back({{-0.5f, 0.5f},{0.0f, 1.0f, 0.0f}});

    indices = {
            //0, 1, 2, 2, 3, 0
            0, 1, 3, 3, 1, 2
    };

    vk::DeviceSize bufferSize = sizeof(vertices[0]) * vertices.size();

    // Staging buffer creation
    //---------------------------------
    Buffer stagingBuffer(
            resources,
            bufferSize,
            vk::BufferUsageFlagBits::eVertexBuffer | vk::BufferUsageFlagBits::eTransferSrc,
            vk::MemoryPropertyFlagBits::eHostVisible | vk::MemoryPropertyFlagBits::eHostCoherent
    );

    stagingBuffer.writeToBuffer(vertices.data(), bufferSize);

    // Device local buffer creation
    //---------------------------------
    vertexBuffer = std::make_unique<Buffer>(
            resources,
            bufferSize,
            vk::BufferUsageFlagBits::eTransferDst | vk::BufferUsageFlagBits::eVertexBuffer,
            vk::MemoryPropertyFlagBits::eDeviceLocal
    );

    vertexBuffer->bind(0);

    // Copying data from staging buffer to local buffer
    //---------------------------------
    Buffer::copyBuffer(stagingBuffer.getBufferHandle(),vertexBuffer->getBufferHandle(),bufferSize);

    createIndexBuffer();
    createUniformBuffers();
    createDescriptors();
}

void Vltava::Model::createIndexBuffer() {
    vk::DeviceSize bufferSize = sizeof(indices[0]) * indices.size();

    // Staging buffer creation
    //---------------------------------
    Buffer stagingBuffer(
            resources,
            bufferSize,
            vk::BufferUsageFlagBits::eTransferSrc,
            vk::MemoryPropertyFlagBits::eHostVisible | vk::MemoryPropertyFlagBits::eHostCoherent
    );

    stagingBuffer.writeToBuffer(indices.data(), bufferSize);

    // Device local buffer creation
    //---------------------------------
    indexBuffer = std::make_unique<Buffer>(
            resources,
            bufferSize,
            vk::BufferUsageFlagBits::eTransferDst | vk::BufferUsageFlagBits::eIndexBuffer,
            vk::MemoryPropertyFlagBits::eDeviceLocal
    );

    indexBuffer->bind(0);

    // Copying data from staging buffer to local buffer
    //---------------------------------
    Buffer::copyBuffer(stagingBuffer.getBufferHandle(),indexBuffer->getBufferHandle(),bufferSize);
}

void Vltava::Model::createDescriptors() {
    std::cout << "createDescriptors" << std::endl;

    // Descriptor pool creation
    //---------------------------------
    std::vector<vk::DescriptorPoolSize> sizes = {
            {vk::DescriptorType::eUniformBuffer, 10}
    };

    vk::DescriptorPoolCreateInfo poolInfo(
            vk::DescriptorPoolCreateFlagBits::eFreeDescriptorSet,
            10,
            sizes
    );

    descriptorPool = std::make_unique<vk::raii::DescriptorPool>(*resources.dev, poolInfo);

    // DescriptorSetLayout creation
    //---------------------------------
    vk::DescriptorSetLayoutBinding binding(
            0,
            vk::DescriptorType::eUniformBuffer,
            1,
            vk::ShaderStageFlagBits::eVertex
    );

    vk::DescriptorSetLayoutCreateInfo layoutInfo(
            {},
            binding
    );

    setLayout = std::make_unique<vk::raii::DescriptorSetLayout>(*resources.dev, layoutInfo);

    // Creating descriptor sets
    //---------------------------------
    vk::DescriptorSetAllocateInfo allocInfo(
            **descriptorPool,
            1,
            &**setLayout
    );

    for (int i = 0; i < resources.FRAMES_IN_FLIGHT; i++) {
        descriptorSets.push_back(std::move(resources.dev->allocateDescriptorSets(allocInfo).front()));

        vk::DescriptorBufferInfo binfo(
                uniformBuffers[i].getBufferHandle(),
                0,
                sizeof(MVP)
        );

        vk::WriteDescriptorSet setWrite(
                *descriptorSets[i],
                0,
                {},
                1,
                vk::DescriptorType::eUniformBuffer,
                nullptr,
                &binfo,
                nullptr
        );

        resources.dev->updateDescriptorSets(setWrite, nullptr);
    }

}

void Vltava::Model::createUniformBuffers() {
    vk::DeviceSize bufferSize = sizeof(MVP);

    uniformBuffers.reserve(resources.FRAMES_IN_FLIGHT);

    for (size_t i = 0; i < resources.FRAMES_IN_FLIGHT; i++) {
        uniformBuffers.emplace_back(
                resources,
                bufferSize,
                vk::BufferUsageFlagBits::eUniformBuffer,
                vk::MemoryPropertyFlagBits::eHostVisible | vk::MemoryPropertyFlagBits::eHostCoherent
        );

        uniformBuffers[i].bind(0);
    }
}

//TODO the way the time measurement here is done is just dumb
void Vltava::Model::updateUniformBuffer(uint32_t currentImage) {
    static auto startTime = std::chrono::high_resolution_clock::now();

    auto currentTime = std::chrono::high_resolution_clock::now();
    float time = std::chrono::duration<float, std::chrono::seconds::period>(currentTime - startTime).count();

    MVP mvp{};
    mvp.model = glm::rotate(glm::mat4(1.0f), time * glm::radians(90.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    mvp.view = glm::lookAt(glm::vec3(2.0f, 2.0f, 2.0f), glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    mvp.proj = glm::perspective(glm::radians(45.0f), aspect, 0.1f, 10.0f);
    mvp.proj[1][1] *= -1;

    uniformBuffers[currentImage].writeToBuffer(&mvp, sizeof(mvp));
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

    vk::Viewport viewport(0.0f, 0.0f, (float) resources.extent.width, (float) resources.extent.height, 0.0f, 1.0f);
    vk::Rect2D scissor({0, 0}, resources.extent);
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
            &**setLayout,
            0,
            nullptr
    );

    pipelineLayout = std::make_unique<vk::raii::PipelineLayout>(*resources.dev, pipelineLayoutInfo);
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
            **pipelineLayout,
            **resources.renderPass,
            0
    );

    graphicsPipeline = std::make_unique<vk::raii::Pipeline>(*resources.dev, nullptr, pipelineInfo);
}

void Vltava::Model::draw(const vk::raii::CommandBuffer& cmdBuffer, uint32_t currentFrame) {
//void Vltava::Model::draw(const vk::raii::CommandBuffer& cmdBuffer) {
    updateUniformBuffer(currentFrame);
    cmdBuffer.bindPipeline(vk::PipelineBindPoint::eGraphics, **graphicsPipeline);

    cmdBuffer.bindDescriptorSets(
            vk::PipelineBindPoint::eGraphics,
            **pipelineLayout,
            0,
            {*descriptorSets[currentFrame]},
            nullptr
    );

    cmdBuffer.bindVertexBuffers(0, vertexBuffer->getBufferHandle(), {0});
    cmdBuffer.bindIndexBuffer(indexBuffer->getBufferHandle(), 0, vk::IndexType::eUint16); // Extra

    //cmdBuffer.draw(static_cast<uint32_t>(vertices.size()), 1, 0, 0); // --> changed to drawIndexed for indexed use-cases
    cmdBuffer.drawIndexed(static_cast<uint32_t>(indices.size()), 1, 0, 0, 0);
}