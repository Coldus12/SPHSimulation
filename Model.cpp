#include "Model.hpp"
#include "VltavaFunctions.hpp"

#include <glm/gtc/matrix_transform.hpp>
#include <chrono>

Vltava::Model::Model(VulkanResources &resources) : resources(resources) {
    aspect = resources.extent.width / (float) resources.extent.height;
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

    aspect = resources.extent.width / (float) resources.extent.height;
    mat->recreatePipeline();
}

void Vltava::Model::loadModel(const std::string &path) {
    createUniformBuffers();
    for (int i = 0; i < uniformBuffers.size(); i++)
        uniformBuffers[i].setSize(sizeof(MVP));

    vertices.reserve(3);
    /*vertices.push_back({{-0.5f, -0.5f},{0.0f, 0.0f, 1.0f}});
    vertices.push_back({{0.5f, -0.5f},{0.0f, 1.0f, 0.0f}});
    vertices.push_back({{0.5f, 0.5f},{0.0f, 0.0f, 1.0f}});
    // For index buffer showcase
    vertices.push_back({{-0.5f, 0.5f},{0.0f, 1.0f, 0.0f}});*/

    vertices.push_back({{-0.2f, -0.2f},{0.0f, 0.0f, 1.0f}});
    vertices.push_back({{0.2f, -0.2f},{0.0f, 1.0f, 0.0f}});
    vertices.push_back({{0.2f, 0.2f},{0.0f, 0.0f, 1.0f}});
    // For index buffer showcase
    vertices.push_back({{-0.2f, 0.2f},{0.0f, 1.0f, 0.0f}});

    /*vertices.push_back({{0.4f, -0.2f},{0.0f, 0.0f, 1.0f}});
    vertices.push_back({{0.8f, -0.2f},{0.0f, 1.0f, 0.0f}});
    vertices.push_back({{0.8f, 0.2f},{0.0f, 0.0f, 1.0f}});
    // For index buffer showcase
    vertices.push_back({{0.4f, 0.2f},{0.0f, 1.0f, 0.0f}});*/

    indices = {
            //0, 1, 2, 2, 3, 0
            0, 1, 3, 1, 2, 3

            //,4, 5, 7, 7, 5, 6
    };

    mat = std::make_unique<Material>(resources, "shaders/vert.spv", "shaders/frag.spv");
    mat->uploadVertexData<Vertex>(vertices);
    mat->uploadIndexData(indices);
    mat->setBuffers(&uniformBuffers, nullptr, vk::ShaderStageFlagBits::eVertex | vk::ShaderStageFlagBits::eFragment);

    std::vector<vk::VertexInputBindingDescription> bindings = {
            vertices[0].getBindingDescription()
    };

    std::vector<vk::VertexInputAttributeDescription> attribs = {
            vertices[0].getAttributeDescriptins()[0],
            vertices[0].getAttributeDescriptins()[1]
    };

    mat->createPipeline(bindings, attribs);
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

    mvp.localPos1 = vertices[0].pos;
    mvp.localPos2 = vertices[1].pos;
    mvp.localPos3 = vertices[2].pos;
    mvp.localPos4 = vertices[3].pos;

    /*mvp.view = glm::mat4(1.0f);
    mvp.proj = glm::mat4(1.0f);
    mvp.model = glm::mat4(1.0f);
    mvp.proj[1][1] *= -1;*/

    uniformBuffers[currentImage].writeToBuffer(&mvp, sizeof(mvp));
}

void Vltava::Model::draw(const vk::raii::CommandBuffer& cmdBuffer, uint32_t currentFrame) {
    updateUniformBuffer(currentFrame);
    mat->draw(cmdBuffer, currentFrame);
}