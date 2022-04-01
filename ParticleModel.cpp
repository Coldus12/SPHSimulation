//
// Created by PCF112021 on 3/31/2022.
//

#include <glm/gtc/matrix_transform.hpp>
#include <chrono>
#include "ParticleModel.hpp"

namespace Vltava {
    // Particle
    //------------------------------------------------------------------------------------------------------------------
    std::vector<vk::VertexInputBindingDescription> Particle::getBindingDescription() {
        std::vector<vk::VertexInputBindingDescription> bindings = {
                {0, sizeof(Particle), vk::VertexInputRate::eVertex}
        };

        return bindings;
    }

    std::vector<vk::VertexInputAttributeDescription> Particle::getAttributeDescription() {
        std::vector<vk::VertexInputAttributeDescription> attribs = {
                {0, 0, vk::Format::eR32G32B32Sfloat, offsetof(Particle, x)},
                {1, 0, vk::Format::eR32Sfloat, offsetof(Particle, h)},
                {2, 0, vk::Format::eR32G32B32Sfloat, offsetof(Particle, v)},
                {3, 0, vk::Format::eR32Sfloat, offsetof(Particle, m)},

                {4, 0, vk::Format::eR32Sfloat, offsetof(Particle, rho)},
                {5, 0, vk::Format::eR32Sfloat, offsetof(Particle, p)},

                {6, 0, vk::Format::eR32Sfloat, offsetof(Particle, padding1)},
                {7, 0, vk::Format::eR32Sfloat, offsetof(Particle, padding2)},
        };

        return attribs;
    }

    // ParticleModel
    //------------------------------------------------------------------------------------------------------------------
    ParticleModel::ParticleModel(VulkanResources& resources, Buffer* simPropsBuffer, std::vector<Buffer>* storageBuffers)
    : Model(resources), simProps(simPropsBuffer), storageBuffers(storageBuffers) {
        loadModel("");
    }

    void ParticleModel::loadModel(const std::string &path) {
        createUniformBuffers();
        mat = std::make_unique<Material>(resources, "shaders/vert_p.spv", "shaders/frag_p.spv");

        std::vector<Particle> vertexBufferData;
        vertexBufferData.resize(64 * sizeof(Particle));
        mat->uploadVertexData<Particle>(vertexBufferData);

        // Uniform buffers scissored together
        uniformBuffersU.reserve(uniformBuffers.size() + resources.FRAMES_IN_FLIGHT /* times provided buffers (1)*/);
        for (auto & uniformBuffer : uniformBuffers) {
            uniformBuffersU.push_back(&uniformBuffer);
        }
        for (int i = 0; i < resources.FRAMES_IN_FLIGHT; i++) {
            uniformBuffersU.push_back(simProps);
        }

        // Updating storageBuffers (vector<Buffer> to vector<Buffer*>)
        storageBuffersU.reserve(storageBuffers->size());
        for (auto & storageBuffer : *storageBuffers) {
            storageBuffersU.push_back(&storageBuffer);
        }

        mat->setBuffers(&uniformBuffersU, &storageBuffersU);
        mat->createPipeline(
                Particle::getBindingDescription(),
                Particle::getAttributeDescription(),
                vk::PrimitiveTopology::ePointList
        );
    }

    void ParticleModel::updateUniformBuffer(uint32_t currentImage) {
        static auto startTime = std::chrono::high_resolution_clock::now();

        auto currentTime = std::chrono::high_resolution_clock::now();
        float time = std::chrono::duration<float, std::chrono::seconds::period>(currentTime - startTime).count();

        MVP mvp{};
        mvp.model = glm::mat4(1.0f);
        mvp.view = glm::lookAt(glm::vec3(2.0f, 2.0f, 2.0f), glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
        mvp.proj = glm::perspective(glm::radians(45.0f), aspect, 0.1f, 10.0f);
        mvp.proj[1][1] *= -1;

        uniformBuffers[currentImage].writeToBuffer(&mvp, sizeof(mvp));
    }

    void ParticleModel::draw(const vk::raii::CommandBuffer &buffer, uint32_t currentFrame) {
        updateUniformBuffer(currentFrame);
        mat->draw(buffer, currentFrame);
    }
}