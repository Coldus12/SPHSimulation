//
// Created by PCF112021 on 3/31/2022.
//

#include <glm/gtc/matrix_transform.hpp>
#include <chrono>
#include "ParticleModel.hpp"
#include "StdWindow.hpp"

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
    ParticleModel::ParticleModel(VulkanResources& resources, std::vector<Buffer>* simPropsBuffers, std::vector<Buffer>* storageBuffers)
    : Model(resources), simPropsBuffers(simPropsBuffers), storageBuffers(storageBuffers) {
        loadModel("");
    }

    void ParticleModel::loadModel(const std::string &path) {
        createUniformBuffers();
        for (int i = 0; i < uniformBuffers.size(); i++)
            uniformBuffers[i].setSize(sizeof(MVP));

        mat = std::make_unique<Material>(resources, "shaders/vert_p.spv", "shaders/frag_p.spv");

        int nr = 64;
        vertexBufferData.reserve(nr * 4);
        indices.clear();
        for (uint16_t i = 0; i < nr; i++) {
            vertexBufferData.emplace_back(-0.2f, -0.2f);
            vertexBufferData.emplace_back(0.2f, -0.2f);
            vertexBufferData.emplace_back(0.2f, 0.2f);
            vertexBufferData.emplace_back(-0.2f, 0.2f);

            uint16_t increment = 4 * i;
            indices.push_back(0 + increment); indices.push_back(1 + increment); indices.push_back(3 + increment);
            indices.push_back(1 + increment); indices.push_back(2 + increment); indices.push_back(3 + increment);
        }
        mat->uploadVertexData<glm::vec2>(vertexBufferData);
        mat->uploadIndexData(indices);

        // Uniform buffers scissored together
        uniformBuffersU.reserve(uniformBuffers.size() + resources.FRAMES_IN_FLIGHT * simPropsBuffers->size());

        for (int f = 0; f < resources.FRAMES_IN_FLIGHT; f++) {
            uniformBuffersU.push_back(&uniformBuffers[f]);
            for (auto & simPropsBuffer : *simPropsBuffers)
                uniformBuffersU.push_back(&simPropsBuffer);
        }

        // Updating storageBuffers (vector<Buffer> to vector<Buffer*>)
        storageBuffersU.reserve(storageBuffers->size());
        for (auto & storageBuffer : *storageBuffers) {
            storageBuffersU.push_back(&storageBuffer);
        }

        mat->setBuffers(&uniformBuffersU, &storageBuffersU);

        std::vector<vk::VertexInputBindingDescription> bindings = {
                {0, sizeof(glm::vec2), vk::VertexInputRate::eVertex}
        };

        std::vector<vk::VertexInputAttributeDescription> attribs = {
                {0, 0, vk::Format::eR32G32Sfloat, 0}
        };

        mat->createPipeline(
                bindings,
                attribs
        );
    }

    void ParticleModel::updateUniformBuffer(uint32_t currentImage) {
        static auto startTime = std::chrono::high_resolution_clock::now();

        auto currentTime = std::chrono::high_resolution_clock::now();
        float time = std::chrono::duration<float, std::chrono::seconds::period>(currentTime - startTime).count();

        MVP mvp{};
        //mvp.model = glm::mat4(1.0f);

        if (StdWindow::rot) {
            mvp.model = glm::rotate(glm::mat4(1.0f), time * glm::radians(90.0f), glm::vec3(0.0f, 0.0f, 1.0f));
        } else {
            mvp.model = glm::identity<glm::mat4>();
        }
        mvp.view = glm::lookAt(glm::vec3(0.0f, 10.0f, 10.0f), glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
        mvp.proj = glm::perspective(glm::radians(45.0f), aspect, 0.1f, 40.0f);
        mvp.proj[1][1] *= -1;

        mvp.localPos1 = vertexBufferData[0];
        mvp.localPos2 = vertexBufferData[1];
        mvp.localPos3 = vertexBufferData[2];
        mvp.localPos4 = vertexBufferData[3];

        uniformBuffers[currentImage].writeToBuffer(&mvp, sizeof(mvp));
    }

    void ParticleModel::draw(const vk::CommandBuffer &buffer, uint32_t currentFrame) {
        updateUniformBuffer(currentFrame);
        mat->draw(buffer, currentFrame);
    }
}