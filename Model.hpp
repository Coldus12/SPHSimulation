#ifndef SPHSIMULATION_MODEL_HPP
#define SPHSIMULATION_MODEL_HPP

#include <iostream>
#include <fstream>
#include <glm/glm.hpp>
#include "vulkan/vulkan.hpp"

#include "Vertex.hpp"
#include "VltavaFunctions.hpp"
#include "Buffer.hpp"
#include "Material.hpp"

namespace Vltava {
    struct MVP {
        glm::mat4 model;
        glm::mat4 view;
        glm::mat4 proj;
        glm::vec2 localPos1;
        glm::vec2 localPos2;
        glm::vec2 localPos3;
        glm::vec2 localPos4;
    };

    class Model {
    public:
        Model() = delete;
        Model(VulkanResources& resources);

        virtual void loadModel(const std::string& path);
        void updateResources(const VulkanResources& res);
        virtual void draw(const vk::CommandBuffer& buffer, uint32_t currentFrame);

        VulkanResources& resources;
    protected:
        std::unique_ptr<Material> mat;
        float aspect;

        std::vector<Vertex> vertices;
        std::vector<uint16_t> indices;

        // Vertex buffer
        std::unique_ptr<Buffer> vertexBuffer;

        // Uniform buffers
        std::vector<Buffer> uniformBuffers;

        void createUniformBuffers();
        void updateUniformBuffer(uint32_t currentImage);
    };
}


#endif //SPHSIMULATION_MODEL_HPP
