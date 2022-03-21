#ifndef SPHSIMULATION_MODEL_HPP
#define SPHSIMULATION_MODEL_HPP

#include <iostream>
#include <fstream>
#include <glm/glm.hpp>
#include "vulkan/vulkan.hpp"
#include "vulkan/vulkan_raii.hpp"

#include "Vertex.hpp"
#include "VltavaFunctions.hpp"
#include "Buffer.hpp"

namespace Vltava {
    struct MVP {
        glm::mat4 model;
        glm::mat4 view;
        glm::mat4 proj;
    };

    class Model {
    public:
        Model() = delete;
        Model(VulkanResources& resources);

        virtual void loadModel(const std::string& path);
        virtual void loadShaders(const std::string& vertShader, const std::string& fragShader);
        virtual void draw(const vk::raii::CommandBuffer& buffer, uint32_t currentFrame);
        //virtual void draw(const vk::raii::CommandBuffer& buffer);
        virtual void updateResources(const VulkanResources& res);

        VulkanResources& resources;
    private:
        float aspect;

        std::string vertShaderPath;
        std::string fragShaderPath;

        std::vector<Vertex> vertices;
        std::vector<uint16_t> indices;

        std::unique_ptr<vk::raii::PipelineLayout> pipelineLayout;
        std::unique_ptr<vk::raii::Pipeline> graphicsPipeline;

        // Vertex buffer
        std::unique_ptr<Buffer> vertexBuffer;

        // Index buffer
        std::unique_ptr<Buffer> indexBuffer;

        // Uniform buffers
        std::vector<Buffer> uniformBuffers;
        std::unique_ptr<vk::raii::DescriptorPool> descriptorPool;
        std::unique_ptr<vk::raii::DescriptorSetLayout> setLayout;
        std::vector<vk::raii::DescriptorSet> descriptorSets;

        void createDescriptors();

        void createIndexBuffer();
        void createUniformBuffers();
        void updateUniformBuffer(uint32_t currentImage);
        virtual void createPipeline();
    };
}


#endif //SPHSIMULATION_MODEL_HPP
