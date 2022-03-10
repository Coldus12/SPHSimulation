#ifndef SPHSIMULATION_MODEL_HPP
#define SPHSIMULATION_MODEL_HPP

#include <iostream>
#include <fstream>
#include <glm/glm.hpp>
#include "vulkan/vulkan.hpp"
#include "vulkan/vulkan_raii.hpp"

#include "Vertex.hpp"

namespace Vltava {

    struct PipelineCreationResources {
        vk::raii::RenderPass* renderPass;
        vk::Extent2D* extent;
        vk::raii::PhysicalDevice* physDev;
        vk::raii::Device* dev;
    };

    class Model {
    public:
        Model() = delete;
        Model(PipelineCreationResources& resources);

        virtual void loadModel(const std::string& path);
        virtual void loadShaders(const std::string& vertShader, const std::string& fragShader);
        virtual void draw(const vk::raii::CommandBuffer& buffer);
        virtual void updateResources(const PipelineCreationResources& res);
    private:
        PipelineCreationResources& resources;

        std::string vertShaderPath;
        std::string fragShaderPath;

        std::vector<Vertex> vertices;

        std::unique_ptr<vk::raii::Pipeline> graphicsPipeline;
        std::unique_ptr<vk::raii::Buffer> buffer;
        std::unique_ptr<vk::raii::DeviceMemory> vertexBuffer;

        std::pair<vk::raii::Buffer, vk::raii::DeviceMemory> createBuffer(vk::DeviceSize bufferSize,
                                                                         vk::BufferUsageFlags usage,
                                                                         vk::MemoryPropertyFlags memFlags);
        virtual void createPipeline();
    };
}


#endif //SPHSIMULATION_MODEL_HPP
