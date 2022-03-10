#ifndef SPHSIMULATION_VERTEX_HPP
#define SPHSIMULATION_VERTEX_HPP

#include <glm/glm.hpp>
#include "vulkan/vulkan.hpp"
#include "vulkan/vulkan_raii.hpp"

namespace Vltava {
    struct Vertex {
        glm::vec2 pos;
        glm::vec3 color;

        static vk::VertexInputBindingDescription getBindingDescription();
        static std::array<vk::VertexInputAttributeDescription, 2> getAttributeDescriptins();
    };
}


#endif //SPHSIMULATION_VERTEX_HPP
