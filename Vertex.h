//
// Created by PCF112021 on 3/8/2022.
//

#ifndef SPHSIMULATION_VERTEX_H
#define SPHSIMULATION_VERTEX_H

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


#endif //SPHSIMULATION_VERTEX_H
