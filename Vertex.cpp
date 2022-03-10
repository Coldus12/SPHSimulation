//
// Created by PCF112021 on 3/8/2022.
//

#include "Vertex.h"
namespace Vltava {
    vk::VertexInputBindingDescription Vertex::getBindingDescription() {
        vk::VertexInputBindingDescription bindingDescription(0, sizeof(Vertex), vk::VertexInputRate::eVertex);
        return bindingDescription;
    }

    std::array<vk::VertexInputAttributeDescription, 2> Vertex::getAttributeDescriptins() {
        std::array<vk::VertexInputAttributeDescription, 2> attributeDescriptions = {
                vk::VertexInputAttributeDescription(0, 0, vk::Format::eR32G32Sfloat, offsetof(Vertex, pos)),
                vk::VertexInputAttributeDescription(1, 0, vk::Format::eR32G32B32Sfloat, offsetof(Vertex, color))
        };

        return attributeDescriptions;
    }
}