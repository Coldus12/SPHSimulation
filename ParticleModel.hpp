//
// Created by PCF112021 on 3/31/2022.
//

#ifndef SPHSIMULATION_PARTICLEMODEL_HPP
#define SPHSIMULATION_PARTICLEMODEL_HPP

#include "Model.hpp"

namespace Vltava {
    // Individual particle
    struct Particle {
        glm::vec3 x;
        float h;
        glm::vec3 v;
        float m;

        float rho;
        float p;

        float padding1 = 0;
        float padding2 = 0;

        static std::vector<vk::VertexInputBindingDescription> getBindingDescription();
        static std::vector<vk::VertexInputAttributeDescription> getAttributeDescription();
    };

    // Properties of the simulation
    struct SimProps {
        float desired_density;
        float k;                    // normalization constant / stiffness constant
        float nr_of_particles;
        float kernelh;
    };

    class ParticleModel : Model {
    public:
        ParticleModel(VulkanResources& resources, Buffer* simPropsBuffer, std::vector<Buffer>* storageBuffers);

        void loadModel(const std::string& path) override;
        void draw(const vk::raii::CommandBuffer& buffer, uint32_t currentFrame) override;
    private:
        Buffer* simProps = nullptr;
        std::vector<Buffer>* storageBuffers = nullptr;
        std::vector<Buffer*> uniformBuffersU/*pdated*/;
        std::vector<Buffer*> storageBuffersU/*pdated*/;

        void updateUniformBuffer(uint32_t currentImage);

        // Helper functions
        //void createVertexBuffer();
    };
}


#endif //SPHSIMULATION_PARTICLEMODEL_HPP
