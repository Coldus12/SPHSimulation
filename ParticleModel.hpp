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

        float staticP = 0;
        float padding = 0;

        static std::vector<vk::VertexInputBindingDescription> getBindingDescription();
        static std::vector<vk::VertexInputAttributeDescription> getAttributeDescription();
    };

    // Properties of the simulation
    struct SimProps {
        float desired_density = 1000.0f;
        float k = 50.0f;                    // normalization constant / stiffness constant
        float nr_of_particles = 512.0f;
        float kernelh = 0.2f;

        float dt = 0.01f;
        float kPCI = 0.0f;
        float padding1 = 0.0f;
        float padding2 = 0.0f;

        glm::vec4 gridA = glm::vec4(1);
        glm::vec4 gridB = glm::vec4(-1);
    };

    class ParticleModel : public Model {
    public:
        ParticleModel(std::vector<Buffer>* simPropsBuffer, std::vector<Buffer>* storageBuffers);

        void loadModel(const std::string& path) override;
        void changeModel(int nrOfParticles, std::vector<Buffer>* simPropsBuffers, std::vector<Buffer>* storageBuffers);
        void draw(const vk::CommandBuffer& buffer, uint32_t currentFrame) override;
    private:
        int nr = 64;
        std::vector<glm::vec2> vertexBufferData;

        //Buffer* simProps = nullptr;
        std::vector<Buffer>* storageBuffers = nullptr;
        std::vector<Buffer>* simPropsBuffers = nullptr;
        std::vector<Buffer*> uniformBuffersU/*pdated*/;
        std::vector<Buffer*> storageBuffersU/*pdated*/;

        void updateUniformBuffer(uint32_t currentImage);

        // Helper functions
        //void createVertexBuffer();
    };
}


#endif //SPHSIMULATION_PARTICLEMODEL_HPP
