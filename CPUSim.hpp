#ifndef SPHSIMULATION_CPUSIM_HPP
#define SPHSIMULATION_CPUSIM_HPP

#include "ParticleModel.hpp"

namespace Vltava {
    class CPUSim {
    public:
        CPUSim();
        void run(int iterNr);
        SimProps simProps;

        void printData();
        void setSimProps();
    private:
        int cellx, celly, cellz;

        bool first = true;
        std::vector<Particle> particles1;
        std::vector<Particle> particles2;

        void calculateRhoAndP();
        void iter();

        float kernel(int i, int j);
        glm::vec3 gradKernel(int i, int j);
    };
}


#endif //SPHSIMULATION_CPUSIM_HPP
