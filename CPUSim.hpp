#ifndef SPHSIMULATION_CPUSIM_HPP
#define SPHSIMULATION_CPUSIM_HPP

#include "ParticleModel.hpp"

namespace Vltava {
    struct Neighbourhood {
        glm::vec3 neighbour[27]; // 3^3
    };

    class CPUSim {
    public:
        CPUSim();
        void run(int iterNr);
        SimProps simProps;

        void printData(int nr = 64);
        void printGridData();
        void setSimProps();
        void setSimProps(SimProps& props);
        void setData(const std::vector<Particle>& p);

        void debugNeighbour();
    private:
        int cellx, celly, cellz;

        bool first = true;
        std::vector<int> grid_data;
        std::vector<Particle> particles1;
        std::vector<Particle> particles2;

        void place();
        void calculateRhoAndP();
        void iter();

        float kernel(int i, int j);
        glm::vec3 gradKernel(int i, int j);

        glm::vec3 determineGridTuple(int particleIdx);
        int getStartIdxOfCell(glm::vec3 tuple);
        Neighbourhood getNeighbouringCells(glm::vec3 cellTuple);
        bool checkBounds(glm::vec3 tuple);
        void placeParticleIntoCell(int particleIdx);
        void printNeighbourhood(Neighbourhood nh);
    };
}


#endif //SPHSIMULATION_CPUSIM_HPP
