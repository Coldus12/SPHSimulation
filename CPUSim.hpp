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
        void runSESPH(int iterNr, bool log = false, bool neighbour = true);
        SimProps simProps;

        void printData(int nr = 64);
        void printGridData();
        void setSimProps();
        void setSimProps(SimProps& props);
        void setData(const std::vector<Particle>& p);
    //private:
        int cellx, celly, cellz;

        bool first = true;
        std::vector<int> grid_data;
        std::vector<Particle> particles1;
        std::vector<Particle> particles2;

        void place();
        void originalCalculateRhoAndP();
        void neighbourCalculateRhoAndP();
        void originalIter();
        void neighbourIter();
        void container(int idx);

        float kernel(int i, int j);
        glm::vec3 gradKernel(int i, int j);

        glm::vec3 determineGridTuple(int particleIdx);
        int getStartIdxOfCell(glm::vec3 tuple);
        Neighbourhood getNeighbouringCells(glm::vec3 cellTuple);
        bool checkBounds(glm::vec3 tuple);
        void placeParticleIntoCell(int particleIdx);
        void printNeighbourhood(Neighbourhood nh);


        // IISPH according to https://cg.informatik.uni-freiburg.de/publications/2013_TVCG_IISPH.pdf
        //--------------------------------------------------------------------------------------------------------------

        // variables
        std::vector<glm::vec3> dii; // displacement
        std::vector<float> aii;
        std::vector<glm::vec3> sumDijPj;
        std::vector<float> rho_adv;
        std::vector<float> rho_pred;
        std::vector<glm::vec3> v_adv;

        // functions
        void calculateRho();
        float calculateAverageError();
        void computeVadvAndDii(float dt);
        void computeRhoadvAndAii(float dt);
        void computeAii(int particleIdx, float dt);
        void predictAdvection(float dt);
        void computeSumDijPj(float dt);
        void updatePressure(float dt);
        void pressureSolve(float dt);
        void integrate(float dt);
        void runIISPH(int iterNr);
    };
}


#endif //SPHSIMULATION_CPUSIM_HPP
