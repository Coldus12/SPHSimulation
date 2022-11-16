#ifndef SPHSIMULATION_SPH_H
#define SPHSIMULATION_SPH_H

#include <vector>
#include "ParticleModel.hpp"

namespace Vltava {
    struct Neighbourhood {
        glm::vec3 neighbour[27]; // 3^3
    };

    // Abstract SPH base class
    class SPH {
    public:
        virtual void gpuTimeStep() = 0;
        virtual void cpuTimeStep() = 0;

        virtual std::string log();

        virtual void setSimProps(SimProps props) {
            this->props = props;
        }

        virtual void setData(std::vector<Particle> particles) {
            this->particles1 = particles;
            this->particles2 = particles;
        }

        // TODO: pointers could theoretically be null, and the vector may be empty. Handle these exceptions!
        virtual void setBuffers(std::vector<Buffer>* uBuffers, std::vector<Buffer>* sBuffers) {
            std::cout << "got here at least?" << std::endl;

            this->uBuffers = uBuffers;
            this->sBuffers = sBuffers;

            // The first 2 storage buffers in the sBuffers
            // should contain the in/out particle buffers.
            this->particles1 = sBuffers->at(0).getData<Particle>();
            this->particles2 = sBuffers->at(1).getData<Particle>();

            // Getting simProps from the uBuffers
            this->props = uBuffers->at(0).getData<SimProps>()[0];
        }

        virtual void setCellSizes(int cellx, int celly, int cellz, int listSize) {
            this->cellx = cellx;
            this->celly = celly;
            this->cellz = cellz;

            this->list_size = listSize;
        }

        // TODO rename first to something more understandable
        bool getFirst() {
            return first;
        }

        std::vector<Particle> getData1() {
            return particles1;
        }

        std::vector<Particle> getData2() {
            return particles2;
        }
    protected:
        std::vector<Buffer>* uBuffers = nullptr;
        std::vector<Buffer>* sBuffers = nullptr;

        SimProps props;
        std::vector<Particle> particles1;
        std::vector<Particle> particles2;

        int cellx, celly, cellz, list_size = 25;
        bool first = true;
        bool neighbour = true;
        std::vector<int> grid_data;

        // false means log will show GPU data instead
        bool logCpuData = true;

        // Functions
        virtual float kernel(int i, int j);
        virtual float kernel(glm::vec3 i, glm::vec3 j);
        virtual glm::vec3 gradKernel(int i, int j);
        virtual glm::vec3 gradKernel(glm::vec3 i, glm::vec3 j);

        void container(int idx);
        void place();
        glm::vec3 determineGridTuple(int particleIdx);
        int getStartIdxOfCell(glm::vec3 tuple);
        Neighbourhood getNeighbouringCells(glm::vec3 cellTuple);
        bool checkBounds(glm::vec3 tuple);
        void placeParticleIntoCell(int particleIdx);
    };

} // Vltava

#endif //SPHSIMULATION_SPH_H
