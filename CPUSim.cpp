#include "CPUSim.hpp"

#define PI 3.1415926538
#define list_size 20
#define ppic false

namespace Vltava {

// Neighbourhood stuff
//----------------------------------------------------------------------------------------------------------------------
    glm::vec3 CPUSim::determineGridTuple(int particleIdx) {
        auto& particles = first ? particles1 : particles2;
        auto& p = particles.at(particleIdx);

        // Checking bounds
        glm::vec3 lower(
                simProps.gridA.x < simProps.gridB.x ? simProps.gridA.x : simProps.gridB.x,
                simProps.gridA.y < simProps.gridB.y ? simProps.gridA.y : simProps.gridB.y,
                simProps.gridA.z < simProps.gridB.z ? simProps.gridA.z : simProps.gridB.z
        );

        glm::vec3 upper(
                simProps.gridA.x > simProps.gridB.x ? simProps.gridA.x : simProps.gridB.x,
                simProps.gridA.y > simProps.gridB.y ? simProps.gridA.y : simProps.gridB.y,
                simProps.gridA.z > simProps.gridB.z ? simProps.gridA.z : simProps.gridB.z
        );

        if ((p.x.x < lower.x) || (p.x.y < lower.y) || (p.x.z < lower.z) ||
            (p.x.x > upper.x) || (p.x.y > upper.y) || (p.x.z > upper.z)) {
            return {-1,-1,-1};
        }

        glm::vec3 diff = glm::vec3(simProps.gridA.x, simProps.gridA.y, simProps.gridA.z) - particles.at(particleIdx).x;
        diff /= simProps.kernelh;

        return floor(abs(diff));
    }

    int CPUSim::getStartIdxOfCell(glm::vec3 tuple) {
        if (!checkBounds(tuple)) return -1;

        return int((tuple.x * cellz * celly + tuple.y * cellz + tuple.z)*list_size);
    }

    Neighbourhood CPUSim::getNeighbouringCells(glm::vec3 cellTuple) {
        Neighbourhood ret;

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    ret.neighbour[i * 9 + j * 3 + k] = glm::vec3(cellTuple.x + (i-1), cellTuple.y + (j-1), cellTuple.z + (k-1));
                }
            }
        }

        return ret;
    }

    bool CPUSim::checkBounds(glm::vec3 tuple) {
        if (tuple.x < 0 || tuple.x >= cellx)
            return false;

        if (tuple.y < 0 || tuple.y >= celly)
            return false;

        if (tuple.z < 0 || tuple.z >= cellz)
            return false;

        return true;
    }

    void CPUSim::placeParticleIntoCell(int particleIdx) {
        glm::vec3 cellTuple = determineGridTuple(particleIdx);
#if ppic
        std::cout << "[ppic] particleIdx" << particleIdx << " placed into (" << cellTuple.x << "," << cellTuple.y << "," << cellTuple.z << ")\n";
#endif
        //if (checkBounds(cellTuple)) {
            int startIdx = getStartIdxOfCell(cellTuple);
            if (startIdx >= 0) {
                int realIdx = /*atomicAdd(*/grid_data.at(startIdx);
                grid_data.at(startIdx) += 1;
                realIdx += startIdx + 1;

                // Placing particleIdx into cell
                grid_data.at(realIdx) = particleIdx;
            }
        //}
    }

    void CPUSim::printNeighbourhood(Neighbourhood nh) {
        for (int i = 0; i < 27; i++) {
            glm::vec3 current = nh.neighbour[i];
            if (!checkBounds(current)) continue;

            int idx = getStartIdxOfCell(current);

            std::cout << "Tuple: " << current.x << " " << current.y << " " << current.z << "\n---------------------" << std::endl;
            int size = grid_data.at(getStartIdxOfCell(current));
            for (int j = 1; j < size+1; j++) {
                std::cout << grid_data.at(idx + j) << " ";
            }

            std::cout << std::endl;
        }
    }

//----------------------------------------------------------------------------------------------------------------------

    void CPUSim::setSimProps(SimProps& props) {
        this->simProps = props;

        cellx = (int) ceil(abs((simProps.gridB.x - simProps.gridA.x)/simProps.kernelh)); // Number of cells in x direction
        celly = (int) ceil(abs((simProps.gridB.y - simProps.gridA.y)/simProps.kernelh)); // Number of cells in y direction
        cellz = (int) ceil(abs((simProps.gridB.z - simProps.gridA.z)/simProps.kernelh)); // Number of cells in z direction
    }

    void CPUSim::setData(const std::vector<Particle>& p) {
        particles1.clear();
        particles2.clear();

        particles1.reserve(p.size());
        particles2.reserve(p.size());

        this->particles1 = p;
        this->particles2 = p;
    }

    CPUSim::CPUSim() {
        // Water rest density = 1 kg/m^3
        //
        // Particle h = 0.05 meter
        // Particle volume = 5.2359 * 10^-4 m^3 = 0.005236 m^3 (4/3 * 0.05^3 * PI)
        // Particle mass = volume * density = 0.005236 kg
        //
        // Smoothing length := 2 * h = 0.1 meter (for now)

        int nrOfP = 64;
        particles1.reserve(nrOfP);
        particles2.reserve(nrOfP);
        float s=0.10; // meter

        int pnrAlongAxis = round(pow(nrOfP, 1.0/3.0));

        nrOfP = pnrAlongAxis * pnrAlongAxis * pnrAlongAxis;

        int r = -1;
        int z = -1;

        float mass = 0.005236f; // kg

        for (int i = 0; i < nrOfP; i++) {
            Particle data;
            if (i % pnrAlongAxis == 0)
                r++;

            if (i % (pnrAlongAxis * pnrAlongAxis) == 0)
                z++;

            data.x = glm::vec3((i % pnrAlongAxis) * s, (r % pnrAlongAxis) * s, z * s);

            data.h = 0.1;
            data.v = glm::vec3(0, 0, 0);
            data.m = mass;
            //data[i].m = 1.0f;

            data.rho = 1;
            data.p = 1;

            particles1.push_back(data);
            particles2.push_back(data);
        }
    }

    void CPUSim::run(int iterNr) {
        for (int i = 0; i < iterNr; i++) {
            place();
            //printGridData();
            calculateRhoAndP();
            iter();
        }
        printData();
    }

    void CPUSim::setSimProps() {
        simProps.desired_density = 1.0f;
        simProps.k = 0.001f;
        simProps.nr_of_particles = 64;
        simProps.kernelh = 0.2f;

        simProps.gridA = glm::vec4(-1,-1,-1,0);
        simProps.gridB = glm::vec4(1,1,1,0);

        cellx = (int) ceil(abs((simProps.gridB.x - simProps.gridA.x)/simProps.kernelh)); // Number of cells in x direction
        celly = (int) ceil(abs((simProps.gridB.y - simProps.gridA.y)/simProps.kernelh)); // Number of cells in y direction
        cellz = (int) ceil(abs((simProps.gridB.z - simProps.gridA.z)/simProps.kernelh)); // Number of cells in z direction

        std::cout << cellx << " " << celly << " " << cellz << std::endl;
    }

    void CPUSim::printGridData() {
        std::cout << "\ngrid_data:\n";
        std::cout << "----------------------\n";
        for (auto& i : grid_data) {
            std::cout << i << std::endl;
        }
    }

    void CPUSim::printData(int nr) {
        auto& particles = first ? particles1 : particles2;

        for (int i = 0; i < nr; i++) {
            auto& p = particles[i];
            std::cout << "Density: " << p.rho <<
                      "; Pressure: " << p.p <<
                      " ; Position: " << p.x.x << " " << p.x.y << " " << p.x.z <<
                      " ; Mass: " << p.m <<
                      " ; Velocity: " << p.v.x << " " << p.v.y << " " << p.v.z <<
                      " ; padding: " << p.padding;
                      if (p.rho - p.padding > 0)
                        std::cout << " ; diff = " << p.rho - p.padding;

                      std::cout << std::endl;
        }
        std::cout << "----------\n";
    }

    void CPUSim::place() {
        auto& particles = first ? particles1 : particles2;

        grid_data.clear();
        grid_data.resize(1, 0);
        grid_data.resize(list_size*cellx*celly*cellz, 0);

        for (int i = 0; i < particles.size(); i++) {
            placeParticleIntoCell(i);
        }
    }

    void CPUSim::calculateRhoAndP() {
        auto& particles = first ? particles1 : particles2;

        // Density calculation
        for (int i = 0; i < particles.size(); i++) {
            float density = 0.0f;
            float ogd = 0.0f;

            // Original
            for (int j = 0; j < particles.size(); j++) {
                if (i == j)
                    continue;

                /*density*/ density += particles[j].m * kernel(i, j);
            }

            // Neighbour
            glm::vec3 tuple = determineGridTuple(i);
            Neighbourhood n = getNeighbouringCells(tuple);
            //for (auto current : n.neighbour) {
            for (int nr = 0; nr < 27; nr++) {
                glm::vec3 current = n.neighbour[nr];
                if (!checkBounds(current)) continue;

                int idx = getStartIdxOfCell(current);
                if (idx >= 0) {
                    int size = grid_data.at(idx);

                    int iterIdx = 0;
                    for (int j = 1; j < size + 1; j++) {
                        iterIdx = grid_data.at(idx + j);
                        if (i == iterIdx) continue;

                        //density += in_data.p[gID].m * kernel(in_data.p[gID].x, in_data.p[iterIdx].x);
                        ogd += particles[i].m * kernel(i, iterIdx);
                    }
                }
            }

            particles[i].rho = density;
            particles[i].padding = ogd;

            // Pressure calculation
            float ratio = density / simProps.desired_density;
            float p = simProps.k * (pow(ratio, 7) - 1.0);
            particles[i].p = p < 0 ? 0 : p;
        }
    }

    void CPUSim::iter() {
        const auto& p1 = first ? particles1 : particles2;
        auto& p2 = first ? particles2 : particles1;
        for (int i = 0; i < p1.size(); i++) {

            if (p1.at(i).staticP == 0) {

                // Calculate pressure
                glm::vec3 pressure = glm::vec3(0);
                glm::vec3 viscosity = glm::vec3(0);

                // Original
                for (int j = 0; j < p1.size(); j++) {
                    if (j == i)
                        continue;

                    // Note to self: as the particles get further and further from each other the density decreases which means rho --> 0
                    // whihc leads to something/0^2, which is either inf or -inf ----> nan or -nan
                    float val = 0;
                    if (p1[j].rho != 0 && p1[i].rho != 0)
                        val = p1[j].m * ((p1[i].p / pow(p1[i].rho, 2)) + (p1[j].p / pow(p1[j].rho, 2)));

                    pressure += val * gradKernel(i, j);

                    glm::vec3 xij = p1[i].x - p1[j].x;

                    float pval = 0;
                    if (p1[j].rho != 0)
                        pval = (p1[j].m / p1[j].rho) *
                               (dot(xij, gradKernel(i, j)) / (dot(xij, xij) + 0.01 * simProps.kernelh));

                    glm::vec3 vij = p1[i].v - p1[j].v;
                    viscosity += pval * vij;
                }

                // Neighbour
                /*glm::vec3 tuple = determineGridTuple(i);
                Neighbourhood n = getNeighbouringCells(tuple);
                for (auto current : n.neighbour) {
                    if (!checkBounds(current)) continue;

                    int idx = getStartIdxOfCell(current);
                    int size = grid_data.at(idx);

                    int iterIdx = 0;
                    for (int j = 1; j < size+1; j++) {
                        iterIdx = grid_data.at(idx + j);
                        if (i == iterIdx) continue;

                        // Note to self: as the particles get further and further from each other the density decreases which means rho --> 0
                        // whihc leads to something/0^2, which is either inf or -inf ----> nan or -nan
                        float val = 0;
                        if (p1[iterIdx].rho != 0 && p1[i].rho != 0)
                            val = p1[iterIdx].m * ((p1[i].p / pow(p1[i].rho, 2)) + (p1[iterIdx].p / pow(p1[iterIdx].rho, 2)));

                        pressure += val * gradKernel(i, iterIdx);

                        glm::vec3 xij = p1[i].x - p1[iterIdx].x;

                        float pval = 0;
                        if (p1[iterIdx].rho != 0)
                            pval = (p1[iterIdx].m / p1[iterIdx].rho) * (dot(xij, gradKernel(i, iterIdx)) /  (dot(xij, xij) + 0.01 * simProps.kernelh));

                        glm::vec3 vij = p1[i].v - p1[iterIdx].v;
                        viscosity += pval * vij;
                    }
                }*/

                pressure *= -p1[i].m;

                float nu = 0.01;
                viscosity *= 2 * nu * p1[i].m;

                glm::vec3 gravity(0, 0, -9.81 * p1[i].m);
                glm::vec3 acc = (pressure + viscosity /*+ gravity*/) / p1[i].m;

                float dt = 0.01;
                glm::vec3 viNext = p1[i].v;
                glm::vec3 xiNext = p1[i].x;
                viNext += acc * (float) dt / 1.0f;
                xiNext += viNext * dt;

                p2[i].x = xiNext;
                p2[i].v = viNext;
                p2[i].rho = p1[i].rho;
                p2[i].p = p1[i].p;
            } else {
                p2.at(i) = p1.at(i);
            }
        }

        first = !first;
    }

    float CPUSim::kernel(int i, int j) {
        auto& particles = first ? particles1 : particles2;

        float q = glm::length(particles[i].x - particles[j].x)/simProps.kernelh;
        float m_k = 8.0 / (PI * pow(simProps.kernelh, 3));
        float m_l = 1.0/(pow(simProps.kernelh, 3)) * 3.0/(2.0*PI);

        float ret = 0;
        if (q <= 1.0) {
            if (q <= 0.5) {
                ret = m_k * (6.0 * (pow(q, 3.0) - pow(q, 2.0)) + 1.0);
            } else {
                ret = m_k * (2.0 * pow(1 - q, 3));
            }
        }

        return ret;
    }

    // One of the errors: vec3 i == vec3 j ---> normalize(i-j) -> nan
    glm::vec3 CPUSim::gradKernel(int i, int j) {
        float m_k = 48.0 / (pow(simProps.kernelh, 3) * PI);
        float m_l = 1.0/(pow(simProps.kernelh, 3)) * 3.0/(2.0*PI);

        auto& particles = first ? particles1 : particles2;

        glm::vec3 r = particles[i].x - particles[j].x;
        float rlength = length(r);
        float q = rlength / simProps.kernelh;
        glm::vec3 ret = glm::vec3(0);

        if (q > 0.0001 && q <= 1.0) {
            glm::vec3 gradq = r / (rlength * simProps.kernelh);

            if (q <= 0.5) {
                gradq *= m_k * q * (3.0 * q - 2);
                ret = gradq;
            } else {
                float factor = 1.0 - q;
                ret = m_k * (-factor * factor) * gradq;
            }
        }

        return ret;
    }
}