#include "CPUSim.hpp"

#define PI 3.1415926538
#define list_size 25
#define ppic false

namespace Vltava {

    void CPUSim::container(int idx) {
        float left = (simProps.gridA.x < simProps.gridB.x ? simProps.gridA.x : simProps.gridB.x) + 0.1;
        float right = (simProps.gridA.x > simProps.gridB.x ? simProps.gridA.x : simProps.gridB.x) - 0.1;

        float front = (simProps.gridA.y > simProps.gridB.y ? simProps.gridA.y : simProps.gridB.y) - 0.1;
        float back = (simProps.gridA.y < simProps.gridB.y ? simProps.gridA.y : simProps.gridB.y) + 0.1;

        float bottom = (simProps.gridA.z < simProps.gridB.z ? simProps.gridA.z : simProps.gridB.z) + 0.1;
        //float bottom = -0.1;
        auto& particles = first ? particles2 : particles1;

        // Bottom
        if (particles[idx].x.z <= bottom) {
            particles[idx].v.z = -particles[idx].v.z;
            particles[idx].x.z = bottom;
        }

        // Left
        if (particles[idx].x.x <= left) {
            particles[idx].v.x = -particles[idx].v.x;
            particles[idx].x.x = left;
        }

        // Right
        if (particles[idx].x.x >= right) {
            particles[idx].v.x = -particles[idx].v.x;
            particles[idx].x.x = right;
        }

        // Front
        if (particles[idx].x.y >= front) {
            particles[idx].v.y = -particles[idx].v.y;
            particles[idx].x.y = front;
        }

        // Back
        if (particles[idx].x.y <= back) {
            particles[idx].v.y = -particles[idx].v.y;
            particles[idx].x.y = back;
        }
    }

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

        diagonals.reserve(p.size());
        sourceTerms.reserve(p.size());
        pressureAccs.reserve(p.size());
        densityError.reserve(p.size());

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

    void CPUSim::runSESPH(int iterNr, bool log, bool neighbour) {
        for (int i = 0; i < iterNr; i++) {
            place();
            //printGridData();
            if (neighbour) {
                neighbourCalculateRhoAndP();
                neighbourIter();
            } else {
                originalCalculateRhoAndP();
                originalIter();
            }
        }
        if (log)
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

    void CPUSim::neighbourCalculateRhoAndP() {
        auto& particles = first ? particles1 : particles2;

        // Density calculation
        for (int i = 0; i < particles.size(); i++) {
            float density = 0.0f;
            float ogd = 0.0f;

            // Original
            /*for (int j = 0; j < particles.size(); j++) {
                if (i == j)
                    continue;

                density += particles[j].m * kernel(i, j);
            }*/

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
                        density += particles[i].m * kernel(i, iterIdx);
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

    void CPUSim::originalCalculateRhoAndP() {
        auto& particles = first ? particles1 : particles2;

        // Density calculation
        for (int i = 0; i < particles.size(); i++) {
            float density = 0.0f;
            float ogd = 0.0f;

            // Original
            for (int j = 0; j < particles.size(); j++) {
                if (i == j)
                    continue;

                density += particles[j].m * kernel(i, j);
            }

            particles[i].rho = density;
            particles[i].padding = ogd;

            // Pressure calculation
            float ratio = density / simProps.desired_density;
            float p = simProps.k * (pow(ratio, 7) - 1.0);
            particles[i].p = p < 0 ? 0 : p;
        }
    }

    void CPUSim::neighbourIter() {
        const auto& p1 = first ? particles1 : particles2;
        auto& p2 = first ? particles2 : particles1;
        for (int i = 0; i < p1.size(); i++) {

            if (p1.at(i).staticP == 0) {

                // Calculate pressure
                glm::vec3 pressure = glm::vec3(0);
                glm::vec3 viscosity = glm::vec3(0);

                // Neighbour
                glm::vec3 tuple = determineGridTuple(i);
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
                }

                pressure *= -p1[i].m;

                float nu = 0.01;
                viscosity *= 2 * nu * p1[i].m;

                glm::vec3 gravity(0, 0, -9.81 * p1[i].m);
                glm::vec3 acc = (pressure + viscosity + gravity) / p1[i].m;

                float dt = 0.01;
                glm::vec3 viNext = p1[i].v;
                glm::vec3 xiNext = p1[i].x;
                viNext += acc * (float) dt / 1.0f;
                xiNext += viNext * dt;

                p2[i].x = xiNext;
                p2[i].v = viNext;
                p2[i].rho = p1[i].rho;
                p2[i].p = p1[i].p;

                container(i);
            } else {
                p2.at(i) = p1.at(i);
            }
        }

        first = !first;
    }

    void CPUSim::originalIter() {
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

                pressure *= -p1[i].m;

                float nu = 0.01;
                viscosity *= 2 * nu * p1[i].m;

                glm::vec3 gravity(0, 0, -9.81 * p1[i].m);
                glm::vec3 acc = (pressure + viscosity + gravity) / p1[i].m;

                float dt = 0.01;
                glm::vec3 viNext = p1[i].v;
                glm::vec3 xiNext = p1[i].x;
                viNext += acc * (float) dt / 1.0f;
                xiNext += viNext * dt;

                p2[i].x = xiNext;
                p2[i].v = viNext;
                p2[i].rho = p1[i].rho;
                p2[i].p = p1[i].p;

                container(i);
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

    // IISPH stuff
    //------------------------------------------------------------------------------------------------------------------

    // According to https://interactivecomputergraphics.github.io/SPH-Tutorial/pdf/SPH_Tutorial.pdf
    //------------------------------------------------------------------------------------------------------------------
    float CPUSim::calculateDiagonal(int particleIdx, float dt) {
        auto& particles = first ? particles1 : particles2;
        float dtmsquared = -1 * dt * dt + particles[0].m * particles[0].m;

        float val = 0;
        glm::vec3 data(0);

        /*for (int j = 0; j < particles.size(); j++) {
            for (int j2 = 0; j2 < particles.size(); j2++) {
                data = gradKernel(particleIdx, j2);
                data /= (particles[j2].rho * particles[j2].rho);

                data += gradKernel(particleIdx, j2) / (particles[particleIdx].rho * particles[particleIdx].rho);
            }
            val += glm::dot(data, gradKernel(particleIdx, j));
        }*/

        for (int j = 0; j < particles.size(); j++) {
            for (int j2 = 0; j2 < particles.size(); j2++) {
                data += gradKernel(particleIdx, j2);
                data /= (particles[j2].rho * particles[j2].rho);
            }
            val += glm::dot(data, gradKernel(particleIdx, j));

            data = gradKernel(particleIdx, j);
            data /= (particles[particleIdx].rho * particles[particleIdx].rho);
            val += glm::dot(data, gradKernel(particleIdx, j));
            data = glm::vec3(0);
        }

        /*for (int j = 0; j < particles.size(); j++) {
            data = gradKernel(particleIdx, j);
            data /= (particles[particleIdx].rho * particles[particleIdx].rho);
            val += glm::dot(data, gradKernel(particleIdx, j));
        }*/

        return val * dtmsquared;
    }

    float CPUSim::calculateSourceTerm(int particleIdx, float dt) {
        auto& particles = first ? particles1 : particles2;
        float val = 0;

        // a^nonp is same for both vi* and vj*, therefore if vi* = vi + dt*anonp
        // and we need vi*-vj* => vi - vj
        /*glm::vec3 anonpdt = glm::vec3(0, 0, -9.81) * dt;

        glm::vec3 vi = particles[particleIdx].v + anonpdt;
        glm::vec3 vj;*/
        glm::vec3 viscosity(0);

        for (int j = 0; j < particles.size(); j++) {
            if (j == particleIdx)
                continue;

            glm::vec3 xij = particles[particleIdx].x - particles[j].x;

            float pval = 0;
            if (particles[j].rho != 0)
                pval = (particles[j].m / particles[j].rho) *
                       (dot(xij, gradKernel(particleIdx, j)) / (dot(xij, xij) + 0.01 * simProps.kernelh));

            glm::vec3 vij = particles[particleIdx].v - particles[j].v;
            viscosity += pval * vij;
        }

        glm::vec3 vi = particles[particleIdx].v + dt * (glm::vec3(0,0,-9.81) + viscosity);
        glm::vec3 vj;

        for (int j = 0; j < particles.size(); j++) {
            vj = particles[particleIdx].v + dt * glm::vec3(0,0,-9.81);
            val += glm::dot(vi - vj, gradKernel(particleIdx, j));
        }

        float retVal = simProps.desired_density - particles[particleIdx].rho - dt * particles[0].m * val;

        return retVal;
    }

    void CPUSim::initPressure() {
        auto& particles = first ? particles1 : particles2;
        for (auto& particle : particles) {
            particle.p = 0;
        }
    }

    glm::vec3 CPUSim::computePressureAcceleration(int particleIdx) {
        auto& particles = first ? particles1 : particles2;
        glm::vec3 api = glm::vec3(0);

        for (int j = 0; j < particles.size(); j++) {
            float val = particles[particleIdx].p / (particles[particleIdx].rho * particles[particleIdx].rho)
                    + particles[j].p / (particles[j].rho * particles[j].rho);

            api += val * gradKernel(particleIdx, j);
        }

        api *= -1 * particles[particleIdx].m;

        return api;
    }

    float CPUSim::computeLaplacian(int particleIdx, float dt) {
        auto& particles = first ? particles1 : particles2;
        float val = 0;
        glm::vec3 ai, aj;

        for (int j = 0; j < particles.size(); j++) {
            //ai = computePressureAcceleration(particleIdx);
            //aj = computePressureAcceleration(j);
            ai = pressureAccs[particleIdx];
            aj = pressureAccs[j];

            val += glm::dot(ai - aj, gradKernel(particleIdx, j));
        }
        val *= particles[0].m * dt * dt;

        return val;
    }

    void CPUSim::updatePressure(int particleIdx, float dt) {
        auto& particles = first ? particles1 : particles2;
        float omega = 0.5f;
        //float si = calculateSourceTerm(particleIdx, dt);
        //float li = computeLaplacian(particleIdx);
        float si = sourceTerms[particleIdx];
        float li = computeLaplacian(particleIdx, dt);

        float val = 0;
        //if (abs(diagonals[particleIdx]) > 0.01)
            val = particles[particleIdx].p + omega/diagonals[particleIdx] * (si - li);
        //else val = particles[particleIdx].p;

        //std::cout << "p " << particles[particleIdx].p << " val: " << val << std::endl;
        assert((simProps.desired_density != 0) && "desired density is 0 -> nan/inf error");
        densityError.push_back((li - si) / simProps.desired_density);

        particles[particleIdx].p = val < 0 ? 0 : val;
    }

    float CPUSim::calculateAverageError() {
        float val = 0;
        for (float i : densityError) {
            val += i > 0 ? i : -i;
        }

        val /= densityError.size();

        return val;
    }

    void CPUSim::calculateRho() {
        auto& particles = first ? particles1 : particles2;
        for (int i = 0; i < particles.size(); i++) {
            float density = 0.0f;

            for (int j = 0; j < particles.size(); j++) {
                if (i == j)
                    continue;

                density += particles[j].m * kernel(i, j);
            }

            particles[i].rho = density;
            //particles[i].p = 0;
        }
    }

    void CPUSim::runIISPH(int iterNr) {
        auto& particles = first ? particles1 : particles2;
        float dt = 0.01;

        for (int iter = 0; iter < iterNr; iter++) {
            std::cout << iter << std::endl;
            // Step 1: Calculate density for all particles
            // Also sets pressure to 0.
            calculateRho();

            // Predict advection
            // Clearing diagonals
            /*diagonals.clear();
            for (int i = 0; i < particles.size(); i++) {
                //particles[i].
            }*/

            // Step 2: Calculate diagonal elements and source terms
            for (int i = 0; i < particles.size(); i++) {
                diagonals.push_back(calculateDiagonal(i, dt));
                sourceTerms.push_back(calculateSourceTerm(i, dt));
            }

            // Step 3: Pressure computation
            float error = 10;
            int nr = 0;
            while (error > 0.1f && nr <= 15000) {
                nr++;
                // Computing pressure accelerations
                for (int i = 0; i < particles.size(); i++) {
                    pressureAccs.push_back(computePressureAcceleration(i));
                }

                // Pressure update
                for (int i = 0; i < particles.size(); i++) {
                    updatePressure(i, dt);
                }

                // Calculating average error
                error = calculateAverageError();

                // Clearing stuff
                pressureAccs.clear();
                densityError.clear();

                //if (error < 0.12)
                    std::cout << "[IISPH] pressure average error: " << error << " nr " << nr << std::endl;
            }

            // Step 4: Iteration
            originalIter();

            // Clearing stuff
            densityError.clear();
            diagonals.clear();
            sourceTerms.clear();
            pressureAccs.clear();

            //printData();
        }
    }

    // According to https://cg.informatik.uni-freiburg.de/publications/2013_TVCG_IISPH.pdf
    //------------------------------------------------------------------------------------------------------------------

    void CPUSim::computeVadvAndDiiAndRho(int particleIdx, float dt) {
        auto& particles = first ? particles1 : particles2;
        glm::vec3 viscosity(0);
        glm::vec3 dii(0);

        //float rho = 0;

        for (int j = 0; j < particles.size(); j++) {
            // Computing dii
            dii += -dt * dt * particles[j].m * gradKernel(particleIdx, j) / (particles[particleIdx].rho * particles[particleIdx].rho);

            if (j == particleIdx)
                continue;

            // Computing rho
            //rho += particles[j].m * kernel(particleIdx, j);

            // Computing viscosity for non-pressure velocity advection
            glm::vec3 xij = particles[particleIdx].x - particles[j].x;

            float pval = 0;
            if (particles[j].rho != 0)
                pval = (particles[j].m / particles[j].rho) *
                       (dot(xij, gradKernel(particleIdx, j)) / (dot(xij, xij) + 0.01 * simProps.kernelh));

            glm::vec3 vij = particles[particleIdx].v - particles[j].v;
            viscosity += pval * vij;
        }

        glm::vec3 vi = particles[particleIdx].v + dt * (glm::vec3(0,0,-9.81) + viscosity);
        /*particles1[particleIdx].v = vi;
        particles2[particleIdx].v = vi;*/

        //rho = (rho - simProps.desired_density) < 0 ? 0 : rho;

        v_advected.push_back(vi);
        this->dii.push_back(dii);
        //particles[particleIdx].rho = rho;
    }

    void CPUSim::computeAdvectedRho(int particleIdx, float dt) {
        auto& particles = first ? particles1 : particles2;
        float val = 0;

        for (int j = 0; j < particles.size(); j++) {
            if (j == particleIdx)
                continue;

            val += particles[j].m * glm::dot(v_advected[particleIdx] - v_advected[j], gradKernel(particleIdx, j));
        }

        val *= dt;

        /*particles1[particleIdx].rho += val;
        particles2[particleIdx].rho += val;*/
        rho_advected.push_back(particles[particleIdx].rho + val);
    }

    void CPUSim::computeAii(int particleIdx, float dt) {
        auto& particles = first ? particles1 : particles2;
        float aii = 0;

        for (int j = 0; j < particles.size(); j++) {
            if (j == particleIdx)
                continue;

            glm::vec3 grad = gradKernel(particleIdx, j);

            glm::vec3 dji = -dt * dt * particles[j].m / (particles[j].rho * particles[j].rho) * grad;
            //glm::vec3 dji = -dt * dt * particles[j].m / (particles[particleIdx].rho * particles[particleIdx].rho) * -grad;

            aii += particles[j].m * glm::dot(dii[particleIdx] - dji, grad);
        }

        if (aii != aii) {
            std::cout << "NAN here in diag calculation; id: " << particleIdx << std::endl;
            std::cout << particles[particleIdx].rho << std::endl;
            std::cout << std::endl;
        }

        diagonals.push_back(aii);
    }

    void CPUSim::predictAdvection(float dt) {
        // Computing rho
        calculateRho();

        dii.clear();
        v_advected.clear();

        // Predict v adv, and compute dii
        for (int i = 0; i < particles1.size(); i++) {
            computeVadvAndDiiAndRho(i, dt);
        }

        diagonals.clear();
        rho_advected.clear();

        // Calculate advected rho
        for (int i = 0; i < particles1.size(); i++) {
            //particles1[i].v = v_advected[i];
            //particles2[i].v = v_advected[i];

            computeAdvectedRho(i, dt);
            particles1[i].p *= 0.5;
            particles2[i].p *= 0.5;
            computeAii(i, dt);
        }
    }

    void CPUSim::computeSumDijPj(int particleIdx, float dt) {
        auto& particles = first ? particles1 : particles2;
        glm::vec3 dijpj(0);

        for (int j = 0; j < particles.size(); j++) {
            float val = -particles[j].m * particles[j].p / (particles[j].rho * particles[j].rho);
            dijpj += gradKernel(particleIdx, j) * val;
        }

        dijpj *= dt*dt;
        sumDijPj.push_back(dijpj);
    }

    void CPUSim::updatePressure2(int particleIdx, float dt) {
        auto& particles = first ? particles1 : particles2;
        glm::vec3 djipi, grad, vecSum;
        float sum = 0;
        //float djivals = -dt*dt * (particles[particleIdx].m * particles[particleIdx].p) / (particles[particleIdx].rho * particles[particleIdx].rho);

        float omega = 0.5;

        for (int j = 0; j < particles.size(); j++) {
            grad = gradKernel(particleIdx, j);
            float djivals = -dt*dt * (particles[j].m * particles[particleIdx].p) / (particles[j].rho * particles[j].rho);
            //float djivals = -dt*dt * (particles[j].m * particles[particleIdx].p) / (particles[particleIdx].rho * particles[particleIdx].rho);
            djipi = djivals * gradKernel(j, particleIdx);
            sum += particles[j].m * glm::dot((sumDijPj[particleIdx] - dii[j] * particles[j].p - (sumDijPj[j] - djipi)), grad);

            /*if (particleIdx == 0) {
                std::cout << "djivals" << djivals << " sumDijPj (" << sumDijPj[particleIdx].x << "," << sumDijPj[particleIdx].y << "," << sumDijPj[particleIdx].z
                << ") dii[j] (" << dii[j].x << "," << dii[j].y << "," << dii[j].z<< ") original p " << particles[j].p
                << " djkpk (" << (sumDijPj[j] - djipi).x << "," << (sumDijPj[j] - djipi).y << "," << (sumDijPj[j] - djipi).z <<")"<< std::endl;
            }*/
        }

        //if(particleIdx == 0) std::cout <<"\n------------------\n";

        //particles[particleIdx].p = (1 - omega) * particles[particleIdx].p + (omega/diagonals[particleIdx]) * (simProps.desired_density - rho_advected[particleIdx] - sum);

        // NOTES
        // p^(l+1) = (1-omega) * p + (omega/diag) * (rho0 - rho_adv - sum)
        // p^(l+1) * diag = (1-omega) * p * diag + omega * (rho0 - rho_adv - sum)
        // p^(l+1) * diag = p * diag - omega * (rho_adv - rho0 + sum + p * diag)
        // rho_pred = rho_adv + sum + p * diag
        // p^(l+1) * diag = p * diag - omega * (rho_pred - rho0)
        // p^(l+1) = p - (omega/diag) * (rho_pred - rho0)

        float previous_p = particles[particleIdx].p;
        float l_rho_pred = rho_advected[particleIdx] + sum + particles[particleIdx].p * diagonals[particleIdx];
        particles[particleIdx].p = particles[particleIdx].p - (omega / diagonals[particleIdx]) * (l_rho_pred - simProps.desired_density);

        //particles[particleIdx].p = particles[particleIdx].p < 0 ? 0 : particles[particleIdx].p;
        std::cout << "[UpdatePressure2] particleIdx: " << particleIdx << " rho_adv: " << rho_advected[particleIdx] << " p: " << particles[particleIdx].p << " diag: " << diagonals[particleIdx] << " sum: " << sum  << " rho_pred " << l_rho_pred << " prevP " << previous_p << std::endl;
        //rho_pred.push_back(rho_advected[particleIdx] + particles[particleIdx].p * diagonals[particleIdx] + sum);
        rho_pred.push_back(l_rho_pred);
    }

    float CPUSim::computeAverageRho() {
        float val = 0;

        for (auto& r : rho_pred) {
            val += r;
        }

        val /= rho_pred.size();

        return val;
    }

    void CPUSim::pressureSolve(float dt) {
        auto& particles = first ? particles1 : particles2;

        int l = 0;
        float error = 10;

        while (error > 0.011 && l < 1000) {
            sumDijPj.clear();
            rho_pred.clear();

            for (int j = 0; j < particles.size(); j++) {
                computeSumDijPj(j, dt);
            }

            for (int j = 0; j < particles.size(); j++) {
                updatePressure2(j, dt);
            }

            float rhoAvg = computeAverageRho();
            error = abs(rhoAvg - simProps.desired_density);
            std::cout << "[IISPH2] error: " << error << " nr: " << l << " rhoavg: " << rhoAvg << " desireddensity: " << simProps.desired_density << " aii: " << std::endl;
            std::cout << "----------------" << std::endl;

            l++;
        }
    }

    void CPUSim::integrate(float dt) {
        const auto& p1 = first ? particles1 : particles2;
        auto& p2 = first ? particles2 : particles1;

        for (int i = 0; i < p1.size(); i++) {

            if (p1.at(i).staticP == 0) {

                // Calculate pressure
                glm::vec3 pressure = glm::vec3(0);

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
                }

                pressure *= -p1[i].m;
                glm::vec3 acc = (pressure) / p1[i].m;

                glm::vec3 viNext = p1[i].v;
                glm::vec3 xiNext = p1[i].x;
                viNext += acc * (float) dt + v_advected[i];
                xiNext += viNext * dt;

                p2[i].x = xiNext;
                p2[i].v = viNext;
                p2[i].rho = p1[i].rho;
                p2[i].p = p1[i].p;

                container(i);
            } else {
                p2.at(i) = p1.at(i);
            }
        }

        first = !first;
    }

    void CPUSim::runIISPH2(int iterNr) {
        float dt = 0.01;

        for (int iter = 0; iter < iterNr; iter++) {
            predictAdvection(dt);
            pressureSolve(dt);
            integrate(dt);
        }
    }
}