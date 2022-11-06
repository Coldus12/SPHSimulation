#include "CPUSim.hpp"

#define PI 3.1415926538
#define list_size 25
#define ppic false

namespace Vltava {
#if 0
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

        dii = std::vector<glm::vec3>(p.size(), glm::vec3(0));
        sumDijPj = std::vector<glm::vec3>(p.size(), glm::vec3(0));
        v_adv = std::vector<glm::vec3>(p.size(), glm::vec3(0));

        aii = std::vector<float>(p.size(), 0.0f);
        rho_adv = std::vector<float>(p.size(), 0.0f);
        rho_pred = std::vector<float>(p.size(), 0.0f);

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

    // IISPH time step
    //----------------------------------------------------
    void CPUSim::runIISPH(int iterNr) {
        float dt = 0.01;

        for (int iter = 0; iter < iterNr; iter++) {
            predictAdvection(dt);
            pressureSolve(dt);
            integrate(dt);
        }
    }

    // Predict advection
    //----------------------------------------------------
    void CPUSim::predictAdvection(float dt) {
        calculateRho();
        computeVadvAndDii(dt);
        computeRhoadvAndAii(dt);
    }

    void CPUSim::calculateRho() {
        auto& particles = first ? particles1 : particles2;

        for (int i = 0; i < particles.size(); i++) {

            // Calculating rho
            float rho = 0;

            for (int j = 0; j < particles.size(); j++) {
                if (i == j)
                    continue;

                rho += particles[j].m * kernel(i, j);
            }

            particles[i].rho = rho;
        }
    }

    void CPUSim::computeVadvAndDii(float dt) {
        auto& particles = first ? particles1 : particles2;

        // Going through all particles
        for (int i = 0; i < particles.size(); i++) {
            // Forces
            glm::vec3 f(0);

            // dii
            glm::vec3 dii(0);

            // Calculating viscosity
            glm::vec3 viscosity = glm::vec3(0);

            for (int j = 0; j < particles.size(); j++) {
                if (j == i)
                    continue;

                glm::vec3 xij = particles[i].x - particles[j].x;

                float pval = 0;
                if (particles[j].rho != 0)
                    pval = (particles[j].m / particles[j].rho) * (dot(xij, gradKernel(i, j)) / (dot(xij, xij) + 0.01 * simProps.kernelh));

                glm::vec3 vij = particles[i].v - particles[j].v;
                viscosity += pval * vij;

                dii += -particles[j].m * gradKernel(i,j);
            }

            float nu = 0.01;
            viscosity *= 2 * nu * particles[i].m;

            if (particles[i].staticP == 0)
                f += viscosity;

            f += glm::vec3(0, 0, -9.81 * particles[i].m);

            // Calculating v_adv
            v_adv[i] = particles[i].v + dt * f / particles[i].m;

            // Calculating dii
            dii *= dt * dt / (particles[i].rho * particles[i].rho); // TODO, what happens if rho == 0?
            this->dii[i] = dii;
        }
    }

    void CPUSim::computeRhoadvAndAii(float dt) {
        auto& particles = first ? particles1 : particles2;

        for (int i = 0; i < particles.size(); i++) {
            float rho_adv = 0;
            float aii = 0;

            for (int j = 0; j < particles.size(); j++) {
                if (i == j)
                    continue;

                // dji = ??????
                // dji = mi / rho_i^2 * grad(j,i)?
                //TODO rho == 0?
                glm::vec3 dji = -dt * dt * particles[i].m / (particles[i].rho * particles[i].rho) * gradKernel(j, i);
                aii += particles[j].m * glm::dot(dii[i] - dji, gradKernel(i,j));

                glm::vec3 vij_adv = v_adv[i] - v_adv[j];
                rho_adv += particles[j].m * glm::dot(vij_adv, gradKernel(i,j));
            }

            rho_adv *= dt;
            rho_adv += particles[i].rho;
            this->rho_adv[i] = rho_adv;
            this->aii[i] = aii;

            // pressure init
            particles[i].p *= 0.5;
        }
    }

    // Pressure solve
    //----------------------------------------------------
    void CPUSim::pressureSolve(float dt) {
        int nr = 0;
        float error = 10;

        while (error > 0.01 && nr < 100) {
            computeSumDijPj(dt);
            updatePressure(dt);

            error = abs(calculateAverageError());
            std::cout << "[PressureSolve] nr: " << nr << " error: " << error << std::endl;
            nr++;
        }
        std::cout << "----" << std::endl;
        auto& particles = first ? particles1 : particles2;
        for (int i = 0; i < particles.size(); i++) {
            std::cout << i << " pressure: " << particles[i].p << " rho_pred " << rho_pred[i] << " actual rho " << particles[i].rho << " rho_adv " << rho_adv[i] << " aii " << aii[i] << std::endl;
        }
        std::cout << "--------------" << std::endl;
    }

    void CPUSim::computeSumDijPj(float dt) {
        auto& particles = first ? particles1 : particles2;

        for (int i = 0; i < particles.size(); i++) {
            glm::vec3 dijpj(0);

            for (int j = 0; j < particles.size(); j++) {
                if (i == j)
                    continue;

                // TODO rho == 0?
                dijpj += particles[j].m * particles[j].p / (particles[j].rho * particles[j].rho) * gradKernel(i,j);
            }

            dijpj *= -dt*dt;

            sumDijPj[i] = dijpj;
        }
    }

    void CPUSim::updatePressure(float dt) {
        auto& particles = first ? particles1 : particles2;

        float omega = 0.5f;

        for (int i = 0; i < particles.size(); i++) {
            float sum = 0;

            for (int j = 0; j < particles.size(); j++) {
                if (i == j)
                    continue;

                glm::vec3 dji = -dt*dt*particles[i].m / (particles[i].rho * particles[i].rho) * gradKernel(j, i);
                sum += particles[j].m * glm::dot(sumDijPj[i] - dii[j] * particles[j].p - (sumDijPj[j] - dji * particles[i].p),
                                                gradKernel(i,j));
            }

            rho_pred[i] = rho_adv[i] + particles[i].p * aii[i] + sum;
            particles[i].p = (1 - omega) * particles[i].p + (omega/aii[i]) * (simProps.desired_density - rho_adv[i] - sum);
        }
    }

    float CPUSim::calculateAverageError() {
        float rho_avg = 0;

        for (auto& rho : rho_pred) {
            rho_avg += rho;
        }

        rho_avg /= rho_pred.size();
        return rho_avg - simProps.desired_density;
    }

    // Integrate
    //----------------------------------------------------
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

                    pressure += -val * gradKernel(i, j);
                }

                /*pressure *= -p1[i].m;
                glm::vec3 acc = (pressure) / p1[i].m;*/

                glm::vec3 viNext = v_adv[i];
                glm::vec3 xiNext = p1[i].x;
                viNext += pressure * dt;
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
#endif
}